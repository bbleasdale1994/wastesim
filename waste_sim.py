"""
waste_sim.py - Coupled closed-loop waste recovery model (SimPy, hourly)
=========================================================================

What this models
----------------
Crew -> urine line:  urine hold -> nitrification (dual tanks) --------------\
Crew -> feces line:  feed prep (dilution) -> AD1 acid stage (dual) ->        > recombination
                     AD2 methane stage (dual) -> solids separation ---------/        |
                                                                          pasteurizer
                                                                                |
     cabin CO2 + biogas CO2 -> CO2 store -> photobioreactor (LED lit) <---------/
                                               |  O2 -> O2 store -> crew
                                               |  harvest -> edible algae -> crew food
                                               |           -> inedible + surplus -> back to AD1
                                               |  water recovered -> water store -> AD dilution

SimPy runs the clock and the discrete events (each toilet use, each batch-tank
hold timer, each daily meal). The biology (hydrolysis, methanogenesis, urea
hydrolysis, nitrification, algae growth) is integrated every hour.

Every tank tracks what is IN it (water, solids, methane potential, nitrogen
species, salts), so outputs now depend on tank sizes, residence times, light
and nutrients, rather than being computed straight from crew inputs.

Units: time in hours inside the sim; masses kg; volumes L (density 1 kg/L);
gases kg unless labelled m3 (STP).

Run:  python waste_sim.py          (prints a summary for the default design)
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field, asdict, replace

import numpy as np
import pandas as pd
import simpy

# ---------------------------------------------------------------------------
# Composition vector: what every vessel holds
# ---------------------------------------------------------------------------
W, VSD, VSI, ASH, SALT, CH4P, VFA, UREA, TAN, NO3, ORGN, CHAR = range(12)
# degradable solids + methane potential kept in separate pools, because each breaks down at its own rate:
# VSD/CH4P = feces (and anything else), VSDA/CH4PA = algae, VSDC/CH4PC = cellulose (paper, wipes),
# VSDP/CH4PP = PHA containers.
VSDA, CH4PA, VSDC, CH4PC, VSDP, CH4PP = range(12, 18)
NC = 18
COMP_NAMES = ["water", "vs_degradable", "vs_inert", "ash", "salts",
              "ch4_potential_m3", "vfa_as_ch4_m3", "urea_N", "ammonium_N",
              "nitrate_N", "organic_N", "biochar",
              "vs_degradable_algae", "ch4_potential_algae_m3", "vs_degradable_cellulose",
              "ch4_potential_cellulose_m3", "vs_degradable_pha", "ch4_potential_pha_m3"]
VSD_POOLS = (VSD, VSDA, VSDC, VSDP)
CH4P_POOLS = (CH4P, CH4PA, CH4PC, CH4PP)


CH4_KG_PER_M3 = 0.717       # methane gas density at 0 C, 1 atm (volumes here are at 0 C, 1 atm)
LCH4_KG_PER_L = 0.4226      # liquid methane at its boiling point (-162 C, 1 atm)


def vol_of(c: np.ndarray) -> float:
    """Liquid volume in L (water + suspended/dissolved solids, density ~1)."""
    return float(c[W] + c[VSD] + c[VSDA] + c[VSDC] + c[VSDP] + c[VSI] + c[ASH] + c[SALT] + c[CHAR])


def dissolved_kg(c: np.ndarray) -> float:
    """Approximate dissolved salts + N species mass (for TDS limits)."""
    return float(c[SALT] + c[TAN] * 18 / 14 + c[NO3] * 62 / 14 + c[UREA] * 60 / 28)


def n_total(c: np.ndarray) -> float:
    return float(c[UREA] + c[TAN] + c[NO3] + c[ORGN])


from params import SPECIES, Params   # all assumptions live in params.py


# ---------------------------------------------------------------------------
# Vessels
# ---------------------------------------------------------------------------
class Vessel:
    def __init__(self, name, cap_L, high=0.9, low=0.0, start_frac=0.0):
        self.name, self.cap, self.high, self.low = name, cap_L, high, low
        self.c = np.zeros(NC)
        self.c[W] = start_frac * cap_L

    @property
    def vol(self):
        return vol_of(self.c)

    def room(self):
        return max(0.0, self.high * self.cap - self.vol)

    def add(self, comp):
        self.c += comp

    def take(self, volume):
        v = self.vol
        if v <= 1e-12 or volume <= 0:
            return np.zeros(NC)
        out = self.c * min(1.0, volume / v)
        self.c -= out
        return out


class DualTank:
    """Two batch tanks: fill one to 'high', then it holds for the HRT while the
    other fills; after the hold it drains to 'low' (only if downstream has room)."""

    def __init__(self, name, cap, hrt_h, high, low, drain_L_h, start_frac=0.5):
        self.name = name
        self.t = [Vessel(f"{name}_0", cap, high, low, start_frac),
                  Vessel(f"{name}_1", cap, high, low, start_frac)]
        self.state = ["filling", "idle"]
        self.hold_end = [None, None]
        self.active = 0
        self.hrt_h, self.high, self.low, self.cap = hrt_h, high, low, cap
        self.drain_L_h = drain_L_h
        self.blocked_h = 0.0          # hours a finished batch waited on downstream
        self.residence_h = []         # actual hold+wait time of each batch

    def vessels(self):
        return self.t

    def room(self):
        return sum(v.room() for i, v in enumerate(self.t)
                   if self.state[i] in ("filling", "idle"))

    def _close(self, i, now):
        self.state[i] = "holding"
        self.hold_end[i] = now + self.hrt_h
        self._closed_at = getattr(self, "_closed_at", [None, None])
        self._closed_at[i] = now
        self.active = 1 - i

    def add(self, comp, now):
        """Returns the part that could not be accepted."""
        rem = comp.copy()
        for i in (self.active, 1 - self.active):
            if vol_of(rem) <= 1e-9:
                break
            if self.state[i] not in ("filling", "idle"):
                continue
            self.state[i] = "filling"
            self.active = i
            room = self.t[i].room()
            v = vol_of(rem)
            if room > 1e-9:
                part = rem * min(1.0, room / v)
                self.t[i].add(part)
                rem = rem - part
            if self.t[i].room() <= 1e-6:
                self._close(i, now)
        return rem

    def drain_offer(self, now):
        """Volume this tank wants to send downstream this hour."""
        offer = 0.0
        for i in range(2):
            if self.state[i] == "holding" and now >= self.hold_end[i]:
                self.state[i] = "draining"
            if self.state[i] == "draining":
                offer += min(self.drain_L_h, max(0.0, self.t[i].vol - self.low * self.cap))
        return offer

    def drain(self, now, amount):
        out = np.zeros(NC)
        want = self.drain_offer(now)
        if want > 1e-9 and amount < want - 1e-6:
            self.blocked_h += 1.0 * (1 - amount / want)
        left = amount
        for i in range(2):
            if self.state[i] != "draining" or left <= 1e-9:
                continue
            avail = max(0.0, self.t[i].vol - self.low * self.cap)
            amt = min(self.drain_L_h, avail, left)
            out += self.t[i].take(amt)
            left -= amt
            if self.t[i].vol <= self.low * self.cap + 1e-6:
                self.state[i] = "idle"
                ca = getattr(self, "_closed_at", [None, None])[i]
                if ca is not None:
                    self.residence_h.append(now - ca)
        return out


# ---------------------------------------------------------------------------
# Crew gas profile (OCHMO-TB-004 rates, scaled to the published daily totals)
# ---------------------------------------------------------------------------
def crew_gas_profile(o2_day, co2_day):
    sleep_h = [23, 0, 1, 2, 3, 4, 5, 6]
    ex_h = 17
    o2 = np.zeros(24); co2 = np.zeros(24)
    o2_sleep, co2_sleep = 3.6e-4 * 60, 4.55e-4 * 60          # kg/h
    o2_ex, co2_ex = 39.4e-4 * 60, 49.85e-4 * 60
    for h in sleep_h:
        o2[h], co2[h] = o2_sleep, co2_sleep
    rest = [h for h in range(24) if h not in sleep_h]
    # exercise hour: 0.5 h at 75 % VO2max + 0.5 h nominal
    fixed_o2 = o2_sleep * 8 + 0.5 * o2_ex
    fixed_co2 = co2_sleep * 8 + 0.5 * co2_ex
    nom_o2 = (o2_day - fixed_o2) / (len(rest) - 0.5)
    nom_co2 = (co2_day - fixed_co2) / (len(rest) - 0.5)
    for h in rest:
        o2[h], co2[h] = nom_o2, nom_co2
    o2[ex_h] = 0.5 * o2_ex + 0.5 * nom_o2
    co2[ex_h] = 0.5 * co2_ex + 0.5 * nom_co2
    return o2, co2


# ---------------------------------------------------------------------------
# The system
# ---------------------------------------------------------------------------
class System:
    def __init__(self, env, p: Params):
        self.env, self.p = env, p
        self.sp = p.sp()
        self.rng = np.random.default_rng(p.seed)
        dr = lambda vol: vol * p.drain_frac_per_day / 24.0

        self.urine = Vessel("urine", p.urine_vol, p.fill_high, p.urine_low, 0.3)
        self.urine_draining = False
        self.nit = DualTank("nit", p.nit_vol, p.nit_hrt * 24, p.fill_high, p.tank_low, dr(p.nit_vol), 0.5)
        self.ad1 = DualTank("ad1", p.ad1_vol, p.ad1_hrt * 24, p.fill_high, p.tank_low, dr(p.ad1_vol), 0.5)
        self.ad2 = DualTank("ad2", p.ad2_vol, p.ad2_hrt * 24, p.fill_high, p.tank_low, dr(p.ad2_vol), 0.5)
        self.rec = Vessel("rec", p.rec_vol, p.fill_high, p.rec_low, 0.3)
        self.feed_buf = Vessel("feed_buffer", p.feed_buffer_L, 1.0, 0.0, 0.0)
        self.pbr = Vessel("pbr", p.pbr_vol_L, 1.0, 0.0, 1.0)
        self.pbr.c[TAN] = p.pbr_n_set * p.pbr_vol_L / 1000.0
        self.biomass = p.x_set * p.pbr_vol_L / 1000.0      # kg dry algae in culture

        # seed nitrification tanks with partly nitrified urine, AD with digestate
        for v in self.nit.vessels():
            v.c[TAN] = v.c[NO3] = 0.0025 * v.vol
            v.c[SALT] = 0.012 * v.vol
        self.rec.c[TAN] = p.startup_tan_g_L / 1000 * self.rec.vol
        for v in self.ad1.vessels() + self.ad2.vessels():
            v.c[TAN] = p.startup_tan_g_L / 1000 * v.vol

        # stores and accumulators
        self.o2 = p.o2_store_init_kg
        self.co2 = p.co2_store_init_kg
        self.water = max(p.water_store_init_L, p.water_reserve_L)   # start with the reserve full
        self.food = 0.0                    # kg dry edible algae
        self.char = 0.0                    # kg biochar in store
        self.ad_feed_queue = np.zeros(NC)  # algae waiting to go to feed prep
        self.inbox_urine = np.zeros(NC)
        self.inbox_feces = np.zeros(NC)
        self.n_exc_g = p.diet_protein_frac * p.diet_kcal / 4 / 6.25   # g N per CM-day
        self.algae_intake_g = 0.0
        self.o2_h, self.co2_h = crew_gas_profile(p.crew_o2_kg_cmd, p.crew_co2_kg_cmd)
        # consumables: today's AD-bound share is released into the feed over 24 hourly portions
        self.cons_today = np.zeros(NC)
        self.cn_prev = 0.0                       # digester-feed C/N over the previous day
        self.day_c = self.day_n = self.day_c_cons = 0.0
        self.prev_c_other = self.prev_n = 0.0
        self.cum_ch4p = {"feces": 0.0, "algae": 0.0, "cellulose": 0.0, "pha": 0.0}

        self.cum = dict(ch4_m3=0.0, o2_algae=0.0, o2_demand=0.0, o2_deficit=0.0,
                        o2_vented=0.0, co2_vented=0.0, co2_fixed=0.0, biomass=0.0,
                        food_eaten=0.0, food_to_ad=0.0, overflow=0.0, urine_overflow=0.0,
                        n_lost=0.0, sludge=0.0, base_kg=0.0, water_short=0.0,
                        energy_kWh=0.0, fert_L=0.0, fert_N=0.0, makeup_water=0.0,
                        n_excreted=0.0, n_overflow=0.0, n_eaten=0.0,
                        flow_ad1=0.0, flow_ad2=0.0, flow_nit=0.0, flow_urine=0.0, flow_rec=0.0,
                        water_out=0.0, rxn_water=0.0, cake_water_lost=0.0, brine_water_lost=0.0,
                        brine_salt=0.0, cake_dry=0.0, char_cake=0.0, char_waste=0.0, char_dosed=0.0,
                        pyro_volatiles=0.0, waste_water_in=0.0, crew_water_in=0.0,
                        fert_water_lost=0.0, fert_solid=0.0,
                        cons_ad_kg=0.0, cons_pyro_kg=0.0, char_cons=0.0, cons_water_in=0.0, te_kg=0.0,
                        cons_dry_kg=0.0, cn_capped_kg=0.0, food_dryer_lost=0.0,
                        ch4p_in_feces=0.0, ch4p_in_char=0.0, overflow_water=0.0)
        self.log = []

    # ---- helpers that build streams --------------------------------------
    def urine_event(self, mass):
        p = self.p
        c = np.zeros(NC)
        un_kg = p.urine_n_frac * self.n_exc_g / 1000.0 * mass / p.urine_kg_cmd
        c[UREA] = un_kg
        solids = p.urine_solids_kg_cmd * mass / p.urine_kg_cmd
        c[SALT] = max(solids - un_kg * 60 / 28, 0.25 * solids)
        c[W] = mass - c[SALT]
        return c

    def feces_event(self, mass):
        p = self.p
        c = np.zeros(NC)
        ts = mass * p.feces_ts_kg_cmd / p.feces_wet_kg_cmd
        vs = ts * p.feces_vs_of_ts
        c[VSD] = vs * p.feces_vs_degradable
        c[VSI] = vs - c[VSD]
        c[ASH] = ts - vs
        c[CH4P] = p.b0_feces * vs
        c[ORGN] = (1 - p.urine_n_frac) * self.n_exc_g / 1000.0 * mass / p.feces_wet_kg_cmd
        c[W] = mass - ts
        return c

    def algae_stream(self, dry_kg, water_kg):
        p = self.p
        c = np.zeros(NC)
        vs = dry_kg * p.algae_vs_of_ts
        c[VSDA] = vs * p.algae_vs_degradable
        c[VSI] = vs - c[VSDA]
        c[ASH] = dry_kg - vs
        c[CH4PA] = p.b0_algae * vs
        c[ORGN] = dry_kg * self.sp["n_frac"]
        c[W] = water_kg
        return c

    def consumable_stream(self, kind, dry_kg):
        """Toilet paper ('tp'), wet wipes ('wipe') or PHA containers ('pha') as AD feed.
        Modelled on algae_stream: water, degradable/inert VS, ash and methane potential; no nitrogen."""
        p = self.p
        c = np.zeros(NC)
        if dry_kg <= 0:
            return c
        wf = getattr(p, f"{kind}_water_frac")
        vs = dry_kg * getattr(p, f"{kind}_vs_of_ts")
        deg = vs * getattr(p, f"{kind}_vs_degradable")
        vsd, ch4p = (VSDP, CH4PP) if kind == "pha" else (VSDC, CH4PC)
        c[vsd] = deg
        c[VSI] = vs - deg
        c[ASH] = dry_kg - vs
        c[ch4p] = (p.b0_pha if kind == "pha" else p.b0_cellulose) * vs
        c[W] = dry_kg * wf / (1 - wf)
        return c

    def c_of(self, c):
        """Carbon (kg) in a stream, from each solids pool's carbon content (for the feed C/N)."""
        p = self.p
        return float(c[VSD] * p.feces_c_frac / p.feces_vs_of_ts
                     + c[VSDA] * self.sp.get("c_frac", 0.46) / p.algae_vs_of_ts
                     + c[VSDC] * p.tp_c_frac / p.tp_vs_of_ts + c[VSDP] * p.pha_c_frac / p.pha_vs_of_ts
                     + c[VSI] * 0.5 + c[UREA] * 12 / 28)

    def daily_consumables(self, e):
        """Once a day: the crew's paper, wipes and PHA containers. PHA always goes to AD. Cellulose goes
        to AD up to the feed C/N cap (judged on yesterday's feed); the rest is dried and charred."""
        p = self.p
        cr = p.crew
        dry = {k: getattr(p, f"{k}_kg_cmd") * cr for k in ("tp", "wipe", "pha")}
        to_ad = {k: dry[k] * getattr(p, f"{k}_to_ad") for k in dry}
        if self.prev_n > 0 and p.ad_cn_max > 0:
            room_c = p.ad_cn_max * self.prev_n - self.prev_c_other - to_ad["pha"] * p.pha_c_frac
            cell_c = to_ad["tp"] * p.tp_c_frac + to_ad["wipe"] * p.wipe_c_frac
            if cell_c > max(room_c, 0.0):
                keep = max(room_c, 0.0) / cell_c
                for k in ("tp", "wipe"):
                    self.cum["cn_capped_kg"] += to_ad[k] * (1 - keep)
                    to_ad[k] *= keep
        s = np.zeros(NC)
        for k in dry:
            s += self.consumable_stream(k, to_ad[k])
            self.cum["cons_ad_kg"] += to_ad[k]
            e["shred"] += to_ad[k] * (p.shred_kWh_per_kg_pha if k == "pha" else p.shred_kWh_per_kg_cell)
            self.cum["cons_dry_kg"] += dry[k]
            wf = getattr(p, f"{k}_water_frac")
            self.cum["cons_water_in"] += dry[k] * wf / (1 - wf)
            # the rest: dried (water back to the loop) and charred
            rest = dry[k] - to_ad[k]
            if rest > 1e-12:
                ww = rest * wf / (1 - wf)
                rec_w = ww * p.cake_dryer_recovery
                self.water += rec_w
                self.cum["cake_water_lost"] += ww - rec_w
                e["cake"] += ww * p.dryer_kWh_per_kg_water
                ch = rest * (p.pha_char_yield if k == "pha" else p.tp_char_yield)
                self.char += ch
                self.cum["char_cons"] += ch
                self.cum["cons_pyro_kg"] += rest
                self.cum["pyro_volatiles"] += rest - ch
                e["pyro"] += rest * p.pyro_kWh_per_kg_dry
        self.cons_today = s

    # ---- SimPy processes --------------------------------------------------
    def crew_member(self, i):
        """Toilet events for one crew member (random times, fixed daily mass)."""
        p, rng = self.p, self.rng
        wake = [7, 10, 13, 16, 19, 22][: p.urinations_per_day]
        while True:
            day0 = self.env.now - (self.env.now % 24)
            events = []
            stagger = 0.25 * (i % 4)          # crew don't all go at once (repeats every 4 people)
            for h in wake:
                events.append((day0 + h + rng.uniform(-0.5, 0.5) + stagger, "u",
                               p.urine_kg_cmd / len(wake)))
            n_def = 1 + (rng.random() < (p.defecations_per_day - 1))
            for k in range(n_def):
                events.append((day0 + 8 + 10 * k + rng.uniform(0, 2), "f",
                               p.feces_wet_kg_cmd / p.defecations_per_day))
            # keep every event inside today so the wait until tomorrow is never negative
            events = [(min(t, day0 + 23.95), kind, m) for t, kind, m in events]
            for t, kind, m in sorted(events):
                if t > self.env.now:
                    yield self.env.timeout(t - self.env.now)
                ev = self.urine_event(m) if kind == "u" else self.feces_event(m)
                self.cum["n_excreted"] += n_total(ev)
                self.cum["crew_water_in"] += ev[W]
                if kind == "u":
                    self.inbox_urine += ev
                else:
                    self.inbox_feces += ev
            yield self.env.timeout(max(0.0, day0 + 24 - self.env.now))

    def meals(self):
        """Once a day the crew eats algae up to the safety caps; surplus goes to AD."""
        p, sp = self.p, self.sp
        yield self.env.timeout(12.5)
        while True:
            cap_na = p.algae_na_cap_g / sp["na_frac"]
            base_prot = p.diet_protein_frac * p.diet_kcal / 4
            cap_prot = ((p.protein_cap_g_per_kg * p.body_mass_kg - base_prot)
                        / max(1e-6, sp["protein"] - p.diet_protein_frac * sp["kcal_per_g"] / 4))
            want_g = max(0.0, min(cap_na, cap_prot))
            have_g = self.food * 1000 / p.crew
            a = min(want_g, have_g)
            self.food -= a * p.crew / 1000
            self.cum["food_eaten"] += a * p.crew / 1000
            self.cum["n_eaten"] += a * p.crew / 1000 * sp["n_frac"]
            self.algae_intake_g = a
            protein = a * sp["protein"] + p.diet_protein_frac * (p.diet_kcal - a * sp["kcal_per_g"]) / 4
            self.n_exc_g = protein / 6.25
            # surplus beyond the buffer -> AD as rehydrated paste
            buffer_kg = p.food_buffer_days * p.crew * want_g / 1000
            if p.algae_surplus_to_ad and self.food > buffer_kg:
                extra = self.food - buffer_kg
                self.food = buffer_kg
                self.ad_feed_queue += self.algae_stream(extra, 0.0)
                self.cum["food_to_ad"] += extra
            yield self.env.timeout(24)

    def plant(self):
        while True:
            self.step()
            yield self.env.timeout(1)

    # ---- one hour of plant operation ---------------------------------------
    def step(self):
        p, sp = self.p, self.sp
        now = self.env.now
        hod = int(now % 24)
        e = dict(led=0.0, aer=0.0, ad_heat=0.0, past=0.0, memb=0.0, dryer=0.0, mix=0.0, comp=0.0, strip=0.0,
                 cake=0.0, pyro=0.0, shred=0.0)
        if hod == 0:                     # day boundary: close yesterday's feed C/N, release today's consumables
            if self.day_n > 0:
                self.cn_prev = self.day_c / self.day_n
            self.prev_c_other, self.prev_n = self.day_c - self.day_c_cons, self.day_n
            self.day_c = self.day_n = self.day_c_cons = 0.0
            self.daily_consumables(e)

        # 1. toilet inputs -------------------------------------------------
        u = self.inbox_urine; self.inbox_urine = np.zeros(NC)
        uv = vol_of(u)
        if uv > 0:
            room = max(0.0, self.urine.cap - self.urine.vol)   # can fill to 100 % in emergency
            acc = min(1.0, room / uv)
            self.urine.add(u * acc)
            self.cum["flow_urine"] += uv * acc
            if acc < 1:
                self.cum["n_overflow"] += n_total(u) * (1 - acc)
                self.cum["urine_overflow"] += uv * (1 - acc)
                self.cum["overflow"] += uv * (1 - acc)
                self.cum["overflow_water"] += u[W] * (1 - acc)

        # feed prep: feces + queued algae, diluted to target TS with stored water
        cons = self.cons_today / 24.0
        feed = self.inbox_feces + self.ad_feed_queue + cons
        self.inbox_feces = np.zeros(NC); self.ad_feed_queue = np.zeros(NC)
        c_cons = self.c_of(cons)
        self.cum["ch4p_in_feces"] += feed[CH4P]          # only feces carry this pool into the feed
        self.day_c += self.c_of(feed); self.day_c_cons += c_cons; self.day_n += n_total(feed)
        ts = feed[VSD] + feed[VSDA] + feed[VSDC] + feed[VSDP] + feed[VSI] + feed[ASH] + feed[SALT]
        if ts > 0:
            need = ts / p.ad_target_ts - vol_of(feed)
            if need > 0:
                give = min(need, max(0.0, self.water))
                feed[W] += give; self.water -= give
                self.cum["water_short"] += need - give
            self.feed_buf.add(feed)
            te = p.te_fe_mg_L + p.te_ni_mg_L + p.te_co_mg_L + p.te_se_mo_w_mg_L
            self.cum["te_kg"] += te * vol_of(feed) * 1e-6
        if self.feed_buf.vol > 0:
            v_in = self.feed_buf.vol
            rej = self.ad1.add(self.feed_buf.c.copy(), now)
            self.cum["flow_ad1"] += v_in - vol_of(rej)
            self.feed_buf.c = rej
            if self.feed_buf.vol > p.feed_buffer_L:
                ov = self.feed_buf.take(self.feed_buf.vol - p.feed_buffer_L)
                self.cum["overflow"] += vol_of(ov)
                self.cum["overflow_water"] += ov[W]
                self.cum["n_overflow"] += n_total(ov)

        # 2. urine hold -> nitrification ---------------------------------------
        if self.urine.vol >= p.fill_high * self.urine.cap:
            self.urine_draining = True
        if self.urine_draining:
            amt = min(self.urine.cap * p.drain_frac_per_day / 24,
                      self.urine.vol - p.urine_low * self.urine.cap, self.nit.room())
            if amt > 0:
                rej = self.nit.add(self.urine.take(amt), now)
                self.urine.add(rej)
                self.cum["flow_nit"] += amt - vol_of(rej)
            if self.urine.vol <= p.urine_low * self.urine.cap + 1e-6:
                self.urine_draining = False

        # 3. biology in tanks ---------------------------------------------------
        fu = 1 - math.exp(-p.k_urea_d / 24)
        for v in [self.urine, self.rec] + self.nit.vessels():
            d = v.c[UREA] * fu; v.c[UREA] -= d; v.c[TAN] += d
            v.c[W] -= 0.643 * d; self.cum["rxn_water"] -= 0.643 * d     # urea + H2O -> 2 NH3 + CO2
        o2_nit = 0.0
        for v in self.nit.vessels():
            if v.vol <= 0:
                continue
            allowed = p.nit_max_frac * (v.c[TAN] + v.c[NO3]) - v.c[NO3]
            r = max(0.0, min(v.c[TAN], allowed, p.nit_rate_gN_L_d / 1000 / 24 * v.vol))
            v.c[TAN] -= r
            v.c[NO3] += r * (1 - p.n2o_frac)
            self.cum["n_lost"] += r * p.n2o_frac
            o2_nit += 4.57 * r
            v.c[W] += 1.286 * r; self.cum["rxn_water"] += 1.286 * r       # NH4+ + 2 O2 -> NO3- + H2O + 2 H+
            if p.nit_max_frac > 0.5:
                self.cum["base_kg"] += 12.0 * r * (p.nit_max_frac - 0.5) / p.nit_max_frac
        e["aer"] = o2_nit / p.aeration_kgO2_per_kWh   # kWh this hour

        # keep the biochar dose in every digester vessel topped up from the char store
        if p.char_dose_g_L > 0:
            for v in self.ad1.vessels() + self.ad2.vessels():
                want = p.char_dose_g_L / 1000 * v.vol - v.c[CHAR]
                if want > 0 and self.char > 0:
                    add = min(want, self.char)
                    v.c[CHAR] += add; self.char -= add; self.cum["char_dosed"] += add
        ch4 = 0.0
        self._ad_state = []
        for stage, tank in ((1, self.ad1), (2, self.ad2)):
            for v in tank.vessels():
                if v.vol <= 0:
                    continue
                cg = v.c[CHAR] * 1000 / max(v.vol, 1e-6)                    # biochar g/L
                f = p.char_effect * cg / (cg + p.char_half_g_L) if cg > 0 else 0.0
                hb = 1 + p.char_max_hyd_boost * f
                codig = p.codig_kh_boost * p.codig_effect * min(1.0, max(0.0, (self.cn_prev - 6) / 14))
                kh_pha = p.kh_pha_d if self.env.now >= p.pha_lag_d * 24 else 0.0
                fa0 = v.c[VSD] + v.c[VSDA]
                dv_fa = 0.0
                # (pool, potential, rate, water taken up per kg VS hydrolysed (Buswell), source)
                for vsd, ch4p, kh, wtr, src in ((VSD, CH4P, p.kh_d, 0.35, "feces"),
                                                (VSDA, CH4PA, p.kh_d * (1 + codig), 0.35, "algae"),
                                                (VSDC, CH4PC, p.kh_cellulose_d, 0.111, "cellulose"),
                                                (VSDP, CH4PP, kh_pha, 0.314, "pha")):
                    if v.c[vsd] <= 0 and v.c[ch4p] <= 0:
                        continue
                    fh = 1 - math.exp(-kh * hb / 24)
                    dv = v.c[vsd] * fh; dp = v.c[ch4p] * fh
                    v.c[vsd] -= dv; v.c[ch4p] -= dp; v.c[VFA] += dp
                    v.c[W] -= wtr * dv; self.cum["rxn_water"] -= wtr * dv    # hydrolysis takes up water
                    self.cum_ch4p[src] += dp
                    if vsd in (VSD, VSDA):
                        dv_fa += dv
                # organic N (feces + algae protein) is released as their solids break down
                dn = v.c[ORGN] * (dv_fa / fa0) if fa0 > 1e-12 else 0.0
                v.c[ORGN] -= dn; v.c[TAN] += dn
                if f > 0 and v.c[VSI] > 0:                                   # biochar opens up 'inert' VS
                    di = v.c[VSI] * (1 - math.exp(-p.kh_d * p.char_inert_degrade * f / 24))
                    v.c[VSI] -= di; v.c[VSD] += di; v.c[CH4P] += p.b0_algae * di
                    self.cum["ch4p_in_char"] += p.b0_algae * di      # biochar bonus shares the feces pool
                free_tan = max(0.0, v.c[TAN] - p.char_n_sorb * v.c[CHAR])   # ammonium held on char
                tan_gl = free_tan * 1000 / max(v.vol, 1e-6)
                ki = p.ki_tan * (1 + p.char_max_ki_boost * f)
                inhib = 1.0 / (1.0 + tan_gl / ki)
                k = (p.k_meth_d * (1 + p.char_max_meth_boost * f) / 24 * inhib
                     * (p.acid_stage_meth_frac if stage == 1 else 1.0))
                dm = v.c[VFA] * (1 - math.exp(-k))
                v.c[VFA] -= dm
                ch4 += dm
                if stage == 2:
                    self._ad_state.append((tan_gl, ki, cg))
        # side-stream ammonia stripping: AD2 ammonium -> recombination (stays in the N loop)
        if p.ad_strip_per_d > 0:
            fs = 1 - math.exp(-p.ad_strip_per_d / 24)
            stripped = 0.0
            for v in self.ad2.vessels():
                dn = v.c[TAN] * fs; v.c[TAN] -= dn; stripped += dn
            self.rec.c[TAN] += stripped
            e["strip"] = stripped * p.strip_kWh_per_kgN
        biogas_co2 = ch4 * (1 - p.ch4_frac_biogas) / p.ch4_frac_biogas * 1.964
        self.cum["ch4_m3"] += ch4
        e["comp"] = ch4 * p.ch4_compress_kWh_per_m3
        ad_vol = sum(v.vol for v in self.ad1.vessels() + self.ad2.vessels())
        e["ad_heat"] = ad_vol * p.ad_loss_W_per_L / 1000
        e["mix"] = (ad_vol + sum(v.vol for v in self.nit.vessels()) + self.rec.vol) * p.mixing_W_per_L / 1000

        # 4. transfers between stages -------------------------------------------
        if vol_of(feed) > 0:   # heat the new feed to digester temperature
            e["ad_heat"] += vol_of(feed) * 4.186 * (p.ad_temp_c - p.cabin_temp_c) * (1 - p.ad_heat_recovery) / 3600
        out = self.ad1.drain(now, min(self.ad1.drain_offer(now), self.ad2.room()))
        if vol_of(out) > 0:
            self.cum["flow_ad2"] += vol_of(out)
            rej = self.ad2.add(out, now)
            self.ad1.t[0].add(rej * 0.5); self.ad1.t[1].add(rej * 0.5)

        off_n, off_a = self.nit.drain_offer(now), self.ad2.drain_offer(now)
        room = self.rec.room()
        scale = min(1.0, room / max(off_n + off_a, 1e-9))
        nit_out = self.nit.drain(now, off_n * scale)
        ad_out = self.ad2.drain(now, off_a * scale)
        if vol_of(ad_out) > 0:       # dewater digestate: solids to sludge cake
            cap = p.solids_capture
            solids = ad_out[[VSD, VSDA, VSDC, VSDP, VSI, ASH]].sum() * cap
            cake = np.zeros(NC)
            for idx in (VSD, VSI, ASH, CH4P, ORGN, CHAR) + (VSDA, CH4PA, VSDC, CH4PC, VSDP, CH4PP):
                cake[idx] = ad_out[idx] * cap
            solids += ad_out[CHAR] * cap
            cake_w = min(ad_out[W], solids * (1 / p.cake_ts - 1))
            cake[W] = cake_w
            ad_out = ad_out - cake
            self.cum["sludge"] += vol_of(cake)
            self.cum["n_lost"] += cake[ORGN]
            # dry the cake (water back to the internal loop), then char it
            rec_w = cake_w * p.cake_dryer_recovery
            self.water += rec_w
            self.cum["cake_water_lost"] += cake_w - rec_w
            e["cake"] += cake_w * p.dryer_kWh_per_kg_water
            dry = cake[VSD] + cake[VSDA] + cake[VSDC] + cake[VSDP] + cake[VSI] + cake[ASH] + cake[CHAR]
            new_char = ((cake[VSD] + cake[VSDA] + cake[VSI]) * p.vs_char_yield + cake[VSDC] * p.tp_char_yield
                        + cake[VSDP] * p.pha_char_yield + cake[ASH])
            old_char = cake[CHAR] * p.char_recycle_yield
            self.char += new_char + old_char
            self.cum["char_cake"] += new_char
            self.cum["cake_dry"] += dry
            self.cum["pyro_volatiles"] += dry - new_char - old_char
            e["pyro"] += dry * p.pyro_kWh_per_kg_dry
        self.rec.add(nit_out + ad_out)
        self.cum["flow_rec"] += vol_of(nit_out + ad_out)
        # recover clean water from the recombined liquor (keeps nutrients, returns water)
        # recover water down to rec_min_frac (limited by how concentrated the membrane can go);
        # the store keeps its reserve and passes any surplus to the crew
        floor_L = max(p.rec_min_frac * self.rec.cap, dissolved_kg(self.rec.c) * 1000 / p.rec_tds_max)
        if self.rec.vol > floor_L:
            xs = min(self.rec.vol - floor_L, self.rec.c[W])
            self.rec.c[W] -= xs
            self.water += xs
            e["memb"] += xs / 1000 * p.membrane_kWh_per_m3
        if self.rec.vol > p.rec_export_frac * self.rec.cap:      # liquid fertilizer product
            fert = self.rec.take(self.rec.vol - (p.rec_export_frac - 0.05) * self.rec.cap)
            self.cum["fert_L"] += vol_of(fert)
            self.cum["fert_N"] += n_total(fert)
            # dry it to solid fertilizer salts; water goes back to the loop
            rw = fert[W] * p.cake_dryer_recovery
            self.water += rw
            self.cum["fert_water_lost"] += fert[W] - rw
            self.cum["fert_solid"] += vol_of(fert) - fert[W] + n_total(fert) * 2.5
            e["cake"] += fert[W] * p.dryer_kWh_per_kg_water

        # extra crew waste: dried (water back to the loop) and charred
        if p.extra_waste_kg_cmd > 0:
            wet = p.extra_waste_kg_cmd * p.crew / 24
            ww = wet * p.extra_waste_water_frac
            rec_w = ww * p.cake_dryer_recovery
            self.water += rec_w
            self.cum["waste_water_in"] += ww
            self.cum["cake_water_lost"] += ww - rec_w
            e["cake"] += ww * p.dryer_kWh_per_kg_water
            dry = wet - ww
            ch = dry * p.extra_waste_char_yield
            self.char += ch
            self.cum["char_waste"] += ch
            self.cum["pyro_volatiles"] += dry - ch
            e["pyro"] += dry * p.pyro_kWh_per_kg_dry

        # 5. photobioreactor ---------------------------------------------------
        light_on = hod < p.photoperiod_h
        B = self.biomass
        V = max(self.pbr.vol, 1e-6)
        over_lit = 0.0
        if light_on:
            ppf = p.led_kw * 1000 * p.led_ppe * p.delivery()              # umol/s into the culture
            geo = pbr_geometry(p)
            x_gm3 = B * 1e6 / V
            if p.pbr_shape == "tube":       # outside wall and inside rods each have their own light limit
                fin = p.inside_share()
                d_in = p.internal_delivery * (1 - p.internal_fouling)
                ppf_in = p.led_kw * 1000 * p.led_ppe * fin * d_in
                ppf_out = ppf - ppf_in
                e_out = min(ppf_out, p.i_max * geo["area_out"])
                e_in = min(ppf_in, p.i_max * geo["area_in"])
                ppf_eff = e_out + e_in
                abs_out = 1 - math.exp(-sp["k_abs"] * x_gm3 * geo["path_out"])
                abs_in = 1 - math.exp(-sp["k_abs"] * x_gm3 * geo["path_in"])
                absf = (e_out * abs_out + e_in * abs_in) / max(ppf_eff, 1e-9)
                # photons used, after light saturation at each surface's own flux
                q_out = e_out / max(geo["area_out"], 1e-9)
                q_in = e_in / max(geo["area_in"], 1e-9)
                used = (e_out * abs_out * p.photon_efficiency(q_out) + e_in * abs_in * p.photon_efficiency(q_in))
            else:
                ppf_eff = min(ppf, p.i_max * geo["area"])
                absf = 1 - math.exp(-sp["k_abs"] * x_gm3 * p.light_path_m)
                used = ppf_eff * absf * p.photon_efficiency(ppf_eff / max(geo["area"], 1e-9))
            over_lit = ppf - ppf_eff
            mol = ppf_eff * absf * 3600e-6
            y0, K = p.light_curve()
            if K is None:                                                    # old linear model (Pirt maintenance)
                gross = y0 * (used * 3600e-6 - sp["maint"] * B * 1000)      # g/h
            else:
                # tier yields are net of light-phase upkeep (fitted to measured productivity), so only
                # culture outside the lit layer (the dark core) respires
                d90 = math.log(10) / (sp["k_abs"] * max(x_gm3, 1e-9))
                dark = 1 - min(1.0, geo["a_light"] * d90)
                gross = y0 * used * 3600e-6 - p.r_dark * B * 1000 * dark    # g/h
            e["led"] = p.led_kw
            radiant_kw = p.led_kw * p.led_ppe * p.photon_J_per_umol            # light leaving the LEDs
            absorbed_kw = ppf_eff * absf * p.photon_J_per_umol / 1000         # light absorbed by culture
            inside_loss_kw = p.inside_share() * (p.led_kw - radiant_kw)        # inside-LED waste heat
        else:
            absorbed_kw = inside_loss_kw = 0.0
            gross = -p.r_dark * B * 1000
        limit = "light" if light_on else "dark"
        stored_kw = 0.0
        o2_alg = 0.0
        if gross > 0:
            n_av = self.pbr.c[TAN] + self.pbr.c[UREA] + self.pbr.c[NO3]
            n_need = gross / 1000 * sp["n_frac"]
            co2_need = gross / 1000 * sp["co2_per_g"]
            s = 1.0
            if n_need > n_av:
                s = n_av / n_need; limit = "N"
            if co2_need * s > self.co2:
                s = self.co2 / co2_need; limit = "CO2"
            g = gross * s / 1000                                            # kg
            nn = g * sp["n_frac"]
            take_nh = min(nn, self.pbr.c[TAN] + self.pbr.c[UREA])
            fr = take_nh / max(self.pbr.c[TAN] + self.pbr.c[UREA], 1e-12)
            self.pbr.c[TAN] -= self.pbr.c[TAN] * fr; self.pbr.c[UREA] -= self.pbr.c[UREA] * fr
            take_no3 = min(nn - take_nh, self.pbr.c[NO3])
            self.pbr.c[NO3] -= take_no3
            pq = (take_nh * p.pq_nh4 + take_no3 * p.pq_no3) / max(nn, 1e-12)
            co2_f = g * sp["co2_per_g"]
            self.co2 -= co2_f
            o2_alg = co2_f * 32 / 44 * pq
            # photosynthesis takes up water (CO2 + H2O + N -> biomass + O2)
            nh_share = take_nh / max(nn, 1e-12)
            wu = co2_f * (0.245 * nh_share + 0.327 * (1 - nh_share))
            self.pbr.c[W] -= wu; self.cum["rxn_water"] -= wu
            stored_kw = g * 1000 * p.biomass_kJ_per_g / 3600
            self.biomass += g
            self.cum["co2_fixed"] += co2_f
            self.cum["biomass"] += g
        elif gross < 0:
            g = -gross / 1000
            g = min(g, self.biomass)
            self.biomass -= g
            o2_alg = -g * sp["co2_per_g"] * 32 / 44
            self.co2 += g * sp["co2_per_g"]
            self.pbr.c[TAN] += g * sp["n_frac"]
            wp = g * sp["co2_per_g"] * 0.245                                  # respiration gives water back
            self.pbr.c[W] += wp; self.cum["rxn_water"] += wp
            stored_kw = -g * 1000 * p.biomass_kJ_per_g / 3600

        # harvest to hold density at x_set
        X = self.biomass * 1000 / V
        if X > p.x_set:
            Vh = V * (X - p.x_set) / X
            hb = self.biomass * Vh / V
            self.biomass -= hb
            liq = self.pbr.take(Vh)
            paste_w = min(liq[W], hb * (1 / p.paste_solids - 1))
            paste_frac = paste_w / max(liq[W], 1e-12)
            paste_part = liq * paste_frac           # dissolved stuff carried in paste
            centrate = liq - paste_part
            self.pbr.add(centrate)                   # spent medium recycled to the culture
            edible = hb * p.harvest_edible_frac
            self.food += edible
            w_dry = paste_w * p.harvest_edible_frac
            self.water += w_dry * p.dryer_recovery
            self.cum["food_dryer_lost"] += w_dry * (1 - p.dryer_recovery)   # vented by the food dryer
            e["dryer"] += w_dry * p.dryer_kWh_per_kg_water
            paste_part[W] -= w_dry                  # edible paste water went to the dryer
            inedible = self.algae_stream(hb - edible, 0.0) + paste_part
            if p.algae_surplus_to_ad:
                self.ad_feed_queue += inedible
            else:
                self.cum["n_lost"] += n_total(inedible)

        # make-up feed from recombination (volume + nitrogen demand)
        V = self.pbr.vol
        deficit = max(0.0, p.pbr_vol_L - V)
        n_now = self.pbr.c[TAN] + self.pbr.c[UREA] + self.pbr.c[NO3]
        extra = 0.0
        rec_n = (self.rec.c[TAN] + self.rec.c[UREA] + self.rec.c[NO3]) / max(self.rec.vol, 1e-9)
        if n_now < p.pbr_n_set * p.pbr_vol_L / 1000 and rec_n > 1e-9:
            extra = (p.pbr_n_set * p.pbr_vol_L / 1000 - n_now) / rec_n
            extra = min(extra, 0.05 * p.pbr_vol_L)
        avail = max(0.0, self.rec.vol - p.rec_low * self.rec.cap)
        n_room = max(0.0, p.pbr_n_max * p.pbr_vol_L / 1000 - n_now)
        fv = min(deficit + extra, avail, n_room / rec_n if rec_n > 1e-9 else avail)
        clean = max(0.0, deficit - fv)                 # rest of make-up from clean water
        clean = min(clean, max(0.0, self.water))
        if clean > 0:
            self.water -= clean
            self.pbr.c[W] += clean
            self.cum["makeup_water"] += clean
        o2_org = 0.0
        if fv > 0:
            f = self.rec.take(fv)
            e["past"] += fv * 4.186 * (p.past_temp_c - p.cabin_temp_c) * (1 - p.past_heat_recovery) / 3600
            # residual organics oxidised by bacteria in the PBR (uses algae O2)
            o2_org = (f[VFA] * 2.857 + (f[VSD] + f[VSDA]) * 1.42 + f[VSDC] * 1.185 + f[VSDP] * 1.67
                      + f[CH4P] * 0)
            self.co2 += o2_org * 44 / 32
            f[W] += 0.5625 * o2_org; self.cum["rxn_water"] += 0.5625 * o2_org   # oxidation makes water
            for k in (VFA,) + VSD_POOLS + CH4P_POOLS:
                f[k] = 0
            n_min = f[ORGN]; f[ORGN] = 0; f[TAN] += n_min
            self.pbr.add(f)
        # membrane extraction of excess water back to the water store
        if self.pbr.vol > p.pbr_vol_L:
            xs = self.pbr.vol - p.pbr_vol_L
            self.pbr.c[W] -= xs
            self.water += xs
            e["memb"] += xs / 1000 * p.membrane_kWh_per_m3
        # inert solids that reached the culture leave with the harvest -> AD
        inert = self.pbr.c[[VSI, ASH, CHAR]].sum()
        if inert > 0:
            s = np.zeros(NC)
            for idx in (VSI, ASH, CHAR):
                s[idx] = self.pbr.c[idx]; self.pbr.c[idx] = 0
            self.ad_feed_queue += s
        # salt blowdown through brine processor
        tds = self.pbr.c[SALT] * 1000 / max(self.pbr.vol, 1e-9)
        if tds > p.tds_max:
            vb = self.pbr.vol * (tds - p.tds_max) / tds
            br = self.pbr.take(vb)
            self.water += br[W] * p.brine_recovery
            self.cum["brine_water_lost"] += br[W] * (1 - p.brine_recovery)
            self.cum["brine_salt"] += br[SALT]
            e["memb"] += br[W] * p.brine_recovery / 1000 * p.membrane_kWh_per_m3 * 3
            self.cum["n_lost"] += br[UREA] + br[TAN] + br[NO3]

        # 6. gases ----------------------------------------------------------------
        crew_o2 = p.crew * self.o2_h[hod]
        crew_co2 = p.crew * self.co2_h[hod]
        self.co2 += crew_co2 + biogas_co2 * p.biogas_co2_capture
        if self.co2 > p.co2_store_max_kg:
            self.cum["co2_vented"] += self.co2 - p.co2_store_max_kg
            self.co2 = p.co2_store_max_kg
        o2_net = o2_alg - o2_org
        demand = crew_o2 + o2_nit
        self.o2 += o2_net - demand
        if self.o2 < 0:
            self.cum["o2_deficit"] += -self.o2
            self.o2 = 0.0
        if self.o2 > p.o2_store_max_kg:
            self.cum["o2_vented"] += self.o2 - p.o2_store_max_kg
            self.o2 = p.o2_store_max_kg
        self.cum["o2_algae"] += o2_net
        self.cum["o2_demand"] += demand

        # purified water: keep the internal reserve, send the surplus to the crew
        if self.water > p.water_reserve_L:
            self.cum["water_out"] += self.water - p.water_reserve_L
            self.water = p.water_reserve_L

        kwh = sum(e.values())
        self.cum["energy_kWh"] += kwh
        # heat to reject: all electricity ends up as heat except what algae store as chemical energy
        cool_total = max(0.0, kwh - stored_kw)
        cool_culture = max(0.0, absorbed_kw - stored_kw + inside_loss_kw)
        ad_tan, ad_ki, ad_char = (np.mean([t[i] for t in self._ad_state]) if self._ad_state else 0.0
                                  for i in range(3))

        # 7. log --------------------------------------------------------------------
        X = self.biomass * 1000 / max(self.pbr.vol, 1e-9)
        self.log.append((
            now / 24.0,
            self.urine.vol, self.nit.t[0].vol, self.nit.t[1].vol,
            self.ad1.t[0].vol, self.ad1.t[1].vol, self.ad2.t[0].vol, self.ad2.t[1].vol,
            self.rec.vol, self.pbr.vol, self.feed_buf.vol,
            X, (self.pbr.c[TAN] + self.pbr.c[UREA] + self.pbr.c[NO3]) * 1000 / max(self.pbr.vol, 1e-9),
            tds, limit, over_lit,
            self.o2, self.co2, self.water, self.food,
            ch4, o2_net, demand, crew_co2 + biogas_co2 * p.biogas_co2_capture,
            kwh, e["led"], e["aer"], e["ad_heat"], e["past"], e["memb"], e["dryer"], e["mix"], e["comp"], e["strip"],
            e["cake"], e["pyro"], cool_total, cool_culture,
            self.algae_intake_g, self.char, ad_tan, ad_ki, ad_char,
            self.cum["water_out"], self.cum["water_short"], self.cum["char_cake"], self.cum["char_waste"],
            self.cum["rxn_water"],
            self.cum["ch4_m3"], self.cum["o2_algae"], self.cum["o2_demand"], self.cum["o2_deficit"],
            self.cum["biomass"], self.cum["food_eaten"], self.cum["overflow"], self.cum["urine_overflow"],
            self.cn_prev, self.cum_ch4p["feces"], self.cum_ch4p["algae"], self.cum_ch4p["cellulose"],
            self.cum_ch4p["pha"], self.cum["cons_ad_kg"], self.cum["cons_pyro_kg"], self.cum["char_cons"],
            e["shred"], self.cum["te_kg"],
        ))


LOG_COLS = ["day", "urine", "nit_0", "nit_1", "ad1_0", "ad1_1", "ad2_0", "ad2_1",
            "rec", "pbr", "feed_buffer",
            "algae_g_L", "pbr_N_g_L", "pbr_tds_g_L", "limit", "overlit_umol_s",
            "o2_store_kg", "co2_store_kg", "water_store_L", "food_store_kg",
            "ch4_m3_h", "o2_algae_kg_h", "o2_demand_kg_h", "co2_in_kg_h",
            "kW", "kW_led", "kW_aeration", "kW_ad_heat", "kW_pasteurize", "kW_membrane",
            "kW_dryer", "kW_mixing", "kW_compress", "kW_strip", "kW_cake_dryer", "kW_pyrolysis",
            "cool_kW", "cool_culture_kW",
            "algae_intake_g_cmd", "char_store_kg", "ad2_free_tan_g_L", "ad2_ki_g_L", "ad2_char_g_L",
            "cum_water_out_L", "cum_water_short_L", "cum_char_cake_kg", "cum_char_waste_kg", "cum_rxn_water_kg",
            "cum_ch4_m3", "cum_o2_algae_kg", "cum_o2_demand_kg", "cum_o2_deficit_kg",
            "cum_biomass_kg", "cum_food_eaten_kg", "cum_overflow_L", "cum_urine_overflow_L",
            "ad_feed_cn", "cum_ch4p_feces", "cum_ch4p_algae", "cum_ch4p_cellulose", "cum_ch4p_pha",
            "cum_cons_ad_kg", "cum_cons_pyro_kg", "cum_char_cons_kg", "kW_shred", "cum_te_kg"]


# ---------------------------------------------------------------------------
# Run + summarise
# ---------------------------------------------------------------------------
def run_sim(p: Params, return_system=False):
    env = simpy.Environment()
    s = System(env, p)
    n0 = n_inventory(s)
    for i in range(p.crew):
        env.process(s.crew_member(i))
    env.process(s.meals())
    env.process(s.plant())
    env.run(until=p.days * 24)
    df = pd.DataFrame(s.log, columns=LOG_COLS)
    res = summarize(p, s, df)
    res.update(n_audit(s, n0))
    return (df, res, s) if return_system else (df, res)


def daily(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()
    d["d"] = d["day"].astype(int)
    agg = {c: "last" for c in df.columns if c not in ("day", "limit")}
    for c in ["ch4_m3_h", "o2_algae_kg_h", "o2_demand_kg_h", "co2_in_kg_h"]:
        agg[c] = "sum"
    for c in [c for c in df.columns if c.startswith("kW") or c.startswith("cool")]:
        agg[c] = "mean"
    out = d.groupby("d").agg(agg).reset_index().rename(columns={
        "d": "day", "ch4_m3_h": "ch4_m3_d", "o2_algae_kg_h": "o2_algae_kg_d",
        "o2_demand_kg_h": "o2_demand_kg_d", "co2_in_kg_h": "co2_in_kg_d"})
    out["o2_coverage"] = out["o2_algae_kg_d"] / out["o2_demand_kg_d"]
    lim = d.groupby("d")["limit"].agg(lambda s: s[s != "dark"].mode().iat[0] if (s != "dark").any() else "dark")
    out["limit"] = lim.values
    return out


def summarize(p: Params, s: System, df: pd.DataFrame) -> dict:
    n = p.steady_days * 24
    w = df.iloc[-n:]
    days = p.steady_days
    first, last = df.iloc[-n - 1] if len(df) > n else df.iloc[0], df.iloc[-1]
    o2_prod = w["o2_algae_kg_h"].sum()
    o2_dem = w["o2_demand_kg_h"].sum()
    deficit = last["cum_o2_deficit_kg"] - first["cum_o2_deficit_kg"]
    kw = {c: w[c].mean() for c in df.columns if c.startswith("kW")}
    res = dict(
        o2_coverage=o2_prod / o2_dem,
        o2_deficit_kg=deficit,
        o2_algae_kg_d=o2_prod / days,
        o2_demand_kg_d=o2_dem / days,
        ch4_m3_d=w["ch4_m3_h"].sum() / days,
        ch4_kg_d=w["ch4_m3_h"].sum() / days * CH4_KG_PER_M3,
        # running total over the whole run (from start-up), and at the end of the mission
        ch4_total_m3=float(df["cum_ch4_m3"].iloc[-1]),
        ch4_total_kg=float(df["cum_ch4_m3"].iloc[-1]) * CH4_KG_PER_M3,
        ch4_total_liquid_L=float(df["cum_ch4_m3"].iloc[-1]) * CH4_KG_PER_M3 / LCH4_KG_PER_L,
        ch4_mission_m3=(float(df.loc[df["day"] <= p.mission_days, "cum_ch4_m3"].iloc[-1])
                        if p.days >= p.mission_days else float("nan")),
        biomass_kg_d=(last["cum_biomass_kg"] - first["cum_biomass_kg"]) / days,
        algae_eaten_g_cmd=(last["cum_food_eaten_kg"] - first["cum_food_eaten_kg"]) / days / p.crew * 1000,
        overflow_L=df["cum_overflow_L"].iloc[-1],
        overflow_steady_L=last["cum_overflow_L"] - first["cum_overflow_L"],
        **{f"{k}_L_d": v / p.days for k, v in s.cum.items() if k.startswith("flow_")},
        urine_overflow_L=s.cum["urine_overflow"],
        urine_overflow_steady_L=last["cum_urine_overflow_L"] - first["cum_urine_overflow_L"],
        light_limited_frac=(w["limit"] == "light").sum() / max(1, (w["limit"] != "dark").sum()),
        n_limited_frac=(w["limit"] == "N").sum() / max(1, (w["limit"] != "dark").sum()),
        co2_limited_frac=(w["limit"] == "CO2").sum() / max(1, (w["limit"] != "dark").sum()),
        overlit=w["overlit_umol_s"].mean() > 1,
        water_net_L_d=(last["water_store_L"] - first["water_store_L"]) / days,
        water_short_L=s.cum["water_short"],
        co2_vented_kg=s.cum["co2_vented"],
        n_lost_g_d=s.cum["n_lost"] / p.days * 1000,
        fert_L_d=s.cum["fert_L"] / p.days,
        fert_N_g_d=s.cum["fert_N"] / p.days * 1000,
        fert_solid_g_d=s.cum["fert_solid"] / p.days * 1000,
        fert_water_lost_L_d=s.cum["fert_water_lost"] / p.days,
        sludge_kg_d=s.cum["sludge"] / p.days,
        base_kg=s.cum["base_kg"],
        ad1_blocked_h=s.ad1.blocked_h, ad2_blocked_h=s.ad2.blocked_h, nit_blocked_h=s.nit.blocked_h,
        ad2_mean_residence_d=(np.mean(s.ad2.residence_h) / 24) if s.ad2.residence_h else float("nan"),
        pbr_tan_ad2_g_L=float(np.mean([v.c[TAN] * 1000 / max(v.vol, 1e-9) for v in s.ad2.vessels()])),
        ad2_free_tan_g_L=float(w["ad2_free_tan_g_L"].mean()),
        ad2_ki_g_L=float(w["ad2_ki_g_L"].mean()),
        ad2_char_g_L=float(w["ad2_char_g_L"].mean()),
        **{f"avg_{k}": v for k, v in kw.items()},
        avg_cool_kW=float(w["cool_kW"].mean()),
        avg_cool_culture_kW=float(w["cool_culture_kW"].mean()),
        # water
        water_out_L_d=(last["cum_water_out_L"] - first["cum_water_out_L"]) / days,
        water_short_steady_L=last["cum_water_short_L"] - first["cum_water_short_L"],
        rxn_water_L_d=(last["cum_rxn_water_kg"] - first["cum_rxn_water_kg"]) / days,
        crew_water_in_L_d=s.cum["crew_water_in"] / p.days,
        waste_water_in_L_d=s.cum["waste_water_in"] / p.days,
        cake_water_lost_L_d=s.cum["cake_water_lost"] / p.days,
        brine_water_lost_L_d=s.cum["brine_water_lost"] / p.days,
        food_dryer_water_lost_L_d=s.cum["food_dryer_lost"] / p.days,
        brine_salt_g_d=s.cum["brine_salt"] / p.days * 1000,
        # biochar
        char_cake_g_d=(last["cum_char_cake_kg"] - first["cum_char_cake_kg"]) / days * 1000,
        char_waste_g_d=(last["cum_char_waste_kg"] - first["cum_char_waste_kg"]) / days * 1000,
        char_store_kg=float(last["char_store_kg"]),
        char_store_growth_g_d=(last["char_store_kg"] - first["char_store_kg"]) / days * 1000,
        cake_dry_kg_d=s.cum["cake_dry"] / p.days,
        pyro_volatiles_g_d=s.cum["pyro_volatiles"] / p.days * 1000,
        fecal_solids_g_d=p.feces_ts_kg_cmd * p.crew * 1000,
        # crew consumables (paper, wipes, PHA containers) and the digester feed
        ad_feed_cn=float(w["ad_feed_cn"].mean()),
        cons_ad_kg_d=(last["cum_cons_ad_kg"] - first["cum_cons_ad_kg"]) / days,
        cons_pyro_kg_d=(last["cum_cons_pyro_kg"] - first["cum_cons_pyro_kg"]) / days,
        char_cons_g_d=(last["cum_char_cons_kg"] - first["cum_char_cons_kg"]) / days * 1000,
        cons_water_in_L_d=s.cum["cons_water_in"] / p.days,
        cn_capped_kg_d=s.cum["cn_capped_kg"] / p.days,
        te_g_d=(last["cum_te_kg"] - first["cum_te_kg"]) / days * 1000,
        **ch4_by_source(w, first, last, days, s),
        # peak liquid on board (all tanks + culture + water store) over the steady window
        peak_liquid_L=float((w[["urine", "nit_0", "nit_1", "ad1_0", "ad1_1", "ad2_0", "ad2_1", "rec", "pbr",
                               "feed_buffer", "water_store_L"]].sum(axis=1)).max()),
    )
    res.update(esm(p, s, res))
    return res


def ch4_by_source(w, first, last, days, s=None) -> dict:
    """Methane split by where it came from: total CH4 shared out by the methane potential each pool
    released in the steady window (feces, algae, cellulose = paper + wipes, PHA containers, and the
    extra potential biochar opens up in 'inert' solids). Biochar's extra potential is added to the
    feces pool, so that pool's release is split by how much potential each put in over the run."""
    rel = {k: last[f"cum_ch4p_{k}"] - first[f"cum_ch4p_{k}"] for k in ("feces", "algae", "cellulose", "pha")}
    tot = sum(rel.values())
    ch4 = w["ch4_m3_h"].sum() / days
    out = {f"ch4_{k}_m3_d": ch4 * v / tot if tot > 0 else 0.0 for k, v in rel.items()}
    fe, ch = (s.cum["ch4p_in_feces"], s.cum["ch4p_in_char"]) if s is not None else (1.0, 0.0)
    char_share = ch / (fe + ch) if fe + ch > 0 else 0.0
    out["ch4_biochar_m3_d"] = out["ch4_feces_m3_d"] * char_share
    out["ch4_feces_m3_d"] *= 1 - char_share
    return out


def n_inventory(s: System) -> float:
    vs = [s.urine, s.rec, s.feed_buf, s.pbr] + s.nit.vessels() + s.ad1.vessels() + s.ad2.vessels()
    inv = sum(n_total(v.c) for v in vs)
    inv += n_total(s.ad_feed_queue) + n_total(s.inbox_urine) + n_total(s.inbox_feces)
    inv += (s.biomass + s.food) * s.sp["n_frac"]
    return inv


def n_audit(s: System, n0: float) -> dict:
    """Nitrogen bookkeeping: start + excreted = end + eaten + fertilizer + losses."""
    c = s.cum
    lhs = n0 + c["n_excreted"]
    rhs = n_inventory(s) + c["n_eaten"] + c["fert_N"] + c["n_lost"] + c["n_overflow"]
    return dict(n_start=n0, n_in=c["n_excreted"], n_end=n_inventory(s), n_eaten=c["n_eaten"],
                n_fert=c["fert_N"], n_lost=c["n_lost"], n_overflow=c["n_overflow"],
                n_error_frac=(lhs - rhs) / max(lhs, 1e-9))


def esm(p: Params, s: System, res: dict) -> dict:
    tanks_L = p.urine_vol + 2 * p.nit_vol + 2 * p.ad1_vol + 2 * p.ad2_vol + p.rec_vol + p.feed_buffer_L
    g = pbr_geometry(p)
    area, wall = g["area"], g["wall"]
    hw = (p.tank_kg_per_L * tanks_L + p.pbr_kg_per_m2 * wall + p.led_kg_per_kw * p.led_kw
          + p.internal_kg_per_m2 * g["area_in"])
    inventory = res["peak_liquid_L"]            # the water you have to launch = peak liquid on board
    vol_m3 = p.packing_factor * tanks_L / 1000 + g["bay_m3"]
    kw, cool = res["avg_kW"], res["avg_cool_kW"]
    # O2 beyond the share other systems are meant to supply counts as a consumable
    o2_allowed = (1 - p.o2_target) * res["o2_demand_kg_d"]
    o2_extra = max(0.0, res["o2_deficit_kg"] / p.steady_days - o2_allowed)
    # crew paper, wipes and PHA containers replace items launched anyway unless told otherwise
    cons_items = 0.0 if p.consumables_replace_existing else s.cum["cons_dry_kg"] / p.days
    consumables = (res["base_kg"] / p.days + o2_extra + res["te_g_d"] / 1000 + cons_items
                   + res["water_short_steady_L"] / p.steady_days) * p.mission_days
    power_cool = kw * p.peq_kg_per_kw + cool * p.ceq_kg_per_kw
    total = hw + inventory + vol_m3 * p.veq_kg_per_m3 + power_cool + consumables
    return dict(esm_kg=total, esm_hardware_kg=hw, esm_inventory_kg=inventory,
                esm_volume_kg=vol_m3 * p.veq_kg_per_m3,
                esm_power_cooling_kg=power_cool,
                esm_power_kg=kw * p.peq_kg_per_kw, esm_cooling_kg=cool * p.ceq_kg_per_kw,
                esm_consumables_kg=consumables, tanks_L=tanks_L,
                pbr_area_m2=area, pbr_wall_area_m2=wall, pressurized_m3=vol_m3,
                pbr_bay_m3=g["bay_m3"], n_tubes=g["n_tubes"], area_per_L_cm2=area / max(p.pbr_vol_L, 1e-9) * 1e4,
                **light_depth(p))


def pbr_geometry(p: Params) -> dict:
    """Lit surfaces, light paths and space taken by the culture, for panels or vertical tubes.
    area_out / area_in: lit surface (m2) facing the outside LEDs / the inside LEDs.
    wall: reactor wall area that has to be built (m2). bay_m3: pressurized space of the reactor bay."""
    V = p.pbr_vol_L / 1000
    if p.pbr_shape == "tube":
        D, H, d, k = p.tube_d_m, p.tube_h_m, p.rod_d_m, p.rods_per_tube
        check_rods(p)
        a_free = math.pi / 4 * (D ** 2 - k * d ** 2)                   # culture cross-section, m2
        v_tube = a_free * H
        n = V / v_tube
        wall = n * math.pi * D * H
        a_out = wall * p.tube_light_coverage
        a_in = n * k * math.pi * d * H
        # mean distance light crosses before reaching another lit surface (mean chord, Cauchy):
        # pi * culture cross-section / lit perimeter
        lit_perim = math.pi * D * p.tube_light_coverage + k * math.pi * d
        path = math.pi * a_free / max(lit_perim, 1e-9)
        return dict(area=a_out + a_in, area_out=a_out, area_in=a_in, wall=wall,
                    path_out=path, path_in=path if k else 0.0,
                    bay_m3=n * (D + p.tube_gap_m) ** 2 * H, n_tubes=n, v_tube_L=v_tube * 1000,
                    a_light=(a_out + a_in) / max(V, 1e-12))
    wall = V / p.light_path_m
    area = wall * p.lit_area_factor
    f = p.internal_light_frac
    return dict(area=area, area_out=area * (1 - f), area_in=area * f, wall=wall,
                path_out=p.light_path_m, path_in=p.light_path_m,
                bay_m3=wall * p.panel_pitch_m, n_tubes=0.0, v_tube_L=0.0, a_light=area / max(V, 1e-12))


def max_rods(p: Params) -> int:
    """Most rods of diameter rod_d_m that fit in the tube on a hexagonal grid with rod_gap_m of culture
    between neighbouring rods and between the rods and the tube wall."""
    D, d, g = p.tube_d_m, p.rod_d_m, p.rod_gap_m
    pitch = d + g
    usable = D - 2 * g - d                     # circle the rod centres must stay inside
    if usable < 0:
        return 0
    if usable < pitch:
        return 1
    return int(math.pi / 4 * usable ** 2 / (math.sqrt(3) / 2 * pitch ** 2)) + 1


def check_rods(p: Params):
    """Refuse tube layouts that cannot be built."""
    if p.tube_d_m <= 0 or p.tube_h_m <= 0:
        raise ValueError("tube diameter and height must be positive")
    if p.rods_per_tube < 0:
        raise ValueError("rods_per_tube cannot be negative")
    if p.rods_per_tube > 0:
        if p.rod_d_m <= 0:
            raise ValueError("rod diameter must be positive")
        k_max = max_rods(p)
        if p.rods_per_tube > k_max:
            raise ValueError(f"{p.rods_per_tube} rods of {p.rod_d_m*100:.1f} cm do not fit in a "
                             f"{p.tube_d_m*100:.1f} cm tube with {p.rod_gap_m*100:.1f} cm gaps "
                             f"(at most {k_max})")


def pbr_wall_area_m2(p: Params) -> float:
    return pbr_geometry(p)["wall"]


def pbr_area_m2(p: Params) -> float:
    """Lit culture surface (m2)."""
    return pbr_geometry(p)["area"]


def light_depth(p: Params) -> dict:
    """How deep light gets into the culture (90 % absorbed) and the share of culture inside that lit
    layer at any moment. The rest is a dark core: algae there only grow when mixing carries them
    to the wall, but they still respire (the model charges maintenance on all biomass)."""
    d90 = math.log(10) / (p.sp()["k_abs"] * p.x_set * 1000)          # m
    if p.pbr_shape == "tube":
        g = pbr_geometry(p)
        frac = min(1.0, g["a_light"] * d90)          # lit layer d90 deep on every lit surface
    else:
        frac = min(1.0, p.lit_area_factor * d90 / p.light_path_m)
    return dict(light_depth_mm=d90 * 1000, lit_volume_frac=frac)


def pbr_vol_for_led(led_kw: float, p: Params, margin: float = 1.05) -> float:
    """Culture volume whose lit surface takes the LED light at the design flux (design_pfd, capped at i_max)."""
    ppf = led_kw * 1000 * p.led_ppe * p.delivery()
    if p.pbr_shape == "tube":
        g1 = pbr_geometry(replace(p, pbr_vol_L=1000.0))                # surfaces per m3 of culture
        fin = p.inside_share()
        ppf_in = led_kw * 1000 * p.led_ppe * fin * p.internal_delivery * (1 - p.internal_fouling)
        q = min(p.design_pfd, p.i_max)
        need = (ppf - ppf_in) / (q * g1["area_out"])
        if g1["area_in"] > 0:
            need = max(need, ppf_in / (q * g1["area_in"]))
        return need * margin * 1000
    return ppf / min(p.design_pfd, p.i_max) * margin * p.light_path_m * 1000 / p.lit_area_factor


if __name__ == "__main__":
    import time
    t0 = time.time()
    df, r = run_sim(Params())
    print(f"ran in {time.time() - t0:.1f} s")
    for k, v in r.items():
        print(f"{k:28s} {v:.4g}" if isinstance(v, float) else f"{k:28s} {v}")
