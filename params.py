"""
params.py - every assumption in the closed-loop model, in one place
====================================================================
Change a number here (or with a slider in the app) to change the model.
Comments give the source or say when a value is a rough placeholder.
"""
from dataclasses import dataclass


# ---------------------------------------------------------------------------
# Algae species data
# ---------------------------------------------------------------------------
SPECIES = {
    # yx: biomass yield on absorbed light (g dry / mol photons). Wageningen panel-PBR
    #     work: 0.75-0.78 g/mol (Chlorella / Dunaliella, Zijffers 2010) and 1.0 g/mol
    #     (Chlorella, Cuaresma 2009). Arthrospira reported higher; 0.8 is a design value.
    # maint: Pirt maintenance (mol photons / g / h), 0.0068 for C. sorokiniana.
    # n_frac: N in dry biomass = protein / 6.25.
    # na_frac: nucleic acid content of dry biomass (Spirulina 4.8-5.7 %).
    # c_frac: carbon in dry biomass (typical elemental analysis), used for the digester-feed C/N.
    "Spirulina": dict(yx=0.80, maint=0.0068, n_frac=0.575 / 6.25, co2_per_g=1.85,
                      protein=0.575, na_frac=0.05, kcal_per_g=2.90, k_abs=0.15, c_frac=0.46),
    "Chlorella": dict(yx=0.78, maint=0.0068, n_frac=0.50 / 6.25, co2_per_g=1.85,
                      protein=0.50, na_frac=0.05, kcal_per_g=3.50, k_abs=0.15, c_frac=0.50),
}
SPECIES["50/50 mix"] = {k: 0.5 * (SPECIES["Spirulina"][k] + SPECIES["Chlorella"][k])
                        for k in SPECIES["Spirulina"]}


# ---------------------------------------------------------------------------
# Light response: how much algae each photon makes, by technology level
# ---------------------------------------------------------------------------
# Biomass made per m2 of lit surface = y0 * K * ln(1 + q/K), q = photon flux on that surface
# (Cornet's light-limited growth law, arXiv:2011.05320). At low light every photon gives y0 g/mol;
# as light rises the yield per photon falls (the culture saturates). K sets how much light a strain
# can use well. Yield per photon at flux q = y0 * (K/q) * ln(1 + q/K). Respiration (maint) comes off
# on top, so dark culture and high density still cost biomass.
#   y0 in g dry biomass per mol photons absorbed at low light; K in umol/m2/s.
LIGHT_TIERS = {
    # old model: constant yield up to i_max (no saturation). Kept so older designs can be compared.
    "linear": dict(y0=None, k_sat=None,
                   label="Old model: no saturation (optimistic at high light)"),
    # Arthrospira, fitted to 31 published reactor runs (+/-15 %): rho_M 0.8, phi 2.4e-9 kg/umol on NH4,
    # alpha 0.9, K 90 (Cornet 2020).
    "spirulina_cornet": dict(y0=1.82, k_sat=90.0,
                             label="Spirulina, published fit (saturates early)"),
    # same on nitrate (phi 1.85e-9 kg/umol): the conditions of the published runs; used to validate the code
    "spirulina_cornet_no3": dict(y0=1.40, k_sat=90.0,
                                 label="Spirulina on nitrate, published fit (validation only)"),
    # Chlorella vulgaris as run long-term in life-support hardware: ISS-candidate membrane raceway
    # (Helisch et al. 2020, 188 days, up to 12.2 g/L, max 1.29 g/L/d) and MELiSSA-type airlifts.
    # y0 1.0 g/mol is ~55 % of the theoretical 1.8; K 300 = C. vulgaris saturating at moderate light.
    # Gives ~1-1.5 g/L/d in a 1.5-2 cm deep culture at 300-500 umol/m2/s (est. calibration).
    "iss": dict(y0=1.0, k_sat=300.0,
                label="ISS-type: Chlorella as run for months in life-support hardware"),
    # well-engineered middle: 0.6 g/mol at 1,000 umol/m2/s
    "engineered": dict(y0=1.2, k_sat=450.0,
                       label="Engineered: tuned strain and reactor (middle of published range)"),
    # best lab result: Chlorella sorokiniana, 14 mm panel, 2,100 umol/m2/s red LED, 1.0 g/mol,
    # 0.5 g/L/h (Cuaresma et al. 2009); y0 1.8 = theoretical maximum (8 photons per O2).
    "best_lab": dict(y0=1.8, k_sat=1080.0,
                     label="Best lab record (C. sorokiniana, short runs)"),
}


# ---------------------------------------------------------------------------
# Parameters (every number you might want to change lives here)
# ---------------------------------------------------------------------------
@dataclass
class Params:
    # --- mission -----------------------------------------------------------
    crew: int = 4
    o2_target: float = 1.0          # share of crew O2 this system must supply (rest from other ECLSS)
    days: int = 365
    steady_days: int = 120          # results averaged over the last N days
    mission_days: int = 250         # mission length; scales consumables in ESM
    seed: int = 1

    # --- crew metabolic interface (NASA BVAD 2022 / HIDH, see spreadsheet) ---
    urine_kg_cmd: float = 1.562
    urine_solids_kg_cmd: float = 0.061
    feces_wet_kg_cmd: float = 0.123
    feces_ts_kg_cmd: float = 0.032
    feces_vs_of_ts: float = 0.84
    feces_vs_degradable: float = 0.75   # share of fecal VS that can digest
    crew_o2_kg_cmd: float = 0.82        # with exercise
    crew_co2_kg_cmd: float = 1.04
    urinations_per_day: int = 6
    defecations_per_day: float = 1.2    # Rose et al. 2015

    # --- diet & the "crew eats the algae" loop ------------------------------
    diet_kcal: float = 3037.0           # BVAD exploration energy target
    diet_protein_frac: float = 0.15     # of energy, for the stored-food part
    urine_n_frac: float = 0.86          # share of excreted N in urine
    body_mass_kg: float = 82.0
    algae_na_cap_g: float = 2.0         # max nucleic acid from microbial food, g/CM-d (UN PAG)
    protein_cap_g_per_kg: float = 2.0   # upper protein intake, g/kg body mass/day
    food_buffer_days: float = 14.0      # algae food kept in store; the rest goes to AD
    algae_surplus_to_ad: bool = True    # recycle surplus/inedible algae through AD

    # --- algae / photobioreactor -------------------------------------------
    species: str = "Chlorella"
    light_tier: str = "iss"             # light response, see LIGHT_TIERS
    design_pfd: float = 1500.0          # photon flux the reactor is sized for on its lit surfaces, umol/m2/s
                                        # (<= i_max; lower = more culture and lit area per LED kW)
    pbr_vol_L: float = 1400.0
    light_path_m: float = 0.03          # panel thickness (short light path)
    led_kw: float = 12.0                # electrical LED power
    led_ppe: float = 3.0                # fixture efficacy, umol/J (3.4 limit -10 % for waterproofing)
    light_delivery: float = 0.85        # outside LEDs: fraction of photons that get through the wall
    i_max: float = 1500.0               # max useful photon flux per m2 of lit culture surface, umol/m2/s
    # lighting layout: 1 = panels lit from one side, 2 = both sides, >2 = LEDs inside the culture too
    lit_area_factor: float = 1.0        # lit culture surface / panel face area
    internal_light_frac: float = 0.0    # share of LED power mounted inside the culture
    internal_delivery: float = 0.95     # inside LEDs: no wall / air gap, small sleeve loss
    internal_fouling: float = 0.05      # light lost to biofilm on inside emitters (needs cleaning)
    internal_kg_per_m2: float = 3.0     # rough: sealed LED rods / sleeves per m2 of lit surface they serve
    # reactor shape: "panel" = flat panels (light_path_m thick); "tube" = vertical tubes floor to ceiling
    pbr_shape: str = "panel"
    tube_d_m: float = 0.305             # tube inside diameter (12 in)
    tube_h_m: float = 2.4               # tube height, floor to ceiling (assumed 8 ft deck)
    tube_gap_m: float = 0.10            # space between tubes for LED strips and access
    tube_light_coverage: float = 1.0    # share of each tube's outside surface that LEDs light (1 = all round)
    rods_per_tube: int = 0              # LED rods running down inside each tube
    rod_d_m: float = 0.05               # rod (sealed sleeve) outside diameter
    rod_gap_m: float = 0.015            # smallest culture gap between rods and to the tube wall (mixing, cleaning)
    photoperiod_h: int = 24
    r_dark: float = 0.002               # dark respiration, 1/h
    x_set: float = 2.0                  # culture density kept by harvesting, g/L
    pbr_n_set: float = 0.10             # dissolved N kept in culture (floor), g N/L
    pbr_n_max: float = 0.30             # ceiling: above this make-up uses clean water, g N/L
    tds_max: float = 15.0               # salt limit in culture, g/L (blowdown above)
    pq_nh4: float = 1.05                # O2/CO2 molar ratio on ammonium
    pq_no3: float = 1.30                # O2/CO2 molar ratio on nitrate
    harvest_edible_frac: float = 0.70
    paste_solids: float = 0.15          # centrifuge paste dry fraction
    dryer_recovery: float = 0.95        # water recovered when drying food paste
    brine_recovery: float = 0.90        # water recovered from salt blowdown

    # --- anaerobic digestion -------------------------------------------------
    ad_target_ts: float = 0.06          # feed diluted to this total-solids fraction
    kh_d: float = 0.20                  # first-order hydrolysis, 1/d (mesophilic)
    k_meth_d: float = 1.0               # VFA -> CH4 rate in AD2, 1/d
    acid_stage_meth_frac: float = 0.15  # AD1 methanogenesis relative to AD2
    b0_feces: float = 0.40              # m3 CH4 / kg VS added (ultimate)
    b0_algae: float = 0.35              # Spirulina continuous AD 0.30-0.35
    algae_vs_of_ts: float = 0.92
    algae_vs_degradable: float = 0.80
    ch4_frac_biogas: float = 0.62
    ki_tan: float = 3.5                 # ammonium-N giving 50 % methanogen inhibition, g/L
    solids_capture: float = 0.90        # digestate solids removed before recombination
    cake_ts: float = 0.20
    biogas_co2_capture: float = 0.95
    ad_strip_per_d: float = 0.0         # side-stream NH3 stripping from AD2, fraction of TAN/day
    strip_kWh_per_kgN: float = 5.0      # rough: membrane-contactor / air stripping incl. heat
    feed_buffer_L: float = 40.0         # space to hold feed if AD1 is full

    # --- digestate cake -> dryer -> pyrolysis -> biochar -------------------------
    cake_dryer_recovery: float = 0.90   # water recovered when drying cake / waste before charring
    pyro_kWh_per_kg_dry: float = 0.61   # 2.18 MJ/kg dry sewage sludge to 570 C (Mendoza-Geney 2022)
    vs_char_yield: float = 0.23         # char from organic matter; with ash gives ~35 % char from raw
                                        # feces, as measured (Yacob et al. 2018: 35.1-35.8 % dry basis)
    char_recycle_yield: float = 0.95    # biochar that comes round again survives re-charring
    char_dose_g_L: float = 0.0          # biochar kept in the digesters (literature optimum ~10 g/L)
    char_effect: float = 1.0            # 0 = biochar does nothing; 1 = mid-literature effects
    char_half_g_L: float = 5.0          # dose giving half of the maximum effect
    char_max_meth_boost: float = 1.0    # methanogen rate up to +100 % (studies: +30 to +130 %, shorter lag)
    char_max_hyd_boost: float = 0.15    # hydrolysis up to +15 %
    char_max_ki_boost: float = 0.5      # ammonia tolerance up to +50 %
    char_n_sorb: float = 0.05           # kg ammonium-N held per kg char (Cai et al. 2022 avg ~50 mg/g)
    char_inert_degrade: float = 0.20    # 'inert' VS opened up to digestion -> extra methane yield
    extra_waste_kg_cmd: float = 0.0     # OTHER crew trash charred directly, wet kg/CM-d: cotton clothing,
                                        # towels, food scraps. Not toilet paper, wipes or PHA containers
                                        # (those are the crew consumables below and go to the digesters)
    extra_waste_water_frac: float = 0.50
    extra_waste_char_yield: float = 0.25   # char per kg dry mixed waste

    # --- urine / nitrification ----------------------------------------------
    k_urea_d: float = 1.0               # urea hydrolysis, 1/d
    nit_rate_gN_L_d: float = 0.42       # pilot MBBR on urine (Etter et al. 2013)
    nit_max_frac: float = 0.50          # 0.5 without base dosing (urine alkalinity limit)
    n2o_frac: float = 0.005
    aeration_kgO2_per_kWh: float = 1.5

    # --- tanks (L) and residence times (days) ------------------------------
    urine_vol: float = 200.0
    nit_vol: float = 175.0
    ad1_vol: float = 50.0
    ad2_vol: float = 150.0
    rec_vol: float = 1000.0
    nit_hrt: float = 7.0
    ad1_hrt: float = 4.0
    ad2_hrt: float = 15.0
    fill_high: float = 0.90
    tank_low: float = 0.30              # dual tanks drain down to this
    urine_low: float = 0.30
    rec_low: float = 0.05
    rec_min_frac: float = 0.15          # membrane recovers water from the liquor down to this level
    rec_target_frac: float = 0.50       # (no longer used; kept so older design files still load)
    rec_tds_max: float = 60.0           # membrane can concentrate liquor up to this, g/L
    rec_export_frac: float = 0.85       # above this level liquor is drawn off as fertilizer
    startup_tan_g_L: float = 5.0        # nutrient charge in AD inoculum / recombination at start
    drain_frac_per_day: float = 0.50    # pump rate, fraction of tank volume per day

    # --- stores ---------------------------------------------------------------
    co2_store_max_kg: float = 30.0
    o2_store_max_kg: float = 30.0
    o2_store_init_kg: float = 15.0
    co2_store_init_kg: float = 3.0
    water_store_init_L: float = 200.0
    water_reserve_L: float = 300.0      # internal water kept for the loop; surplus goes to the crew

    # --- energy ---------------------------------------------------------------
    ad_temp_c: float = 35.0
    cabin_temp_c: float = 22.0
    ad_heat_recovery: float = 0.70
    ad_loss_W_per_L: float = 0.05
    past_temp_c: float = 70.0
    past_heat_recovery: float = 0.80
    membrane_kWh_per_m3: float = 3.0
    dryer_kWh_per_kg_water: float = 0.21   # heat-pump dryer, COP ~3
    mixing_W_per_L: float = 0.02
    ch4_compress_kWh_per_m3: float = 0.15
    photon_J_per_umol: float = 0.20     # energy of red/blue PAR photons
    biomass_kJ_per_g: float = 21.0      # chemical energy stored in algae (protein-rich)

    # --- ESM (equivalent system mass) factors ----------------------------------
    # Lunar surface habitat values from NASA ICES-2022-417 Table 1 (citing BVAD Rev 2).
    veq_kg_per_m3: float = 133.1
    # Power and cooling are counted free: the ship has large surplus solar and radiator arrays.
    # (Lunar-habitat values if they had to be built for this system: 72.9 and 102 kg/kW.)
    peq_kg_per_kw: float = 0.0
    ceq_kg_per_kw: float = 0.0
    tank_kg_per_L: float = 0.25          # rough: tank shell + fittings
    pbr_kg_per_m2: float = 10.0          # rough: panel / tube walls, frame, manifolds per m2 of wall
    panel_pitch_m: float = 0.10          # spacing of panels incl. LED boards (sets PBR bay volume)
    led_kg_per_kw: float = 15.0          # rough: LED fixtures + drivers
    packing_factor: float = 2.0          # pressurized volume per litre of tank

    # --- digestible crew consumables -> AD feed ---------------------------------------------------
    # Three streams join the digester feed: cellulose toilet paper, cellulose wet wipes and PHA food
    # containers. (est.) = engineering estimate. Sources: "Digestible crew consumables" spec.
    # Dry mass used (NASA trash 1.607 kg/CM-d; hygiene towels ~9.5 % split paper/wipes (est.);
    # packaging ~27 %, about half becomes PHA trays/cups (est.)). Set to 0 to switch a stream off.
    tp_kg_cmd: float = 0.035
    wipe_kg_cmd: float = 0.115
    pha_kg_cmd: float = 0.20            # range 0.10-0.43
    tp_water_frac: float = 0.07         # kg water / kg wet at disposal
    wipe_water_frac: float = 0.70       # wipes hold 2-3 g liquid per g fabric
    pha_water_frac: float = 0.05        # food residue (est.)
    tp_vs_of_ts: float = 0.98           # pure cellulose / unfilled PHB: under 2 % ash
    wipe_vs_of_ts: float = 0.98
    pha_vs_of_ts: float = 0.99
    tp_vs_degradable: float = 0.85      # lignin-free cellulose digests almost fully
    wipe_vs_degradable: float = 0.85
    pha_vs_degradable: float = 0.70     # range 0.50-0.90 (49-52 % VS removal in 20 d; PHB ~90 % in 10 d)
    b0_cellulose: float = 0.35          # m3 CH4 / kg total VS (theoretical 0.415; tissue 0.33-0.38)
    b0_pha: float = 0.30                # m3 CH4 / kg total VS (PHB theoretical 0.586; PHBV 0.271 at 35 C)
    kh_cellulose_d: float = 0.15        # 1/d
    kh_pha_d: float = 0.15              # 1/d (0.23 if shredded to powder)
    pha_lag_d: float = 5.0              # PHA degraders build up after inoculation (start-up only, est.)
    tp_c_frac: float = 0.444            # kg C / kg dry (C6H10O5)
    wipe_c_frac: float = 0.444
    pha_c_frac: float = 0.558           # PHB C4H6O2
    tp_char_yield: float = 0.20         # kg char / kg dry if pyrolysed (cellulose 15-25 % at ~500 C)
    pha_char_yield: float = 0.03        # PHB unzips to volatiles: poor char feed
    tp_to_ad: float = 1.0               # share sent to AD (rest to pyrolysis)
    wipe_to_ad: float = 1.0
    pha_to_ad: float = 1.0              # PHA always goes to AD
    shred_kWh_per_kg_cell: float = 0.02  # paper and wipes only need pulping in water (est.)
    shred_kWh_per_kg_pha: float = 0.05  # PHA shredded under 2 mm (est.)
    feces_c_frac: float = 0.45          # kg C / kg fecal TS (est.)
    ad_cn_target: float = 20.0          # optimum 20-30; 18-25 for algae + paper (shown as a reference)
    ad_cn_max: float = 30.0             # cellulose carbon above this feed C/N goes to pyrolysis instead
    codig_kh_boost: float = 1.0         # algae hydrolysis x (1 + boost*effect*clip((C/N-6)/14, 0, 1))
    codig_effect: float = 0.5           # 0 = no co-digestion effect; 1 = full; wide spread in studies (est.)
    te_fe_mg_L: float = 5.0             # trace-element dosing into AD feed, mg/L (paper and PHA bring none)
    te_ni_mg_L: float = 0.5
    te_co_mg_L: float = 0.2
    te_se_mo_w_mg_L: float = 0.05
    consumables_replace_existing: bool = True   # they replace TP, wipes and packaging launched anyway

    def sp(self) -> dict:
        return SPECIES[self.species]

    def light_curve(self) -> tuple:
        """(y0 g/mol, K umol/m2/s) of the light response; K None = old linear model with the species yx."""
        t = LIGHT_TIERS[self.light_tier]
        if t["k_sat"] is None:
            return self.sp()["yx"], None
        return t["y0"], t["k_sat"]

    def photon_efficiency(self, q: float) -> float:
        """Share of the low-light yield kept at photon flux q (1 at low light, falling as light rises)."""
        _, K = self.light_curve()
        if K is None or q <= 1e-9:
            return 1.0
        return K / q * __import__("math").log1p(q / K)

    def inside_share(self) -> float:
        """Share of LED power mounted inside the culture. For tubes it follows the rod surface, so
        the rods and the outside wall reach their light limit together."""
        if self.pbr_shape == "tube":
            out = self.tube_light_coverage * self.tube_d_m
            ins = self.rods_per_tube * self.rod_d_m
            return ins / (out + ins) if ins > 0 else 0.0
        return self.internal_light_frac

    def delivery(self) -> float:
        """Fraction of LED photons that reach the culture, blending outside and inside LEDs."""
        f = self.inside_share()
        return (1 - f) * self.light_delivery + f * self.internal_delivery * (1 - self.internal_fouling)
