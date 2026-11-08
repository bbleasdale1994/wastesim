# wastesim

Simulation of a semi-closed loop waste recycling system.

An hourly SimPy model of a bioregenerative life-support loop for spaceflight: crew urine and feces →
nitrification and two-stage anaerobic digestion → methane, recovered water and fertilizer → LED algae
photobioreactor → oxygen and food for the crew. Crew inputs follow NASA's BVAD (Baseline Values and
Assumptions Document).

## Run it

- **Online:** open the Streamlit app link.
- **Locally:** `pip install -r requirements.txt`, then `streamlit run app.py`.

## Files

- `app.py`: interactive Streamlit app (sliders, charts, checks) with the system diagram built in.
- `waste_sim.py`: the simulation model.
- `params.py`: every assumption as a named parameter, with sources.
- `best_design.json`: the current optimized design (10 crew, algae supplying 60 % of O₂).
