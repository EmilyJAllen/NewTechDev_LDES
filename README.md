# Thermal LDES annual charge/discharge model

Run `ldes_main.m` from this folder. All inputs are in the INPUTS block at the top. Each one cites its source tag ([LOAD], [SOM], [eGRID], [EIA], [ANT], [TEA], [OMB]), and the full references are in the file header.

## Site and data

- **Location:** western Pennsylvania (PJM market, eGRID subregion RFCW).
- **Load:** `Load_Profile_1000_MWitcapacity_1.2_PUE.csv`, column `total_MW` (IT × PUE 1.2, mean ~907 MW). The CSV covers one week (Jan 1–7 2026), which is tiled to the full year.
- **Price:** a 15-min series calibrated to PJM's 2025 real-time LMP statistics. It reproduces the load-weighted mean ($50.73/MWh), median ($38.11) and std ($58.74). The daily, seasonal and weekend *timing* is stylized; the State of the Market report gives the distribution, not the hourly series. To use real timing, export PJM real-time hourly LMPs (Data Miner 2, feed `rt_hrl_lmps`) and set `in.price.mode = 'csv'`, `in.price.file`, `in.price.column`.
- **Pmax:** the 2025 median LMP, $38.11/MWh.
- **Outages:** Pennsylvania 2025 SAIFI and SAIDI including major events: 1.408 per year, 5.42 h average. The count is drawn as Poisson and start times are uniform, over 100 random years.
- **Emissions:** eGRID2023 RFCW, 911.4 lb CO2/MWh (4.134e-4 t/kWh).
- **Cost:** $50/kWh-e total system, from the graphite-TPV TEA at GWh scale. The kWh-e is the electricity from one full block discharge, taken from the 3D calibration run. Financing is 30 yr at 7 % real.
- **Battery:** one block exactly as in `thermal_block_3D_full.m` (35.3 MWh_th between 1000 and 2000 K, 3.87 MW heater).

## Files

| File | What it does |
|---|---|
| `ldes_main.m` | Inputs, runs everything, prints results, plots |
| `thermal_block_params.m` / `thermal_block_run.m` | Your 3D model as callable functions (same physics) |
| `ldes_calibrate.m` | One-time 3D discharge run, which builds the surrogate stored in `ldes_surrogate.mat` |
| `ldes_cycle.m` | Efficiency model: kWh in + discharge window gives electricity, coolant heat, insulation loss and energy left |
| `ldes_dispatch.m` | 15-min annual dispatch: outage backup, discharge when price > Pmax, charge when price < Pmax |
| `ldes_load_profile.m`, `ldes_read_csv_column.m` | Reads the CSV and handles quoted timestamps that contain commas |
| `ldes_outages.m` | Outage mask |
| `ldes_price.m` | PJM-2025-calibrated series, or a CSV, flat price or vector |

## Still yours to own

The 3D model's material, insulation and TPV parameters (`thermal_block_3D_full.m` and `thermal_block_params.m`) are the values in your model, and that file labels them as placeholders. They set the TPV power, the round trip and the self-discharge (1.35 %/day here, vs Antora's claimed <1 %/day). After changing them, set `in.sur.recalibrate = true`.
