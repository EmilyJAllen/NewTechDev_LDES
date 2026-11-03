# Thermal LDES annual charge/discharge model

Run `ldes_main.m` from this folder. All inputs are in the INPUTS block at the top. Each one cites its source tag ([LOAD], [SOM], [eGRID], [EIA], [ANT], [TEA], [OMB], [TEAM]), and the full references are in the file header.

The run sizes the battery to the data center, then runs two dispatch cases against all three outage scenarios:

- **Arbitrage:** the stated price rule. Discharge when price > Pmax, charge when price < Pmax, nothing held back for outages.
- **Backup-only:** the blocks are held full for outages. They never discharge for price, and are topped up only when price < Pmax.

It prints one results column per case × scenario (6 columns) and writes the same table to `ldes_results_by_case.csv`.

## Plots

Saved as PNGs in `plots/` (set `in.plotDir = ''` to skip saving):

| File | What it shows |
|---|---|
| `case_arbitrage.png`, `case_backup_only.png` | One figure per dispatch case, one row per outage scenario: the year's state of charge with outages shaded; a zoom on that year's worst outage (load, TPV output, unserved load, grid charging); and state of charge through the same window. The year shown is the Monte Carlo year whose unserved load is closest to the scenario's mean. |
| `summary_requested_outputs.png` | The three requested outputs (charging cost, unserved load, charging CO2), case vs case for each scenario, with values labelled |
| `price.png` | The year's electricity price against Pmax |

## Files

| File | What it does |
|---|---|
| `ldes_main.m` | Inputs, sizing, runs everything, prints results, plots |
| `thermal_block_params.m` / `thermal_block_run.m` | The 3D block model as callable functions (same physics as `thermal_block_3D_full.m`) |
| `ldes_calibrate.m` | One-time 3D discharge run, which builds the one-block surrogate stored in `ldes_surrogate.mat` |
| `ldes_size_blocks.m` | Fewest blocks that carry a given load through a given outage with zero unserved load |
| `ldes_cycle.m` | Efficiency model: kWh in + discharge window gives electricity, coolant heat, insulation loss and energy left |
| `ldes_dispatch.m` | 15-min annual dispatch: outage backup, discharge when price > Pmax, charge when price < Pmax |
| `ldes_load_profile.m`, `ldes_read_csv_column.m` | Reads the CSV and handles quoted timestamps that contain commas |
| `ldes_outages.m` | Outage time series: the fraction of each 15-min step the grid is out |
| `ldes_price.m` | PJM-2025-calibrated series, or a CSV, flat price or vector |

## Assumptions

### Site and load
1. **Location:** the data center is in western Pennsylvania: PJM market, eGRID subregion RFCW. Prices, emissions and outage frequency are all sourced for that location.
2. **Load covers the whole facility:** the battery must serve `total_MW` (IT × PUE 1.2, mean 907 MW, peak 970 MW), not IT load alone.
3. **One week stands for the year:** the CSV covers Jan 1–7 2026 and is tiled 52×. There is no seasonal change in load or PUE (no summer cooling peak).
4. **No load shedding:** during an outage the full facility load must be served.

### Electricity price
5. **Price statistics are real, timing is stylized:** the 15-min price series matches PJM's 2025 real-time load-weighted mean ($50.73/MWh), median ($38.11) and std ($58.74). The hour-to-hour timing is synthetic, because the State of the Market report publishes the distribution, not the series.
6. **Energy price only:** no capacity, transmission, demand or retail charges. Charging cost is LMP × kWh.
7. **Price taker:** the battery's charging does not move the price, even at GW scale.
8. **Single threshold:** Pmax is the 2025 median LMP. The rule reacts to the current price only, with no forecasting.

### Outages
9. **One frequency for all scenarios:** every scenario uses 1.408 outages/yr, Pennsylvania's 2025 SAIFI including major events. That is a distribution-level statistic; a ~1 GW transmission-connected site would likely see fewer outages. It is also much too frequent for the 70.3 h Catastrophic case.
10. **Fixed duration per scenario:** Median 0.473 h, Bad 0.605 h, Catastrophic 70.3 h [TEAM].
11. **Random timing:** the outage count is Poisson and start times are uniform over the year. There is no correlation with weather, load or price. The same counts and start times are used for every scenario.
12. **Total grid loss:** during an outage the grid supplies nothing, so there is no charging. In a 15-min step only partly covered by an outage, the rest of the step is standby.
13. **Instant response:** TPV output is available immediately. The UPS that would bridge the first seconds is not modelled, and neither are parasitic loads (coolant pumps, controls).

### Battery physics (one block)
14. **3D model parameters are the team's placeholders:** material, insulation, gap and TPV values come from `thermal_block_3D_full.m`, which labels them as placeholders. They set the TPV power, round trip and self-discharge (1.35 %/day at full charge, vs Antora's claimed <1 %/day).
15. **Temperature range:** "full" is a uniform 2000 K (1727 °C) and "empty" is 1000 K, where the model's TPV efficiency reaches zero.
16. **One-state surrogate:** each block is reduced to a single state, its stored energy. Power and losses at a given energy follow the full 3D discharge run, so differences in internal temperature distribution after partial charges are ignored. Validated against the 3D model to within about 2.5 % on electricity out.
17. **Soaked insulation:** standby loss uses the steady-state solution through the insulation stack, which assumes the insulation has soaked. A freshly heated block loses more.
18. **Lossless heater:** all charging electricity becomes heat in the block (heater efficiency 1). Transformer and rectifier losses are not modelled.
19. **Coolant heat is wasted:** heat into the TPV coolant counts as a loss. No heat recovery.
20. **No degradation:** no TPV cell or block degradation over the life of the system.

### Scaling to the data center
21. **Identical, independent blocks:** the battery is N copies of the modelled block. Capacity, power and losses scale linearly with N. That is conservative on losses, since bigger blocks or shared insulation would lose less per kWh.
22. **Blocks discharge in unison:** all blocks share the same state of charge. Staggered discharge (keeping some blocks full) is not modelled; it would need fewer blocks for short outages.
23. **Sizing basis:** N is the fewest blocks that carry the peak facility load (970 MW) through one Median outage (0.473 h) from 100 % charge with zero unserved load. TPV power falls steeply as blocks drain (3.0 MW at full, 0.65 MW at 80 % charge), so the count is set by power, not energy.
24. **Charging cap:** total charging power is capped at the facility peak, 970 MW. This cap is on charging alone; grid import can reach load + charging, roughly 2× peak.

### Dispatch
25. **Two dispatch cases:** Arbitrage (`reserveFrac = 0`) lets price-driven discharge drain the battery before an outage, so the sizing guarantee (which assumes full charge) doesn't hold. Backup-only (`reserveFrac = 1`) keeps the blocks full and only charges when price < Pmax, so recharging after an outage waits for a cheap period.
26. **No export:** price-driven discharge only offsets the data center's own load.
27. **Warm start:** the year starts from the state of charge at the end of an outage-free year.

### Economics and emissions
28. **Capex:** $50/kWh-e total system, from the graphite-TPV TEA at >1 GWh scale. kWh-e is the electricity from one full discharge per block, taken from the 3D run. No O&M or replacements; 30 years at a 7 % real discount rate.
29. **Average emissions factor:** charging CO2 uses the fixed average factor (911.4 lb/MWh). The RFCW non-baseload (marginal) rate is 1,757 lb/MWh and would roughly double charging CO2.
30. **Monte Carlo:** 100 random outage years per scenario, seed 1.
