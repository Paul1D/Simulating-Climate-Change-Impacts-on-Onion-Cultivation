# Onion Climate Impact Model — NRW (Köln/Bonn)

Monte Carlo simulation model for assessing climate change impacts on onion (*Allium cepa* L.) production in North Rhine-Westphalia, Germany. The model estimates yield under historical and future climate scenarios (SSP1-2.6 to SSP5-8.5) for both rainfed and irrigated conditions.

## Repository Structure

```
├── recalculate_weather_HUMID_CLIMATE.R   # Weather parameter derivation
├── 01_weather_preprocessing.R            # Splits weather into scenario structure
├── 02_helper_functions.R                 # Sub-models: biomass, stress, utilities
├── 03_model.R                            # Main model function + Monte Carlo
├── input_values.csv                      # Parameter distributions for MC simulation
├── SCIENTIFIC_DOCUMENTATION.md           # Full documentation of weather derivation methods
│
├── weather_koeln-bonn_raw.RDS            # Raw DWD weather data (input)
├── weather_koeln-bonn_corrected1.rds     # Processed weather data (output of weather script)
├── final25.RDS                           # Saved MC simulation results
│
└── README.md                             # This file
```

## How to Run the Model

### Step 1: Process raw weather data

If you have updated raw weather data (daily Tmin, Tmax, Prec from DWD), run the weather processing script first:

```r
source("weather_param_derivation.R")
```

**Input:** `weather_koeln-bonn_raw.RDS` — raw weather data containing daily Tmin, Tmax, and Prec for historical observations and CMIP6 climate projections.

**Output:** `weather_koeln-bonn_parameterised.rds` — same data enriched with derived parameters: Tavg, relative humidity (RH_mean, RH_max, RH_min), soil temperature (Ts_5cm), reference evapotranspiration (ET0_mm, FAO-56 Penman-Monteith), photosynthetically active radiation (PAR), growing degree days (GDD_daily), and consecutive wet/dry day counts.

> **Note:** If you are using the already-processed file `weather_koeln-bonn_corrected1.rds`, you can skip this step.

### Step 2: Run the model

Run all three scripts in order:

```r
source("01_weather_preprocessing.R")   # Creates weather_precomputed (nested list)
source("02_helper_functions.R")        # Loads stress & biomass functions
source("03_model.R")                   # Runs Monte Carlo simulation
```

The scripts source each other automatically. `03_model.R` calls `mcSimulation()` from the `decisionSupport` package and produces the `onion_mc_simulation` object.

### Step 3: Analyze results

Results are stored in `onion_mc_simulation$y` with columns named `{scenario}.{variable}`, e.g.:

```r
# Mean historical rainfed yield (t/ha fresh weight)
mean(onion_mc_simulation$y$historical.final_yield_per_ha_rainfed)

# Mean irrigated yield under SSP5-8.5
mean(onion_mc_simulation$y$ssp585.final_yield_per_ha_irrigated)

# Yield gap between irrigated and rainfed
mean(onion_mc_simulation$y$historical.yield_gap_rainfed_irrigated)
```

## Model Overview

### Biomass Accumulation

Daily biomass increment (g m⁻²) is calculated as:

```
Biomass = RUE × f(CO₂) × PAR × (1 - exp(-k × LAI)) × f(T) × f(W)
```

- **RUE**: Radiation use efficiency, phase-specific (g MJ⁻¹ iPAR)
- **f(CO₂)**: CO₂ fertilization factor (SIMPLE model, Zhao et al. 2019)
- **PAR**: Photosynthetically active radiation (MJ m⁻² day⁻¹)
- **f(IPAR)**: Beer-Lambert light interception
- **f(T)**: PAR-weighted beta temperature response (Wang & Engel 1998), integrated hourly using sinusoidal diurnal temperature (Parton & Logan 1981) and PAR (Goudriaan 1986)
- **f(W)**: Linear water stress function based on relative soil moisture (theta_rel)

### Water Balance

Irrigation scheduling follows the Geisenheimer Bewässerungssteuerung (Schmidt & Zinkernagel 2017). The model runs two parallel water balances per season:
1. **Rainfed** — no irrigation applied (max_irrig_mm = 0)
2. **Irrigated** — irrigation triggered when cumulative deficit exceeds threshold

Crop coefficients (kc) and rooting depths are phase-specific, following Geisenheim field trials with onion.

### Biotic Stress

Each pathogen is scored daily using a beta temperature response gated by a moisture/humidity condition. The seasonal risk index is the mean daily score across the growth phase.

| Pathogen | Temperature variable | Moisture gate | Reference |
|---|---|---|---|
| Downy mildew (*P. destructor*) | Tmin (night proxy) | RH_max ≥ 93-97% | Jesperson & Sutton (1987); Gilles et al. (2004) |
| Botrytis leaf blight (*B. squamosa*) | Tmin (night proxy) | Prec ≥ 3-5 mm/day | Sutton et al. (1986) |
| Fusarium basal rot (*F. oxysporum*) | Ts_5cm (soil) | theta_rel > 0.7-0.8 (waterlogging) | Abawi & Lorbeer (1972) |

Extreme rainfall and hail are handled as abiotic stressors.

### Monte Carlo Parameters

All model parameters are defined in `input_values.csv` with lower, median, upper bounds and distribution type. The `decisionSupport` package draws random values from these distributions for each MC iteration.

## Required R Packages

- `data.table` (≥ 1.14.0) — fast data manipulation
- `lubridate` (≥ 1.8.0) — date handling (weather script only)
- `decisionSupport` — Monte Carlo simulation framework
- `compiler` — byte-code compilation for speed
- `zoo` — time series tools
- `RcppRoll` — fast rolling operations

Install all:
```r
install.packages(c("data.table", "lubridate", "decisionSupport", "compiler", "zoo", "RcppRoll"))
```

## Data Sources

- **Weather data:** Deutscher Wetterdienst (DWD), Köln/Bonn Airport station (ID 2667). Basic data publicly available from [DWD Climate Data Center](https://opendata.dwd.de/climate_environment/CDC/).
- **Climate projections:** CMIP6 models under SSP1-2.6, SSP2-4.5, SSP3-7.0, SSP5-8.5.
- **CO₂ concentrations:** Meinshausen et al. (2020), MAGICC7.0 projections per SSP pathway.

## Key References

- Allen, R.G., Pereira, L.S., Raes, D., Smith, M. (1998). *Crop evapotranspiration.* FAO Irrigation and Drainage Paper 56.
- Schmidt, N., Zinkernagel, J. (2017). Model and growth stage based variability of the irrigation demand of onion crops with predicted climate change. *Water* 9:693.
- Zhao, C. et al. (2019). A SIMPLE crop model. *Eur. J. Agron.* 104:97-106.
- Sutton, J.C., James, T.D.W., Rowell, P.M. (1986). BOTCAST. *Agric. Ecosyst. Environ.* 18:123-143.
- Jesperson, G.D., Sutton, J.C. (1987). Evaluation of a forecaster for downy mildew of onion. *Crop Prot.* 6:95-103.
- Gilles, T. et al. (2004). MILIONCAST. *Plant Dis.* 88:695-702.
- Abawi, G.S., Lorbeer, J.W. (1972). Several aspects of the ecology and pathology of *Fusarium oxysporum* f. sp. *cepae*. *Phytopathology* 62:870-876.


## Data Availability

- all Data utilized in the model is uploaded in this repository

Alternatively, the processed data files are available from the corresponding author upon request.
