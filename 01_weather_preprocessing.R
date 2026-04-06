  ################################################################################
# 01 — Weather Data Preprocessing
# ──────────────────────────────────
# Purpose:
#   Loads the corrected weather file (output of recalculate_weather_HUMID_CLIMATE.R)
#   and splits it into a nested list structure for the Monte Carlo model.
#
# Input:
#   weather_koeln-bonn_corrected1.rds
#     Must already contain derived columns (Tavg, PAR, ET0_mm, GDD_daily, etc.)
#     produced by recalculate_weather_HUMID_CLIMATE.R.
#
# Output:
#   weather_precomputed — a nested list:
#     weather_precomputed[[scenario]][[id_season]] -> data.table of daily weather
#
#   Scenarios: historical, ssp126, ssp245, ssp370, ssp585
#
# Required packages: zoo, RcppRoll, decisionSupport, compiler, data.table
################################################################################

# ── Load packages ─────────────────────────────────────────────────────────────

load_if_needed <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
    }
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE)
    )
  }
}

load_if_needed(c(
  "zoo",             # Time series tools
  "RcppRoll",        # Fast rolling operations
  "decisionSupport", # Monte Carlo simulation framework
  "compiler",        # Byte-code compilation
  "data.table"       # Fast data manipulation
))


# ── Process weather data ──────────────────────────────────────────────────────

process_weather_data <- function(
    file_path = "weather_koeln-bonn_parameterised.rds",
    scenarios = c("historical", "ssp126", "ssp245", "ssp370", "ssp585")
) {
  # Read weather file
  weather_combined <- readRDS(file_path)
  data.table::setDT(weather_combined)
  
  # Sort: by X column if present (preserves original row order), else by id_season + DATE
  if ("X" %in% names(weather_combined)) {
    data.table::setorder(weather_combined, X)
  } else if ("DATE" %in% names(weather_combined)) {
    data.table::setorder(weather_combined, id_season, DATE)
  } else {
    data.table::setkey(weather_combined, id_season, yday)
  }
  
  # Extract scenario name from id_season (format: "scenario--gcm--year")
  weather_combined[
    ,
    scenario := data.table::tstrsplit(id_season, "--", fixed = TRUE, keep = 1L)
  ]
  
  # Split into nested list: scenario -> id_season -> data.table
  weather_precomputed <- split(
    weather_combined,
    by      = "scenario",
    keep.by = FALSE
  )
  
  weather_precomputed <- lapply(weather_precomputed, function(dt) {
    split(dt, by = "id_season", keep.by = FALSE)
  })
  
  return(weather_precomputed)
}


# ── Run preprocessing ─────────────────────────────────────────────────────────

weather_precomputed <- process_weather_data(
  file_path = "weather_koeln-bonn_corrected1.rds"
)
