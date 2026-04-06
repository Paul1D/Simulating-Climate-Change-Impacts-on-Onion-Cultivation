################################################################################
# Weather Parameter Derivation Function
# ──────────────────────────────────────
# Station:  Köln/Bonn Airport (DWD Station ID 2667)
#           50.87°N, 7.16°E, 91 m AMSL
# Climate:  Humid oceanic (Köppen Cfb)
#
# Purpose:
#   Derives meteorological parameters from daily Tmin, Tmax, and Prec.
#   These derived variables are required by the onion climate impact model
#   (03_model.R).
#
# Derived parameters:
#   Tavg            Mean air temperature [°C]
#   RH_mean/max/min Relative humidity [%]
#   Ts_5cm          Soil temperature at 5 cm depth [°C]
#   Ts_5cm_smooth   7-day smoothed soil temperature [°C]
#   Ra              Extraterrestrial radiation [MJ m⁻² day⁻¹]
#   PAR             Photosynthetically active radiation [MJ m⁻² day⁻¹]
#   ET0_mm          Reference evapotranspiration (FAO-56 PM) [mm day⁻¹]
#   GDD_daily       Growing degree days [°C·day]
#   day_consec_dry  Consecutive dry days (Prec < 1 mm)
#   day_consec_wet  Consecutive wet days (Prec >= 1 mm)
#
# Key references:
#   Allen et al. (1998) FAO Irrigation & Drainage Paper 56.
#   Hargreaves & Samani (1985) Appl. Eng. Agric. 1(2):96-99.
#   Van Wijk & De Vries (1963) in: Physics of Plant Environment, pp. 102-143.
#   Howell et al. (1983) Agric. Meteorol. 28(2):157-175.
#   McMaster & Wilhelm (1997) Agric. For. Meteorol. 87(4):291-300.
#
# For full scientific documentation see: SCIENTIFIC_DOCUMENTATION.md
################################################################################

recalculate_weather_parameters <- function(met, base_temp_gdd = 0) {
  
  require(data.table)
  require(lubridate)
  
  dt <- as.data.table(met)
  
  # ── Station & physical constants ──────────────────────────────────────────
  LATITUDE     <- 50.87       # degrees N
  ELEVATION    <- 91          # m above mean sea level
  GSC          <- 0.0820      # solar constant [MJ m⁻² min⁻¹]
  KRS          <- 0.16        # Hargreaves radiation coefficient (inland)
  SOIL_DEPTH   <- 0.05        # target soil depth [m]
  OMEGA        <- 2 * pi / 86400  # angular frequency for daily cycle [s⁻¹]
  PAR_FRACTION <- 0.48        # PAR as fraction of Rs
                              # (Howell et al. 1983; Tsubo & Walker 2005)
  
  # ── Data cleaning ────────────────────────────────────────────────────────
  # Fix known typo in raw DWD data column name
  if ("id_seaon" %in% names(dt) && !("id_season" %in% names(dt))) {
    setnames(dt, "id_seaon", "id_season")
  }
  
  # Remove CSV row-index column (artifact from write.csv)
  if ("X" %in% names(dt)) {
    setorder(dt, X)
    dt[, X := NULL]
  }
  
  if (!inherits(dt$DATE, "Date")) {
    dt[, DATE := as.Date(DATE)]
  }
  
  if (!("id_season" %in% names(dt))) {
    dt[, id_season := "single_season"]
  }
  
  setorder(dt, id_season, DATE)
  
  # Remove old calculated columns to avoid duplicates on re-runs
  old_cols <- c("Tavg", "RH_mean", "RH_max", "RH_min", "Ts_5cm", 
                "Ts_5cm_smooth", "Ra", "Rs", "PAR", "ET0_mm", "GDD_daily",
                "day_consec_dry", "day_consec_wet", "yday")
  old_cols <- intersect(old_cols, names(dt))
  if (length(old_cols) > 0) {
    dt[, (old_cols) := NULL]
  }
  
  ##############################################################################
  # 1. MEAN AIR TEMPERATURE
  #    Standard WMO method: arithmetic mean of daily extremes.
  #    Ref: WMO (2008) Guide to Meteorological Instruments, WMO-No. 8.
  ##############################################################################
  
  dt[, Tavg := (Tmax + Tmin) / 2]
  dt[, yday := yday(DATE)]
  
  ##############################################################################
  # 2. RELATIVE HUMIDITY
  #    Estimated from temperature using the Tetens formula for saturation
  #    vapor pressure (FAO-56 Eq. 11) and dewpoint depression.
  #
  #    Dewpoint depression is adjusted for Köln's humid oceanic climate.
  #    In humid climates, dewpoint remains close to Tmin.
  #    Ref: Allen et al. (1998) FAO-56, Eq. 11-12, 48.
  #         Tetens (1930) Z. Geophysik 6:297-309.
  ##############################################################################
  
  # Saturation vapor pressure [kPa] (FAO-56 Eq. 11)
  dt[, es_Tmax := 0.6108 * exp((17.27 * Tmax) / (Tmax + 237.3))]
  dt[, es_Tmin := 0.6108 * exp((17.27 * Tmin) / (Tmin + 237.3))]
  dt[, es_mean := (es_Tmax + es_Tmin) / 2]   # FAO-56 Eq. 12
  
  dt[, temp_range := Tmax - Tmin]
  
  # Dewpoint depression [°C] — calibrated for Köln's humid climate (Cfb).
  # Precipitation indicates high atmospheric moisture; small diurnal range
  # indicates cloudy/humid conditions.
  dt[, dewpoint_depression := fcase(
    Prec > 5,          0.4,    # Heavy rain: RH_max ~98-99%
    Prec > 1,          0.5,    # Light rain: RH_max ~96-98%
    Prec > 0,          0.7,    # Trace rain: RH_max ~94-96%
    temp_range < 5,    0.6,    # Cloudy (small diurnal range): RH_max ~95-97%
    temp_range < 8,    0.9,    # Moderate range: RH_max ~92-95%
    temp_range < 12,   1.2,    # Larger range: RH_max ~88-92%
    temp_range >= 12,  1.6,    # Very clear (rare in Köln): RH_max ~85-88%
    default = 0.8
  )]
  
  # Actual vapor pressure from estimated dewpoint
  dt[, Tdew := Tmin - dewpoint_depression]
  dt[, ea := 0.6108 * exp((17.27 * Tdew) / (Tdew + 237.3))]
  
  # Relative humidity, constrained to 0-100%
  dt[, RH_max  := pmin(100, pmax(0, 100 * ea / es_Tmin))]
  dt[, RH_min  := pmin(100, pmax(0, 100 * ea / es_Tmax))]
  dt[, RH_mean := pmin(100, pmax(0, 100 * ea / es_mean))]
  
  ##############################################################################
  # 3. SOIL TEMPERATURE AT 5 cm DEPTH
  #    Harmonic heat diffusion model.
  #    Thermal diffusivity varies with soil moisture (estimated from Prec).
  #    Ref: Van Wijk & De Vries (1963) in: Physics of Plant Environment.
  #         Zheng et al. (1993) Climate Research 2(3):183-191.
  ##############################################################################
  
  dt[, Tmean7 := frollmean(Tavg, n = 7, align = "right", fill = NA_real_),
     by = id_season]
  dt[, Tmean7 := fifelse(is.na(Tmean7), Tavg, Tmean7)]
  dt[, Amp := (Tmax - Tmin) / 2]
  
  # Thermal diffusivity [m² s⁻¹] — depends on soil moisture
  dt[, alpha := fcase(
    Prec == 0,             0.3e-6,    # Dry soil
    Prec > 0 & Prec <= 5,  0.6e-6,    # Moist soil
    Prec > 5,              1.0e-6     # Wet soil
  )]
  dt[, delta := sqrt(2 * alpha / OMEGA)]
  dt[, Ts_5cm := Tmean7 + Amp * exp(-SOIL_DEPTH / delta) * 
       sin(2 * pi * (yday %% 1) - SOIL_DEPTH / delta)]
  dt[, Ts_5cm_smooth := frollmean(Ts_5cm, n = 7, align = "right", 
                                  fill = NA_real_), by = id_season]
  
  ##############################################################################
  # 4. EXTRATERRESTRIAL RADIATION (Ra)
  #    Based on solar geometry.
  #    Ref: Allen et al. (1998) FAO-56, Eq. 21-25.
  ##############################################################################
  
  phi <- LATITUDE * pi / 180    # latitude in radians
  dt[, dr := 1 + 0.033 * cos(2 * pi * yday / 365)]               # Eq. 23
  dt[, delta_sun := 0.409 * sin((2 * pi * yday / 365) - 1.39)]   # Eq. 24
  dt[, omega_s := acos(pmax(-1, pmin(1, -tan(phi) * tan(delta_sun))))]  # Eq. 25
  dt[, Ra := (24 * 60 / pi) * GSC * dr * 
       (omega_s * sin(phi) * sin(delta_sun) + 
          cos(phi) * cos(delta_sun) * sin(omega_s))]              # Eq. 21
  
  ##############################################################################
  # 5. SOLAR RADIATION (Rs) — Hargreaves-Samani method
  #    Rs = kRs × sqrt(Tmax - Tmin) × Ra
  #    Large diurnal range = clear sky = high Rs.
  #    Ref: Hargreaves & Samani (1985) Appl. Eng. Agric. 1(2):96-99.
  #         Allen et al. (1998) FAO-56, Eq. 50.
  ##############################################################################
  
  dt[, Rs := KRS * sqrt(pmax(temp_range, 0)) * Ra]
  
  ##############################################################################
  # 6. PHOTOSYNTHETICALLY ACTIVE RADIATION (PAR)
  #    PAR = 48% of total solar radiation (wavelengths 400-700 nm).
  #    Ref: Howell et al. (1983) Agric. Meteorol. 28(2):157-175.
  ##############################################################################
  
  dt[, PAR := PAR_FRACTION * Rs]
  
  ##############################################################################
  # 7. REFERENCE EVAPOTRANSPIRATION (ET0) — FAO-56 Penman-Monteith
  #    Ref: Allen et al. (1998) FAO-56, Eq. 6.
  #    Assumptions:
  #      - Wind speed u2 = 2.0 m/s (FAO default when not measured)
  #      - Soil heat flux G = 0 for daily time step (FAO-56 recommendation)
  ##############################################################################
  
  # Atmospheric pressure at station elevation [kPa]
  P     <- 101.3 * ((293 - 0.0065 * ELEVATION) / 293)^5.26
  gamma <- 0.000665 * P    # Psychrometric constant [kPa °C⁻¹]
  
  # Slope of saturation vapor pressure curve (FAO-56 Eq. 13)
  dt[, Delta := 4098 * (0.6108 * exp((17.27 * Tavg) / (Tavg + 237.3))) / 
       (Tavg + 237.3)^2]
  
  # Clear-sky radiation (FAO-56 Eq. 37)
  dt[, Rso := (0.75 + 2e-5 * ELEVATION) * Ra]
  
  # Net shortwave radiation (FAO-56 Eq. 38)
  dt[, Rns := 0.77 * Rs]
  
  # Net longwave radiation (FAO-56 Eq. 39)
  sigma <- 4.903e-9    # Stefan-Boltzmann [MJ K⁻⁴ m⁻² day⁻¹]
  dt[, ea_kPa := 0.6108 * exp((17.27 * (Tmin - dewpoint_depression)) / 
                                ((Tmin - dewpoint_depression) + 237.3))]
  dt[, Tmax_K := Tmax + 273.16]
  dt[, Tmin_K := Tmin + 273.16]
  dt[, temp_term  := (Tmax_K^4 + Tmin_K^4) / 2]
  dt[, ea_term    := 0.34 - 0.14 * sqrt(ea_kPa)]
  dt[, Rs_Rso     := pmin(Rs / Rso, 1.0)]
  dt[, cloud_term := pmax(0.05, pmin(1.35 * Rs_Rso - 0.35, 1.0))]
  dt[, Rnl := sigma * temp_term * ea_term * cloud_term]
  
  # Net radiation
  dt[, Rn := pmax(Rns - Rnl, 0)]
  
  G  <- 0      # Soil heat flux [MJ m⁻² day⁻¹] — negligible for daily step
  u2 <- 2.0    # Wind speed at 2 m height [m s⁻¹]
  
  # Vapor pressure deficit [kPa]
  dt[, VPD := pmax(es_mean - ea_kPa, 0)]
  
  # FAO-56 Penman-Monteith equation (Eq. 6)
  dt[, ET0_mm := (0.408 * Delta * (Rn - G) + 
                    gamma * (900 / (Tavg + 273)) * u2 * VPD) /
       (Delta + gamma * (1 + 0.34 * u2))]
  dt[, ET0_mm := pmax(0, pmin(ET0_mm, 12))]    # Physically reasonable range
  
  # Remove intermediate ET0 columns
  dt[, c("Rso", "Rns", "ea_kPa", "Tmax_K", "Tmin_K", "temp_term", 
         "ea_term", "Rs_Rso", "cloud_term", "Rnl") := NULL]
  
  ##############################################################################
  # 8. GROWING DEGREE DAYS (GDD)
  #    GDD = max(0, Tavg - Tbase)
  #    Default Tbase = 0°C (adjustable via function argument).
  #    Ref: McMaster & Wilhelm (1997) Agric. For. Meteorol. 87(4):291-300.
  ##############################################################################
  
  dt[, GDD_daily := pmax(0, Tavg - base_temp_gdd)]
  
  ##############################################################################
  # 9. CONSECUTIVE WET/DRY DAYS
  #    Wet day: Prec >= 1.0 mm (WMO standard threshold).
  #    Uses loop-based counting to avoid memory issues with large datasets.
  #    Ref: WMO (2017) Guidelines on Calculation of Climate Normals, WMO-1203.
  ##############################################################################
  
  dt[, is_wet := as.integer(Prec >= 1)]
  
  dt[, day_consec_wet := {
    n <- .N
    result <- integer(n)
    if (n > 0) {
      result[1] <- is_wet[1]
      if (n > 1) {
        for (i in 2:n) {
          if (is_wet[i] == 1) {
            result[i] <- if (is_wet[i-1] == 1) result[i-1] + 1 else 1
          } else {
            result[i] <- 0
          }
        }
      }
    }
    result
  }, by = id_season]
  
  dt[, day_consec_dry := {
    n <- .N
    result <- integer(n)
    if (n > 0) {
      result[1] <- 1 - is_wet[1]
      if (n > 1) {
        for (i in 2:n) {
          if (is_wet[i] == 0) {
            result[i] <- if (is_wet[i-1] == 0) result[i-1] + 1 else 1
          } else {
            result[i] <- 0
          }
        }
      }
    }
    result
  }, by = id_season]
  
  dt[, is_wet := NULL]
  
  ##############################################################################
  # CLEANUP — remove intermediate columns
  ##############################################################################
  
  cols_to_remove <- c("es_Tmax", "es_Tmin", "es_mean", "ea", "Tmean7", "Amp",
                      "alpha", "delta", "dr", "delta_sun", "omega_s", 
                      "temp_range", "Delta", "Rn", "VPD", "Tdew", 
                      "dewpoint_depression", "Rs")
  cols_to_remove <- intersect(cols_to_remove, names(dt))
  if (length(cols_to_remove) > 0) {
    dt[, (cols_to_remove) := NULL]
  }
  
  ##############################################################################
  # ROUND VALUES
  ##############################################################################
  
  dt[, Tavg          := round(Tavg, 2)]
  dt[, RH_mean       := round(RH_mean, 1)]
  dt[, RH_max        := round(RH_max, 1)]
  dt[, RH_min        := round(RH_min, 1)]
  dt[, Ts_5cm        := round(Ts_5cm, 2)]
  dt[, Ts_5cm_smooth := round(Ts_5cm_smooth, 2)]
  dt[, Ra            := round(Ra, 2)]
  dt[, PAR           := round(PAR, 2)]
  dt[, ET0_mm        := round(ET0_mm, 2)]
  dt[, GDD_daily     := round(GDD_daily, 2)]
  
  # Reorder columns
  target_order <- c("DATE", "Year", "Month", "Day", "Tmin", "Tmax", "Tavg", 
                    "Prec", "RH_mean", "Ts_5cm", "Ra", "day_consec_dry", 
                    "day_consec_wet", "ssp", "gcm", "scenario_year", "season",
                    "id", "id_season", "yday", "RH_max", "RH_min", "ET0_mm",
                    "Ts_5cm_smooth", "PAR", "GDD_daily")
  target_order <- intersect(target_order, names(dt))
  setcolorder(dt, target_order)
  
  return(dt)
}

################################################################################
# USAGE
################################################################################

met <- readRDS("weather_koeln-bonn_raw.RDS")
met <- recalculate_weather_parameters(met, base_temp_gdd = 0)

# Expected results for Köln's humid climate:
# summary(met$RH_max)   # 92-100% (often 98-100% in winter)
# summary(met$RH_min)   # 55-80%
# summary(met$RH_mean)  # 75-90%

saveRDS(met, "weather_koeln-bonn_corrected1.rds", compress = "xz")
