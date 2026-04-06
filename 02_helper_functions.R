################################################################################
# 02 — Helper Functions for Onion Climate Impact Model
# ──────────────────────────────────────────────────────
# Purpose:
#   Defines all sub-models used by the main model function (03_model.R):
#     - Biomass accumulation   (calc_bio_vectorized)
#     - Biotic stress scoring  (botrytis, downy mildew, fusarium)
#     - Abiotic stress scoring (extreme rain)
#     - Utility functions      (sat, consec_counter, etc.)
#
# Biotic stress approach:
#   Each day receives a favorability score (0–1) based on:
#     score = Beta_T(daily temperature) × I(moisture/humidity gate)
#   The seasonal risk index = mean(daily scores) over the growth phase.
#   This represents "infection pressure" — the fraction of days with
#   conditions favorable to the pathogen, weighted by how favorable.
#   Ref: Madden et al. (2007) The Study of Plant Disease Epidemics.
#
# Disease model references:
#   Botrytis leaf blight:  Sutton et al. (1986) BOTCAST.
#                          Agric. Ecosyst. Environ. 18:123-143.
#   Downy mildew:          Jesperson & Sutton (1987) DOWNCAST.
#                          Crop Prot. 6:95-103.
#                          Gilles et al. (2004) MILIONCAST.
#                          Plant Dis. 88:695-702.
#   Fusarium basal rot:    Abawi & Lorbeer (1972) Phytopathology 62:870.
#                          Cramer (2000) Plant Health Progress.
#
# Biomass model references:
#   Temperature response:  Wang & Engel (1998) beta function.
#   Hourly temperature:    Parton & Logan (1981) sinusoidal diurnal model.
#   Hourly PAR:            Goudriaan (1986) sinusoidal clear-sky arc.
#   CO2 fertilization:     Zhao et al. (2019) SIMPLE model, Eq. 10.
#                          Eur. J. Agron. 104:97-106.
#   Light interception:    Beer-Lambert law with extinction coefficient k.
#
# Required packages: zoo, RcppRoll, decisionSupport, compiler, data.table
################################################################################

# ── Install required packages ─────────────────────────────────────────────────

ensure_installed <- function(pkgs) {
  for (pkg in pkgs) {
    if (!base::requireNamespace(pkg, quietly = TRUE)) {
      utils::install.packages(pkg, dependencies = TRUE)
    }
  }
}

ensure_installed(c(
  "zoo",
  "RcppRoll",
  "decisionSupport",
  "compiler",
  "data.table"
))


# ── Helper Function Factory ───────────────────────────────────────────────────
# Returns a list of all helper functions. If attach_to_global = TRUE,
# functions are also attached to the global environment for interactive use.

helper_function <- function(attach_to_global = FALSE) {
  
  # ── Utility functions ──────────────────────────────────────────────────────
  
  # Saturating transform: maps x >= 0 to [0,1) with slope controlled by `impact`.
  # Used to convert counts (e.g. extreme rain days) into a bounded risk score.
  sat <- compiler::cmpfun(function(x, impact) {
    base::pmin(1 - base::exp(-impact * base::pmax(x, 0)), 1)
  })
  
  # Count consecutive TRUEs in a logical vector.
  # Example: c(T,T,F,T,T,T) -> c(1,2,0,1,2,3)
  consec_counter <- compiler::cmpfun(function(x_is_true) {
    x <- base::as.logical(x_is_true)
    n <- base::length(x)
    out <- base::integer(n)
    if (n == 0L) return(out)
    k <- 0L
    for (i in base::seq_len(n)) {
      if (base::isTRUE(x[i])) {
        k <- k + 1L
      } else {
        k <- 0L
      }
      out[i] <- k
    }
    out
  })
  
  # Stochastic gate: returns value_if with probability p, else value_if_not.
  # Used in Monte Carlo to decide if a stress event occurs in a given season.
  fast_chance_event <- compiler::cmpfun(function(p, value_if, value_if_not = 1) {
    if (base::is.na(p)) return(value_if_not)
    if (stats::runif(1) < p) value_if else value_if_not
  })
  
  # Expected loss: risk × loss_fraction, NA-safe, bounded to [0,1]
  expected_loss <- compiler::cmpfun(function(risk, loss_fraction) {
    if (base::is.na(risk) || base::is.na(loss_fraction)) return(0)
    base::pmin(base::pmax(risk, 0), 1) * base::pmin(base::pmax(loss_fraction, 0), 1)
  })
  
  
  # ── Abiotic stress functions ─────────────────────────────────────────────────
  
  # Extreme rain stress: weighted count of extreme rainfall days.
  # Medium events count 1×, high events count 2×.
  # Risk = 1 - exp(-impact × weighted_days), bounded to [0,1].
  get_extreme_rain_stress <- compiler::cmpfun(function(Prec,
                                                       prec_extreme_rain_medium_p,
                                                       prec_extreme_rain_high_p,
                                                       impact_days_extreme_rain_t) {
    
    extreme_days_medium <- base::sum(
      Prec >= prec_extreme_rain_medium_p & Prec < prec_extreme_rain_high_p,
      na.rm = TRUE
    )
    extreme_days_high <- base::sum(
      Prec >= prec_extreme_rain_high_p,
      na.rm = TRUE
    )
    
    weighted_extreme_days <- extreme_days_medium + 2 * extreme_days_high
    risk_extreme_rain <- 1 - base::exp(- (impact_days_extreme_rain_t / 2) * weighted_extreme_days)
    base::pmin(risk_extreme_rain, 1)
  })
  
  # Hail: removed — modeled as pure seasonal chance_event in 3_irrig.R
  # (no reliable daily proxy exists for hail from surface weather data)
  
  # Heat stress: removed — temperature effect on biomass already captured by
  # Beta f_T(Tavg) in calc_bio_vectorized. Onions tolerate heat well (Brewster 2008).
  # Including a separate heat stressor would double-count the temperature penalty.
  
  
  # ── Biotic stress functions ──────────────────────────────────────────────────

  # Risk index aggregation: mean daily score over the growth phase.
  #
  # For each day: score_i = Beta_T(temperature) × I(moisture gate passes)
  # Phase risk  = mean(score_i) across all days in the phase.
  #
  # Why mean score (not cumulative hazard or worst-run)?
  #   (1) Never saturates — preserves discrimination across scenarios.
  #   (2) Directly reflects the climate signal (warmer/wetter futures shift
  #       the temperature distribution toward pathogen optima).
  #   (3) Interpretable: mean=0.30 means 30% of days were favorable.
  #   (4) Zero free parameters beyond the daily scoring functions.
  #
  # Ref: Madden et al. (2007) The Study of Plant Disease Epidemics.
  #      Sutton et al. (1986) BOTCAST; Jesperson & Sutton (1987) DOWNCAST.
  aggregate_infection_periods <- compiler::cmpfun(function(scores,
                                                           min_score_t = 0.05) {
    scores_clean <- base::pmin(base::pmax(
      base::ifelse(base::is.na(scores), 0, scores), 0), 1)
    
    n <- base::length(scores_clean)
    if (n == 0L) return(0)
    
    # Risk index = mean daily score over the phase.
    # Favourable-day gate (min_score_t) retained for consistency: days with
    # score below threshold are treated as zero (no meaningful pathogen activity).
    scores_clean[scores_clean < min_score_t] <- 0
    
    base::mean(scores_clean)
  })
  
  
  # Botrytis leaf blight (B. squamosa) — daily favorable-day scoring.
  #
  # Gate 1: Daily precipitation >= threshold (proxy for leaf surface wetness).
  #         BOTCAST uses wetness duration; precipitation is a proxy.
  # Gate 2: Tmax < 30°C (high temperatures inhibit sporulation).
  # Score:  Beta_T(Tmin) — Tmin as proxy for night temperature during dew.
  #         Optimum 15-20°C, range 1-27°C.
  #
  # Ref: Sutton et al. (1986) Agric. Ecosyst. Environ. 18:123-143 (BOTCAST).
  #      Sutton (1990) Can. J. Plant Pathol. 12:100-110.
  
  get_botrytis_stress <- compiler::cmpfun(function(Tmin_daily,
                                                   Prec_daily,
                                                   prec_botrytis_wet_p,
                                                   Topt_botrytis_p,
                                                   Tmax_botrytis_p,
                                                   Tmin_botrytis_p,
                                                   min_score_botrytis = 0.05) {
    n <- base::length(Tmin_daily)
    if (n == 0L) return(0)
    
    scores <- base::numeric(n)
    for (i in base::seq_len(n)) {
      t_i    <- Tmin_daily[i]
      prec_i <- if (base::is.na(Prec_daily[i])) 0 else Prec_daily[i]
      if (base::is.na(t_i)) { scores[i] <- 0; next }
      
      # Precipitation leaf-wetness gate
      if (prec_i < prec_botrytis_wet_p) { scores[i] <- 0; next }
      
      # Beta_T(Tmin): temperature score on gate-passing days
      if (t_i <= Tmin_botrytis_p || t_i >= Tmax_botrytis_p) {
        scores[i] <- 0
      } else {
        alpha <- (Tmax_botrytis_p - Topt_botrytis_p) /
          (Topt_botrytis_p  - Tmin_botrytis_p)
        scores[i] <- base::max(
          ((t_i - Tmin_botrytis_p) / (Topt_botrytis_p - Tmin_botrytis_p)) *
            ((Tmax_botrytis_p - t_i) / (Tmax_botrytis_p - Topt_botrytis_p))^alpha,
          0)
      }
    }
    
    aggregate_infection_periods(scores,
                                min_score_t = min_score_botrytis)
  })
  
  
  # Downy mildew (P. destructor) — daily favorable-day scoring.
  #
  # Gate:  RH_max >= threshold (93-95%, proxy for nighttime humidity).
  # Score: Beta_T(Tmin) — sporulation occurs at night, Tmin 4-24°C.
  #
  # Ref: Jesperson & Sutton (1987) DOWNCAST
  #      Gilles et al. (2004) MILIONCAST
  get_downy_mildew_stress <- compiler::cmpfun(function(Tmin_daily,
                                                       RH_max_daily,
                                                       rh_mildew_threshold_p,
                                                       Topt_mildew_p,
                                                       Tmax_mildew_p,
                                                       Tmin_mildew_p,
                                                       min_score_mildew = 0.05) {
    n <- base::length(Tmin_daily)
    if (n == 0L) return(0)
    
    scores <- base::numeric(n)
    for (i in base::seq_len(n)) {
      t_i  <- Tmin_daily[i]
      rh_i <- RH_max_daily[i]
      if (base::is.na(t_i) || base::is.na(rh_i)) { scores[i] <- 0; next }
      if (rh_i < rh_mildew_threshold_p) { scores[i] <- 0; next }
      if (t_i <= Tmin_mildew_p || t_i >= Tmax_mildew_p) {
        scores[i] <- 0
      } else {
        alpha <- (Tmax_mildew_p - Topt_mildew_p) / (Topt_mildew_p - Tmin_mildew_p)
        scores[i] <- base::max(
          ((t_i - Tmin_mildew_p) / (Topt_mildew_p - Tmin_mildew_p)) *
            ((Tmax_mildew_p - t_i) / (Tmax_mildew_p - Topt_mildew_p))^alpha,
          0)
      }
    }
    
    aggregate_infection_periods(scores,
                                min_score_t = min_score_mildew)
  })
  
  # Fusarium basal rot (F. oxysporum f. sp. cepae) — daily favorable-day scoring.
  #
  # Gate:  One-sided upper moisture threshold (waterlogging only).
  #        A day is scored only when theta_rel > fusarium_theta_optimal_upper_p.
  #        Saturated soils favor chlamydospore germination and impair root defense.
  #
  #
  # Score: Beta_T(Ts_5cm) — soil temperature drives Fusarium activity.
  #
  # Ref: Abawi & Lorbeer (1972) Phytopathology 62:870-876.
  #      Cramer (2000) Plant Health Progress.
  get_fusarium_stress <- compiler::cmpfun(function(Ts_5cm_daily,
                                                   theta_rel_daily,
                                                   fusarium_theta_optimal_upper_p,
                                                   Topt_fusarium_p,
                                                   Tmax_fusarium_p,
                                                   Tmin_fusarium_p,
                                                   min_score_fusarium = 0.05) {
    n <- base::length(Ts_5cm_daily)
    if (n == 0L) return(0)
    
    scores <- base::numeric(n)
    for (i in base::seq_len(n)) {
      ts_i    <- Ts_5cm_daily[i]
      theta_i <- theta_rel_daily[i]
      if (base::is.na(ts_i))    { scores[i] <- 0; next }
      if (base::is.na(theta_i)) { scores[i] <- 0; next }
      
      # One-sided gate: only waterlogging (theta > upper threshold) triggers risk.
      # Drought stress is NOT modelled for NRW conditions — see rationale above.
      if (theta_i <= fusarium_theta_optimal_upper_p) {
        scores[i] <- 0; next
      }
      
      # Waterlogging confirmed: apply soil temperature response
      if (ts_i <= Tmin_fusarium_p || ts_i >= Tmax_fusarium_p) {
        scores[i] <- 0
      } else {
        alpha <- (Tmax_fusarium_p - Topt_fusarium_p) / (Topt_fusarium_p - Tmin_fusarium_p)
        scores[i] <- base::max(
          ((ts_i - Tmin_fusarium_p) / (Topt_fusarium_p - Tmin_fusarium_p)) *
            ((Tmax_fusarium_p - ts_i) / (Tmax_fusarium_p - Topt_fusarium_p))^alpha,
          0)
      }
    }
    
    aggregate_infection_periods(scores,
                                min_score_t = min_score_fusarium)
  })
  
  
  
  
  
  # ── Biomass accumulation (vectorized) ────────────────────────────────────────
  # Calculates daily biomass increment [g m⁻²] using:
  #   Biomass = LUE × f(CO2) × PAR × f(IPAR) × f(T) × f(W)
  #
  # where:
  #   LUE     = Radiation use efficiency [g MJ⁻¹ iPAR], phase-specific
  #   f(CO2)  = CO2 fertilization factor (Zhao et al. 2019, SIMPLE model)
  #   PAR     = photosynthetically active radiation [MJ m⁻² day⁻¹]
  #   f(IPAR) = 1 - exp(-k × LAI), Beer-Lambert light interception
  #   f(T)    = PAR-weighted daylight beta temperature response
  #   f(W)    = linear water stress function (theta_rel based)
  #
  # Temperature integration uses 24 hourly time steps:
  #   - Hourly temperature: sinusoidal diurnal model (Parton & Logan 1981)
  #   - Hourly PAR: sinusoidal clear-sky arc (Goudriaan 1986)
  #   - f(T) weighted by hourly PAR so night temperatures (PAR=0) are excluded
  #   - Beta function: Wang & Engel (1998) with cardinal temperatures
  #
  # Ref: Zhao et al. (2019) Eur. J. Agron. 104:97-106.
  #      Allen et al. (1998) FAO-56, Eq. 24, 33 (solar geometry).
  #      Parton & Logan (1981) Agric. Meteorol. 23:305-317.
  calc_bio_vectorized <- compiler::cmpfun(function(PAR, LAI, Tmax_day, Tmin_day, yday,
                                                   Windex,
                                                   Tmin_growth, Topt_growth, Tmax_growth,
                                                   f_W_1_lower, f_W_floor,
                                                   LUE_onion, lec_k,
                                                   f_CO2=1) {
    
    # PAR-weighted daylight f(T) integration (fully vectorized, matrix-based).
    # All operations use matrix broadcasting (n_days × 24); no R-level day loops.
    # See file header for scientific references.
    
    # --- Constants -----------------------------------------------------------
    lat_rad <- 51.0 * base::pi / 180   # NRW: Koeln-Bonn area (51 degrees N)
    h_tmax  <- 14.0                    # hour of Tmax: solar noon + 2 h lag
    hh      <- (0:23) + 0.5            # 24 hour-centre values
    
    # --- Scalar pre-computations ---------------------------------------------
    alpha_g <- (Tmax_growth - Topt_growth) / (Topt_growth - Tmin_growth)
    n_days  <- base::length(Tmax_day)
    valid   <- !base::is.na(Tmax_day) & !base::is.na(Tmin_day) & !base::is.na(yday)
    
    # --- Solar geometry: vectorised over days --------------------------------
    decl <- 0.409 * base::sin(2 * base::pi * yday / 365 - 1.39)
    ws   <- base::acos(base::pmax(-1, base::pmin(1,
                                                 -base::tan(lat_rad) * base::tan(decl))))
    DL   <- 24 * ws / base::pi           # daylength [h], length n_days
    t_r  <- 12 - DL / 2                  # sunrise
    t_s  <- 12 + DL / 2                  # sunset
    
    # --- Build n_days x 24 broadcast matrices --------------------------------
    hh_mat   <- base::matrix(hh,      nrow = n_days, ncol = 24L, byrow = TRUE)
    tr_mat   <- base::matrix(t_r,     nrow = n_days, ncol = 24L)
    ts_mat   <- base::matrix(t_s,     nrow = n_days, ncol = 24L)
    tmax_mat <- base::matrix(Tmax_day, nrow = n_days, ncol = 24L)
    tmin_mat <- base::matrix(Tmin_day, nrow = n_days, ncol = 24L)
    DL_mat   <- base::matrix(DL,      nrow = n_days, ncol = 24L)
    
    # --- Hourly temperature matrix [n_days x 24] -----------------------------
    rising  <- hh_mat >= tr_mat & hh_mat <= h_tmax
    falling <- hh_mat >  h_tmax & hh_mat <  ts_mat
    
    T_mat            <- tmin_mat   # night floor = Tmin (irrelevant: PAR=0)
    T_mat[rising]  <- (tmin_mat + (tmax_mat - tmin_mat) *
                         base::sin(base::pi / 2 * (hh_mat - tr_mat) / (h_tmax - tr_mat)))[rising]
    T_mat[falling] <- (tmin_mat + (tmax_mat - tmin_mat) *
                         base::sin(base::pi / 2 * (ts_mat - hh_mat) / (ts_mat - h_tmax)))[falling]
    T_mat[!valid, ] <- NA_real_
    
    # --- Hourly PAR matrix [n_days x 24] -------------------------------------
    daylight <- hh_mat >= tr_mat & hh_mat < ts_mat
    PAR_raw  <- base::matrix(0, nrow = n_days, ncol = 24L)
    PAR_raw[daylight] <- base::pmax(
      base::sin(base::pi * (hh_mat - tr_mat) / DL_mat)[daylight], 0)
    row_sums  <- base::rowSums(PAR_raw)
    safe_sums <- base::ifelse(row_sums > 0 & valid, row_sums, 1)
    PAR_mat   <- PAR_raw * (PAR / safe_sums)   # scaled to daily PAR total
    
    # --- Beta f(T) matrix [n_days x 24] --------------------------------------
    inside  <- T_mat > Tmin_growth & T_mat < Tmax_growth
    fT_mat  <- base::matrix(0, nrow = n_days, ncol = 24L)
    fT_mat[inside] <- base::pmax(
      ((T_mat  - Tmin_growth) / (Topt_growth - Tmin_growth)) *
        ((Tmax_growth - T_mat) / (Tmax_growth  - Topt_growth))^alpha_g,
      0)[inside]
    
    # --- PAR-weighted daily f(T) ---------------------------------------------
    # f_T[d] = sum_h[fT(h) * PAR(h)] / sum_h PAR(h)
    num_fT <- base::rowSums(fT_mat * PAR_mat)
    den_fT <- base::rowSums(PAR_mat)
    f_T    <- base::ifelse(valid & den_fT > 0, num_fT / den_fT, 0)
    
    # --- Water stress function -----------------------------------------------
    # Linear ramp from f_W_floor at wilting point to 1.0 above f_W_1_lower.
    # Floor > 0: retains minimal residual conductance at wilting (Brewster 2008).
    Windex <- base::ifelse(base::is.na(Windex), f_W_1_lower, Windex)
    f_W    <- base::ifelse(
      Windex >= f_W_1_lower,
      1.0,
      f_W_floor + (1.0 - f_W_floor) * (Windex / f_W_1_lower)
    )
    f_W <- base::pmin(base::pmax(f_W, f_W_floor), 1.0)
    
    # --- CO2 fertilisation effect -------------------------------------------
    # f_CO2 scales LUE via the SIMPLE model formulation (Zhao et al., 2019):
    #   f_CO2 = 1 + S_co2 * (CO2_ppm - CO2_ref)   for CO2 < 700 ppm
    #   f_CO2 = 1 + S_co2 * (700 - CO2_ref)        for CO2 >= 700 ppm (saturation)
    # CO2_ref = co2_ppm_historical (~415 ppm): ensures f_CO2 = 1.0 for baseline.
    # f_CO2 is computed externally per scenario and passed in as a scalar.
    # Default f_CO2 = 1.0 (no effect) if scenario not matched.
    # C3 crop (onion): S_co2 ~ 0.0006-0.0008 per ppm (Zhao et al., 2019).
    LUE_onion * f_CO2 * PAR * (1 - base::exp(-lec_k * LAI)) * f_T * f_W
  })
  
  
  # Collect helpers and optionally attach to global environment ----
  
  helpers <- base::list(
    sat                          = sat,
    expected_loss                = expected_loss,
    fast_chance_event            = fast_chance_event,
    consec_counter               = consec_counter,
    aggregate_infection_periods  = aggregate_infection_periods,
    get_extreme_rain_stress      = get_extreme_rain_stress,
    # get_hail_stress removed — hail is a pure chance_event
    # get_heat_stress removed — temperature captured by Beta f_T in biomass
    get_botrytis_stress          = get_botrytis_stress,
    get_downy_mildew_stress      = get_downy_mildew_stress,
    get_fusarium_stress          = get_fusarium_stress,
    calc_bio_vectorized          = calc_bio_vectorized
  )
  
  if (attach_to_global) {
    base::list2env(helpers, envir = base::.GlobalEnv)
    base::message("✅ Helper functions (DAILY SCORING) attached to global environment.")
  }
  
  helpers
}


# Attach helpers for interactive use ----

helper_function(attach_to_global = TRUE)
