################################################################################
# 03 — Onion Climate Impact Model (Main Model Function)
# ──────────────────────────────────────────────────────
# Purpose:
#   Implements the core onion climate impact model for Monte Carlo simulation.
#   For each MC iteration and each climate scenario, the model:
#     1. Randomly samples one growing season from the weather data
#     2. Computes soil water balance (rainfed & irrigated)
#     3. Calculates biomass accumulation via RUE × PAR × f(T) × f(W) × f(CO2)
#     4. Scores biotic and abiotic stress risks
#     5. Applies stochastic yield loss gates
#     6. Returns yield and diagnostic outputs per scenario
#
# Water balance:
#   Geisenheimer Bewässerungssteuerung — a FAO-56 Kc-ET0 approach with
#   nFK (nutzbare Feldkapazität) as plant-available water in the root zone.
#   Ref: Schmidt & Zinkernagel (2017) Water 9:693.
#        Zinkernagel et al. (2022) Hochschule Geisenheim University.
#        Allen et al. (1998) FAO-56.
#
# Growth phases (GDD-based):
#   Emergence → Vegetative → Bulbing → Maturation
#   Phase transitions defined by cumulative GDD thresholds.
#   Ref: Schmidt & Zinkernagel (2017), Table 2.
#
# Dependencies:
#   01_weather_preprocessing.R  (provides weather_precomputed)
#   02_helper_functions.R       (provides stress & biomass functions)
#   input_irrig1.csv            (Monte Carlo parameter distributions)
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

load_if_needed(base::c(
  "zoo",
  "RcppRoll",
  "decisionSupport",
  "compiler",
  "data.table"
))

compiler::enableJIT(3)   # Enable JIT compilation (level 3 = aggressive)


# ── Phase definition via cumulative GDD ──────────────────────────────────────
# Assigns each day to a growth phase based on cumulative GDD after sowing.
# Phase transitions follow Schmidt & Zinkernagel (2017), Table 2.
precompute_phase_data <- compiler::cmpfun(function(df, sow,
                                                   GDD_field_emergence_required_p,
                                                   GDD_vegetative_required_p,
                                                   GDD_bulbing_required_p,
                                                   GDD_maturation_required_p) {
  n <- base::nrow(df)
  if (n == 0L) base::return(NULL)
  
  # Initialize cumulative GDD
  df[, GDD_cum := 0]
  # Only accumulate after sowing
  after_sow <- df$yday >= sow
  df$GDD_cum[after_sow] <- base::cumsum(df$GDD_daily[after_sow])
  
  # Phase thresholds (GDD-based)
  em_idx <- after_sow & df$GDD_cum <= GDD_field_emergence_required_p
  vg_idx <- after_sow & df$GDD_cum > GDD_field_emergence_required_p &
    df$GDD_cum <= GDD_vegetative_required_p
  bl_idx <- after_sow & df$GDD_cum > GDD_vegetative_required_p &
    df$GDD_cum <= GDD_bulbing_required_p
  mt_idx <- after_sow & df$GDD_cum > GDD_bulbing_required_p &
    df$GDD_cum <= GDD_maturation_required_p
  
  base::list(
    em = base::list(idx = em_idx, range = base::which(em_idx)),
    vg = base::list(idx = vg_idx, range = base::which(vg_idx)),
    bl = base::list(idx = bl_idx, range = base::which(bl_idx)),
    mt = base::list(idx = mt_idx, range = base::which(mt_idx))
  )
})


# ── Batch risk computation per phase ──────────────────────────────────────────
# Computes all biotic and abiotic stress scores for each growth phase.
# Returns a list: list(em = risks_em, vg = risks_vg, bl = risks_bl, mt = risks_mt)
compute_all_risks <- compiler::cmpfun(function(phase_data, df, params) {
  
  # Helper: safe subset returning zero-length numeric if idx is empty
  safe_sub <- function(col, idx) {
    if (is.null(col) || !base::any(idx, na.rm = TRUE)) return(base::numeric(0L))
    col[idx]
  }
  
  # Soil temperature: use Ts_5cm if available, otherwise estimate from Tavg
  # Soil temp at 5cm lags and dampens air temp; rough approx: Ts ≈ 0.8*Tavg + 3
  safe_soil_temp <- function(df, idx) {
    if (!is.null(df$Ts_5cm)) return(safe_sub(df$Ts_5cm, idx))
    tavg <- safe_sub(df$Tavg, idx)
    if (length(tavg) == 0L) return(base::numeric(0L))
    0.8 * tavg + 3
  }
  
  # --- EMERGENCE phase risks ---
  em_idx <- phase_data$em$idx
  risks_em <- base::c(
    exrain = get_extreme_rain_stress(
      df$Prec[em_idx],
      params$prec_extreme_rain_medium_p,
      params$prec_extreme_rain_high_p,
      params$impact_days_extreme_rain_t
    )
    # Fusarium basal rot excluded from emergence phase.
    # F. oxysporum f. sp. cepae attacks the stem plate and basal bulb tissue
    # as the bulb develops — it is a bulbing/maturation-phase disease.
    # em-phase theta_rel is structurally elevated (mean ~0.92) due to
    # post-sowing soil moisture and near-zero ET demand, which would
    # spuriously inflate risk_fusarium regardless of pathogen biology.
    # Ref: Abawi & Lorbeer (1972) Phytopathology 62:870-876.
  )
  
  # --- VEGETATIVE phase risks ---
  vg_idx <- phase_data$vg$idx
  risks_vg <- base::c(
    exrain = get_extreme_rain_stress(
      df$Prec[vg_idx],
      params$prec_extreme_rain_medium_p,
      params$prec_extreme_rain_high_p,
      params$impact_days_extreme_rain_t
    ),
    mildew = get_downy_mildew_stress(
      safe_sub(df$Tmin,   vg_idx),
      safe_sub(df$RH_max, vg_idx),
      params$rh_mildew_threshold_p,
      params$Topt_mildew_p,
      params$Tmax_mildew_p,
      params$Tmin_mildew_p
    ),
    fusarium = get_fusarium_stress(
      safe_soil_temp(df, vg_idx),
      safe_sub(df$theta_rel_rf, vg_idx),
      params$fusarium_theta_optimal_upper_p,
      params$Topt_fusarium_p,
      params$Tmax_fusarium_p,
      params$Tmin_fusarium_p
    )
  )
  
  # --- BULBING phase risks ---
  bl_idx <- phase_data$bl$idx
  risks_bl <- base::c(
    exrain = get_extreme_rain_stress(
      df$Prec[bl_idx],
      params$prec_extreme_rain_medium_p,
      params$prec_extreme_rain_high_p,
      params$impact_days_extreme_rain_t
    ),
    botrytis = get_botrytis_stress(
      safe_sub(df$Tmin, bl_idx),
      safe_sub(df$Prec, bl_idx),
      params$prec_botrytis_wet_p,
      params$Topt_botrytis_p,
      params$Tmax_botrytis_p,
      params$Tmin_botrytis_p
    ),
    mildew = get_downy_mildew_stress(
      safe_sub(df$Tmin,   bl_idx),
      safe_sub(df$RH_max, bl_idx),
      params$rh_mildew_threshold_p,
      params$Topt_mildew_p,
      params$Tmax_mildew_p,
      params$Tmin_mildew_p
    ),
    fusarium = get_fusarium_stress(
      safe_soil_temp(df, bl_idx),
      safe_sub(df$theta_rel_rf, bl_idx),
      params$fusarium_theta_optimal_upper_p,
      params$Topt_fusarium_p,
      params$Tmax_fusarium_p,
      params$Tmin_fusarium_p
    )
  )
  
  # --- MATURATION phase risks ---
  mt_idx <- phase_data$mt$idx
  risks_mt <- base::c(
    exrain = get_extreme_rain_stress(
      df$Prec[mt_idx],
      params$prec_extreme_rain_medium_p,
      params$prec_extreme_rain_high_p,
      params$impact_days_extreme_rain_t
    ),
    botrytis = get_botrytis_stress(
      safe_sub(df$Tmin, mt_idx),
      safe_sub(df$Prec, mt_idx),
      params$prec_botrytis_wet_p,
      params$Topt_botrytis_p,
      params$Tmax_botrytis_p,
      params$Tmin_botrytis_p
    ),
    fusarium = get_fusarium_stress(
      safe_soil_temp(df, mt_idx),
      safe_sub(df$theta_rel_rf, mt_idx),
      params$fusarium_theta_optimal_upper_p,
      params$Topt_fusarium_p,
      params$Tmax_fusarium_p,
      params$Tmin_fusarium_p
    ),
    mildew = get_downy_mildew_stress(
      safe_sub(df$Tmin,   mt_idx),
      safe_sub(df$RH_max, mt_idx),
      params$rh_mildew_threshold_p,
      params$Topt_mildew_p,
      params$Tmax_mildew_p,
      params$Tmin_mildew_p
    )
  )
  
  base::list(em = risks_em, vg = risks_vg, bl = risks_bl, mt = risks_mt)
})


# ── Load helper functions (stress + biomass) from 02_helper_functions.R ───────
helper_function(attach_to_global = TRUE)


# ── Geisenheimer Bewässerungssteuerung ───────────────────────────────────────
# Irrigation scheduling algorithm based on climatic water balance.
#
# Algorithm:
#   1. Daily deficit = ET0 × kc × Ks − Precipitation
#      Ks = FAO-56 stress coefficient (reduces ETc when soil dries below RAW)
#      RAW = p × TAW (p = 0.25 for onion, FAO-56 Table 22)
#   2. Cumulative deficit (Gesamtdefizit) tracked daily
#   3. Irrigation triggered when Gesamtdefizit >= Einzelgabe + Puffer
#      Einzelgabe = 30% × TAW (refills from 60% to 90% nFK)
#      Puffer     = 10% × TAW (buffer zone, not refilled)
#   4. theta_rel = 1 − Gesamtdefizit / TAW (for water stress function f(W))
#
# Ref: Schmidt & Zinkernagel (2017) Water 9:693.
#      Zinkernagel et al. (2022) Geisenheimer Bewässerungssteuerung.
#      Allen et al. (1998) FAO-56.
compute_irrigation_geisenheim <- compiler::cmpfun(function(df, sow_yday, params, phase_data) {
  n <- base::nrow(df)
  irrig     <- numeric(n)
  theta_rel <- rep(NA_real_, n)
  
  # map each day to a numeric phase: 1=em, 2=vg, 3=bl, 4=mt
  phase_id <- integer(n)
  phase_id[phase_data$em$idx] <- 1L
  phase_id[phase_data$vg$idx] <- 2L
  phase_id[phase_data$bl$idx] <- 3L
  phase_id[phase_data$mt$idx] <- 4L
  
  # per-phase kc and rooting depth (cm)
  kc_phase   <- base::c(params$kc_em, params$kc_vg, params$kc_bl, params$kc_mt)
  root_phase <- base::c(params$root_em_cm, params$root_vg_cm, params$root_bl_cm, params$root_mt_cm)
  
  # soil parameters
  nFK_mm_per_cm <- params$nFK_mm_per_cm
  max_irrig_mm  <- params$max_irrig_mm
  
  # state variables
  Gesamtdefizit <- 0    # accumulated climatic water deficit (mm)
  prev_TAW      <- NA_real_
  
  for (i in base::seq_len(n)) {
    # before sowing: skip
    if (df$yday[i] < sow_yday) next
    
    pid <- phase_id[i]
    if (pid == 0L) next  # outside em/vg/bl/mt
    
    kc   <- kc_phase[pid]
    root <- root_phase[pid]
    
    # Handle maturation phase: kc=0 means no water demand, no irrigation
    if (base::is.na(kc) || base::is.na(root) || root <= 0) next
    
    # Total available water (TAW) in current root zone
    TAW <- root * nFK_mm_per_cm
    
    # Rescale deficit when rooting depth changes (phase transition)
    if (!base::is.na(prev_TAW) && prev_TAW > 0 && TAW > 0 && 
        !base::isTRUE(all.equal(TAW, prev_TAW))) {
      Gesamtdefizit <- Gesamtdefizit * (TAW / prev_TAW)
    }
    prev_TAW <- TAW
    
    # --- Geisenheimer thresholds ---
    Einzelgabe <- 0.30 * TAW   # refill amount: 60% -> 90% nFK
    Puffer     <- 0.10 * TAW   # buffer zone:   90% -> 100% nFK (not refilled)
    
    # --- Step 1: Daily deficit accumulation ---
    # FAO-56 Ks: reduce actual ETc when soil dries below RAW (p=0.25 for onion)
    # This prevents unrealistic deficit accumulation under rainfed conditions.
    # When plants are stressed they close stomata → actual ET drops.
    RAW <- 0.25 * TAW   # readily available water (FAO-56 Table 22: p=0.25 onion)
    Dr  <- Gesamtdefizit # current root zone depletion
    if (Dr > RAW && TAW > RAW) {
      Ks <- base::max(0, (TAW - Dr) / (TAW - RAW))
    } else {
      Ks <- 1.0
    }
    ETc  <- df$ET0_mm[i] * kc * Ks
    Peff <- df$Prec[i]
    if (base::is.na(Peff)) Peff <- 0
    
    Tagesdefizit  <- ETc - Peff
    Gesamtdefizit <- Gesamtdefizit + Tagesdefizit
    
    # --- Step 2: Irrigation decision ---
    irrigation_today <- 0
    if (kc > 0 && Gesamtdefizit >= (Einzelgabe + Puffer)) {
      # Apply one Einzelgabe (capped at max_irrig_mm)
      gift <- base::min(Einzelgabe, max_irrig_mm)
      irrigation_today <- gift
      Gesamtdefizit    <- Gesamtdefizit - gift
    }
    irrig[i] <- irrigation_today
    
    # --- Step 3: Constrain deficit ---
    # Heavy rain can push deficit negative => excess drains, reset to 0
    if (Gesamtdefizit < 0) Gesamtdefizit <- 0
    # Cannot lose more water than TAW
    if (Gesamtdefizit > TAW) Gesamtdefizit <- TAW
    
    # --- Step 4: Relative soil moisture for stress function ---
    theta_rel[i] <- 1 - Gesamtdefizit / TAW
  }
  
  base::list(irrig_mm = irrig, theta_rel = theta_rel)
})


# ── MAIN MODEL FUNCTION ──────────────────────────────────────────────────────
# Called by decisionSupport::mcSimulation() for each Monte Carlo iteration.
# Parameters are drawn from input_irrig1.csv distributions and injected
# into the function scope by decisionSupport (plainNames syntax).

onion_climate_impact <- compiler::cmpfun(function() {
  
  # 1) For each scenario, randomly sample one season (one id_season) ----
  weather_scenario_list <- base::lapply(weather_precomputed, function(list_per_scenario) {
    sample(list_per_scenario, 1)[[1]]
  })
  
  # 2) Bundle parameters into a single list for passing around ----
  params <- base::list(
    # irrigation (soil water & kc / root by phase)
    nFK_mm_per_cm      = nfk_mm_per_cm_p,
    max_irrig_mm       = max_irrig_mm_p,
    
    root_em_cm         = root_em_cm_p,
    root_vg_cm         = root_vg_cm_p,
    root_bl_cm         = root_bl_cm_p,
    root_mt_cm         = root_mt_cm_p,
    
    kc_em              = kc_em_p,
    kc_vg              = kc_vg_p,
    kc_bl              = kc_bl_p,
    kc_mt              = kc_mt_p,
    
    LAI_emergence  = LAI_emergence_p,
    LAI_veg        = LAI_veg_p,
    LAI_bulbing    = LAI_bulbing_p,
    LAI_maturation = LAI_maturation_p,
    
    
    # Phase-specific LUE values
    LUE_em_p           = LUE_em_p,
    LUE_vg_p           = LUE_vg_p,
    LUE_bl_p           = LUE_bl_p,
    LUE_mt_p           = LUE_mt_p,
    
    # Temperature response (Beta curve)
    Tmin_growth_p   = Tmin_growth_p,
    Topt_growth_p   = Topt_growth_p,
    Tmax_growth_p   = Tmax_growth_p,
    
    # Water response
    f_W_1_lower_t   = f_W_1_lower_t,
    f_W_floor_t     = f_W_floor_t,
    
    lec_k_c            = lec_k_c,
    HI_onions_t        = HI_onions_t,
    
    # extreme rain
    prec_extreme_rain_medium_p = prec_extreme_rain_medium_p,
    prec_extreme_rain_high_p   = prec_extreme_rain_high_p,
    impact_days_extreme_rain_t = impact_days_extreme_rain_t,
    
    # hail: pure chance_event — no proxy params needed (uses chance_hail_* directly)
    
    # heat: removed — temperature effect already captured by Beta f_T in biomass module
    
    # botrytis
    prec_botrytis_wet_p        = prec_botrytis_wet_p,
    Topt_botrytis_p            = Topt_botrytis_p,
    Tmin_botrytis_p            = Tmin_botrytis_p,
    Tmax_botrytis_p            = Tmax_botrytis_p,
    
    # mildew
    rh_mildew_threshold_p      = rh_mildew_threshold_p,
    Topt_mildew_p              = Topt_mildew_p,
    Tmin_mildew_p              = Tmin_mildew_p,
    Tmax_mildew_p              = Tmax_mildew_p,
    
    # fusarium
    fusarium_theta_optimal_upper_p = fusarium_theta_optimal_upper_p,  # BUG FIX: was incorrectly assigned lower_p
    Topt_fusarium_p                = Topt_fusarium_p,
    Tmin_fusarium_p                = Tmin_fusarium_p,
    Tmax_fusarium_p                = Tmax_fusarium_p,
    
    # onion fly excluded — see 02_new.R changelog for scientific rationale
    
    # yield losses
    yield_reduction_extreme_rain_t = yield_reduction_extreme_rain_t,
    yield_reduction_fusarium_t     = yield_reduction_fusarium_t,
    yield_reduction_botrytis_t     = yield_reduction_botrytis_t,
    yield_reduction_downy_mildew_t = yield_reduction_downy_mildew_t
  )
  
  results <- vector("list", base::length(weather_scenario_list))
  base::names(results) <- base::names(weather_scenario_list)
  
  # 4) Loop over scenarios ----
  scenario_count <- 0
  for (sc in base::names(weather_scenario_list)) {
    
    scenario_count <- scenario_count + 1
    
    df <- weather_scenario_list[[sc]]
    data.table::setDT(df)
    
    # Derived streak variables
    # (No pre-computed streak variables needed — all stressors use daily scoring)
    
    # ── CO2 fertilization factor ─────────────────────────────────────────────
    # SIMPLE model (Zhao et al. 2019, Eq. 10):
    #   f_CO2 = 1 + S_co2 × (CO2 - CO2_ref)  for CO2 < 700 ppm
    #   f_CO2 = 1 + S_co2 × (700 - CO2_ref)  for CO2 >= 700 ppm (saturation)
    # S_co2 = 0.0006-0.0008 per ppm for C3 crops (onion is C3).
    # CO2_ref = ~415 ppm (2020 baseline).
    # Ref: Zhao et al. (2019) Eur. J. Agron. 104:97-106.
    #      Ainsworth & Long (2005) New Phytol. 165:351-372.
    co2_ppm_scenario <- if (grepl("ssp119|ssp126|ssp1", sc, ignore.case = TRUE)) {
      co2_ppm_ssp126
    } else if (grepl("ssp245|ssp2", sc, ignore.case = TRUE)) {
      co2_ppm_ssp245
    } else if (grepl("ssp370|ssp3", sc, ignore.case = TRUE)) {
      co2_ppm_ssp370
    } else if (grepl("ssp585|ssp5", sc, ignore.case = TRUE)) {
      co2_ppm_ssp585
    } else {
      co2_ppm_historical   # historical ~415 ppm -> f_CO2 ~ 1.045 (small effect)
    }
    co2_effect <- base::min(co2_ppm_scenario, 700) - co2_ppm_historical
    f_CO2 <- 1 + S_co2_p * base::max(co2_effect, 0)
    f_CO2 <- base::max(f_CO2, 1.0)   # floor: no negative CO2 effect below reference
    
    # Sowing date
    sow <- pmax(base::round(sowing_day_p), 45)
    
    # Phase data
    phase_data <- precompute_phase_data(
      df, sow,
      GDD_field_emergence_required_p,
      GDD_vegetative_required_p,
      GDD_bulbing_required_p,
      GDD_maturation_required_p
    )
    
    # Water balance: RAINFED
    params_rf <- params
    params_rf$max_irrig_mm <- 0
    res_rf <- compute_irrigation_geisenheim(df, sow, params_rf, phase_data)
    df$theta_rel_rf <- res_rf$theta_rel
    
    # Water balance: IRRIGATED
    res_ir <- compute_irrigation_geisenheim(df, sow, params, phase_data)
    df$Irrig_mm    <- res_ir$irrig_mm
    df$theta_rel_ir <- res_ir$theta_rel
    
    
    
    # ================================================================
    # BIOTIC RISK SCORING — separate rainfed / irrigated theta
    #
    # Stressors that depend on soil moisture (Fusarium) must use the
    # theta from the matching water balance run. All other stressors
    # (mildew, botrytis, exrain) are weather-driven and do not
    # depend on theta, so a single call suffices for those.
    # ================================================================
    
    # Rainfed risks: uses df$theta_rel_rf (set above)
    all_risks_rf <- compute_all_risks(phase_data, df, params)
    
    # Irrigated risks: temporarily swap theta column, score, then restore
    df$theta_rel_rf <- df$theta_rel_ir
    all_risks_ir    <- compute_all_risks(phase_data, df, params)
    df$theta_rel_rf <- res_rf$theta_rel   # restore rainfed theta
    
    # ── Seasonal risk = max daily score across worst-exposed phase ───────
    max_phase_risk <- function(risks, phases, name) {
      vals <- sapply(phases, function(ph) {
        v <- risks[[ph]][name]
        if (is.null(v) || is.na(v)) -Inf else v
      })
      base::min(1, base::max(0, base::max(vals, na.rm = TRUE)))
    }
    
    # Shared weather-driven risks (identical for rainfed and irrigated)
    risk_mildew   <- max_phase_risk(all_risks_rf, c("vg","bl","mt"),       "mildew")
    risk_botrytis <- max_phase_risk(all_risks_rf, c("bl","mt"),            "botrytis")
    risk_exrain   <- max_phase_risk(all_risks_rf, c("em","vg","bl","mt"),  "exrain")
    
    # Fusarium: soil-moisture dependent — separate score per water balance
    risk_fusarium_rf <- max_phase_risk(all_risks_rf, c("vg","bl","mt"), "fusarium")  # em excluded: see risks_em
    risk_fusarium_ir <- max_phase_risk(all_risks_ir, c("vg","bl","mt"), "fusarium")  # em excluded
    
    # Headline risk_fusarium = rainfed (more conservative; used in output columns)
    risk_fusarium <- risk_fusarium_rf
    
    # ── Unified yield loss pool ─────────────────────────────────────────────
    # All stressors follow the same principle:
    #   1. Risk score (0-1) from daily environmental scoring
    #   2. Stochastic gate: chance_event(risk, yield_loss_fraction)
    #   3. Loss = loss_fraction if event fires, else 0
    
    # Stochastic gate: does each stressor cause loss this season?
    # yield_reduction values = fraction of yield LOST (e.g. 0.30 = 30% loss)
    # Shared stressors (weather-driven, same for both water balance scenarios)
    loss_mildew   <- chance_event(chance = risk_mildew,       value_if = params$yield_reduction_downy_mildew_t, value_if_not = 0)
    loss_botrytis <- chance_event(chance = risk_botrytis,     value_if = params$yield_reduction_botrytis_t,     value_if_not = 0)
    loss_exrain   <- chance_event(chance = risk_exrain,       value_if = params$yield_reduction_extreme_rain_t, value_if_not = 0)
    
    # Fusarium: separate loss per water balance (irrigated fields have lower risk)
    loss_fusarium_rf <- chance_event(chance = risk_fusarium_rf, value_if = params$yield_reduction_fusarium_t, value_if_not = 0)
    loss_fusarium_ir <- chance_event(chance = risk_fusarium_ir, value_if = params$yield_reduction_fusarium_t, value_if_not = 0)
    
    # Headline loss_fusarium = rainfed (used in output columns for consistency)
    loss_fusarium <- loss_fusarium_rf
    
    # Hail: pure chance_event (no daily risk scoring, just annual probability)
    # yield_reduction_hail = fraction of yield LOST if hail occurs
    if (grepl("historical", sc, ignore.case = TRUE)) {
      loss_hail <- chance_event(
        chance       = chance_hail_historical,
        value_if     = yield_reduction_hail_historical,
        value_if_not = 0
      )
    } else {
      loss_hail <- chance_event(
        chance       = chance_hail_future,
        value_if     = yield_reduction_hail_future,
        value_if_not = 0
      )
    }
    
    # Total loss: ALL stressors additive, capped
    # Rainfed uses rainfed Fusarium loss; irrigated uses lower irrigated Fusarium loss
    shared_loss <- loss_mildew + loss_botrytis + loss_exrain + loss_hail
    
    total_stress_loss_rf <- base::min(
      shared_loss + loss_fusarium_rf
    )
    total_stress_loss_ir <- base::min(
      shared_loss + loss_fusarium_ir
    )
    
    # Headline total_stress_loss = rainfed (used in output columns)
    total_stress_loss <- total_stress_loss_rf
    
    # Separate yield multipliers per water balance scenario
    stress_yield_multiplier_rf <- 1 - total_stress_loss_rf
    stress_yield_multiplier_ir <- 1 - total_stress_loss_ir
    
    # Headline multiplier = rainfed (for output column backward compatibility)
    stress_yield_multiplier <- stress_yield_multiplier_rf
    
    # ── Biomass calculation (abiotic only: T + W + PAR + CO2) ───────────────
    
    # RAINFED biomass calculation
    water_rf <- df$theta_rel_rf
    
    biomass_em_pot_rf <- base::sum(calc_bio_vectorized(
      df$PAR[phase_data$em$idx], params$LAI_emergence,
      df$Tmax[phase_data$em$idx], df$Tmin[phase_data$em$idx], df$yday[phase_data$em$idx], water_rf[phase_data$em$idx],
      params$Tmin_growth_p, params$Topt_growth_p, params$Tmax_growth_p,
      params$f_W_1_lower_t, params$f_W_floor_t,
      params$LUE_em_p, params$lec_k_c, f_CO2
    ))
    
    biomass_vg_pot_rf <- base::sum(calc_bio_vectorized(
      df$PAR[phase_data$vg$idx], params$LAI_veg,
      df$Tmax[phase_data$vg$idx], df$Tmin[phase_data$vg$idx], df$yday[phase_data$vg$idx], water_rf[phase_data$vg$idx],
      params$Tmin_growth_p, params$Topt_growth_p, params$Tmax_growth_p,
      params$f_W_1_lower_t, params$f_W_floor_t,
      params$LUE_vg_p, params$lec_k_c, f_CO2
    ))
    
    biomass_bl_pot_rf <- base::sum(calc_bio_vectorized(
      df$PAR[phase_data$bl$idx], params$LAI_bulbing,
      df$Tmax[phase_data$bl$idx], df$Tmin[phase_data$bl$idx], df$yday[phase_data$bl$idx], water_rf[phase_data$bl$idx],
      params$Tmin_growth_p, params$Topt_growth_p, params$Tmax_growth_p,
      params$f_W_1_lower_t, params$f_W_floor_t,
      params$LUE_bl_p, params$lec_k_c, f_CO2
    ))
    
    biomass_mt_pot_rf <- base::sum(calc_bio_vectorized(
      df$PAR[phase_data$mt$idx], params$LAI_maturation,
      df$Tmax[phase_data$mt$idx], df$Tmin[phase_data$mt$idx], df$yday[phase_data$mt$idx], water_rf[phase_data$mt$idx],
      params$Tmin_growth_p, params$Topt_growth_p, params$Tmax_growth_p,
      params$f_W_1_lower_t, params$f_W_floor_t,
      params$LUE_mt_p, params$lec_k_c, f_CO2
    ))
    
    # Attainable biomass (abiotic only: T + W + PAR)
    attainable_biomass_g_m2_rf <- biomass_em_pot_rf + biomass_vg_pot_rf +
      biomass_bl_pot_rf + biomass_mt_pot_rf
    
    # Convert to yield, then apply biotic loss ONCE.
    # Hard cap at 150 t/ha fresh weight — the absolute physiological ceiling
    # for onion. Values above this are MC tail artefacts from stochastic
    # parameter combinations at their upper bounds (high RUE + high CO2 + long season).
    MAX_YIELD_FW_T_HA        <- 150
    attainable_biomass_t_ha_rf  <- attainable_biomass_g_m2_rf * 0.01
    attainable_yield_DM_t_ha_rf <- base::min(
      attainable_biomass_t_ha_rf * params$HI_onions_t,
      MAX_YIELD_FW_T_HA * dry_matter_content_t   # cap in DM space
    )
    final_yield_DM_t_ha_rf <- attainable_yield_DM_t_ha_rf * stress_yield_multiplier_rf
    
    # IRRIGATED biomass calculation
    water_ir <- df$theta_rel_ir
    
    biomass_em_pot_ir <- base::sum(calc_bio_vectorized(
      df$PAR[phase_data$em$idx], params$LAI_emergence,
      df$Tmax[phase_data$em$idx], df$Tmin[phase_data$em$idx], df$yday[phase_data$em$idx], water_ir[phase_data$em$idx],
      params$Tmin_growth_p, params$Topt_growth_p, params$Tmax_growth_p,
      params$f_W_1_lower_t, params$f_W_floor_t,
      params$LUE_em_p, params$lec_k_c, f_CO2
    ))
    
    biomass_vg_pot_ir <- base::sum(calc_bio_vectorized(
      df$PAR[phase_data$vg$idx], params$LAI_veg,
      df$Tmax[phase_data$vg$idx], df$Tmin[phase_data$vg$idx], df$yday[phase_data$vg$idx], water_ir[phase_data$vg$idx],
      params$Tmin_growth_p, params$Topt_growth_p, params$Tmax_growth_p,
      params$f_W_1_lower_t, params$f_W_floor_t,
      params$LUE_vg_p, params$lec_k_c, f_CO2
    ))
    
    biomass_bl_pot_ir <- base::sum(calc_bio_vectorized(
      df$PAR[phase_data$bl$idx], params$LAI_bulbing,
      df$Tmax[phase_data$bl$idx], df$Tmin[phase_data$bl$idx], df$yday[phase_data$bl$idx], water_ir[phase_data$bl$idx],
      params$Tmin_growth_p, params$Topt_growth_p, params$Tmax_growth_p,
      params$f_W_1_lower_t, params$f_W_floor_t,
      params$LUE_bl_p, params$lec_k_c, f_CO2
    ))
    
    biomass_mt_pot_ir <- base::sum(calc_bio_vectorized(
      df$PAR[phase_data$mt$idx], params$LAI_maturation,
      df$Tmax[phase_data$mt$idx], df$Tmin[phase_data$mt$idx], df$yday[phase_data$mt$idx], water_ir[phase_data$mt$idx],
      params$Tmin_growth_p, params$Topt_growth_p, params$Tmax_growth_p,
      params$f_W_1_lower_t, params$f_W_floor_t,
      params$LUE_mt_p, params$lec_k_c, f_CO2
    ))
    
    # Attainable biomass (abiotic only)
    attainable_biomass_g_m2_ir <- biomass_em_pot_ir + biomass_vg_pot_ir +
      biomass_bl_pot_ir + biomass_mt_pot_ir
    
    # Convert to yield, then apply biotic loss ONCE (same cap as rainfed)
    attainable_biomass_t_ha_ir  <- attainable_biomass_g_m2_ir * 0.01
    attainable_yield_DM_t_ha_ir <- base::min(
      attainable_biomass_t_ha_ir * params$HI_onions_t,
      MAX_YIELD_FW_T_HA * dry_matter_content_t   # cap in DM space
    )
    final_yield_DM_t_ha_ir <- attainable_yield_DM_t_ha_ir * stress_yield_multiplier_ir
    
    
    # ── Diagnostic: heat stress and drought stress ─────────────────────────
    # These diagnostics do NOT feed back into yield. They quantify the
    # fraction of growing-season days under stress and biomass foregone.
    # Heat stress:    Tmax >= Tmax_growth (f(T) = 0, growth halted)
    # Drought stress: theta_rel < 0.75 (FAO-56: p = 0.25 for onion)
    
    # Growing-season day indices (all phases combined)
    gs_idx   <- phase_data$em$idx | phase_data$vg$idx |
      phase_data$bl$idx | phase_data$mt$idx
    gs_Tmax  <- df$Tmax[gs_idx]
    n_gs     <- base::max(base::sum(gs_idx, na.rm = TRUE), 1L)
    
    # ── HEAT STRESS ──────────────────────────────────────────────────────────
    # Heat stress day: Tmax >= Tmax_growth_p — the daily peak temperature
    # reaches or exceeds the biological ceiling, meaning the crop experiences
    # at least some hours with f(T) = 0 in the sinusoidal diurnal integration.
    # Using Tmax (not Tavg) is now consistent with the biomass calculation,
    # which integrates f(T) across the full diurnal cycle.
    # Ref: Brewster (2008); McMaster et al. (2008).
    
    heat_days_mask_gs <- !base::is.na(gs_Tmax) & gs_Tmax >= params$Tmax_growth_p
    heat_stress_day_fraction <- base::sum(heat_days_mask_gs) / n_gs
    # (identical for ir and rf — temperature does not depend on irrigation)
    
    # Biomass lost on heat-stress days: counterfactual increment if f(T) = 1
    # on those days (holding actual f(W) in place so only heat is attributed).
    heat_biomass_lost_ir <- 0
    heat_biomass_lost_rf <- 0
    
    for (ph in base::list(
      base::list(idx = phase_data$em$idx, LAI = params$LAI_emergence, LUE = params$LUE_em_p),
      base::list(idx = phase_data$vg$idx, LAI = params$LAI_veg,       LUE = params$LUE_vg_p),
      base::list(idx = phase_data$bl$idx, LAI = params$LAI_bulbing,   LUE = params$LUE_bl_p),
      base::list(idx = phase_data$mt$idx, LAI = params$LAI_maturation,LUE = params$LUE_mt_p)
    )) {
      if (!base::any(ph$idx, na.rm = TRUE)) next
      
      ph_Tmax <- df$Tmax[ph$idx]
      ph_PAR  <- df$PAR[ph$idx]
      ph_heat <- !base::is.na(ph_Tmax) & ph_Tmax >= params$Tmax_growth_p
      if (!base::any(ph_heat)) next
      
      # Light interception term (same for both water scenarios)
      light_int <- ph$LUE * ph_PAR * (1 - base::exp(-params$lec_k_c * ph$LAI))
      
      # f(W) irrigated on heat-stress days
      theta_ir_ph <- water_ir[ph$idx]
      w_ir_safe <- base::ifelse(base::is.na(theta_ir_ph), params$f_W_1_lower_t, theta_ir_ph)
      f_W_ir_ph <- base::pmin(base::pmax(
        base::ifelse(w_ir_safe >= params$f_W_1_lower_t, 1.0,
                     params$f_W_floor_t + (1 - params$f_W_floor_t) *
                       (w_ir_safe / params$f_W_1_lower_t)),
        params$f_W_floor_t), 1.0)
      heat_biomass_lost_ir <- heat_biomass_lost_ir +
        base::sum(light_int[ph_heat] * f_W_ir_ph[ph_heat])
      
      # f(W) rainfed on heat-stress days
      theta_rf_ph <- water_rf[ph$idx]
      w_rf_safe <- base::ifelse(base::is.na(theta_rf_ph), params$f_W_1_lower_t, theta_rf_ph)
      f_W_rf_ph <- base::pmin(base::pmax(
        base::ifelse(w_rf_safe >= params$f_W_1_lower_t, 1.0,
                     params$f_W_floor_t + (1 - params$f_W_floor_t) *
                       (w_rf_safe / params$f_W_1_lower_t)),
        params$f_W_floor_t), 1.0)
      heat_biomass_lost_rf <- heat_biomass_lost_rf +
        base::sum(light_int[ph_heat] * f_W_rf_ph[ph_heat])
    }
    
    # ── DROUGHT STRESS ────────────────────────────────────────────────────────
    # Drought stress day: theta_rel < RAW threshold.
    # RAW = p * TAW with p = 0.25 for onion (FAO-56 Table 22).
    # theta_rel = 1 - Dr/TAW, so the stress threshold is theta_rel < 1 - 0.25 = 0.75.
    # Below this point Ks < 1 in the FAO-56 framework, meaning actual ET
    # is reduced and f(W) begins to suppress biomass accumulation.
    # Ref: Allen et al. (1998) FAO-56.
    
    RAW_threshold <- 0.75
    
    # Fraction of growing-season days below RAW threshold
    gs_theta_rf <- df$theta_rel_rf[gs_idx]
    drought_stress_day_fraction_rf <- base::sum(
      !base::is.na(gs_theta_rf) & gs_theta_rf < RAW_threshold
    ) / n_gs
    
    gs_theta_ir <- df$theta_rel_ir[gs_idx]
    drought_stress_day_fraction_ir <- base::sum(
      !base::is.na(gs_theta_ir) & gs_theta_ir < RAW_threshold
    ) / n_gs
    
    # Biomass lost on drought-stress days:
    # The daily biomass cost of water deficit = (f_W=1 - actual f_W) * base_flux,
    # where base_flux = LUE * PAR * light_interception * f_T (actual temperature).
    # f_T is held at its actual value so only the water deficit is attributed.
    
    drought_biomass_lost_rf <- 0
    drought_biomass_lost_ir <- 0
    
    for (ph in base::list(
      base::list(idx = phase_data$em$idx, LAI = params$LAI_emergence, LUE = params$LUE_em_p),
      base::list(idx = phase_data$vg$idx, LAI = params$LAI_veg,       LUE = params$LUE_vg_p),
      base::list(idx = phase_data$bl$idx, LAI = params$LAI_bulbing,   LUE = params$LUE_bl_p),
      base::list(idx = phase_data$mt$idx, LAI = params$LAI_maturation,LUE = params$LUE_mt_p)
    )) {
      if (!base::any(ph$idx, na.rm = TRUE)) next
      
      ph_PAR    <- df$PAR[ph$idx]
      ph_Tmax_d <- df$Tmax[ph$idx]
      ph_Tmin_d <- df$Tmin[ph$idx]
      ph_yday_d <- df$yday[ph$idx]
      
      # f_T: PAR-weighted daylight integration using same matrix method as
      # calc_bio_vectorized -- no day-level loops.
      {
        n_ph    <- base::length(ph_PAR)
        hh_d    <- (0:23) + 0.5
        ag_d    <- (params$Tmax_growth_p - params$Topt_growth_p) /
          (params$Topt_growth_p - params$Tmin_growth_p)
        valid_d <- !base::is.na(ph_Tmax_d) & !base::is.na(ph_Tmin_d)
        
        decl_d  <- 0.409 * base::sin(2 * base::pi * ph_yday_d / 365 - 1.39)
        ws_d    <- base::acos(base::pmax(-1, base::pmin(1,
                                                        -base::tan(51 * base::pi / 180) * base::tan(decl_d))))
        DLd     <- 24 * ws_d / base::pi
        trd     <- 12 - DLd / 2
        tsd     <- 12 + DLd / 2
        
        hm_d  <- base::matrix(hh_d,      nrow = n_ph, ncol = 24L, byrow = TRUE)
        trm_d <- base::matrix(trd,        nrow = n_ph, ncol = 24L)
        tsm_d <- base::matrix(tsd,        nrow = n_ph, ncol = 24L)
        txm_d <- base::matrix(ph_Tmax_d,  nrow = n_ph, ncol = 24L)
        tnm_d <- base::matrix(ph_Tmin_d,  nrow = n_ph, ncol = 24L)
        dlm_d <- base::matrix(DLd,        nrow = n_ph, ncol = 24L)
        
        rise_d <- hm_d >= trm_d & hm_d <= 14
        fall_d <- hm_d >  14    & hm_d <  tsm_d
        Tm_d   <- tnm_d
        Tm_d[rise_d] <- (tnm_d + (txm_d - tnm_d) *
                           base::sin(base::pi/2*(hm_d-trm_d)/(14-trm_d)))[rise_d]
        Tm_d[fall_d] <- (tnm_d + (txm_d - tnm_d) *
                           base::sin(base::pi/2*(tsm_d-hm_d)/(tsm_d-14)))[fall_d]
        Tm_d[!valid_d, ] <- NA_real_
        
        dld   <- hm_d >= trm_d & hm_d < tsm_d
        PR_d  <- base::matrix(0, nrow = n_ph, ncol = 24L)
        PR_d[dld] <- base::pmax(
          base::sin(base::pi*(hm_d-trm_d)/dlm_d)[dld], 0)
        rs_d  <- base::rowSums(PR_d)
        ss_d  <- base::ifelse(rs_d > 0 & valid_d, rs_d, 1)
        PM_d  <- PR_d * (ph_PAR / ss_d)   # scale to daily PAR total (matches calc_bio_vectorized)
        
        insd  <- Tm_d > params$Tmin_growth_p & Tm_d < params$Tmax_growth_p
        fTm_d <- base::matrix(0, nrow = n_ph, ncol = 24L)
        fTm_d[insd] <- base::pmax(
          ((Tm_d  - params$Tmin_growth_p) / (params$Topt_growth_p - params$Tmin_growth_p)) *
            ((params$Tmax_growth_p - Tm_d)  / (params$Tmax_growth_p - params$Topt_growth_p))^ag_d,
          0)[insd]
        
        num_d <- base::rowSums(fTm_d * PM_d)
        den_d <- base::rowSums(PM_d)
        f_T_ph <- base::ifelse(valid_d & den_d > 0, num_d / den_d, 0)
      }
      
      # Base flux: light × temperature (water = 1 counterfactual)
      base_flux <- ph$LUE * ph_PAR *
        (1 - base::exp(-params$lec_k_c * ph$LAI)) * f_T_ph
      
      # Rainfed: biomass lost = (1 - f_W_rf) * base_flux on drought days
      theta_rf_ph <- water_rf[ph$idx]
      w_rf_safe <- base::ifelse(base::is.na(theta_rf_ph), params$f_W_1_lower_t, theta_rf_ph)
      f_W_rf <- base::pmin(base::pmax(
        base::ifelse(w_rf_safe >= params$f_W_1_lower_t, 1.0,
                     params$f_W_floor_t + (1 - params$f_W_floor_t) *
                       (w_rf_safe / params$f_W_1_lower_t)),
        params$f_W_floor_t), 1.0)
      drought_rf <- !base::is.na(theta_rf_ph) & theta_rf_ph < RAW_threshold
      drought_biomass_lost_rf <- drought_biomass_lost_rf +
        base::sum((1.0 - f_W_rf[drought_rf]) * base_flux[drought_rf])
      
      # Irrigated: same logic — typically near-zero by design
      theta_ir_ph <- water_ir[ph$idx]
      w_ir_safe <- base::ifelse(base::is.na(theta_ir_ph), params$f_W_1_lower_t, theta_ir_ph)
      f_W_ir <- base::pmin(base::pmax(
        base::ifelse(w_ir_safe >= params$f_W_1_lower_t, 1.0,
                     params$f_W_floor_t + (1 - params$f_W_floor_t) *
                       (w_ir_safe / params$f_W_1_lower_t)),
        params$f_W_floor_t), 1.0)
      drought_ir <- !base::is.na(theta_ir_ph) & theta_ir_ph < RAW_threshold
      drought_biomass_lost_ir <- drought_biomass_lost_ir +
        base::sum((1.0 - f_W_ir[drought_ir]) * base_flux[drought_ir])
    }
    
    # Fresh-weight yield gap (t/ha): absolute benefit of irrigation each season
    yield_gap_rainfed_irrigated <-
      (attainable_yield_DM_t_ha_ir - attainable_yield_DM_t_ha_rf) /
      dry_matter_content_t
    
    # Harvest day
    after_sow_idx <- df$yday >= sow
    harvest_yday <- if (base::any(phase_data$mt$idx)) {
      base::max(df$yday[phase_data$mt$idx])
    } else {
      base::max(df$yday[after_sow_idx], na.rm = TRUE)
    }
    
    # Store results
    results[[sc]] <- base::list(
      sowing_yday  = sow,
      harvest_yday = harvest_yday,
      
      raw_yield_per_ha_rainfed   = attainable_yield_DM_t_ha_rf/dry_matter_content_t,
      final_yield_per_ha_rainfed = final_yield_DM_t_ha_rf/dry_matter_content_t,
      
      
      raw_yield_per_ha_irrigated   = attainable_yield_DM_t_ha_ir/dry_matter_content_t,
      final_yield_per_ha_irrigated = final_yield_DM_t_ha_ir/dry_matter_content_t,
      
      # Seasonal stress diagnostics (unified pool)
      stress_yield_multiplier = stress_yield_multiplier,
      total_stress_loss       = total_stress_loss,
      
      # Phase-specific risk scores
      risk_mildew_vg  = all_risks_rf$vg["mildew"],
      risk_mildew_bl  = all_risks_rf$bl["mildew"],
      risk_mildew_mt  = all_risks_rf$mt["mildew"],
      
      risk_botrytis_bl = all_risks_rf$bl["botrytis"],
      risk_botrytis_mt = all_risks_rf$mt["botrytis"],
      
      risk_fusarium_vg = all_risks_rf$vg["fusarium"],
      risk_fusarium_bl = all_risks_rf$bl["fusarium"],
      risk_fusarium_mt = all_risks_rf$mt["fusarium"],
      
      risk_exrain_em  = all_risks_rf$em["exrain"],
      risk_exrain_vg  = all_risks_rf$vg["exrain"],
      risk_exrain_bl  = all_risks_rf$bl["exrain"],
      risk_exrain_mt  = all_risks_rf$mt["exrain"],
      
      # Individual seasonal risk scores (0-1, from daily scoring)
      risk_mildew    = risk_mildew,
      risk_botrytis  = risk_botrytis,
      risk_fusarium  = risk_fusarium,
      risk_exrain    = risk_exrain,
      
      # Individual seasonal losses (0-1, after chance_event gate)
      loss_mildew    = loss_mildew,
      loss_botrytis  = loss_botrytis,
      loss_fusarium  = loss_fusarium,
      loss_exrain    = loss_exrain,
      loss_hail      = loss_hail,
      
      total_irrigation_mm = base::sum(df$Irrig_mm, na.rm = TRUE),
      
      # Per-phase irrigation totals (mm) — for validation against practice
      irrig_mm_em = base::sum(df$Irrig_mm[phase_data$em$idx], na.rm = TRUE),
      irrig_mm_vg = base::sum(df$Irrig_mm[phase_data$vg$idx], na.rm = TRUE),
      irrig_mm_bl = base::sum(df$Irrig_mm[phase_data$bl$idx], na.rm = TRUE),
      irrig_mm_mt = base::sum(df$Irrig_mm[phase_data$mt$idx], na.rm = TRUE),
      
      # Number of irrigation events per phase
      irrig_events_total = base::sum(df$Irrig_mm > 0, na.rm = TRUE),
      
      # Mean theta_rel per phase — rainfed (diagnostic)
      theta_rel_rf_em = base::mean(df$theta_rel_rf[phase_data$em$idx], na.rm = TRUE),
      theta_rel_rf_vg = base::mean(df$theta_rel_rf[phase_data$vg$idx], na.rm = TRUE),
      theta_rel_rf_bl = base::mean(df$theta_rel_rf[phase_data$bl$idx], na.rm = TRUE),
      theta_rel_rf_mt = base::mean(df$theta_rel_rf[phase_data$mt$idx], na.rm = TRUE),
      
      # Mean theta_rel per phase — irrigated (diagnostic)
      theta_rel_ir_em = base::mean(df$theta_rel_ir[phase_data$em$idx], na.rm = TRUE),
      theta_rel_ir_vg = base::mean(df$theta_rel_ir[phase_data$vg$idx], na.rm = TRUE),
      theta_rel_ir_bl = base::mean(df$theta_rel_ir[phase_data$bl$idx], na.rm = TRUE),
      theta_rel_ir_mt = base::mean(df$theta_rel_ir[phase_data$mt$idx], na.rm = TRUE),
      
      # ── DIAGNOSTIC: Heat stress (Tavg >= Tmax_growth) ─────────────────────
      # heat_stress_day_fraction: fraction of growing-season days with complete
      #   growth halt. Identical for ir/rf since temperature is independent of
      #   irrigation. Access: onion_mc_simulation$y$historical.heat_stress_day_fraction
      heat_stress_day_fraction    = heat_stress_day_fraction,
      
      # heat_biomass_lost_ir/rf: biomass foregone on heat-stress days (g m-2)
      #   under irrigated and rainfed water conditions respectively.
      #   Access: onion_mc_simulation$y$ssp585.heat_biomass_lost_ir
      heat_biomass_lost_ir        = heat_biomass_lost_ir,
      heat_biomass_lost_rf        = heat_biomass_lost_rf,
      
      # ── DIAGNOSTIC: Drought stress (theta_rel < RAW = 0.75) ───────────────
      # drought_stress_day_fraction_rf: fraction of days below FAO-56 RAW
      #   threshold in the rainfed scenario.
      # drought_stress_day_fraction_ir: same for irrigated (sanity check, ~0).
      #   Access: onion_mc_simulation$y$historical.drought_stress_day_fraction_rf
      drought_stress_day_fraction_rf = drought_stress_day_fraction_rf,
      drought_stress_day_fraction_ir = drought_stress_day_fraction_ir,
      
      # drought_biomass_lost_rf/ir: biomass suppressed by f(W) < 1 on drought
      #   days (g m-2), holding f(T) constant so only water deficit is attributed.
      #   Access: onion_mc_simulation$y$ssp585.drought_biomass_lost_rf
      drought_biomass_lost_rf     = drought_biomass_lost_rf,
      drought_biomass_lost_ir     = drought_biomass_lost_ir,
      
      # ── Rainfed–irrigated yield gap ────────────────────────────────────────
      #   Access: onion_mc_simulation$y$historical.yield_gap_rainfed_irrigated
      yield_gap_rainfed_irrigated = yield_gap_rainfed_irrigated
    )
  }
  
  base::return(results)
})


# ── Monte Carlo simulation ────────────────────────────────────────────────────

input_variables <- read.csv("input_irrig1.csv", header = TRUE, sep = ";")

onion_mc_simulation <- mcSimulation(
  estimate          = as.estimate(input_variables),
  model_function    = onion_climate_impact,
  numberOfModelRuns = 10,       # Set to your preference
  functionSyntax    = "plainNames"
)

# Quick check: mean historical yields
mean(onion_mc_simulation$y$historical.final_yield_per_ha_rainfed)
mean(onion_mc_simulation$y$historical.final_yield_per_ha_irrigated)

# Save / load results
# saveRDS(onion_mc_simulation, "results.RDS")
#onion_mc_simulation <- readRDS("results.RDS")

