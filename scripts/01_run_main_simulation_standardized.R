# scripts/11_run_rescue_simulations.R
# ------------------------------------------------------------
# Paper B runner: pressure-indexed robustness benchmark
# - Reuses Paper A DGP + SPD(t) code
# - Legacy PP-IPCW Cox log-HR retained as a diagnostic consistency check
# - Primary risk comparison: baseline-standardized IPCW versus DR-AIPW
#   for the same marginal fixed-horizon per-protocol risk estimand
# - Optional: ML outcome regression for DR (secondary)
# - Optional: simple parametric TMLE for risk (secondary)
#
# Outputs (OUT_DIR):
#   raw/replicate_results_rescue.csv
#   raw/spd_curves_rescue.csv
#   raw/weight_diagnostics_rescue.csv
#   raw/ipcw_risk_rescue.csv
#   raw/dr_risk_rescue.csv
#   raw/tmle_risk_rescue.csv (optional)
#   perf_summary_rescue_hr.csv
#   perf_summary_rescue_ipcw_risk.csv
#   perf_summary_rescue_dr.csv
#   perf_summary_rescue_tmle.csv (optional)
# ------------------------------------------------------------

source("scripts/00_utils.R")

cfg <- init_project(
  n_cores = as.integer(Sys.getenv("N_CORES", "3")),
  seed    = 2026L,
  out_dir = Sys.getenv("OUT_DIR", "output_jci_final")
)
cfg$require_pkgs(c("data.table", "survival", "future.apply", "future"))

# --- knobs ---
B     <- as.integer(Sys.getenv("B", "200"))
N     <- as.integer(Sys.getenv("N", "2000"))
t_max <- as.numeric(Sys.getenv("T_MAX", "5"))
dt    <- as.numeric(Sys.getenv("DT", "0.25"))

# Horizon for risk estimators
T_CUT <- as.numeric(Sys.getenv("T_CUT", as.character(t_max)))

# Nuisance-regime levels
MIS_SPEC_LEVELS <- strsplit(Sys.getenv("MIS_SPEC_LEVELS", "0,1"), ",", fixed = TRUE)[[1]]
MIS_SPEC_LEVELS <- as.integer(trimws(MIS_SPEC_LEVELS))
MIS_SPEC_LEVELS <- MIS_SPEC_LEVELS[MIS_SPEC_LEVELS %in% c(0L, 1L)]
if (length(MIS_SPEC_LEVELS) == 0) MIS_SPEC_LEVELS <- c(0L, 1L)

# Truth (default: causal null only; acceptance-risk minimized)
parse_truth_levels <- function() {
  raw <- trimws(Sys.getenv("TRUTH_LEVELS", unset = "null"))
  if (!nzchar(raw) || tolower(raw) %in% c("null", "none")) return(c("null"))

  # Comma-separated list, e.g., "null,non_null"
  out <- trimws(strsplit(raw, ",", fixed = TRUE)[[1]])
  out <- out[!is.na(out)]
  out <- out[out %in% c("null", "non_null")]
  if (length(out) == 0) out <- c("null")
  out
}
truth_levels <- parse_truth_levels()

# Parse numeric environment settings robustly. This allows either a plain
# numeric string (recommended for reproducibility, e.g. -0.22314355131420976)
# or a simple R expression such as log(0.8). The previous version used
# as.numeric("log(0.8)"), which produced NA and then caused simulated hazards
# and event times to become missing.
parse_numeric_env <- function(name, default = "0") {
  raw <- trimws(Sys.getenv(name, unset = default))
  out <- suppressWarnings(as.numeric(raw))
  if (length(out) == 1L && is.finite(out)) return(out)

  # Restricted fallback for simple numeric expressions. This is mainly for
  # backwards compatibility with BETA_TRUE_PP=log(0.8).
  if (!grepl("^[0-9eE+./*()[:space:]a-zA-Z_-]+$", raw)) {
    stop(name, " must be numeric or a simple numeric expression; got: ", raw)
  }
  val <- tryCatch(eval(parse(text = raw), envir = baseenv()), error = function(e) NA_real_)
  if (!is.numeric(val) || length(val) != 1L || !is.finite(val)) {
    stop(name, " could not be parsed as a finite numeric value; got: ", raw)
  }
  as.numeric(val)
}

beta_true_pp <- parse_numeric_env("BETA_TRUE_PP", "0")

if ("non_null" %in% truth_levels && isTRUE(beta_true_pp == 0)) {
  warning("TRUTH_LEVELS includes non_null but BETA_TRUE_PP is 0. Set e.g. BETA_TRUE_PP=-0.22314355131420976 for a non-null RR of 0.8.")
}

# Monte Carlo "truth" for risk-scale estimand (logRR/RD at T_CUT)
# Needed to compute bias/RMSE/coverage under non-null on the risk scale.
MC_TRUTH <- as.integer(Sys.getenv("MC_TRUTH", "1")) == 1L
N_TRUTH  <- as.integer(Sys.getenv("N_TRUTH", "200000"))
B_TRUTH  <- as.integer(Sys.getenv("B_TRUTH", "5"))


# tipping thresholds
c_ess  <- as.numeric(Sys.getenv("C_ESS", "0.25"))
c_tail <- as.numeric(Sys.getenv("C_TAIL", "0.10"))

# SPD(t) piecewise breaks
breaks <- seq(0, t_max, by = 1.0)
if (tail(breaks, 1) < t_max) breaks <- c(breaks, t_max)

# Comparator settings
RUN_DR   <- as.integer(Sys.getenv("RUN_DR", "1")) == 1L
RUN_ML   <- as.integer(Sys.getenv("RUN_ML", "0")) == 1L
ML_METHOD <- Sys.getenv("ML_METHOD", "glmnet")
CROSSFIT <- as.integer(Sys.getenv("CROSSFIT", "0")) == 1L
CF_FOLDS <- as.integer(Sys.getenv("CF_FOLDS", "2"))
RUN_TMLE <- as.integer(Sys.getenv("RUN_TMLE", "0")) == 1L
IPCW_ESTIMATOR <- Sys.getenv("IPCW_ESTIMATOR", "hajek")
if (!IPCW_ESTIMATOR %in% c("hajek", "ht")) IPCW_ESTIMATOR <- "hajek"

# Supplementary reviewer-guardrail sensitivity estimator:
# IPCW with a person-period deviation model that includes observed
# time-updated prognosis z_obs(t). This is not the main estimator, but it
# addresses the concern that the primary baseline/horizon-level nuisance
# models may be too simple for a time-varying deviation DGP.
RUN_TIME_UPDATED_SENS <- as.integer(Sys.getenv("RUN_TIME_UPDATED_SENS", "1")) == 1L


# file tag
TAG <- sprintf("B%d_N%d", B, N)

# scenario grid (reuse Paper A)
grid <- make_scenario_grid()

# expanded grid with mis-spec
scn_grid <- data.table::CJ(
  scenario_id = grid$scenario_id,
  mis_spec = MIS_SPEC_LEVELS,
  unique = TRUE
)
scn_grid <- merge(scn_grid, grid, by = "scenario_id", all.x = TRUE)
data.table::setorder(scn_grid, scenario_num, mis_spec)

# optional scenario subset
subset_ids <- strsplit(Sys.getenv("SCENARIO_SUBSET", ""), ",", fixed = TRUE)[[1]]
subset_ids <- trimws(subset_ids)
subset_ids <- subset_ids[nzchar(subset_ids)]
if (length(subset_ids) > 0) {
  scn_grid <- scn_grid[scenario_id %in% subset_ids]
  data.table::setorder(scn_grid, scenario_num, mis_spec)
}


# ------------------------------------------------------------
# Intervention-based Monte Carlo truth for risk-scale estimand
# ------------------------------------------------------------
# Acceptance-risk minimisation for JCI:
#   The true fixed-horizon risks are computed by intervention on treatment
#   assignment under the known DGP, with artificial deviation disabled. This
#   avoids defining truth by one of the estimators being evaluated.
#
#   Null contrasts are set exactly to zero. Non-null contrasts are computed
#   by common-random-number Monte Carlo for A=1 versus A=0.
# ------------------------------------------------------------
mc_truth_risk <- NULL
if (MC_TRUTH && ("non_null" %in% truth_levels || "null" %in% truth_levels)) {
  truth_path <- file.path(cfg$out_dir, "raw", "mc_truth_risk.csv")
  dir.create(dirname(truth_path), showWarnings = FALSE, recursive = TRUE)

  # Reuse cached truth only if compatible and intervention-based.
  if (file.exists(truth_path)) {
    tmp <- tryCatch(data.table::fread(truth_path, showProgress = FALSE), error = function(e) NULL)
    needed_cols <- c("scenario_id", "truth", "t_cut", "beta_true_pp", "N_truth", "B_truth",
                     "truth_method", "logrr_true", "rd_true")
    if (!is.null(tmp) && all(needed_cols %in% names(tmp))) {
      ok <- isTRUE(all(tmp$t_cut == T_CUT)) &&
        isTRUE(all(tmp$beta_true_pp == beta_true_pp)) &&
        isTRUE(all(tmp$N_truth == N_TRUTH)) &&
        isTRUE(all(tmp$B_truth == B_TRUTH)) &&
        isTRUE(all(grepl("^intervention", tmp$truth_method)))
      if (ok) {
        mc_truth_risk <- tmp
        message("Loaded cached intervention MC truth: ", truth_path)
      }
    }
  }

  if (is.null(mc_truth_risk)) {
    message("Computing intervention-based MC truth for risk logRR/RD (N_TRUTH=", N_TRUTH,
            ", B_TRUTH=", B_TRUTH, ") ...")

    base_grid <- make_scenario_grid()[, .(scenario_id, axisA, axisB, rho_meas)]
    truth_jobs <- data.table::CJ(scenario_id = base_grid$scenario_id, truth = truth_levels, unique = TRUE)

    mc_list <- lapply(seq_len(nrow(truth_jobs)), function(i) {
      sid <- truth_jobs$scenario_id[i]
      tr  <- truth_jobs$truth[i]
      seed_i <- seed_for_job(cfg$seed * 900000L, sid, truth = tr, rep_id = 999001L)
      compute_intervention_truth_risk(
        scenario_id = sid,
        truth = tr,
        N_truth = N_TRUTH,
        B_truth = B_TRUTH,
        t_cut = T_CUT,
        t_max = t_max,
        dt = dt,
        beta_true_pp = beta_true_pp,
        seed_base = seed_i
      )
    })

    mc_truth_risk <- data.table::rbindlist(mc_list, fill = TRUE)
    data.table::fwrite(mc_truth_risk, truth_path)
    message("Saved intervention MC truth: ", truth_path)
  }
}


# logging
log_file <- file.path(cfg$log_dir, sprintf("rescue_progress_%s.log", TAG))
if (file.exists(log_file)) file.remove(log_file)

cfg$log_line(log_file, sprintf(
  "[%s] Rescue simulation started | OUT_DIR=%s | scenarios=%d | mis_spec=%s | B=%d | N=%d | truths=%s | DR=%d | ML=%d (%s) | TMLE=%d | IPCW_ESTIMATOR=%s | TIME_UPDATED_SENS=%d\n",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  cfg$out_dir, nrow(scn_grid), paste(MIS_SPEC_LEVELS, collapse = ","),
  B, N, paste(truth_levels, collapse = ","), as.integer(RUN_DR), as.integer(RUN_ML), ML_METHOD, as.integer(RUN_TMLE), IPCW_ESTIMATOR,
  as.integer(RUN_TIME_UPDATED_SENS)
))

cat("Scenarios:", nrow(scn_grid), " | truths:", paste(truth_levels, collapse = ","), " | B=", B, " | N=", N, "\n")
cat("OUT_DIR:", cfg$out_dir, "\n")

future::plan(future::multisession, workers = cfg$n_cores)

run_one <- function(scn_row, truth, rep_id) {
  axisA <- scn_row$axisA
  axisB <- scn_row$axisB
  rho   <- scn_row$rho_meas
  mis_spec <- scn_row$mis_spec

  truth <- as.character(truth)
  if (length(truth) != 1L || is.na(truth) || !nzchar(truth)) {
    stop(sprintf("Invalid truth level passed to run_one: length=%d value=%s", length(truth), paste(truth, collapse=",")))
  }

  beta <- if (isTRUE(truth == "non_null")) beta_true_pp else 0

  seed <- seed_for_job(cfg$seed * 100000L, scn_row$scenario_id, truth = truth, rep_id = rep_id)

  # simulate dataset
  sim <- simulate_one_dataset(
    N = N, t_max = t_max, dt = dt,
    beta_true = beta,
    axisA = axisA,
    axisB = axisB,
    rho_meas = rho,
    seed = seed
  )

  long_dt <- sim$long_dt
  pp_dt   <- sim$pp_dt

  # SPD(t)
  spd_fit <- estimate_spd_piecewise(
    long_dt, breaks = breaks,
    z_col = "z_obs", dev_col = "dev", id_col = "id", robust = TRUE
  )
  spd_dt <- compute_cum_pressure(spd_fit$spd)

  # IPCW Cox log-HR
  cens_fml <- if (mis_spec == 1L) {
    survival::Surv(time_pp, dev_ind) ~ A0 + age + sex + bmi
  } else {
    survival::Surv(time_pp, dev_ind) ~ A0 + age + sex + bmi + egfr + util + gall
  }
  pp_est <- fit_pp_ipcw(pp_dt, cens_formula = cens_fml)

  # Baseline-standardized IPCW risk at horizon.
  # This is the JCI-primary IPCW comparator. It uses the same baseline
  # covariate set and horizon-level adherence/censoring target as DR-AIPW,
  # so the IPCW and DR-AIPW risk estimators target the same marginal
  # fixed-horizon per-protocol risk contrast.
  ipcw_risk <- fit_ipcw_risk_standardized(
    pp_dt = pp_dt,
    t_cut = T_CUT,
    mis_spec = mis_spec,
    w_floor = 0.05,
    ps_floor = 1e-3,
    estimator = IPCW_ESTIMATOR
  )

  # Time-updated IPCW sensitivity estimator. This is supplementary and should
  # not replace the main standardized-IPCW versus DR-AIPW comparison.
  ipcw_tu <- NULL
  if (RUN_TIME_UPDATED_SENS) {
    ipcw_tu <- fit_ipcw_risk_timeupdated(
      long_dt = long_dt,
      pp_dt = pp_dt,
      t_cut = T_CUT,
      mis_spec = mis_spec,
      w_floor = 0.05,
      ps_floor = 1e-3,
      estimator = IPCW_ESTIMATOR,
      ns_df = 3L
    )
  }

  # weights diagnostics on the time grid use the Cox deviation model from
  # the legacy PP-IPCW Cox analysis. These diagnostics remain useful for
  # time-resolved weight stability, but they are not the primary risk
  # estimator in the JCI package.
  fit_cens <- pp_est$fit_cens
  lp <- stats::predict(fit_cens, type = "lp")
  t_eval <- sort(unique(long_dt$tstart))
  weights_dt <- data.table::CJ(id = pp_dt$id, t = t_eval)
  lp_dt <- data.table::data.table(id = pp_dt$id, lp = lp)
  weights_dt <- merge(weights_dt, lp_dt, by = "id", all.x = TRUE)
  weights_dt[, S := get_S_from_cox(fit_cens, time_vec = t, lp_vec = lp)]
  weights_dt[, w := 1 / pmax(S, 1e-3)]

  diag_dt <- compute_weight_diagnostics(long_dt, weights_dt,
                                        time_col = "tstart", id_col = "id", w_col = "w", q_tail = 0.99)
  tip <- detect_tipping(diag_dt, c_ess = c_ess, c_tail = c_tail)

  # DR risk (logRR/RD) at horizon
  dr <- NULL
  if (RUN_DR) {
    dr <- fit_dr_risk_rescue(
      pp_dt = pp_dt,
      t_cut = T_CUT,
      mis_spec = mis_spec,
      use_ml_Q = RUN_ML,
      ml_method = ML_METHOD,
      crossfit = CROSSFIT,
      cf_folds = CF_FOLDS,
      seed = seed + 17L
    )
  }

  # TMLE risk (supplementary)
  tm <- NULL
  if (RUN_TMLE) {
    tm <- fit_tmle_or_dr(pp_dt = pp_dt, t_cut = T_CUT, mis_spec = mis_spec)
  }

  list(
    rep = data.table::data.table(
      scenario_id = scn_row$scenario_id,
      axisA = axisA, axisB = axisB, rho_meas = rho,
      mis_spec = mis_spec,
      truth = truth,
      replicate = rep_id,
      beta_hat = pp_est$beta,
      se_hat = pp_est$se,
      ess_final = pp_est$ess,
      tip_time = tip$t_star,
      max_spd = max(abs(spd_dt$gamma_hat), na.rm = TRUE),
      Gamma_end = max(spd_dt$Gamma_hat, na.rm = TRUE)
    ),
    spd = spd_dt[, .(
      scenario_id = scn_row$scenario_id,
      mis_spec = mis_spec,
      truth = truth,
      replicate = rep_id,
      interval, tstart, tstop, t_mid,
      gamma_hat, se, HR_1SD, CI_low, CI_high,
      Gamma_hat
    )],
    diag = tip$diag[, .(
      scenario_id = scn_row$scenario_id,
      mis_spec = mis_spec,
      truth = truth,
      replicate = rep_id,
      t, N, ESS, rESS, tail_share, tipped
    )],
    dr = if (!is.null(dr)) data.table::data.table(
      scenario_id = scn_row$scenario_id,
      axisA = axisA, axisB = axisB, rho_meas = rho,
      mis_spec = mis_spec,
      truth = truth,
      replicate = rep_id,
      t_cut = dr$t,
      method_Q = dr$method_Q,
      logrr = dr$logrr,
      se_logrr = dr$se_logrr,
      rd = dr$rd,
      se_rd = dr$se_rd,
      risk1 = dr$risk1,
      risk0 = dr$risk0
    ) else data.table::data.table(),
    ipcw_risk = data.table::data.table(
      scenario_id = scn_row$scenario_id,
      axisA = axisA, axisB = axisB, rho_meas = rho,
      mis_spec = mis_spec,
      truth = truth,
      replicate = rep_id,
      t_cut = ipcw_risk$t,
      method = ipcw_risk$method,
      logrr = ipcw_risk$logrr,
      se_logrr = ipcw_risk$se_logrr,
      rd = ipcw_risk$rd,
      se_rd = ipcw_risk$se_rd,
      risk1 = ipcw_risk$risk1,
      risk0 = ipcw_risk$risk0,
      den1 = ipcw_risk$den1,
      den0 = ipcw_risk$den0,
      min_g1 = ipcw_risk$min_g1,
      max_g1 = ipcw_risk$max_g1,
      min_c_hat = ipcw_risk$min_c_hat
    ),
    ipcw_tu_risk = if (!is.null(ipcw_tu)) data.table::data.table(
      scenario_id = scn_row$scenario_id,
      axisA = axisA, axisB = axisB, rho_meas = rho,
      mis_spec = mis_spec,
      truth = truth,
      replicate = rep_id,
      t_cut = ipcw_tu$t,
      method = ipcw_tu$method,
      dev_model_status = ipcw_tu$dev_model_status,
      logrr = ipcw_tu$logrr,
      se_logrr = ipcw_tu$se_logrr,
      rd = ipcw_tu$rd,
      se_rd = ipcw_tu$se_rd,
      risk1 = ipcw_tu$risk1,
      risk0 = ipcw_tu$risk0,
      den1 = ipcw_tu$den1,
      den0 = ipcw_tu$den0,
      min_g1 = ipcw_tu$min_g1,
      max_g1 = ipcw_tu$max_g1,
      min_c_hat = ipcw_tu$min_c_hat,
      median_c_hat = ipcw_tu$median_c_hat,
      median_intervals = ipcw_tu$median_intervals,
      mean_interval_dev_prob = ipcw_tu$mean_interval_dev_prob
    ) else data.table::data.table(),
    tm = if (!is.null(tm)) data.table::data.table(
      scenario_id = scn_row$scenario_id,
      axisA = axisA, axisB = axisB, rho_meas = rho,
      mis_spec = mis_spec,
      truth = truth,
      replicate = rep_id,
      t_cut = tm$t,
      method_Q = tm$method_Q,
      logrr = tm$logrr,
      se_logrr = tm$se_logrr,
      rd = tm$rd,
      se_rd = tm$se_rd,
      risk1 = tm$risk1,
      risk0 = tm$risk0
    ) else data.table::data.table()
  )
}

# scenario loop
rep_chunks  <- vector("list", nrow(scn_grid))
spd_chunks  <- vector("list", nrow(scn_grid))
diag_chunks <- vector("list", nrow(scn_grid))
dr_chunks   <- vector("list", nrow(scn_grid))
tm_chunks   <- vector("list", nrow(scn_grid))
ipcw_risk_chunks <- vector("list", nrow(scn_grid))
ipcw_tu_risk_chunks <- vector("list", nrow(scn_grid))

for (i in seq_len(nrow(scn_grid))) {
  scn <- scn_grid[i]

  cat(sprintf("\n[Scenario %d/%d] %s | A=%s | B=%s | rho=%.1f | mis_spec=%d\n",
              i, nrow(scn_grid), scn$scenario_id, scn$axisA, scn$axisB, scn$rho_meas, scn$mis_spec))
  flush.console()

  cfg$log_line(log_file, sprintf("[%s] SCENARIO %d/%d start | %s | mis_spec=%d\n",
                                 format(Sys.time(), "%H:%M:%S"), i, nrow(scn_grid), scn$scenario_id, scn$mis_spec))

  jobs <- data.table::CJ(truth = truth_levels, replicate = 1:B)
  jobs <- jobs[order(truth, replicate)]

  res_list <- future.apply::future_lapply(
    X = seq_len(nrow(jobs)),
    FUN = function(j) {
      truth_j <- as.character(jobs$truth[j])[1]
      if (length(truth_j) != 1L || is.na(truth_j) || !nzchar(truth_j)) truth_j <- "null"
      rep_j   <- as.integer(jobs$replicate[j])[1]
      t0 <- proc.time()[3]
      out <- tryCatch(
        run_one(scn, truth_j, rep_j),
        error = function(e) {
          msg <- conditionMessage(e)

          # Ensure we always return a consistent schema so downstream summaries can be written
          rep_dt <- data.table::data.table(
            scenario_id = scn$scenario_id,
            axisA = scn$axisA, axisB = scn$axisB, rho_meas = scn$rho_meas,
            mis_spec = scn$mis_spec,
            truth = truth_j,
            replicate = rep_j,
            beta_hat = NA_real_, se_hat = NA_real_, ess_final = NA_real_,
            tip_time = NA_real_, max_spd = NA_real_, Gamma_end = NA_real_,
            error = msg
          )

          ipcw_risk_dt <- data.table::data.table(
            scenario_id = scn$scenario_id,
            axisA = scn$axisA, axisB = scn$axisB, rho_meas = scn$rho_meas,
            mis_spec = scn$mis_spec,
            truth = truth_j,
            replicate = rep_j,
            t_cut = T_CUT,
            method = paste0("ipcw_std_", IPCW_ESTIMATOR),
            logrr = NA_real_, se_logrr = NA_real_,
            rd = NA_real_, se_rd = NA_real_,
            risk1 = NA_real_, risk0 = NA_real_,
            error = msg
          )

          dr_dt <- if (RUN_DR) data.table::data.table(
            scenario_id = scn$scenario_id,
            axisA = scn$axisA, axisB = scn$axisB, rho_meas = scn$rho_meas,
            mis_spec = scn$mis_spec,
            truth = truth_j,
            replicate = rep_j,
            t_cut = T_CUT,
            method_Q = if (RUN_ML) paste0("ml_", ML_METHOD) else "glm",
            logrr = NA_real_, se_logrr = NA_real_,
            rd = NA_real_, se_rd = NA_real_,
            risk1 = NA_real_, risk0 = NA_real_,
            error = msg
          ) else data.table::data.table()

          ipcw_tu_dt <- if (RUN_TIME_UPDATED_SENS) data.table::data.table(
            scenario_id = scn$scenario_id,
            axisA = scn$axisA, axisB = scn$axisB, rho_meas = scn$rho_meas,
            mis_spec = scn$mis_spec,
            truth = truth_j,
            replicate = rep_j,
            t_cut = T_CUT,
            method = paste0("ipcw_timeupdated_", IPCW_ESTIMATOR),
            dev_model_status = NA_character_,
            logrr = NA_real_, se_logrr = NA_real_,
            rd = NA_real_, se_rd = NA_real_,
            risk1 = NA_real_, risk0 = NA_real_,
            error = msg
          ) else data.table::data.table()

          tm_dt <- if (RUN_TMLE) data.table::data.table(
            scenario_id = scn$scenario_id,
            axisA = scn$axisA, axisB = scn$axisB, rho_meas = scn$rho_meas,
            mis_spec = scn$mis_spec,
            truth = truth_j,
            replicate = rep_j,
            t_cut = T_CUT,
            method_Q = "tmle",
            logrr = NA_real_, se_logrr = NA_real_,
            rd = NA_real_, se_rd = NA_real_,
            risk1 = NA_real_, risk0 = NA_real_,
            error = msg
          ) else data.table::data.table()

          list(
            rep = rep_dt,
            spd = data.table::data.table(
              scenario_id = scn$scenario_id,
              mis_spec = scn$mis_spec,
              truth = truth_j,
              replicate = rep_j,
              interval = NA_integer_, tstart = NA_real_, tstop = NA_real_, t_mid = NA_real_,
              gamma_hat = NA_real_, se = NA_real_, HR_1SD = NA_real_, CI_low = NA_real_, CI_high = NA_real_,
              Gamma_hat = NA_real_,
              error = msg
            ),
            diag = data.table::data.table(
              scenario_id = scn$scenario_id,
              mis_spec = scn$mis_spec,
              truth = truth_j,
              replicate = rep_j,
              t = NA_real_, N = NA_integer_, ESS = NA_real_, rESS = NA_real_, tail_share = NA_real_, tipped = NA_integer_,
              error = msg
            ),
            ipcw_risk = ipcw_risk_dt,
            ipcw_tu_risk = ipcw_tu_dt,
            dr = dr_dt,
            tm = tm_dt
          )
        }
      )
      t1 <- proc.time()[3]
      out$rep[, runtime_sec := (t1 - t0)]
      out
    },
    future.seed = TRUE
  )

  rep_dt_i  <- data.table::rbindlist(lapply(res_list, `[[`, "rep"),  fill = TRUE)
  spd_dt_i  <- data.table::rbindlist(lapply(res_list, `[[`, "spd"),  fill = TRUE)
  diag_dt_i <- data.table::rbindlist(lapply(res_list, `[[`, "diag"), fill = TRUE)
  ipcw_risk_dt_i <- data.table::rbindlist(lapply(res_list, `[[`, "ipcw_risk"), fill = TRUE)
  ipcw_tu_risk_dt_i <- data.table::rbindlist(lapply(res_list, function(x) x[["ipcw_tu_risk"]]), fill = TRUE)
  dr_dt_i   <- data.table::rbindlist(lapply(res_list, `[[`, "dr"),   fill = TRUE)
  tm_dt_i   <- data.table::rbindlist(lapply(res_list, `[[`, "tm"),   fill = TRUE)

  # log errors
  if ("error" %in% names(rep_dt_i)) {
    err_i <- rep_dt_i[!is.na(error) & nzchar(error)]
    if (nrow(err_i) > 0) {
      cfg$log_line(log_file, sprintf("[%s] SCENARIO %s mis_spec=%d had %d replicate errors (kept as NA)\n",
                                     format(Sys.time(), "%H:%M:%S"), scn$scenario_id, scn$mis_spec, nrow(err_i)))
    }
  }

  rep_chunks[[i]]  <- rep_dt_i
  spd_chunks[[i]]  <- spd_dt_i
  diag_chunks[[i]] <- diag_dt_i
  ipcw_risk_chunks[[i]] <- ipcw_risk_dt_i
  ipcw_tu_risk_chunks[[i]] <- ipcw_tu_risk_dt_i
  dr_chunks[[i]]   <- dr_dt_i
  tm_chunks[[i]]   <- tm_dt_i

  cfg$log_line(log_file, sprintf("[%s] SCENARIO %d/%d done | %s | mis_spec=%d\n",
                                 format(Sys.time(), "%H:%M:%S"), i, nrow(scn_grid), scn$scenario_id, scn$mis_spec))
}

rep_dt  <- data.table::rbindlist(rep_chunks,  fill = TRUE)
spd_dt  <- data.table::rbindlist(spd_chunks,  fill = TRUE)
diag_dt <- data.table::rbindlist(diag_chunks, fill = TRUE)
ipcw_risk_dt <- data.table::rbindlist(ipcw_risk_chunks, fill = TRUE)
ipcw_tu_risk_dt <- data.table::rbindlist(ipcw_tu_risk_chunks, fill = TRUE)
dr_dt   <- data.table::rbindlist(dr_chunks,   fill = TRUE)
tm_dt   <- data.table::rbindlist(tm_chunks,   fill = TRUE)
# Guardrail: "truth" should never be NA in Paper B simulations
coerce_truth <- function(dt, dt_name, required = FALSE) {
  if (!("truth" %in% names(dt))) {
    if (required) stop(paste0("Missing 'truth' column in ", dt_name, "."))
    return(invisible(TRUE))
  }
  dt[, truth := as.character(truth)]
  if (anyNA(dt$truth)) {
    bad <- dt[is.na(truth)]
    if (nrow(bad) > 5) bad <- bad[1:5]
    msg <- paste(utils::capture.output(print(bad)), collapse = "\n")
    stop(paste0("truth has NA in ", dt_name, ". TRUTH_LEVELS=", Sys.getenv("TRUTH_LEVELS", unset = ""),
                " unique(truth)=", paste(unique(dt$truth), collapse = ","), "\n", msg))
  }
  invisible(TRUE)
}
coerce_truth(rep_dt, "rep_dt", required = TRUE)
coerce_truth(ipcw_risk_dt, "ipcw_risk_dt", required = TRUE)
if (nrow(ipcw_tu_risk_dt) > 0) coerce_truth(ipcw_tu_risk_dt, "ipcw_tu_risk_dt")
coerce_truth(diag_dt, "diag_dt", required = TRUE)
if (nrow(dr_dt) > 0) coerce_truth(dr_dt, "dr_dt")
if (nrow(tm_dt) > 0) coerce_truth(tm_dt, "tm_dt")



get_true_logrr <- function(scenario_id, truth) {
  truth <- as.character(truth)
  if (is.na(truth) || !nzchar(truth) || truth == "null") return(0)
  if (is.null(mc_truth_risk)) return(NA_real_)
  sid <- scenario_id
  tr  <- truth
  v <- mc_truth_risk[scenario_id == sid & truth == tr, logrr_true]
  if (length(v) > 0 && is.finite(v[1])) return(as.numeric(v[1]))
  NA_real_
}

get_true_rd <- function(scenario_id, truth) {
  truth <- as.character(truth)
  if (is.na(truth) || !nzchar(truth) || truth == "null") return(0)
  if (is.null(mc_truth_risk)) return(NA_real_)
  sid <- scenario_id
  tr  <- truth
  v <- mc_truth_risk[scenario_id == sid & truth == tr, rd_true]
  if (length(v) > 0 && is.finite(v[1])) return(as.numeric(v[1]))
  NA_real_
}


# performance summaries
perf_hr <- rep_dt[, {
  beta_true <- if (isTRUE(.BY$truth == "non_null")) beta_true_pp else 0
  s <- summarize_perf(beta_hat, se_hat, beta_true)
  sign_error <- if (is.na(beta_true) || isTRUE(beta_true == 0)) NA_real_ else mean(sign(beta_hat) != sign(beta_true), na.rm = TRUE)
  .(bias = s$bias, rmse = s$rmse, cover = s$cover, ci_excl0 = s$ci_excl0,
    sign_error = sign_error,
    ess_median = median(ess_final, na.rm = TRUE),
    tip_time_median = stats::median(tip_time, na.rm = TRUE),
    runtime_median_sec = stats::median(runtime_sec, na.rm = TRUE))
}, by = .(truth, scenario_id, axisA, axisB, rho_meas, mis_spec)]

perf_dr <- NULL
if (nrow(dr_dt) > 0) {
  perf_dr <- dr_dt[, {
    beta_true <- get_true_logrr(.BY$scenario_id, .BY$truth)
    s <- summarize_perf(logrr, se_logrr, beta_true)
    sign_error <- if (is.na(beta_true) || isTRUE(beta_true == 0)) NA_real_ else mean(sign(logrr) != sign(beta_true), na.rm = TRUE)
    .(bias = s$bias, rmse = s$rmse, cover = s$cover, ci_excl0 = s$ci_excl0,
      sign_error = sign_error,
      logrr_median = stats::median(logrr, na.rm = TRUE),
      rd_median = stats::median(rd, na.rm = TRUE))
  }, by = .(truth, scenario_id, axisA, axisB, rho_meas, mis_spec, method_Q)]
}

perf_ipcw_risk <- ipcw_risk_dt[, {
  beta_true <- get_true_logrr(.BY$scenario_id, .BY$truth)
  s <- summarize_perf(logrr, se_logrr, beta_true)
  sign_error <- if (is.na(beta_true) || isTRUE(beta_true == 0)) NA_real_ else mean(sign(logrr) != sign(beta_true), na.rm = TRUE)
  .(bias = s$bias, rmse = s$rmse, cover = s$cover, ci_excl0 = s$ci_excl0,
    sign_error = sign_error,
    logrr_median = stats::median(logrr, na.rm = TRUE),
    rd_median = stats::median(rd, na.rm = TRUE))
}, by = .(truth, scenario_id, axisA, axisB, rho_meas, mis_spec)]

perf_ipcw_tu <- NULL
if (nrow(ipcw_tu_risk_dt) > 0) {
  perf_ipcw_tu <- ipcw_tu_risk_dt[, {
    beta_true <- get_true_logrr(.BY$scenario_id, .BY$truth)
    s <- summarize_perf(logrr, se_logrr, beta_true)
    sign_error <- if (is.na(beta_true) || isTRUE(beta_true == 0)) NA_real_ else mean(sign(logrr) != sign(beta_true), na.rm = TRUE)
    .(bias = s$bias, rmse = s$rmse, cover = s$cover, ci_excl0 = s$ci_excl0,
      sign_error = sign_error,
      logrr_median = stats::median(logrr, na.rm = TRUE),
      rd_median = stats::median(rd, na.rm = TRUE),
      dev_model_status = paste(sort(unique(dev_model_status)), collapse = ";"),
      min_c_hat_median = stats::median(min_c_hat, na.rm = TRUE))
  }, by = .(truth, scenario_id, axisA, axisB, rho_meas, mis_spec)]
}

perf_tm <- NULL
if (nrow(tm_dt) > 0) {
  perf_tm <- tm_dt[, {
    beta_true <- get_true_logrr(.BY$scenario_id, .BY$truth)
    s <- summarize_perf(logrr, se_logrr, beta_true)
    sign_error <- if (is.na(beta_true) || isTRUE(beta_true == 0)) NA_real_ else mean(sign(logrr) != sign(beta_true), na.rm = TRUE)
    .(bias = s$bias, rmse = s$rmse, cover = s$cover, ci_excl0 = s$ci_excl0,
      sign_error = sign_error,
      logrr_median = stats::median(logrr, na.rm = TRUE),
      rd_median = stats::median(rd, na.rm = TRUE))
  }, by = .(truth, scenario_id, axisA, axisB, rho_meas, mis_spec, method_Q)]
}

# save outputs
raw_dir <- file.path(cfg$out_dir, "raw")
dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)

data.table::fwrite(rep_dt,  file.path(raw_dir, "replicate_results_rescue.csv"))
data.table::fwrite(spd_dt,  file.path(raw_dir, "spd_curves_rescue.csv"))
data.table::fwrite(diag_dt, file.path(raw_dir, "weight_diagnostics_rescue.csv"))
data.table::fwrite(ipcw_risk_dt, file.path(raw_dir, "ipcw_risk_rescue.csv"))
if (nrow(ipcw_tu_risk_dt) > 0) data.table::fwrite(ipcw_tu_risk_dt, file.path(raw_dir, "ipcw_timeupdated_risk_rescue.csv"))
if (nrow(dr_dt) > 0) data.table::fwrite(dr_dt, file.path(raw_dir, "dr_risk_rescue.csv"))
if (nrow(tm_dt) > 0) data.table::fwrite(tm_dt, file.path(raw_dir, "tmle_risk_rescue.csv"))

data.table::fwrite(perf_hr, file.path(cfg$out_dir, "perf_summary_rescue_hr.csv"))
data.table::fwrite(perf_ipcw_risk, file.path(cfg$out_dir, "perf_summary_rescue_ipcw_risk.csv"))
if (!is.null(perf_ipcw_tu)) data.table::fwrite(perf_ipcw_tu, file.path(cfg$out_dir, "perf_summary_rescue_ipcw_timeupdated.csv"))
if (!is.null(perf_dr)) data.table::fwrite(perf_dr, file.path(cfg$out_dir, "perf_summary_rescue_dr.csv"))
if (!is.null(perf_tm)) data.table::fwrite(perf_tm, file.path(cfg$out_dir, "perf_summary_rescue_tmle.csv"))

cat("\nSaved outputs to:", cfg$out_dir, "\n")
cfg$write_session_info("rescue")
cfg$log_line(log_file, sprintf("[%s] Rescue simulation finished\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
