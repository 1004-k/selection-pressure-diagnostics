#!/usr/bin/env Rscript
# scripts/05_rerun_timeupdated_sensitivity_only.R
# ------------------------------------------------------------
# Rerun ONLY the time-updated IPCW sensitivity estimator.
#
# Use case:
# The full JCI v5 run completed and primary outputs are valid, but the
# time-updated sensitivity raw output used glm_fallback in every replicate.
# This script recreates the same simulated datasets using the same seed scheme
# and overwrites only:
#   output_jci_final/raw/ipcw_timeupdated_risk_rescue.csv
#   output_jci_final/perf_summary_rescue_ipcw_timeupdated.csv
#   time-updated supplementary tables and audit files
#
# It does not overwrite the primary standardized IPCW, DR-AIPW, SPD, or weight
# diagnostic outputs.
# ------------------------------------------------------------

source("scripts/00_utils.R")

cfg <- init_project(
  n_cores = as.integer(Sys.getenv("N_CORES", "3")),
  seed    = 2026L,
  out_dir = Sys.getenv("OUT_DIR", "output_jci_final")
)
cfg$require_pkgs(c("data.table", "future.apply", "future", "splines"))

B     <- as.integer(Sys.getenv("B", "200"))
N     <- as.integer(Sys.getenv("N", "2000"))
t_max <- as.numeric(Sys.getenv("T_MAX", "5"))
dt    <- as.numeric(Sys.getenv("DT", "0.25"))
T_CUT <- as.numeric(Sys.getenv("T_CUT", as.character(t_max)))
IPCW_ESTIMATOR <- Sys.getenv("IPCW_ESTIMATOR", "hajek")
if (!IPCW_ESTIMATOR %in% c("hajek", "ht")) IPCW_ESTIMATOR <- "hajek"

parse_numeric_env <- function(name, default = "0") {
  raw <- trimws(Sys.getenv(name, unset = default))
  out <- suppressWarnings(as.numeric(raw))
  if (length(out) == 1L && is.finite(out)) return(out)
  if (!grepl("^[0-9eE+./*()[:space:]a-zA-Z_-]+$", raw)) {
    stop(name, " must be numeric or a simple numeric expression; got: ", raw)
  }
  val <- tryCatch(eval(parse(text = raw), envir = baseenv()), error = function(e) NA_real_)
  if (!is.numeric(val) || length(val) != 1L || !is.finite(val)) {
    stop(name, " could not be parsed as a finite numeric value; got: ", raw)
  }
  as.numeric(val)
}
beta_true_pp <- parse_numeric_env("BETA_TRUE_PP", "-0.22314355131420976")

truth_levels <- trimws(strsplit(Sys.getenv("TRUTH_LEVELS", "null,non_null"), ",", fixed = TRUE)[[1]])
truth_levels <- truth_levels[truth_levels %in% c("null", "non_null")]
if (length(truth_levels) == 0) truth_levels <- c("null", "non_null")

MIS_SPEC_LEVELS <- trimws(strsplit(Sys.getenv("MIS_SPEC_LEVELS", "0,1"), ",", fixed = TRUE)[[1]])
MIS_SPEC_LEVELS <- as.integer(MIS_SPEC_LEVELS[MIS_SPEC_LEVELS %in% c("0", "1")])
if (length(MIS_SPEC_LEVELS) == 0) MIS_SPEC_LEVELS <- c(0L, 1L)

# Scenario grid, with optional subset for debugging.
grid <- make_scenario_grid()
scn_grid <- data.table::CJ(scenario_id = grid$scenario_id, mis_spec = MIS_SPEC_LEVELS, unique = TRUE)
scn_grid <- merge(scn_grid, grid, by = "scenario_id", all.x = TRUE)
data.table::setorder(scn_grid, scenario_num, mis_spec)
subset_ids <- trimws(strsplit(Sys.getenv("SCENARIO_SUBSET", ""), ",", fixed = TRUE)[[1]])
subset_ids <- subset_ids[nzchar(subset_ids)]
if (length(subset_ids) > 0) {
  scn_grid <- scn_grid[scenario_id %in% subset_ids]
  data.table::setorder(scn_grid, scenario_num, mis_spec)
}

# Load or compute intervention truth. Full reruns should load the existing full-run truth.
truth_path <- file.path(cfg$out_dir, "raw", "mc_truth_risk.csv")
mc_truth_risk <- NULL
if (file.exists(truth_path)) {
  mc_truth_risk <- data.table::fread(truth_path, showProgress = FALSE)
  message("Loaded existing intervention MC truth: ", truth_path)
} else {
  N_TRUTH <- as.integer(Sys.getenv("N_TRUTH", "5000"))
  B_TRUTH <- as.integer(Sys.getenv("B_TRUTH", "1"))
  message("No existing truth file found. Computing intervention MC truth for this output folder (N_TRUTH=", N_TRUTH, ", B_TRUTH=", B_TRUTH, ") ...")
  truth_jobs <- data.table::CJ(scenario_id = grid$scenario_id, truth = truth_levels, unique = TRUE)
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
  dir.create(dirname(truth_path), showWarnings = FALSE, recursive = TRUE)
  data.table::fwrite(mc_truth_risk, truth_path)
  message("Saved intervention MC truth: ", truth_path)
}

get_true_logrr <- function(scenario_id, truth) {
  truth <- as.character(truth)
  if (is.na(truth) || !nzchar(truth) || truth == "null") return(0)
  v <- mc_truth_risk[scenario_id == scenario_id & truth == truth, logrr_true]
  # avoid data.table scoping ambiguity
  v <- mc_truth_risk[mc_truth_risk$scenario_id == scenario_id & mc_truth_risk$truth == truth, "logrr_true"][[1]]
  if (length(v) > 0 && is.finite(v[1])) return(as.numeric(v[1]))
  NA_real_
}

# Corrected helper avoiding data.table variable masking.
get_true_logrr <- function(sid, tr) {
  if (is.na(tr) || !nzchar(tr) || tr == "null") return(0)
  v <- mc_truth_risk[scenario_id == sid & truth == tr, logrr_true]
  if (length(v) > 0 && is.finite(v[1])) return(as.numeric(v[1]))
  NA_real_
}

# Logging.
log_dir <- file.path(cfg$out_dir, "logs")
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(log_dir, sprintf("timeupdated_sensitivity_rerun_B%d_N%d.log", B, N))
if (file.exists(log_file)) file.remove(log_file)
log_line <- function(...) {
  txt <- paste0(...)
  cat(txt)
  cat(txt, file = log_file, append = TRUE)
}

log_line(sprintf("[%s] Time-updated IPCW sensitivity rerun started | OUT_DIR=%s | scenarios=%d | mis_spec=%s | B=%d | N=%d | truths=%s | N_CORES=%d\n",
                 format(Sys.time(), "%Y-%m-%d %H:%M:%S"), cfg$out_dir, nrow(scn_grid),
                 paste(MIS_SPEC_LEVELS, collapse = ","), B, N, paste(truth_levels, collapse = ","), cfg$n_cores))

future::plan(future::multisession, workers = cfg$n_cores)

run_one_tu <- function(scn_row, truth, rep_id) {
  beta <- if (identical(as.character(truth), "non_null")) beta_true_pp else 0
  seed <- seed_for_job(cfg$seed * 100000L, scn_row$scenario_id, truth = truth, rep_id = rep_id)
  sim <- simulate_one_dataset(
    N = N,
    t_max = t_max,
    dt = dt,
    beta_true = beta,
    axisA = scn_row$axisA,
    axisB = scn_row$axisB,
    rho_meas = scn_row$rho_meas,
    seed = seed
  )
  tu <- fit_ipcw_risk_timeupdated(
    long_dt = sim$long_dt,
    pp_dt = sim$pp_dt,
    t_cut = T_CUT,
    mis_spec = scn_row$mis_spec,
    estimator = IPCW_ESTIMATOR,
    w_floor = 0.05,
    ps_floor = 1e-3,
    ns_df = 3L
  )
  data.table::data.table(
    scenario_id = scn_row$scenario_id,
    axisA = scn_row$axisA,
    axisB = scn_row$axisB,
    rho_meas = scn_row$rho_meas,
    mis_spec = scn_row$mis_spec,
    truth = as.character(truth),
    replicate = as.integer(rep_id),
    t_cut = tu$t,
    method = tu$method,
    dev_model_status = tu$dev_model_status,
    logrr = tu$logrr,
    se_logrr = tu$se_logrr,
    rd = tu$rd,
    se_rd = tu$se_rd,
    risk1 = tu$risk1,
    risk0 = tu$risk0,
    den1 = tu$den1,
    den0 = tu$den0,
    min_g1 = tu$min_g1,
    max_g1 = tu$max_g1,
    min_c_hat = tu$min_c_hat,
    median_c_hat = tu$median_c_hat,
    median_intervals = tu$median_intervals,
    mean_interval_dev_prob = tu$mean_interval_dev_prob
  )
}

chunks <- vector("list", nrow(scn_grid))
for (i in seq_len(nrow(scn_grid))) {
  scn <- scn_grid[i]
  msg <- sprintf("[%s] SCENARIO %d/%d start | %s | mis_spec=%d\n",
                 format(Sys.time(), "%H:%M:%S"), i, nrow(scn_grid), scn$scenario_id, scn$mis_spec)
  log_line(msg)
  jobs <- data.table::CJ(truth = truth_levels, replicate = seq_len(B))
  jobs <- jobs[order(truth, replicate)]
  res <- future.apply::future_lapply(
    seq_len(nrow(jobs)),
    function(j) {
      truth_j <- as.character(jobs$truth[j])[1]
      rep_j <- as.integer(jobs$replicate[j])[1]
      tryCatch(
        run_one_tu(scn, truth_j, rep_j),
        error = function(e) data.table::data.table(
          scenario_id = scn$scenario_id,
          axisA = scn$axisA,
          axisB = scn$axisB,
          rho_meas = scn$rho_meas,
          mis_spec = scn$mis_spec,
          truth = truth_j,
          replicate = rep_j,
          t_cut = T_CUT,
          method = paste0("ipcw_timeupdated_", IPCW_ESTIMATOR),
          dev_model_status = "error",
          logrr = NA_real_, se_logrr = NA_real_, rd = NA_real_, se_rd = NA_real_,
          risk1 = NA_real_, risk0 = NA_real_,
          error = conditionMessage(e)
        )
      )
    },
    future.seed = TRUE
  )
  chunks[[i]] <- data.table::rbindlist(res, fill = TRUE)
  log_line(sprintf("[%s] SCENARIO %d/%d done | %s | mis_spec=%d\n",
                   format(Sys.time(), "%H:%M:%S"), i, nrow(scn_grid), scn$scenario_id, scn$mis_spec))
}

ipcw_tu <- data.table::rbindlist(chunks, fill = TRUE)
ipcw_tu[, truth := as.character(truth)]
raw_dir <- file.path(cfg$out_dir, "raw")
dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)
raw_path <- file.path(raw_dir, "ipcw_timeupdated_risk_rescue.csv")
data.table::fwrite(ipcw_tu, raw_path)
message("Saved patched time-updated raw output: ", raw_path)

# Performance summary.
perf_tu <- ipcw_tu[, {
  beta_true <- get_true_logrr(.BY$scenario_id, .BY$truth)
  s <- summarize_perf(logrr, se_logrr, beta_true)
  sign_error <- if (is.na(beta_true) || isTRUE(beta_true == 0)) NA_real_ else mean(sign(logrr) != sign(beta_true), na.rm = TRUE)
  .(
    bias = s$bias,
    rmse = s$rmse,
    cover = s$cover,
    ci_excl0 = s$ci_excl0,
    sign_error = sign_error,
    logrr_median = stats::median(logrr, na.rm = TRUE),
    rd_median = stats::median(rd, na.rm = TRUE),
    dev_model_status = paste(sort(unique(dev_model_status)), collapse = ";"),
    pct_glm_timeupdated = mean(grepl("^glm_timeupdated", dev_model_status), na.rm = TRUE),
    min_c_hat_median = stats::median(min_c_hat, na.rm = TRUE)
  )
}, by = .(truth, scenario_id, axisA, axisB, rho_meas, mis_spec)]
perf_path <- file.path(cfg$out_dir, "perf_summary_rescue_ipcw_timeupdated.csv")
data.table::fwrite(perf_tu, perf_path)
message("Saved patched time-updated performance summary: ", perf_path)

# Supplementary tables.
tab_dir <- file.path(cfg$out_dir, "tables")
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
data.table::fwrite(perf_tu, file.path(tab_dir, "OnlineTable_ipcw_timeupdated_sensitivity_performance.csv"))

tu_compact <- perf_tu[, .(
  TypeI_timeupdated_IPCW = median(ci_excl0[truth == "null"], na.rm = TRUE),
  Coverage_timeupdated_IPCW = median(cover[truth == "null"], na.rm = TRUE),
  Power_timeupdated_IPCW = median(ci_excl0[truth == "non_null"], na.rm = TRUE),
  SignError_timeupdated_IPCW = median(sign_error[truth == "non_null"], na.rm = TRUE),
  Median_min_c_hat = median(min_c_hat_median, na.rm = TRUE),
  Median_pct_glm_timeupdated = median(pct_glm_timeupdated, na.rm = TRUE)
), by = .(nuisance_regime = mis_spec)]
data.table::setorder(tu_compact, nuisance_regime)
data.table::fwrite(tu_compact, file.path(tab_dir, "OnlineTable_timeupdated_sensitivity_compact.csv"))

base_perf_path <- file.path(cfg$out_dir, "perf_summary_rescue_ipcw_risk.csv")
if (file.exists(base_perf_path)) {
  base <- data.table::fread(base_perf_path, showProgress = FALSE)
  base[, truth := as.character(truth)]
  base_vs_tu <- merge(
    base[truth == "non_null", .(scenario_id, mis_spec, axisA, axisB, rho_meas,
                                power_baseline = ci_excl0,
                                sign_error_baseline = sign_error,
                                rmse_baseline = rmse)],
    perf_tu[truth == "non_null", .(scenario_id, mis_spec,
                                   power_timeupdated = ci_excl0,
                                   sign_error_timeupdated = sign_error,
                                   rmse_timeupdated = rmse,
                                   pct_glm_timeupdated)],
    by = c("scenario_id", "mis_spec"), all = FALSE
  )
  base_vs_tu[, `:=`(
    delta_power_timeupdated_minus_baseline = power_timeupdated - power_baseline,
    delta_sign_error_timeupdated_minus_baseline = sign_error_timeupdated - sign_error_baseline,
    delta_rmse_baseline_minus_timeupdated = rmse_baseline - rmse_timeupdated
  )]
  data.table::fwrite(base_vs_tu, file.path(tab_dir, "OnlineTable_baseline_vs_timeupdated_ipcw_sensitivity.csv"))
}

status_tab <- ipcw_tu[, .N, by = .(dev_model_status)]
status_tab[, pct := N / sum(N)]
data.table::fwrite(status_tab, file.path(tab_dir, "OnlineTable_timeupdated_model_status.csv"))

# Dedicated audit.
audit_dir <- file.path(cfg$out_dir, "audit")
dir.create(audit_dir, showWarnings = FALSE, recursive = TRUE)
expected_rows <- nrow(scn_grid) * length(truth_levels) * B
n_bad <- ipcw_tu[, sum(!is.finite(logrr) | !is.finite(se_logrr) | !is.finite(rd) | !is.finite(se_rd))]
pct_tu <- ipcw_tu[, mean(grepl("^glm_timeupdated", dev_model_status), na.rm = TRUE)]
audit_lines <- c(
  sprintf("Time-updated IPCW sensitivity patch rerun at %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("OUT_DIR=%s", cfg$out_dir),
  sprintf("B=%d; N=%d; expected_rows=%d; observed_rows=%d", B, N, expected_rows, nrow(ipcw_tu)),
  sprintf("nonfinite core estimates: %d", n_bad),
  sprintf("proportion glm_timeupdated status: %.4f", pct_tu),
  "",
  "Model status counts:",
  paste(utils::capture.output(print(status_tab)), collapse = "\n")
)
if (nrow(ipcw_tu) != expected_rows) audit_lines <- c(audit_lines, "FAIL row count mismatch")
if (n_bad > 0) audit_lines <- c(audit_lines, "WARN nonfinite core estimates present")
if (!is.finite(pct_tu) || pct_tu < 0.80) audit_lines <- c(audit_lines, "FAIL too few rows used an actual time-updated glm model")
writeLines(audit_lines, file.path(audit_dir, "timeupdated_sensitivity_patch_qa.txt"))
cat(paste(audit_lines, collapse = "\n"), "\n")

log_line(sprintf("[%s] Time-updated IPCW sensitivity rerun finished\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
