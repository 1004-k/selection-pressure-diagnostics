#!/usr/bin/env Rscript
# QA checks for the standardized-IPCW rerun outputs, patched to flag whether
# the time-updated IPCW sensitivity actually used the time-updated glm model.

suppressPackageStartupMessages(library(data.table))
out_dir <- Sys.getenv("OUT_DIR", "output_jci_final")
qa_dir <- file.path(out_dir, "audit")
dir.create(qa_dir, showWarnings = FALSE, recursive = TRUE)

B <- as.integer(Sys.getenv("B", "200"))
truth_levels <- trimws(strsplit(Sys.getenv("TRUTH_LEVELS", "null,non_null"), ",", fixed = TRUE)[[1]])
truth_levels <- truth_levels[truth_levels %in% c("null", "non_null")]
if (length(truth_levels) == 0) truth_levels <- c("null", "non_null")
mis <- trimws(strsplit(Sys.getenv("MIS_SPEC_LEVELS", "0,1"), ",", fixed = TRUE)[[1]])
mis <- as.integer(mis[mis %in% c("0", "1")])
if (length(mis) == 0) mis <- c(0L, 1L)

expected_rows <- 18L * length(mis) * length(truth_levels) * B
checks <- list()
check_file <- function(rel, expected = expected_rows, finite_cols = character()) {
  path <- file.path(out_dir, rel)
  if (!file.exists(path)) {
    checks[[length(checks) + 1L]] <<- sprintf("FAIL missing: %s", rel)
    return(invisible(NULL))
  }
  dt <- fread(path, showProgress = FALSE)
  checks[[length(checks) + 1L]] <<- sprintf("%s rows: observed=%d expected=%d", rel, nrow(dt), expected)
  if (nrow(dt) != expected) checks[[length(checks) + 1L]] <<- sprintf("FAIL row count mismatch: %s", rel)
  for (cc in finite_cols) {
    if (cc %in% names(dt)) {
      n_bad <- dt[, sum(!is.finite(get(cc)))]
      checks[[length(checks) + 1L]] <<- sprintf("%s nonfinite %s: %d", rel, cc, n_bad)
      if (n_bad > 0) checks[[length(checks) + 1L]] <<- sprintf("WARN nonfinite values found in %s:%s", rel, cc)
    }
  }
  invisible(dt)
}

check_file("raw/replicate_results_rescue.csv", finite_cols = c("beta_hat", "se_hat", "max_spd"))
check_file("raw/ipcw_risk_rescue.csv", finite_cols = c("logrr", "se_logrr", "rd", "se_rd", "risk1", "risk0"))
check_file("raw/dr_risk_rescue.csv", finite_cols = c("logrr", "se_logrr", "rd", "se_rd", "risk1", "risk0"))
check_file("raw/spd_curves_rescue.csv", expected = expected_rows * 5L, finite_cols = c("gamma_hat"))
check_file("raw/weight_diagnostics_rescue.csv", expected = expected_rows * 20L, finite_cols = c("rESS", "tail_share"))

summary_expected <- 18L * length(mis) * length(truth_levels)
check_file("perf_summary_rescue_ipcw_risk.csv", expected = summary_expected, finite_cols = c("ci_excl0", "cover"))
check_file("perf_summary_rescue_dr.csv", expected = summary_expected, finite_cols = c("ci_excl0", "cover"))

truth_path <- file.path(out_dir, "raw", "mc_truth_risk.csv")
if ("non_null" %in% truth_levels) {
  if (!file.exists(truth_path)) {
    checks[[length(checks) + 1L]] <- "FAIL missing intervention truth file: raw/mc_truth_risk.csv"
  } else {
    tr <- fread(truth_path, showProgress = FALSE)
    checks[[length(checks) + 1L]] <- sprintf("raw/mc_truth_risk.csv rows: observed=%d expected=%d", nrow(tr), 18L * length(truth_levels))
    if (!("truth_method" %in% names(tr))) {
      checks[[length(checks) + 1L]] <- "FAIL truth_method missing in mc_truth_risk.csv"
    } else if (!all(grepl("^intervention", tr$truth_method))) {
      checks[[length(checks) + 1L]] <- "FAIL mc_truth_risk.csv is not intervention-based"
    }
    for (cc in c("logrr_true", "rd_true")) {
      if (cc %in% names(tr)) {
        n_bad <- tr[, sum(!is.finite(get(cc)))]
        checks[[length(checks) + 1L]] <- sprintf("raw/mc_truth_risk.csv nonfinite %s: %d", cc, n_bad)
        if (n_bad > 0) checks[[length(checks) + 1L]] <- sprintf("FAIL nonfinite truth values in %s", cc)
      }
    }
  }
}

run_tu <- as.integer(Sys.getenv("RUN_TIME_UPDATED_SENS", "1")) == 1L
if (run_tu) {
  tu_raw <- check_file("raw/ipcw_timeupdated_risk_rescue.csv", finite_cols = c("logrr", "se_logrr", "rd", "se_rd", "risk1", "risk0"))
  check_file("perf_summary_rescue_ipcw_timeupdated.csv", expected = summary_expected, finite_cols = c("ci_excl0", "cover"))
  if (!is.null(tu_raw) && "dev_model_status" %in% names(tu_raw)) {
    status_tab <- tu_raw[, .N, by = dev_model_status][order(-N)]
    status_tab[, pct := N / sum(N)]
    status_txt <- paste(utils::capture.output(print(status_tab)), collapse = "\n")
    checks[[length(checks) + 1L]] <- "raw/ipcw_timeupdated_risk_rescue.csv dev_model_status counts:"
    checks[[length(checks) + 1L]] <- status_txt
    pct_tu <- tu_raw[, mean(grepl("^glm_timeupdated", dev_model_status), na.rm = TRUE)]
    checks[[length(checks) + 1L]] <- sprintf("raw/ipcw_timeupdated_risk_rescue.csv proportion glm_timeupdated status: %.4f", pct_tu)
    if (!is.finite(pct_tu) || pct_tu < 0.80) {
      checks[[length(checks) + 1L]] <- "FAIL time-updated sensitivity did not primarily use a time-updated glm model"
    } else if (pct_tu < 0.95) {
      checks[[length(checks) + 1L]] <- "WARN some time-updated sensitivity replicates used fallback/error status"
    }
  }
}

out <- c(
  sprintf("JCI standardized-IPCW QA run at %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("OUT_DIR=%s", out_dir),
  sprintf("B=%d; truth_levels=%s; mis_spec=%s; time_updated_sensitivity=%d", B, paste(truth_levels, collapse = ","), paste(mis, collapse = ","), as.integer(run_tu)),
  "",
  unlist(checks)
)
writeLines(out, con = file.path(qa_dir, "jci_output_qa_report.txt"))
cat(paste(out, collapse = "\n"), "\n")
