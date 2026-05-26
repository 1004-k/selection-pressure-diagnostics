#!/usr/bin/env Rscript
# run_jci_timeupdated_only.R
# ------------------------------------------------------------
# R-only runner for the patched time-updated IPCW sensitivity rerun.
#
# Usage from project root:
#   Rscript run_jci_timeupdated_only.R smoke
#   Rscript run_jci_timeupdated_only.R full
#
# The full mode is intended to be run in the existing JCI v5 project folder
# after the primary full run has completed. It reuses output_jci_final/raw/
# mc_truth_risk.csv and overwrites only the time-updated IPCW sensitivity output.
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) >= 1) tolower(args[1]) else "full"
if (!mode %in% c("smoke", "full")) stop("mode must be 'smoke' or 'full'")

set_default_env <- function(name, value) {
  current <- Sys.getenv(name, unset = NA_character_)
  if (is.na(current) || !nzchar(current)) {
    do.call(Sys.setenv, stats::setNames(as.list(value), name))
  }
}

required_files <- c(
  "scripts/00_utils.R",
  "scripts/05_rerun_timeupdated_sensitivity_only.R",
  "R/12_ipcw_timeupdated_sensitivity.R"
)
missing <- required_files[!file.exists(required_files)]
if (length(missing) > 0) {
  stop("Required files missing. Run from the project root. Missing:\n", paste(missing, collapse = "\n"))
}

if (mode == "smoke") {
  set_default_env("OUT_DIR", "output_jci_tu_smoke")
  set_default_env("B", "3")
  set_default_env("N", "300")
  set_default_env("TRUTH_LEVELS", "null,non_null")
  set_default_env("MIS_SPEC_LEVELS", "0,1")
  set_default_env("BETA_TRUE_PP", "-0.22314355131420976")
  set_default_env("N_TRUTH", "5000")
  set_default_env("B_TRUTH", "1")
  set_default_env("IPCW_ESTIMATOR", "hajek")
  set_default_env("N_CORES", "4")
}

if (mode == "full") {
  set_default_env("OUT_DIR", "output_jci_final")
  set_default_env("B", "200")
  set_default_env("N", "2000")
  set_default_env("TRUTH_LEVELS", "null,non_null")
  set_default_env("MIS_SPEC_LEVELS", "0,1")
  set_default_env("BETA_TRUE_PP", "-0.22314355131420976")
  set_default_env("IPCW_ESTIMATOR", "hajek")
  set_default_env("RUN_TIME_UPDATED_SENS", "1")
  # Do not overwrite N_CORES if user already exported it.
  detected <- tryCatch(parallel::detectCores(logical = TRUE), error = function(e) 4L)
  default_cores <- max(1L, min(12L, detected - 2L))
  set_default_env("N_CORES", as.character(default_cores))
}

cat("\n============================================================\n")
cat("JCI patched time-updated IPCW sensitivity rerun\n")
cat("============================================================\n")
cat("Mode: ", mode, "\n", sep = "")
cat("OUT_DIR: ", Sys.getenv("OUT_DIR"), "\n", sep = "")
cat("B: ", Sys.getenv("B"), "\n", sep = "")
cat("N: ", Sys.getenv("N"), "\n", sep = "")
cat("N_CORES: ", Sys.getenv("N_CORES"), "\n", sep = "")
cat("============================================================\n\n")

rscript <- file.path(R.home("bin"), "Rscript")
if (.Platform$OS.type == "windows") rscript <- paste0(rscript, ".exe")
if (!file.exists(rscript)) rscript <- Sys.which("Rscript")
if (!nzchar(rscript)) stop("Could not find Rscript")

run_step <- function(label, script) {
  cat("\n------------------------------------------------------------\n")
  cat("Running: ", label, "\n", sep = "")
  cat("Script: ", script, "\n", sep = "")
  cat("------------------------------------------------------------\n")
  status <- system2(rscript, script, stdout = "", stderr = "")
  if (!identical(status, 0L)) stop("Step failed: ", label, " with exit status ", status)
  cat("Finished: ", label, "\n", sep = "")
}

run_step("patched time-updated IPCW sensitivity only", "scripts/05_rerun_timeupdated_sensitivity_only.R")

# In full mode, refresh the main table pack and standard QA report so the final
# output folder contains consistent supplementary time-updated tables. Smoke mode
# intentionally skips these because primary outputs are absent in the smoke OUT_DIR.
if (mode == "full") {
  run_step("refresh JCI tables", "scripts/02_make_tables_jci.R")
  run_step("refresh JCI QA report", "scripts/04_check_outputs_jci.R")
}

cat("\n============================================================\n")
cat("Patched time-updated sensitivity rerun completed.\n")
cat("Check: ", Sys.getenv("OUT_DIR"), "/audit/timeupdated_sensitivity_patch_qa.txt\n", sep = "")
cat("============================================================\n")
