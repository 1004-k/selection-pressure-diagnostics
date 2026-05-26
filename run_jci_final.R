#!/usr/bin/env Rscript

# run_jci_final.R
# ------------------------------------------------------------
# R-only runner for the JCI final rerun package.
#
# Usage from a terminal:
#   Rscript run_jci_final.R smoke
#   Rscript run_jci_final.R full
#
# Usage from the RStudio Console:
#   system2("Rscript", c("run_jci_final.R", "smoke"))
#   system2("Rscript", c("run_jci_final.R", "full"))
#
# Run this file from the project root, where R/ and scripts/ exist.
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) >= 1) tolower(args[1]) else "smoke"

if (!mode %in% c("smoke", "full")) {
  stop("mode must be either 'smoke' or 'full'. Example: Rscript run_jci_final.R smoke")
}

project_root <- getwd()

required_files <- c(
  "scripts/01_run_main_simulation_standardized.R",
  "scripts/02_make_tables_jci.R",
  "scripts/03_make_figures_jci.R",
  "scripts/04_check_outputs_jci.R",
  "R/02_dgp_simulate.R",
  "R/10_ipcw_risk_standardized.R",
  "R/11_mc_truth_intervention.R",
  "R/12_ipcw_timeupdated_sensitivity.R"
)

missing_files <- required_files[!file.exists(file.path(project_root, required_files))]

if (length(missing_files) > 0) {
  stop(
    "Required files are missing. Are you running this from the project root?\nMissing:\n",
    paste(missing_files, collapse = "\n")
  )
}

set_default_env <- function(name, value) {
  current <- Sys.getenv(name, unset = NA_character_)
  if (is.na(current) || !nzchar(current)) {
    do.call(Sys.setenv, stats::setNames(as.list(value), name))
  }
}

if (mode == "smoke") {
  # Fast smoke test. Do not use these outputs for the manuscript.
  set_default_env("OUT_DIR", "output_jci_smoke_final")
  set_default_env("B", "3")
  set_default_env("N", "300")
  set_default_env("TRUTH_LEVELS", "null,non_null")
  set_default_env("MIS_SPEC_LEVELS", "0,1")
  set_default_env("BETA_TRUE_PP", "-0.22314355131420976")
  set_default_env("MC_TRUTH", "1")
  set_default_env("N_TRUTH", "5000")
  set_default_env("B_TRUTH", "1")
  set_default_env("RUN_DR", "1")
  set_default_env("RUN_ML", "0")
  set_default_env("RUN_TMLE", "0")
  set_default_env("RUN_TIME_UPDATED_SENS", "1")
  set_default_env("IPCW_ESTIMATOR", "hajek")
  set_default_env("N_CORES", "2")
}

if (mode == "full") {
  # Full JCI rerun. Adjust N_CORES if your machine or GCP instance allows more.
  set_default_env("OUT_DIR", "output_jci_final")
  set_default_env("B", "200")
  set_default_env("N", "2000")
  set_default_env("TRUTH_LEVELS", "null,non_null")
  set_default_env("MIS_SPEC_LEVELS", "0,1")
  set_default_env("BETA_TRUE_PP", "-0.22314355131420976")
  set_default_env("MC_TRUTH", "1")
  set_default_env("N_TRUTH", "200000")
  set_default_env("B_TRUTH", "5")
  set_default_env("RUN_DR", "1")
  set_default_env("RUN_ML", "0")
  set_default_env("RUN_TMLE", "0")
  set_default_env("RUN_TIME_UPDATED_SENS", "1")
  set_default_env("IPCW_ESTIMATOR", "hajek")
  set_default_env("N_CORES", "3")
}

cat("\n")
cat("============================================================\n")
cat("JCI final rerun: standardized IPCW + intervention truth\n")
cat("============================================================\n")
cat("Mode: ", mode, "\n", sep = "")
cat("Project root: ", project_root, "\n", sep = "")
for (nm in c("OUT_DIR", "B", "N", "TRUTH_LEVELS", "MIS_SPEC_LEVELS", "BETA_TRUE_PP",
             "MC_TRUTH", "N_TRUTH", "B_TRUTH", "RUN_DR", "RUN_ML", "RUN_TMLE",
             "RUN_TIME_UPDATED_SENS", "IPCW_ESTIMATOR", "N_CORES")) {
  cat(nm, ": ", Sys.getenv(nm), "\n", sep = "")
}
cat("============================================================\n\n")

rscript <- file.path(R.home("bin"), "Rscript")
if (.Platform$OS.type == "windows") rscript <- paste0(rscript, ".exe")
if (!file.exists(rscript)) rscript <- Sys.which("Rscript")
if (!nzchar(rscript)) stop("Could not find Rscript.")

run_step <- function(label, script, extra_args = character()) {
  cat("\n")
  cat("------------------------------------------------------------\n")
  cat("Running: ", label, "\n", sep = "")
  cat("Script: ", script, "\n", sep = "")
  cat("------------------------------------------------------------\n")

  status <- system2(
    command = rscript,
    args = c(script, extra_args),
    stdout = "",
    stderr = ""
  )

  if (!identical(status, 0L)) {
    stop("Step failed: ", label, " with exit status ", status)
  }

  cat("\nFinished: ", label, "\n", sep = "")
}

run_step("main simulation with standardized IPCW, intervention truth, and time-updated sensitivity",
         "scripts/01_run_main_simulation_standardized.R")
run_step("JCI tables", "scripts/02_make_tables_jci.R")
run_step("JCI figures", "scripts/03_make_figures_jci.R", extra_args = Sys.getenv("OUT_DIR"))
run_step("JCI output QA checks", "scripts/04_check_outputs_jci.R")

cat("\n")
cat("============================================================\n")
cat("All steps completed successfully.\n")
cat("Output folder: ", Sys.getenv("OUT_DIR"), "\n", sep = "")
cat("Use these outputs only if the QA report is clean.\n")
cat("============================================================\n")
