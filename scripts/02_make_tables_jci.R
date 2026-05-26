#!/usr/bin/env Rscript
# JCI table builder for the standardized-IPCW rerun.

source("scripts/00_utils.R")

cfg <- init_project(n_cores = 1L, seed = 2026L, out_dir = Sys.getenv("OUT_DIR", "output_jci_final"))
cfg$require_pkgs(c("data.table"))

out_tab <- file.path(cfg$out_dir, "tables")
dir.create(out_tab, showWarnings = FALSE, recursive = TRUE)

grid <- make_scenario_grid()
mis <- strsplit(Sys.getenv("MIS_SPEC_LEVELS", "0,1"), ",", fixed = TRUE)[[1]]
mis <- as.integer(trimws(mis)); mis <- mis[mis %in% c(0L,1L)]
if (length(mis) == 0) mis <- c(0L,1L)

scn <- data.table::CJ(scenario_id = grid$scenario_id, mis_spec = mis)
scn <- merge(scn, grid, by = "scenario_id", all.x = TRUE)
data.table::setorder(scn, scenario_num, mis_spec)

tab1 <- scn[, .(scenario_id, nuisance_regime = mis_spec, axisA, axisB, rho_meas, panel_title)]
data.table::fwrite(tab1, file.path(out_tab, "Table1_scenario_grid.csv"))

ipcw_file <- file.path(cfg$out_dir, "perf_summary_rescue_ipcw_risk.csv")
dr_file   <- file.path(cfg$out_dir, "perf_summary_rescue_dr.csv")
hr_file   <- file.path(cfg$out_dir, "perf_summary_rescue_hr.csv")
truth_file <- file.path(cfg$out_dir, "raw", "mc_truth_risk.csv")

if (file.exists(ipcw_file)) data.table::fwrite(data.table::fread(ipcw_file), file.path(out_tab, "OnlineTable_ipcw_standardized_performance.csv"))
if (file.exists(dr_file)) data.table::fwrite(data.table::fread(dr_file), file.path(out_tab, "OnlineTable_dr_aipw_performance.csv"))
if (file.exists(hr_file)) data.table::fwrite(data.table::fread(hr_file), file.path(out_tab, "OnlineTable_legacy_ipcw_cox_hr_performance.csv"))
if (file.exists(truth_file)) data.table::fwrite(data.table::fread(truth_file), file.path(out_tab, "OnlineTable_mc_truth_risk.csv"))

stopifnot(file.exists(ipcw_file), file.exists(dr_file))
ipcw <- data.table::fread(ipcw_file)
dr   <- data.table::fread(dr_file)
ipcw[, truth := as.character(truth)]
dr[, truth := as.character(truth)]

dr_glm <- dr[method_Q == "glm"]
if (nrow(dr_glm) == 0) dr_glm <- dr

# Main compact operating characteristics table.
main_ipcw <- ipcw[, .(
  TypeI_IPCW = median(ci_excl0[truth == "null"], na.rm = TRUE),
  Coverage_IPCW = median(cover[truth == "null"], na.rm = TRUE),
  Power_IPCW = median(ci_excl0[truth == "non_null"], na.rm = TRUE),
  SignError_IPCW = median(sign_error[truth == "non_null"], na.rm = TRUE)
), by = .(nuisance_regime = mis_spec)]

main_dr <- dr_glm[, .(
  TypeI_DR = median(ci_excl0[truth == "null"], na.rm = TRUE),
  Coverage_DR = median(cover[truth == "null"], na.rm = TRUE),
  Power_DR = median(ci_excl0[truth == "non_null"], na.rm = TRUE),
  SignError_DR = median(sign_error[truth == "non_null"], na.rm = TRUE)
), by = .(nuisance_regime = mis_spec)]

tab2 <- merge(main_ipcw, main_dr, by = "nuisance_regime", all = TRUE)
data.table::setorder(tab2, nuisance_regime)
data.table::fwrite(tab2, file.path(out_tab, "Table2_operating_characteristics.csv"))

# Scenario-level comparison for method-choice map.
nonnull_ip <- ipcw[truth == "non_null", .(
  scenario_id, mis_spec, axisA, axisB, rho_meas,
  power_ipcw = ci_excl0,
  cover_ipcw = cover,
  bias_ipcw = bias,
  rmse_ipcw = rmse,
  sign_error_ipcw = sign_error,
  logrr_median_ipcw = logrr_median,
  rd_median_ipcw = rd_median
)]
nonnull_dr <- dr_glm[truth == "non_null", .(
  scenario_id, mis_spec, axisA, axisB, rho_meas,
  power_dr = ci_excl0,
  cover_dr = cover,
  bias_dr = bias,
  rmse_dr = rmse,
  sign_error_dr = sign_error,
  logrr_median_dr = logrr_median,
  rd_median_dr = rd_median
)]

choice <- merge(nonnull_ip, nonnull_dr,
                by = c("scenario_id", "mis_spec", "axisA", "axisB", "rho_meas"), all = FALSE)
choice[, `:=`(
  delta_power = power_dr - power_ipcw,
  delta_sign_error = sign_error_dr - sign_error_ipcw,
  delta_rmse = rmse_ipcw - rmse_dr
)]
choice[, method_region := data.table::fifelse(
  delta_power <= -0.05 | delta_sign_error >= 0.05, "No rescue",
  data.table::fifelse(delta_power < 0 | delta_sign_error > 0, "DR with caution", "DR default")
)]
data.table::setorder(choice, mis_spec, scenario_id)
data.table::fwrite(choice, file.path(out_tab, "Table3_method_choice_regions.csv"))

# Threshold sensitivity for the operating-map heuristic. This is important for
# JCI because the categories are a decision-support summary, not a universal rule.
thr_grid <- data.table::data.table(threshold = c(0.02, 0.05, 0.10))
thr_sens <- thr_grid[, {
  tt <- threshold
  tmp <- data.table::copy(choice[mis_spec == 1])
  tmp[, region := data.table::fifelse(
    delta_power <= -tt | delta_sign_error >= tt, "No rescue",
    data.table::fifelse(delta_power < 0 | delta_sign_error > 0, "DR with caution", "DR default")
  )]
  tmp[, .N, by = region]
}, by = threshold]
thr_sens <- data.table::dcast(thr_sens, threshold ~ region, value.var = "N", fill = 0)
data.table::fwrite(thr_sens, file.path(out_tab, "OnlineTable_threshold_sensitivity.csv"))

message("Saved JCI tables to: ", out_tab)

# Supplementary time-updated IPCW sensitivity table.
# This is generated only when RUN_TIME_UPDATED_SENS=1 and the corresponding
# full-run output exists. It should be reported as a sensitivity analysis, not
# used to define the main operating-map regions.
tu_file <- file.path(cfg$out_dir, "perf_summary_rescue_ipcw_timeupdated.csv")
if (file.exists(tu_file)) {
  tu <- data.table::fread(tu_file)
  data.table::fwrite(tu, file.path(out_tab, "OnlineTable_ipcw_timeupdated_sensitivity_performance.csv"))

  tu_main <- tu[, .(
    TypeI_timeupdated_IPCW = median(ci_excl0[truth == "null"], na.rm = TRUE),
    Coverage_timeupdated_IPCW = median(cover[truth == "null"], na.rm = TRUE),
    Power_timeupdated_IPCW = median(ci_excl0[truth == "non_null"], na.rm = TRUE),
    SignError_timeupdated_IPCW = median(sign_error[truth == "non_null"], na.rm = TRUE),
    Median_min_c_hat = median(min_c_hat_median, na.rm = TRUE)
  ), by = .(nuisance_regime = mis_spec)]
  data.table::setorder(tu_main, nuisance_regime)
  data.table::fwrite(tu_main, file.path(out_tab, "OnlineTable_timeupdated_sensitivity_compact.csv"))

  if (exists("ipcw") && nrow(ipcw) > 0) {
    base_vs_tu <- merge(
      ipcw[truth == "non_null", .(scenario_id, mis_spec, axisA, axisB, rho_meas,
                                  power_baseline = ci_excl0,
                                  sign_error_baseline = sign_error,
                                  rmse_baseline = rmse)],
      tu[truth == "non_null", .(scenario_id, mis_spec,
                                power_timeupdated = ci_excl0,
                                sign_error_timeupdated = sign_error,
                                rmse_timeupdated = rmse)],
      by = c("scenario_id", "mis_spec"), all = FALSE
    )
    base_vs_tu[, `:=`(
      delta_power_timeupdated_minus_baseline = power_timeupdated - power_baseline,
      delta_sign_error_timeupdated_minus_baseline = sign_error_timeupdated - sign_error_baseline,
      delta_rmse_baseline_minus_timeupdated = rmse_baseline - rmse_timeupdated
    )]
    data.table::fwrite(base_vs_tu, file.path(out_tab, "OnlineTable_baseline_vs_timeupdated_ipcw_sensitivity.csv"))
  }
}
