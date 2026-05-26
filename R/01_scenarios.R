# R/01_scenarios.R
# ------------------------------------------------------------
# Scenario grid for the 3-axis simulation design (18 scenarios)
# Axis A: SPD(t) path = flat / increasing / late-surge
# Axis B: nonlinear selection = none / threshold-jump
# Axis C: imperfect prognosis = Corr(z_obs, z_true) = 0.4 / 0.7 / 1.0
#
# We assign stable IDs S01–S18 in a Table-1-friendly order:
#   Axis A (flat → increasing → late-surge) ×
#   Axis B (none → threshold-jump) ×
#   rho (0.4 → 0.7 → 1.0)
# ------------------------------------------------------------

make_scenario_grid <- function() {
  cfgA <- data.table::data.table(axisA = c("flat", "increasing", "late_surge"))
  cfgB <- data.table::data.table(axisB = c("none", "threshold_jump"))
  cfgC <- data.table::data.table(rho_meas = c(0.4, 0.7, 1.0))

  grid <- data.table::CJ(
    axisA = cfgA$axisA,
    axisB = cfgB$axisB,
    rho_meas = cfgC$rho_meas,
    unique = TRUE
  )

  # deterministic ordering
  grid[, axisA_ord := factor(axisA, levels = c("flat","increasing","late_surge"))]
  grid[, axisB_ord := factor(axisB, levels = c("none","threshold_jump"))]
  data.table::setorder(grid, axisA_ord, axisB_ord, rho_meas)
  grid[, `:=`(axisA_ord = NULL, axisB_ord = NULL)]

  grid[, scenario_num := .I]
  grid[, scenario_id  := sprintf("S%02d", scenario_num)]

  # pretty labels for printing/plotting
  grid[, axisA_pretty := ifelse(axisA == "late_surge", "late-surge", axisA)]
  grid[, axisB_pretty := ifelse(axisB == "threshold_jump", "threshold-jump", axisB)]
  grid[, panel_title  := sprintf("A: %s; B: %s; ρ=%.1f", axisA_pretty, axisB_pretty, rho_meas)]

  grid[]
}

# Helper: map legacy IDs (A_flat__B_none__C_rho1.0) to Sxx
# This keeps backwards-compatibility if users have older output files.
normalize_scenario_id <- function(dt, grid = NULL, id_col = "scenario_id") {
  if (is.null(grid)) grid <- make_scenario_grid()
  if (!id_col %in% names(dt)) return(dt)

  x <- dt[[id_col]]

  # if already Sxx, keep
  is_s <- grepl("^S\\d\\d$", x)

  # legacy parse: A_<axisA>__B_<axisB>__C_rho<rho>
  axisA <- ifelse(is_s, NA_character_, sub("^A_([^_]+).*", "\\1", x))
  axisB <- ifelse(is_s, NA_character_, sub("^A_[^_]+__B_([^_]+).*", "\\1", x))
  rho   <- ifelse(is_s, NA_character_, sub(".*__C_rho([0-9.]+)$", "\\1", x))

  tmp <- data.table::data.table(axisA = axisA, axisB = axisB, rho_meas = as.numeric(rho))
  tmp[, old := x]

  map <- grid[, .(axisA, axisB, rho_meas, new = scenario_id)]
  tmp <- merge(tmp, map, by = c("axisA","axisB","rho_meas"), all.x = TRUE)

  x2 <- ifelse(is_s, x, ifelse(is.na(tmp$new), x, tmp$new))
  dt[[id_col]] <- x2
  dt[]
}


# Table 1 helper (CSV-friendly)
scenario_table <- function(grid) {
  g <- data.table::copy(grid)
  data.table::setorder(g, scenario_num)
  g[, .(
    scenario_id,
    axisA,
    axisB,
    rho_meas,
    panel_title
  )]
}
