#!/usr/bin/env Rscript
# JCI figure builder for the standardized-IPCW rerun.

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else Sys.getenv("OUT_DIR", "output_jci_final")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

ipcw_file <- file.path(out_dir, "perf_summary_rescue_ipcw_risk.csv")
dr_file   <- file.path(out_dir, "perf_summary_rescue_dr.csv")
spd_file  <- file.path(out_dir, "raw", "spd_curves_rescue.csv")

stopifnot(file.exists(ipcw_file), file.exists(dr_file), file.exists(spd_file))

ipcw <- fread(ipcw_file)
dr   <- fread(dr_file)
spd  <- fread(spd_file)

ipcw[, truth := fifelse(is.na(truth), "null", as.character(truth))]
dr[,   truth := fifelse(is.na(truth), "null", as.character(truth))]
spd[,  truth := fifelse(is.na(truth), "null", as.character(truth))]

# Pressure summary under the causal null: median over replicates of max_t |SPD(t)|.
spd_null <- spd[truth == "null"]
spd_null[, abs_gamma := abs(gamma_hat)]
max_by_rep <- spd_null[, .(max_abs_gamma = suppressWarnings(max(abs_gamma, na.rm = TRUE))),
                       by = .(scenario_id, replicate)]
max_by_rep[!is.finite(max_abs_gamma), max_abs_gamma := NA_real_]
pressure <- max_by_rep[, .(pressure = suppressWarnings(median(max_abs_gamma, na.rm = TRUE))),
                       by = .(scenario_id)]
pressure[!is.finite(pressure), pressure := NA_real_]

# Prefer the parametric GLM DR-AIPW row for the main JCI comparison.
dr_glm <- dr[method_Q == "glm"]
if (nrow(dr_glm) == 0) dr_glm <- dr

mk_perf <- function(dt, method_label) {
  dt[, .(scenario_id, mis_spec, truth, rho_meas, ci_excl0, sign_error)][, method := method_label]
}
perf <- rbindlist(list(mk_perf(ipcw, "standardized IPCW"), mk_perf(dr_glm, "DR-AIPW")),
                  use.names = TRUE, fill = TRUE)
perf <- merge(perf, pressure, by = "scenario_id", all.x = TRUE)
perf <- perf[truth %in% c("null", "non_null")]
perf[, metric := fifelse(truth == "null", "Type I error under null", "Power under non-null")]
perf[, metric := factor(metric, levels = c("Type I error under null", "Power under non-null"))]
perf[, mis_spec_f := factor(mis_spec, levels = c(0, 1), labels = c("Richer baseline nuisance set", "Reduced baseline nuisance set"))]
perf[, method := factor(method, levels = c("standardized IPCW", "DR-AIPW"))]
perf[, x_plot := pressure + fifelse(method == "DR-AIPW", -0.02, 0.02)]

hline <- unique(perf[metric == "Type I error under null", .(metric)])
hline[, y := 0.05]

p1 <- ggplot(perf, aes(x = x_plot, y = ci_excl0, shape = method)) +
  geom_point(size = 2.4, alpha = 0.9) +
  geom_hline(data = hline, aes(yintercept = y), linetype = "dashed") +
  facet_grid(metric ~ mis_spec_f, scales = "free_y") +
  scale_shape_manual(values = c("standardized IPCW" = 2, "DR-AIPW" = 16)) +
  labs(
    x = "Selection pressure under null: median of max |SPD(t)|",
    y = NULL,
    shape = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(fig_dir, "Figure1_pressure_type1_power_JCI.pdf"),
       p1, width = 6.85, height = 4.9, units = "in")
ggsave(file.path(fig_dir, "Figure1_pressure_type1_power_JCI.png"),
       p1, width = 6.85, height = 4.9, units = "in", dpi = 300)

# Fig 2: operating map under non-null.
nn_ip <- ipcw[truth == "non_null", .(scenario_id, mis_spec,
                                      power_ipcw = ci_excl0,
                                      sign_ipcw = sign_error,
                                      rho_meas)]
nn_dr <- dr_glm[truth == "non_null", .(scenario_id, mis_spec,
                                        power_dr = ci_excl0,
                                        sign_dr = sign_error)]

m <- merge(nn_ip, nn_dr, by = c("scenario_id", "mis_spec"), all = FALSE)
m <- merge(m, pressure, by = "scenario_id", all.x = TRUE)
m[, delta_power := power_dr - power_ipcw]
m[, delta_sign  := sign_dr  - sign_ipcw]
m[, region := fifelse(delta_power <= -0.05 | delta_sign >= 0.05, "No rescue",
               fifelse(delta_power < 0 | delta_sign > 0, "DR with caution", "DR default"))]
m[, mis_spec_f := factor(mis_spec, levels = c(0, 1), labels = c("Richer baseline nuisance set", "Reduced baseline nuisance set"))]
m[, region := factor(region, levels = c("DR default", "DR with caution", "No rescue"))]

# Deterministic offsets, not random jitter, to reduce overplotting.
m[, x_base := pressure + fifelse(region == "DR with caution", -0.012,
                          fifelse(region == "No rescue", 0.012, 0.0))]
m[, x_plot := x_base]
m[region != "DR default", `:=`(rk = seq_len(.N), n = .N), by = .(mis_spec, rho_meas, region)]
m[region != "DR default", x_plot := x_base + (rk - (n + 1) / 2) * 0.02]
m[, c("rk", "n") := NULL]

fwrite(m[order(mis_spec, pressure)], file.path(fig_dir, "Figure2_operating_map_source_data.csv"))

p2 <- ggplot(m, aes(x = x_plot, y = rho_meas)) +
  geom_point(
    data = m[region == "DR default"],
    aes(shape = region, fill = region),
    size = 2.2, alpha = 0.25, stroke = 0.6, color = "black"
  ) +
  geom_point(
    data = m[region != "DR default"],
    aes(shape = region, fill = region),
    size = 3.4, alpha = 0.95, stroke = 0.6, color = "black"
  ) +
  facet_wrap(~ mis_spec_f, nrow = 1) +
  scale_shape_manual(values = c("DR default" = 21, "DR with caution" = 24, "No rescue" = 22)) +
  scale_fill_manual(values = c("DR default" = "white", "DR with caution" = "black", "No rescue" = "black")) +
  scale_y_continuous(breaks = c(0.4, 0.7, 1.0)) +
  coord_cartesian(ylim = c(0.38, 1.08), clip = "off") +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.06))) +
  labs(
    x = "Selection pressure under null: median of max |SPD(t)|",
    y = "Prognostic information quality (rho)",
    shape = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.margin = margin(5.5, 18, 5.5, 5.5)
  ) +
  guides(fill = "none")

if (requireNamespace("ggrepel", quietly = TRUE)) {
  p2 <- p2 +
    ggrepel::geom_text_repel(
      data = m[region != "DR default"],
      aes(label = scenario_id),
      size = 3,
      seed = 1,
      max.overlaps = Inf,
      force = 3,
      box.padding = 0.25,
      point.padding = 0.18,
      min.segment.length = 0,
      segment.size = 0.2
    )
} else {
  p2 <- p2 +
    geom_text(
      data = m[region != "DR default"],
      aes(label = scenario_id),
      vjust = -0.8, size = 3,
      check_overlap = FALSE
    )
}

ggsave(file.path(fig_dir, "Figure2_operating_map_JCI.pdf"),
       p2, width = 6.85, height = 3.6, units = "in")
ggsave(file.path(fig_dir, "Figure2_operating_map_JCI.png"),
       p2, width = 6.85, height = 3.6, units = "in", dpi = 300)

cat("Saved JCI figures to: ", fig_dir, "\n")
