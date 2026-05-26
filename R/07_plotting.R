# R/07_plotting.R
# ------------------------------------------------------------
# Base-R plotting helpers (black & white friendly)
# Goal: produce figures that match the manuscript PDFs (fonts/margins/lwd).
# ------------------------------------------------------------

# ---- constants (tuned to match the provided manuscript PDFs) ----
SPD_BAND_COL   <- grDevices::adjustcolor("grey80", alpha.f = 0.85)
SPD_ZERO_COL   <- "grey70"
THRESH_COL     <- "grey70"
MED_LWD        <- 3.0
BAND_BORDER    <- NA
AXIS_LWD       <- 1.0
ZERO_LWD       <- 1.6
THRESH_LWD     <- 1.6
THRESH_LTY     <- 3

# ---- summarizers ----
summarize_band <- function(dt, x_col, y_col, by_cols) {
  stopifnot(x_col %in% names(dt), y_col %in% names(dt))
  data.table::as.data.table(dt)[, .(
    med = stats::median(get(y_col), na.rm = TRUE),
    q25 = stats::quantile(get(y_col), 0.25, na.rm = TRUE, type = 7),
    q75 = stats::quantile(get(y_col), 0.75, na.rm = TRUE, type = 7)
  ), by = c(by_cols, x_col)]
}

# ---- device helpers ----
open_pdf <- function(file, width, height, pointsize = 12, family = "Helvetica") {
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  grDevices::pdf(file, width = width, height = height, pointsize = pointsize, family = family)
}

close_device <- function() {
  grDevices::dev.off()
}

# ---- low-level panel drawer ----
draw_band_panel <- function(x, med, q25, q75,
                            xlim, ylim,
                            xlab = "", ylab = "",
                            main = "",
                            show_y_labels = TRUE,
                            show_x_labels = TRUE,
                            x_at = NULL, x_labels = NULL,
                            y_at = NULL, y_labels = NULL,
                            add_zero = FALSE,
                            zero_y = 0,
                            add_thresh = NULL) {

  # frame
  graphics::plot(x, med, type = "n", xlim = xlim, ylim = ylim,
                 xlab = xlab, ylab = ylab, main = main, axes = FALSE)

  # ribbon
  graphics::polygon(c(x, rev(x)), c(q25, rev(q75)),
                    col = SPD_BAND_COL, border = BAND_BORDER)

  # reference lines
  if (isTRUE(add_zero)) {
    graphics::abline(h = zero_y, col = SPD_ZERO_COL, lwd = ZERO_LWD)
  }
  if (!is.null(add_thresh)) {
    for (thr in add_thresh) {
      graphics::abline(h = thr, col = THRESH_COL, lwd = THRESH_LWD, lty = THRESH_LTY)
    }
  }

  # median
  graphics::lines(x, med, lwd = MED_LWD)

  # axes
  if (is.null(x_at)) x_at <- pretty(x)
  if (is.null(x_labels)) x_labels <- x_at
  if (show_x_labels) {
    graphics::axis(1, at = x_at, labels = x_labels, lwd = AXIS_LWD)
  } else {
    graphics::axis(1, at = x_at, labels = FALSE, lwd = AXIS_LWD)
  }

  if (is.null(y_at)) y_at <- pretty(ylim)
  if (is.null(y_labels)) y_labels <- y_at
  if (show_y_labels) {
    graphics::axis(2, at = y_at, labels = y_labels, las = 1, lwd = AXIS_LWD)
  } else {
    graphics::axis(2, at = y_at, labels = FALSE, lwd = AXIS_LWD)
  }

  graphics::box(lwd = AXIS_LWD)
  invisible(TRUE)
}

# ---- figure-specific par settings (tuned to match PDFs) ----
par_tripanel <- function() {
  graphics::par(
    mfrow = c(1,3),
    oma = c(0.6, 0.6, 4.2, 0.4),
    mar = c(4.2, 4.4, 2.6, 1.0),
    mgp = c(2.3, 0.8, 0),
    tcl = -0.3,
    cex.axis = 1.15,
    cex.lab  = 1.25,
    cex.main = 1.25
  )
}

par_tripanel_two_rows <- function() {
  # default; individual panels will adjust mar for top/bottom rows
  graphics::par(
    mfrow = c(2,3),
    oma = c(0.6, 0.6, 4.8, 0.4),
    mar = c(3.2, 4.4, 2.6, 1.0),
    mgp = c(2.3, 0.8, 0),
    tcl = -0.3,
    cex.axis = 1.15,
    cex.lab  = 1.25,
    cex.main = 1.25
  )
}

