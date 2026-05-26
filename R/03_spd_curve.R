# R/03_spd_curve.R
# ------------------------------------------------------------
# SPD(t) estimation (piecewise slope per time interval)
# SPD(t) is the log-HR for deviation per 1 SD higher prognostic score at time t.
#
# Patch v3:
# - If pooled piecewise Cox fails, fall back to one Cox model per interval.
# - Cumulative pressure uses 0 for missing interval estimates (but returns NA if all are missing).
# ------------------------------------------------------------

standardize_by_interval <- function(dt, z_col = "z_obs", interval_col = "interval") {
  dt[, z_std := {
    zz <- get(z_col)
    m <- mean(zz, na.rm = TRUE)
    s <- stats::sd(zz, na.rm = TRUE)
    if (!is.finite(s) || s == 0) rep(0, .N) else (zz - m) / s
  }, by = interval_col]
  dt[]
}

.assemble_out <- function(intervals, breaks) {
  out <- data.table::data.table(
    interval = as.integer(intervals),
    gamma_hat = NA_real_,
    se = NA_real_
  )
  out[, HR_1SD := NA_real_]
  out[, `:=`(CI_low = NA_real_, CI_high = NA_real_)]
  out <- out[order(interval)]
  out[, tstart := breaks[pmax(interval, 1L)]]
  out[, tstop  := breaks[pmin(interval + 1L, length(breaks))]]
  out[, t_mid  := (tstart + tstop) / 2]
  out
}

.fit_one_interval <- function(dt_k, dev_col = "dev", id_col = "id", robust = TRUE) {
  if (sum(dt_k[[dev_col]] == 1L, na.rm = TRUE) == 0L) {
    return(list(beta = NA_real_, se = NA_real_))
  }

  fml_k <- stats::as.formula(paste0(
    "survival::Surv(tstart, tstop, ", dev_col, ") ~ z_std",
    if (robust) paste0(" + cluster(", id_col, ")") else ""
  ))

  fit_k <- tryCatch(survival::coxph(fml_k, data = dt_k), error = function(e) e)
  if (inherits(fit_k, "error")) return(list(beta = NA_real_, se = NA_real_))

  b  <- tryCatch(as.numeric(stats::coef(fit_k)[1]), error = function(e) NA_real_)
  se <- tryCatch(sqrt(diag(stats::vcov(fit_k)))[1], error = function(e) NA_real_)
  list(beta = b, se = se)
}

estimate_spd_piecewise <- function(long_dt,
                                   breaks,
                                   z_col = "z_obs",
                                   dev_col = "dev",
                                   id_col = "id",
                                   robust = TRUE) {
  dt <- data.table::copy(long_dt)
  dt[, tstart := as.numeric(tstart)]
  dt[, tstop  := as.numeric(tstop)]

  eps_dt <- 1e-6
  dt <- dt[is.finite(tstart) & is.finite(tstop) & (tstop - tstart) > eps_dt]

  dt[, interval := findInterval(tstart, vec = breaks, rightmost.closed = TRUE)]
  dt[interval < 1, interval := 1L]
  dt[, interval := as.integer(interval)]
  dt[, interval_f := factor(interval)]

  dt <- standardize_by_interval(dt, z_col = z_col, interval_col = "interval")

  intervals <- sort(unique(dt$interval))
  if (length(intervals) == 0) {
    out <- .assemble_out(integer(0), breaks)
    return(list(fit = NULL, spd = out, dt_used = dt))
  }

  # pooled piecewise model
  fml <- stats::as.formula(paste0(
    "survival::Surv(tstart, tstop, ", dev_col, ") ~ 0 + z_std:interval_f + strata(interval_f)",
    if (robust) paste0(" + cluster(", id_col, ")") else ""
  ))

  fit <- tryCatch(survival::coxph(fml, data = dt), error = function(e) e)

  if (!inherits(fit, "error")) {
    coefs <- stats::coef(fit)
    vc    <- stats::vcov(fit)
    nm    <- names(coefs)

    k  <- as.integer(gsub(".*interval_f", "", nm))
    se <- sqrt(diag(vc))

    out <- data.table::data.table(
      interval = k,
      gamma_hat = as.numeric(coefs),
      se = as.numeric(se)
    )
    out[, HR_1SD := exp(gamma_hat)]
    z <- stats::qnorm(0.975)
    out[, `:=`(
      CI_low = exp(gamma_hat - z * se),
      CI_high= exp(gamma_hat + z * se)
    )]

    out <- out[order(interval)]
    out[, tstart := breaks[pmax(interval, 1L)]]
    out[, tstop  := breaks[pmin(interval + 1L, length(breaks))]]
    out[, t_mid  := (tstart + tstop) / 2]

    return(list(fit = fit, spd = out, dt_used = dt))
  }

  # fallback: fit per interval
  out <- .assemble_out(intervals, breaks)
  for (kk in intervals) {
    dt_k <- dt[interval == kk]
    est <- .fit_one_interval(dt_k, dev_col = dev_col, id_col = id_col, robust = robust)
    out[interval == kk, `:=`(gamma_hat = est$beta, se = est$se)]
  }

  out[, HR_1SD := ifelse(is.na(gamma_hat), NA_real_, exp(gamma_hat))]
  z <- stats::qnorm(0.975)
  out[, `:=`(
    CI_low  = ifelse(is.na(gamma_hat) | is.na(se), NA_real_, exp(gamma_hat - z * se)),
    CI_high = ifelse(is.na(gamma_hat) | is.na(se), NA_real_, exp(gamma_hat + z * se))
  )]

  list(fit = NULL, spd = out, dt_used = dt)
}

compute_cum_pressure <- function(spd_dt) {
  dt <- data.table::copy(spd_dt)
  dt[, dt_len := (tstop - tstart)]
  g <- ifelse(is.na(dt$gamma_hat), 0, dt$gamma_hat)
  Gamma <- cumsum(g * dt$dt_len)
  if (all(is.na(dt$gamma_hat))) Gamma <- rep(NA_real_, length(Gamma))
  dt[, Gamma_hat := Gamma]
  dt[]
}
