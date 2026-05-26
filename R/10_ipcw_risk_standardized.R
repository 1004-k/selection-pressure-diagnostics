# R/10_ipcw_risk_standardized.R
# ------------------------------------------------------------
# Baseline-standardized IPCW estimator for the fixed-horizon
# per-protocol risk estimand.
#
# Why this file exists:
# The legacy IPCW risk code estimated a weighted risk contrast within
# observed treatment groups. That is useful descriptively, but it is not the
# cleanest comparator for the DR-AIPW estimator because DR-AIPW standardizes
# over the baseline covariate distribution through a treatment model and an
# outcome regression.
#
# This replacement estimates the same marginal fixed-horizon risk target as
# the DR-AIPW estimator by combining:
#   1) treatment inverse-probability weights, g(A | W), and
#   2) horizon-level adherence/censoring weights, P(U = 1 | A, W).
#
# The default uses a stabilized Hajek form because it is bounded and usually
# more stable in finite samples. Set estimator = "ht" for the unnormalised
# Horvitz-Thompson form that corresponds directly to the first-order IPW
# component of the AIPW estimating equation.
# ------------------------------------------------------------

.clip_probability <- function(p, eps = 1e-6) {
  pmin(pmax(as.numeric(p), eps), 1 - eps)
}

.make_binary_horizon_for_ipcw <- function(pp_dt, t_cut) {
  d <- data.table::copy(pp_dt)
  d[, Y := as.integer(delta_pp == 1L & time_pp <= t_cut)]
  # U = 1 if the fixed-horizon outcome is observed before artificial censoring.
  # This is true for outcome events by t_cut and for individuals still under
  # per-protocol follow-up at t_cut.
  d[, U := as.integer(Y == 1L | time_pp >= t_cut)]
  d[]
}

.get_baseline_W <- function(mis_spec = 0L) {
  W_full <- c("age", "sex", "bmi", "egfr", "util", "gall")
  W_miss <- c("age", "sex", "bmi")
  if (as.integer(mis_spec) == 1L) W_miss else W_full
}

fit_ipcw_risk_standardized <- function(pp_dt,
                                       t_cut,
                                       mis_spec = 0L,
                                       w_floor = 0.05,
                                       ps_floor = 1e-3,
                                       estimator = c("hajek", "ht"),
                                       return_nuisance = FALSE) {
  estimator <- match.arg(estimator)
  t_cut <- as.numeric(t_cut)
  if (!is.finite(t_cut) || t_cut <= 0) stop("t_cut must be positive.")

  d <- .make_binary_horizon_for_ipcw(pp_dt, t_cut = t_cut)
  W <- .get_baseline_W(mis_spec)
  missing_W <- setdiff(W, names(d))
  if (length(missing_W) > 0) {
    stop("Missing baseline covariates in pp_dt: ", paste(missing_W, collapse = ", "))
  }

  # Censoring/adherence model at the fixed horizon.
  # This is deliberately aligned with fit_dr_risk_rescue() in R/09.
  f_c <- stats::as.formula(paste0("U ~ A0 + ", paste(W, collapse = " + ")))
  fit_c <- stats::glm(f_c, data = d, family = stats::binomial())
  p_unc <- .clip_probability(stats::predict(fit_c, type = "response"), eps = 1e-6)
  c_hat <- pmax(p_unc, w_floor)

  # Baseline treatment model for marginal standardization.
  f_g <- stats::as.formula(paste0("A0 ~ ", paste(W, collapse = " + ")))
  fit_g <- stats::glm(f_g, data = d, family = stats::binomial())
  g1 <- .clip_probability(stats::predict(fit_g, type = "response"), eps = ps_floor)

  A <- d$A0
  Y <- d$Y
  U <- d$U

  H1 <- as.numeric(A == 1L) * U / (g1 * c_hat)
  H0 <- as.numeric(A == 0L) * U / ((1 - g1) * c_hat)

  if (estimator == "ht") {
    risk1 <- mean(H1 * Y, na.rm = TRUE)
    risk0 <- mean(H0 * Y, na.rm = TRUE)
    ic1 <- H1 * Y - risk1
    ic0 <- H0 * Y - risk0
    den1 <- 1
    den0 <- 1
  } else {
    den1 <- mean(H1, na.rm = TRUE)
    den0 <- mean(H0, na.rm = TRUE)
    risk1 <- mean(H1 * Y, na.rm = TRUE) / pmax(den1, 1e-12)
    risk0 <- mean(H0 * Y, na.rm = TRUE) / pmax(den0, 1e-12)
    ic1 <- H1 * (Y - risk1) / pmax(den1, 1e-12)
    ic0 <- H0 * (Y - risk0) / pmax(den0, 1e-12)
  }

  rd <- risk1 - risk0
  ic_rd <- ic1 - ic0
  se_rd <- stats::sd(ic_rd, na.rm = TRUE) / sqrt(nrow(d))

  rr <- risk1 / pmax(risk0, 1e-9)
  logrr <- log(rr)
  ic_logrr <- (ic1 / pmax(risk1, 1e-9)) - (ic0 / pmax(risk0, 1e-9))
  se_logrr <- stats::sd(ic_logrr, na.rm = TRUE) / sqrt(nrow(d))

  out <- list(
    t = t_cut,
    method = paste0("ipcw_std_", estimator),
    mis_spec = as.integer(mis_spec),
    risk1 = as.numeric(risk1),
    risk0 = as.numeric(risk0),
    rd = as.numeric(rd),
    se_rd = as.numeric(se_rd),
    logrr = as.numeric(logrr),
    se_logrr = as.numeric(se_logrr),
    den1 = as.numeric(den1),
    den0 = as.numeric(den0),
    mean_c_hat = mean(c_hat, na.rm = TRUE),
    min_c_hat = min(c_hat, na.rm = TRUE),
    min_g1 = min(g1, na.rm = TRUE),
    max_g1 = max(g1, na.rm = TRUE)
  )

  if (isTRUE(return_nuisance)) {
    out$fit_c <- fit_c
    out$fit_g <- fit_g
    out$c_hat <- c_hat
    out$g1 <- g1
    out$H1 <- H1
    out$H0 <- H0
  }

  out
}
