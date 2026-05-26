# R/08_generality.R
# ------------------------------------------------------------
# Minimal "generality check" estimators for Supplementary Material:
#  1) CCW-IPCW (clone-censor-weight) Cox model for per-protocol HR
#     - Uses two baseline regimes (A0=0 vs A0=1) via cloning.
#     - Censors clones immediately at baseline if they did not initiate the regime.
#     - Uses IPCW for deviation/censoring in the clone dataset and robust SE via cluster(id).
#
#  2) TMLE-like doubly-robust (AIPW) estimator for risk at a fixed horizon t_cut
#     - Uses baseline models for treatment assignment and outcome regression
#     - Uses a simple censoring weight for deviation before t_cut
#     - If `tmle` package is available, it can be used; otherwise falls back to AIPW.
#
# These checks are intentionally light-touch: they are NOT a full bake-off.
# ------------------------------------------------------------

ess_kish <- function(w) (sum(w)^2) / sum(w^2)

# ---- CCW (clone-censor-weight) ----
make_ccw_clones <- function(pp_dt, eps_time = 1e-6) {
  d <- data.table::copy(pp_dt)

  # Two clones per subject: regime 0 and 1
  c0 <- data.table::copy(d); c0[, regime := 0L]
  c1 <- data.table::copy(d); c1[, regime := 1L]
  ccw <- data.table::rbindlist(list(c0, c1), use.names = TRUE, fill = TRUE)

  # Baseline mismatch => immediate deviation (censor at ~0)
  ccw[, follow := as.integer(A0 == regime)]
  ccw[follow == 0L, `:=`(
    time_pp  = eps_time,
    delta_pp = 0L,
    dev_ind  = 1L
  )]

  ccw[]
}

fit_ccw_ipcw <- function(pp_dt,
                         cens_formula = survival::Surv(time_pp, dev_ind) ~ regime + age + sex + bmi + egfr + util + gall,
                         outcome_formula = survival::Surv(time_pp, delta_pp) ~ regime,
                         w_floor = 1e-3) {
  ccw <- make_ccw_clones(pp_dt)

  # Censor model in clone data
  fit_cens <- survival::coxph(cens_formula, data = ccw)
  lp <- stats::predict(fit_cens, type = "lp")

  S <- get_S_from_cox(fit_cens, time_vec = ccw$time_pp, lp_vec = lp)
  w <- 1 / pmax(S, w_floor)
  ccw[, w_ipcw := w]

  # Outcome Cox with robust SE by original subject id
  fit_out <- survival::coxph(outcome_formula, data = ccw, weights = w_ipcw, robust = TRUE, cluster = id)

  b  <- unname(stats::coef(fit_out)[["regime"]])
  se <- sqrt(stats::vcov(fit_out)[["regime","regime"]])

  list(beta = b, se = se, ess = ess_kish(w), fit_cens = fit_cens, fit_out = fit_out)
}

fit_ccw_over_time <- function(pp_dt,
                              t_vec,
                              cens_formula = survival::Surv(time_pp, dev_ind) ~ regime + age + sex + bmi + egfr + util + gall,
                              w_floor = 1e-3) {
  ccw0 <- make_ccw_clones(pp_dt)
  fit_cens <- survival::coxph(cens_formula, data = ccw0)
  lp <- stats::predict(fit_cens, type = "lp")

  res <- data.table::data.table(t = t_vec, beta = NA_real_, se = NA_real_, ess = NA_real_)

  for (j in seq_along(t_vec)) {
    t_cut <- t_vec[j]
    ccw <- data.table::copy(ccw0)

    ccw[, time_t  := pmin(time_pp, t_cut)]
    ccw[, delta_t := as.integer(delta_pp == 1L & time_pp <= t_cut)]

    S  <- get_S_from_cox(fit_cens, time_vec = ccw$time_t, lp_vec = lp)
    w  <- 1 / pmax(S, w_floor)
    ccw[, w_ipcw := w]

    fit_out <- tryCatch({
      survival::coxph(
        survival::Surv(time_t, delta_t) ~ regime,
        data = ccw,
        weights = w_ipcw,
        robust = TRUE,
        cluster = id
      )
    }, error = function(e) NULL)

    if (!is.null(fit_out)) {
      b  <- unname(stats::coef(fit_out)[["regime"]])
      se <- sqrt(stats::vcov(fit_out)[["regime","regime"]])
      res[j, `:=`(beta = b, se = se, ess = ess_kish(w))]
    }
  }

  list(fit_cens = fit_cens, curve = res)
}

# ---- TMLE-like DR risk at fixed horizon ----
make_binary_by_horizon <- function(pp_dt, t_cut) {
  d <- data.table::copy(pp_dt)
  d[, Y := as.integer(delta_pp == 1L & time_pp <= t_cut)]
  # uncensored for Y at t_cut if: event happened by t_cut OR still under follow-up at t_cut (no deviation before)
  d[, uncens := as.integer(Y == 1L | time_pp >= t_cut)]
  d[]
}

fit_dr_risk <- function(pp_dt,
                        t_cut,
                        w_floor = 0.05) {
  d <- make_binary_by_horizon(pp_dt, t_cut = t_cut)

  # censoring model for being uncensored at t_cut
  fit_c <- stats::glm(uncens ~ A0 + age + sex + bmi + egfr + util + gall,
                      data = d, family = stats::binomial())
  p_unc <- stats::predict(fit_c, type = "response")
  w_c   <- 1 / pmax(p_unc, w_floor)

  # treatment PS model (baseline)
  fit_g <- stats::glm(A0 ~ age + sex + bmi + egfr + util + gall,
                      data = d, family = stats::binomial())
  g1 <- stats::predict(fit_g, type = "response")
  g1 <- pmin(pmax(g1, 1e-3), 1 - 1e-3)

  # outcome regression using uncensored subjects, weighted by censoring weights
  d_obs <- d[uncens == 1L]
  fit_Q <- stats::glm(Y ~ A0 + age + sex + bmi + egfr + util + gall,
                      data = d_obs, family = stats::quasibinomial(), weights = w_c[uncens == 1L])

  # predict Q(a, W) for all
  d1 <- data.table::copy(d); d1[, A0 := 1L]
  d0 <- data.table::copy(d); d0[, A0 := 0L]
  Q1 <- stats::predict(fit_Q, newdata = d1, type = "response")
  Q0 <- stats::predict(fit_Q, newdata = d0, type = "response")

  # AIPW with censoring weights
  A <- d$A0
  Y <- d$Y
  U <- d$uncens

  term1 <- U * w_c * (A / g1) * (Y - Q1)
  term0 <- U * w_c * ((1 - A) / (1 - g1)) * (Y - Q0)

  psi1 <- mean(Q1 + term1)
  psi0 <- mean(Q0 + term0)

  # influence curve for RD
  ic1 <- (Q1 + term1) - psi1
  ic0 <- (Q0 + term0) - psi0
  ic_rd <- ic1 - ic0

  rd <- psi1 - psi0
  se_rd <- stats::sd(ic_rd) / sqrt(nrow(d))

  # log risk ratio (delta method)
  rr <- psi1 / pmax(psi0, 1e-9)
  logrr <- log(rr)
  # var(logrr) approx via IC
  ic_logrr <- (ic1 / pmax(psi1, 1e-9)) - (ic0 / pmax(psi0, 1e-9))
  se_logrr <- stats::sd(ic_logrr) / sqrt(nrow(d))

  list(
    t = t_cut,
    risk1 = psi1, risk0 = psi0,
    rd = rd, se_rd = se_rd,
    logrr = logrr, se_logrr = se_logrr
  )
}

fit_tmle_or_dr <- function(pp_dt, t_cut) {
  d <- make_binary_by_horizon(pp_dt, t_cut = t_cut)

  # censoring weights
  fit_c <- stats::glm(uncens ~ A0 + age + sex + bmi + egfr + util + gall,
                      data = d, family = stats::binomial())
  p_unc <- stats::predict(fit_c, type = "response")
  w_c   <- 1 / pmax(p_unc, 0.05)

  if (requireNamespace("tmle", quietly = TRUE)) {
    W <- d[, .(age, sex, bmi, egfr, util, gall)]
    out <- tmle::tmle(
      Y = d$Y,
      A = d$A0,
      W = W,
      family = "binomial",
      obsWeights = w_c
    )
    # Return RD and logRR-like summary when available
    # tmle object typically contains estimates in out$estimates$ATE / out$estimates$RR etc
    # We keep this light and robust: if fields missing, fall back to DR.
    rd <- tryCatch(out$estimates$ATE$psi, error = function(e) NA_real_)
    se_rd <- tryCatch(sqrt(out$estimates$ATE$var.psi), error = function(e) NA_real_)
    rr <- tryCatch(out$estimates$RR$psi, error = function(e) NA_real_)
    se_rr <- tryCatch(sqrt(out$estimates$RR$var.psi), error = function(e) NA_real_)
    list(method = "tmle", t = t_cut, rd = rd, se_rd = se_rd, rr = rr, se_rr = se_rr)
  } else {
    dr <- fit_dr_risk(pp_dt, t_cut = t_cut)
    list(method = "dr", t = dr$t, rd = dr$rd, se_rd = dr$se_rd, logrr = dr$logrr, se_logrr = dr$se_logrr,
         risk1 = dr$risk1, risk0 = dr$risk0)
  }
}
