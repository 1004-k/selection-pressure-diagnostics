# R/04_pp_ipcw.R
# ------------------------------------------------------------
# Per-protocol IPCW (primary estimator)
# - Fit censoring model for deviation time
# - Compute weights w(t)=1/S_cens(t | X)
# - Fit weighted Cox for outcome
# ------------------------------------------------------------


# ---- helpers ----
ess_kish <- function(w) {
  (sum(w)^2) / sum(w^2)
}

# basehaz -> predicted survival at time t for given linear predictor
# basehaz -> predicted survival at time t for given linear predictor
get_S_from_cox <- function(fit, time_vec, lp_vec) {
  # IMPORTANT: basehaz(fit) can fail later if the model was fit with a local
  # data symbol (e.g., data = d inside a function) and that symbol no longer exists.
  # To make downstream calls robust, we store basehaz at fit time as fit$bh0 and use it here.
  bh <- fit$bh0
  if (is.null(bh)) {
    bh_try <- tryCatch(survival::basehaz(fit, centered = FALSE), error = function(e) e)
    if (inherits(bh_try, "error")) {
      stop("Failed to compute basehaz() for censoring model: ", conditionMessage(bh_try), call. = FALSE)
    }
    bh <- bh_try
  }
  bh <- bh[order(bh$time), ]
  idx <- findInterval(time_vec, bh$time)
  H0_t <- ifelse(idx == 0, 0, bh$hazard[idx])
  S <- exp(- H0_t * exp(lp_vec))
  pmax(S, 1e-9)
}

fit_pp_ipcw <- function(pp_dt,
                        cens_formula = survival::Surv(time_pp, dev_ind) ~ A0 + age + sex + bmi + egfr + util + gall,
                        outcome_formula = survival::Surv(time_pp, delta_pp) ~ A0,
                        w_floor = 1e-3) {
  d <- data.table::copy(pp_dt)

  # Censoring model: time to deviation (dev_ind)
  fit_cens <- survival::coxph(cens_formula, data = d, model = TRUE, x = TRUE)

  # Store basehaz now so later helpers do not depend on the local data symbol 'd'
  fit_cens$bh0 <- survival::basehaz(fit_cens, centered = FALSE)

  lp <- stats::predict(fit_cens, type = "lp")
  S  <- get_S_from_cox(fit_cens, time_vec = d$time_pp, lp_vec = lp)

  w <- 1 / pmax(S, w_floor)
  d[, w_ipcw := w]

  fit_out <- survival::coxph(outcome_formula, data = d, weights = w_ipcw)

  b  <- unname(stats::coef(fit_out)[1])
  se <- sqrt(stats::vcov(fit_out)[1,1])

  list(
    beta = b,
    se   = se,
    ess  = ess_kish(w),
    fit_cens = fit_cens,
    fit_out  = fit_out,
    weights  = d$w_ipcw
  )
}

# ------------------------------------------------------------
# IPCW risk at fixed horizon (t_cut)
# - Uses the same censoring model as IPCW-Cox but reports
#   risk-based estimands (RD, logRR) at t_cut.
# - Intended for Paper B so that IPCW and DR/TMLE share the
#   same estimand family.
# ------------------------------------------------------------

.extract_survfit_at <- function(sf, t_cut) {
  ss <- summary(sf, times = t_cut, extend = TRUE)
  strata <- as.character(ss$strata)
  data.table::data.table(strata = strata, surv = ss$surv, se_surv = ss$std.err)
}

fit_ipcw_risk <- function(pp_dt,
                          t_cut,
                          cens_formula = survival::Surv(time_pp, dev_ind) ~ A0 + age + sex + bmi + egfr + util + gall,
                          w_floor = 1e-3,
                          fit_cens = NULL,
                          lp = NULL) {
  d0 <- data.table::copy(pp_dt)

  if (is.null(fit_cens)) {
    fit_cens <- survival::coxph(cens_formula, data = d0, model = TRUE, x = TRUE)
  # Store basehaz now so later helpers do not depend on local data symbol 'd0'
  fit_cens$bh0 <- survival::basehaz(fit_cens, centered = FALSE)

  # Store basehaz for downstream S(t) calculations
  fit_cens$bh0 <- survival::basehaz(fit_cens, centered = FALSE)
  }

  # If the censoring model was passed in from elsewhere and does not carry bh0,
  # try to compute it; if that fails (common when the original data symbol is gone), refit locally.
  if (is.null(fit_cens$bh0)) {
    bh_try <- tryCatch(survival::basehaz(fit_cens, centered = FALSE), error = function(e) NULL)
    if (!is.null(bh_try)) {
      fit_cens$bh0 <- bh_try
    } else {
    fit_cens <- survival::coxph(cens_formula, data = d0, model = TRUE, x = TRUE)
  # Store basehaz now so later helpers do not depend on local data symbol 'd0'
  fit_cens$bh0 <- survival::basehaz(fit_cens, centered = FALSE)
      fit_cens$bh0 <- survival::basehaz(fit_cens, centered = FALSE)
      lp <- NULL
    }
  }
  if (is.null(lp)) {
    lp <- stats::predict(fit_cens, type = "lp")
  }

  d0[, time_t := pmin(time_pp, t_cut)]
  d0[, delta_t := as.integer(delta_pp == 1L & time_pp <= t_cut)]

  S  <- get_S_from_cox(fit_cens, time_vec = d0$time_t, lp_vec = lp)
  w  <- 1 / pmax(S, w_floor)
  d0[, w_ipcw_t := w]

  sf <- survival::survfit(survival::Surv(time_t, delta_t) ~ A0, data = d0, weights = w_ipcw_t)
  at <- .extract_survfit_at(sf, t_cut)

  surv1 <- at[grepl("A0=1", strata), surv]
  surv0 <- at[grepl("A0=0", strata), surv]
  se1_s <- at[grepl("A0=1", strata), se_surv]
  se0_s <- at[grepl("A0=0", strata), se_surv]

  if (length(surv1) == 0 || length(surv0) == 0) {
    if (nrow(at) >= 2) {
      surv0 <- at$surv[1]; se0_s <- at$se_surv[1]
      surv1 <- at$surv[2]; se1_s <- at$se_surv[2]
    } else {
      return(list(
        t = t_cut,
        risk1 = NA_real_, risk0 = NA_real_,
        rd = NA_real_, se_rd = NA_real_,
        logrr = NA_real_, se_logrr = NA_real_
      ))
    }
  }

  risk1 <- 1 - as.numeric(surv1)
  risk0 <- 1 - as.numeric(surv0)
  se1_r <- as.numeric(se1_s)
  se0_r <- as.numeric(se0_s)

  rd <- risk1 - risk0
  se_rd <- sqrt(se1_r^2 + se0_r^2)

  rr <- risk1 / pmax(risk0, 1e-9)
  logrr <- log(rr)
  se_logrr <- sqrt((se1_r / pmax(risk1, 1e-9))^2 + (se0_r / pmax(risk0, 1e-9))^2)

  list(
    t = t_cut,
    risk1 = risk1,
    risk0 = risk0,
    rd = rd,
    se_rd = se_rd,
    logrr = logrr,
    se_logrr = se_logrr,
    fit_cens = fit_cens
  )
}

# Time-specific PP-IPCW estimates (for inference-tipping curves in simulation)
# - Fit censor model once; evaluate weights at each truncation time
fit_pp_ipcw_over_time <- function(pp_dt,
                                  t_vec,
                                  cens_formula = survival::Surv(time_pp, dev_ind) ~ A0 + age + sex + bmi + egfr + util + gall,
                                  outcome_formula = survival::Surv(time_pp, delta_pp) ~ A0,
                                  w_floor = 1e-3) {
  d0 <- data.table::copy(pp_dt)

    fit_cens <- survival::coxph(cens_formula, data = d0, model = TRUE, x = TRUE)
  # Store basehaz now so later helpers do not depend on local data symbol 'd0'
  fit_cens$bh0 <- survival::basehaz(fit_cens, centered = FALSE)
  lp <- stats::predict(fit_cens, type = "lp")

  # Precompute S(t) for each t in t_vec and each subject
  res <- data.table::data.table(t = t_vec, beta = NA_real_, se = NA_real_, ess = NA_real_)

  for (j in seq_along(t_vec)) {
    t_cut <- t_vec[j]

    d <- data.table::copy(d0)
    d[, time_t  := pmin(time_pp, t_cut)]
    d[, delta_t := as.integer(delta_pp == 1 & time_pp <= t_cut)]

    S  <- get_S_from_cox(fit_cens, time_vec = d$time_t, lp_vec = lp)
    w  <- 1 / pmax(S, w_floor)
    d[, w_ipcw := w]

    fit_out <- tryCatch({
      survival::coxph(
        survival::Surv(time_t, delta_t) ~ A0,
        data = d,
        weights = w_ipcw
      )
    }, error = function(e) NULL)

    if (!is.null(fit_out)) {
      b  <- unname(stats::coef(fit_out)[["A0"]])
      se <- sqrt(stats::vcov(fit_out)[["A0","A0"]])
      res[j, `:=`(beta = b, se = se, ess = ess_kish(w))]
    }
  }

  list(fit_cens = fit_cens, curve = res)
}
