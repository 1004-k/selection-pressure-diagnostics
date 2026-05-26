# R/09_dr_aipw_ml.R
# ------------------------------------------------------------
# Doubly robust AIPW estimator for risk at fixed horizon (t_cut)
# with optional flexible outcome regression (ML) as a secondary analysis.
#
# Notes:
# - This module is designed for robustness benchmarking, not for proposing
#   a new estimator.
# - For acceptance-risk minimization, we keep:
#   * censoring model: logistic regression
#   * treatment model: logistic regression
#   * outcome regression: glm (default) or simple ML (optional)
# - ML is used only for the outcome regression (Q) to reduce tuning degrees
#   of freedom and reviewer questions.
# ------------------------------------------------------------

.clip01 <- function(p, eps = 1e-6) pmin(pmax(p, eps), 1 - eps)

.make_binary_by_horizon <- function(pp_dt, t_cut) {
  d <- data.table::copy(pp_dt)
  d[, Y := as.integer(delta_pp == 1L & time_pp <= t_cut)]
  # uncensored at t_cut if event by t_cut OR still under follow-up at t_cut
  d[, uncens := as.integer(Y == 1L | time_pp >= t_cut)]
  d[]
}

# simple fold assignment
.make_folds <- function(n, K = 2L, seed = 1L) {
  set.seed(seed)
  id <- sample.int(n)
  fold <- rep(seq_len(K), length.out = n)
  out <- integer(n)
  out[id] <- fold
  out
}

# ---- ML backends (binomial) ----
.fit_predict_binomial_glmnet <- function(train, test, y_col, w_col = NULL) {
  if (!requireNamespace("glmnet", quietly = TRUE)) stop("Package glmnet is required for ml_method='glmnet'.")
  y <- train[[y_col]]
  w <- if (!is.null(w_col)) train[[w_col]] else rep(1, nrow(train))
  x <- stats::model.matrix(~ . -1, data = train[, setdiff(names(train), c(y_col, w_col)), with = FALSE])
  xt <- stats::model.matrix(~ . -1, data = test[, setdiff(names(test), c(y_col, w_col)), with = FALSE])
  # small CV to choose lambda (kept modest for speed)
  cv <- glmnet::cv.glmnet(x = x, y = y, family = "binomial", weights = w, nfolds = 3)
  p <- as.numeric(stats::predict(cv, newx = xt, s = "lambda.min", type = "response"))
  .clip01(p)
}

.fit_predict_binomial_ranger <- function(train, test, y_col, w_col = NULL) {
  if (!requireNamespace("ranger", quietly = TRUE)) stop("Package ranger is required for ml_method='ranger'.")
  df <- data.table::copy(train)
  df[[y_col]] <- factor(df[[y_col]], levels = c(0, 1))
  form <- stats::as.formula(paste0(y_col, " ~ ."))
  w <- if (!is.null(w_col)) train[[w_col]] else NULL
  fit <- ranger::ranger(form, data = df, probability = TRUE, case.weights = w, num.trees = 200)
  pred <- stats::predict(fit, data = test)$predictions
  # second column corresponds to class '1'
  p <- as.numeric(pred[, 2])
  .clip01(p)
}

.fit_predict_binomial_mgcv <- function(train, test, y_col, w_col = NULL) {
  if (!requireNamespace("mgcv", quietly = TRUE)) stop("Package mgcv is required for ml_method='mgcv'.")
  # Use simple smooths for numeric covariates; factors enter linearly.
  # Build formula programmatically.
  cols <- setdiff(names(train), c(y_col, w_col))
  terms <- character(0)
  for (nm in cols) {
    if (is.numeric(train[[nm]]) || is.integer(train[[nm]])) {
      terms <- c(terms, paste0("s(", nm, ")"))
    } else {
      terms <- c(terms, nm)
    }
  }
  fml <- stats::as.formula(paste0(y_col, " ~ ", paste(terms, collapse = " + ")))
  w <- if (!is.null(w_col)) train[[w_col]] else NULL
  fit <- mgcv::gam(fml, data = train, family = stats::quasibinomial(), weights = w, method = "REML")
  p <- as.numeric(stats::predict(fit, newdata = test, type = "response"))
  .clip01(p)
}

.fit_predict_Q <- function(train, test, y_col = "Y", w_col = "w_c", ml_method = c("glm", "glmnet", "ranger", "mgcv")) {
  ml_method <- match.arg(ml_method)
  if (ml_method == "glm") {
    fit <- stats::glm(stats::as.formula(paste0(y_col, " ~ .")),
                     data = train, family = stats::quasibinomial(), weights = train[[w_col]])
    p <- as.numeric(stats::predict(fit, newdata = test, type = "response"))
    return(.clip01(p))
  }
  if (ml_method == "glmnet") return(.fit_predict_binomial_glmnet(train, test, y_col = y_col, w_col = w_col))
  if (ml_method == "ranger") return(.fit_predict_binomial_ranger(train, test, y_col = y_col, w_col = w_col))
  if (ml_method == "mgcv") return(.fit_predict_binomial_mgcv(train, test, y_col = y_col, w_col = w_col))
  stop("Unknown ml_method")
}

# ---- Main DR function ----
fit_dr_risk_rescue <- function(pp_dt,
                              t_cut,
                              mis_spec = 0L,
                              w_floor = 0.05,
                              use_ml_Q = FALSE,
                              ml_method = c("glm", "glmnet", "ranger", "mgcv"),
                              crossfit = FALSE,
                              cf_folds = 2L,
                              seed = 1L) {
  ml_method <- match.arg(ml_method)
  d <- .make_binary_by_horizon(pp_dt, t_cut = t_cut)

  # ---- Nuisance-regime knobs (keep simple and transparent) ----
  # mis_spec=0: richer baseline covariate set
  # mis_spec=1: reduced baseline covariate set. This is a practical
  #             under-adjustment regime, not a claim about the true model class.
  W_full <- c("age", "sex", "bmi", "egfr", "util", "gall")
  W_miss <- c("age", "sex", "bmi")
  W <- if (as.integer(mis_spec) == 1L) W_miss else W_full

  # censoring model for being uncensored at t_cut
  f_c <- stats::as.formula(paste0("uncens ~ A0 + ", paste(W, collapse = " + ")))
  fit_c <- stats::glm(f_c, data = d, family = stats::binomial())
  p_unc <- .clip01(stats::predict(fit_c, type = "response"))
  w_c <- 1 / pmax(p_unc, w_floor)
  d[, w_c := w_c]

  # treatment PS model
  f_g <- stats::as.formula(paste0("A0 ~ ", paste(W, collapse = " + ")))
  fit_g <- stats::glm(f_g, data = d, family = stats::binomial())
  g1 <- .clip01(stats::predict(fit_g, type = "response"), eps = 1e-3)

  # outcome regression using uncensored subjects, weighted by w_c
  d_obs <- d[uncens == 1L]

  # design for Q includes A0 and W
  Q_cols <- c("A0", W)
  d_obs_Q <- d_obs[, c("Y", "w_c", Q_cols), with = FALSE]
  d_Q_all <- d[, c(Q_cols), with = FALSE]

  # Q1/Q0 predictions
  if (!isTRUE(use_ml_Q) || ml_method == "glm") {
    # standard glm
    fit_Q <- stats::glm(stats::as.formula(paste0("Y ~ ", paste(Q_cols, collapse = " + "))),
                        data = d_obs, family = stats::quasibinomial(), weights = d_obs$w_c)

    d1 <- data.table::copy(d); d1[, A0 := 1L]
    d0 <- data.table::copy(d); d0[, A0 := 0L]
    Q1 <- .clip01(stats::predict(fit_Q, newdata = d1, type = "response"))
    Q0 <- .clip01(stats::predict(fit_Q, newdata = d0, type = "response"))
  } else {
    # ML for Q, optionally cross-fitted
    if (!isTRUE(crossfit)) {
      # fit once, predict counterfactuals by modifying A0
      train <- d_obs[, c("Y", "w_c", Q_cols), with = FALSE]
      test1 <- data.table::copy(d)[, c(Q_cols), with = FALSE]; test1[, A0 := 1L]
      test0 <- data.table::copy(d)[, c(Q_cols), with = FALSE]; test0[, A0 := 0L]
      Q1 <- .fit_predict_Q(train = train, test = test1, y_col = "Y", w_col = "w_c", ml_method = ml_method)
      Q0 <- .fit_predict_Q(train = train, test = test0, y_col = "Y", w_col = "w_c", ml_method = ml_method)
    } else {
      # cross-fitted Q: fold assignment on full data; fit on uncensored in training folds
      fold <- .make_folds(nrow(d), K = as.integer(cf_folds), seed = seed)
      Q1 <- rep(NA_real_, nrow(d))
      Q0 <- rep(NA_real_, nrow(d))
      for (k in seq_len(as.integer(cf_folds))) {
        idx_te <- which(fold == k)
        idx_tr <- which(fold != k)

        tr_all <- d[idx_tr]
        tr <- tr_all[uncens == 1L]
        if (nrow(tr) < 10) next
        train <- tr[, c("Y", "w_c", Q_cols), with = FALSE]

        te1 <- d[idx_te, c(Q_cols), with = FALSE]; te1[, A0 := 1L]
        te0 <- d[idx_te, c(Q_cols), with = FALSE]; te0[, A0 := 0L]

        Q1[idx_te] <- .fit_predict_Q(train = train, test = te1, y_col = "Y", w_col = "w_c", ml_method = ml_method)
        Q0[idx_te] <- .fit_predict_Q(train = train, test = te0, y_col = "Y", w_col = "w_c", ml_method = ml_method)
      }
      # any NAs fallback to non-crossfit
      if (anyNA(Q1) || anyNA(Q0)) {
        train <- d_obs[, c("Y", "w_c", Q_cols), with = FALSE]
        test1 <- data.table::copy(d)[, c(Q_cols), with = FALSE]; test1[, A0 := 1L]
        test0 <- data.table::copy(d)[, c(Q_cols), with = FALSE]; test0[, A0 := 0L]
        Q1 <- ifelse(is.na(Q1), .fit_predict_Q(train, test1, "Y", "w_c", ml_method), Q1)
        Q0 <- ifelse(is.na(Q0), .fit_predict_Q(train, test0, "Y", "w_c", ml_method), Q0)
      }
      Q1 <- .clip01(Q1)
      Q0 <- .clip01(Q0)
    }
  }

  A <- d$A0
  Y <- d$Y
  U <- d$uncens

  term1 <- U * d$w_c * (A / g1) * (Y - Q1)
  term0 <- U * d$w_c * ((1 - A) / (1 - g1)) * (Y - Q0)

  psi1 <- mean(Q1 + term1)
  psi0 <- mean(Q0 + term0)

  ic1 <- (Q1 + term1) - psi1
  ic0 <- (Q0 + term0) - psi0

  rd <- psi1 - psi0
  ic_rd <- ic1 - ic0
  se_rd <- stats::sd(ic_rd) / sqrt(nrow(d))

  rr <- psi1 / pmax(psi0, 1e-9)
  logrr <- log(rr)
  ic_logrr <- (ic1 / pmax(psi1, 1e-9)) - (ic0 / pmax(psi0, 1e-9))
  se_logrr <- stats::sd(ic_logrr) / sqrt(nrow(d))

  list(
    t = t_cut,
    mis_spec = as.integer(mis_spec),
    method_Q = ifelse(isTRUE(use_ml_Q), paste0("ML_", ml_method), "glm"),
    risk1 = psi1,
    risk0 = psi0,
    rd = rd,
    se_rd = se_rd,
    logrr = logrr,
    se_logrr = se_logrr
  )
}

# ---- Simple TMLE (parametric) for risk at horizon with censoring as missingness ----
# This is included as a supplementary comparator (acceptance-risk minimized):
# - All nuisance models are logistic regressions (no SuperLearner dependency).
# - Uses baseline covariates only (as in fit_dr_risk_rescue).

.logit <- function(p) log(p / (1 - p))
.expit <- function(x) 1 / (1 + exp(-x))

fit_tmle_risk_rescue <- function(pp_dt,
                                t_cut,
                                mis_spec = 0L,
                                w_floor = 0.05) {
  d <- .make_binary_by_horizon(pp_dt, t_cut = t_cut)

  W_full <- c("age", "sex", "bmi", "egfr", "util", "gall")
  W_miss <- c("age", "sex", "bmi")
  W <- if (as.integer(mis_spec) == 1L) W_miss else W_full

  # Treatment model g(A|W)
  f_g <- stats::as.formula(paste0("A0 ~ ", paste(W, collapse = " + ")))
  fit_g <- stats::glm(f_g, data = d, family = stats::binomial())
  g1 <- .clip01(stats::predict(fit_g, type = "response"), eps = 1e-3)

  # Missingness model c(Delta=1 | A,W)
  f_c <- stats::as.formula(paste0("uncens ~ A0 + ", paste(W, collapse = " + ")))
  fit_c <- stats::glm(f_c, data = d, family = stats::binomial())
  c_hat <- .clip01(stats::predict(fit_c, type = "response"))
  c_hat <- pmax(c_hat, w_floor)

  # Initial outcome model Q0(Y|A,W) on observed outcomes
  d_obs <- d[uncens == 1L]
  f_Q <- stats::as.formula(paste0("Y ~ A0 + ", paste(W, collapse = " + ")))
  fit_Q <- stats::glm(f_Q, data = d_obs, family = stats::binomial())

  QAW <- .clip01(stats::predict(fit_Q, newdata = d, type = "response"))
  d1 <- data.table::copy(d); d1[, A0 := 1L]
  d0 <- data.table::copy(d); d0[, A0 := 0L]
  Q1W <- .clip01(stats::predict(fit_Q, newdata = d1, type = "response"))
  Q0W <- .clip01(stats::predict(fit_Q, newdata = d0, type = "response"))

  # Clever covariates (combine treatment + censoring)
  H1 <- (d$A0 == 1L) / (g1 * c_hat)
  H0 <- (d$A0 == 0L) / ((1 - g1) * c_hat)

  # Targeting step: logistic fluctuation with offset logit(QAW)
  # Fit on observed outcomes only.
  dat_fluc <- data.table::data.table(
    Y = d_obs$Y,
    off = .logit(QAW[d$uncens == 1L]),
    H1 = H1[d$uncens == 1L],
    H0 = H0[d$uncens == 1L]
  )
  fit_eps <- stats::glm(Y ~ -1 + offset(off) + H1 + H0,
                        data = dat_fluc,
                        family = stats::binomial())
  eps <- stats::coef(fit_eps)
  eps[is.na(eps)] <- 0

  # Update counterfactual predictions
  Q1_star <- .clip01(.expit(.logit(Q1W) + eps["H1"] * (1 / (g1 * c_hat))))
  Q0_star <- .clip01(.expit(.logit(Q0W) + eps["H0"] * (1 / ((1 - g1) * c_hat))))

  psi1 <- mean(Q1_star)
  psi0 <- mean(Q0_star)

  # Influence curves
  A <- d$A0
  Y <- d$Y
  Delta <- d$uncens

  H1_full <- (A == 1L) / (g1 * c_hat)
  H0_full <- (A == 0L) / ((1 - g1) * c_hat)

  # Use updated Q for observed data at A,W
  QAW_star <- .clip01(.expit(.logit(QAW) + eps["H1"] * H1_full + eps["H0"] * H0_full))

  IC1 <- H1_full * Delta * (Y - QAW_star) + Q1_star - psi1
  IC0 <- H0_full * Delta * (Y - QAW_star) + Q0_star - psi0

  rd <- psi1 - psi0
  IC_rd <- IC1 - IC0
  se_rd <- stats::sd(IC_rd) / sqrt(nrow(d))

  rr <- psi1 / pmax(psi0, 1e-9)
  logrr <- log(rr)
  IC_logrr <- (IC1 / pmax(psi1, 1e-9)) - (IC0 / pmax(psi0, 1e-9))
  se_logrr <- stats::sd(IC_logrr) / sqrt(nrow(d))

  list(
    t = t_cut,
    mis_spec = as.integer(mis_spec),
    method_Q = "tmle_glm",
    risk1 = psi1,
    risk0 = psi0,
    rd = rd,
    se_rd = se_rd,
    logrr = logrr,
    se_logrr = se_logrr
  )
}

# Convenience wrapper used by the runner
fit_tmle_or_dr <- function(pp_dt, t_cut, mis_spec = 0L) {
  fit_tmle_risk_rescue(pp_dt = pp_dt, t_cut = t_cut, mis_spec = mis_spec)
}
