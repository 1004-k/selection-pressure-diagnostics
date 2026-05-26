# R/12_ipcw_timeupdated_sensitivity.R
# ------------------------------------------------------------
# Fixed time-updated IPCW sensitivity estimator for the fixed-horizon
# per-protocol risk estimand.
#
# v7 patch rationale:
# The previous v5 implementation created interaction terms directly in the
# formula using expressions such as z_obs:splines::ns(tstart, df = 3). On some
# R/platform combinations this caused glm construction/prediction failures and
# all replicates fell back to a constant deviation probability. This version
# precomputes the spline basis/interactions and fits glm through an explicit numeric design matrix, so
# the time-updated deviation model is actually estimated rather than silently
# replaced by a fallback model.
# ------------------------------------------------------------

.tu_clip_probability <- function(p, eps = 1e-6) {
  pmin(pmax(as.numeric(p), eps), 1 - eps)
}

.tu_make_binary_horizon <- function(pp_dt, t_cut) {
  d <- data.table::copy(pp_dt)
  d[, Y := as.integer(delta_pp == 1L & time_pp <= t_cut)]
  # U = 1 if the fixed-horizon outcome contribution is observed before
  # artificial censoring. Events before t_cut remain observed outcomes.
  d[, U := as.integer(Y == 1L | time_pp >= t_cut)]
  d[]
}

.tu_get_baseline_W <- function(mis_spec = 0L) {
  W_full <- c("age", "sex", "bmi", "egfr", "util", "gall")
  W_red  <- c("age", "sex", "bmi")
  if (as.integer(mis_spec) == 1L) W_red else W_full
}

.tu_safe_glm_binomial <- function(formula, data, weights = NULL, maxit = 50L) {
  fit <- tryCatch(
    suppressWarnings(stats::glm(
      formula = formula,
      data = data,
      family = stats::binomial(),
      weights = weights,
      control = stats::glm.control(maxit = maxit)
    )),
    error = function(e) {
      attr(e, "tu_error") <- conditionMessage(e)
      NULL
    }
  )
  fit
}

.tu_safe_predict_response <- function(fit, newdata, fallback_prob) {
  if (is.null(fit)) return(rep(fallback_prob, nrow(newdata)))
  p <- tryCatch(
    suppressWarnings(as.numeric(stats::predict(fit, newdata = newdata, type = "response"))),
    error = function(e) rep(fallback_prob, nrow(newdata))
  )
  if (length(p) != nrow(newdata)) p <- rep(fallback_prob, nrow(newdata))
  p[!is.finite(p)] <- fallback_prob
  .tu_clip_probability(p)
}



.tu_numeric_design_matrix <- function(dt, rhs_cols) {
  rhs_cols <- unique(rhs_cols[rhs_cols %in% names(dt)])
  n <- nrow(dt)
  if (length(rhs_cols) == 0L) {
    X <- matrix(1, nrow = n, ncol = 1)
    colnames(X) <- "Intercept"
    return(X)
  }
  mats <- lapply(rhs_cols, function(nm) {
    v <- dt[[nm]]
    if (is.factor(v) || is.character(v)) {
      mm <- stats::model.matrix(~ v)[, -1, drop = FALSE]
      colnames(mm) <- paste0(nm, "_", make.names(colnames(mm)))
      return(mm)
    }
    v <- suppressWarnings(as.numeric(v))
    v[!is.finite(v)] <- 0
    matrix(v, ncol = 1, dimnames = list(NULL, nm))
  })
  X <- do.call(cbind, c(list(Intercept = rep(1, n)), mats))
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  X[!is.finite(X)] <- 0

  # Drop zero-variance and duplicate columns, keeping the intercept.
  keep <- rep(TRUE, ncol(X))
  if (ncol(X) > 1L) {
    for (j in 2:ncol(X)) {
      sj <- stats::sd(X[, j], na.rm = TRUE)
      if (!is.finite(sj) || sj == 0) keep[j] <- FALSE
    }
  }
  X <- X[, keep, drop = FALSE]
  X
}

.tu_fit_binomial_matrix <- function(y, X, maxit = 75L) {
  y <- as.integer(y == 1L)
  ok <- is.finite(y) & rowSums(!is.finite(X)) == 0
  X2 <- X[ok, , drop = FALSE]
  y2 <- y[ok]
  if (nrow(X2) < 20L || length(unique(y2)) < 2L) {
    return(NULL)
  }
  fit <- tryCatch(
    suppressWarnings(stats::glm.fit(
      x = X2,
      y = y2,
      family = stats::binomial(),
      control = stats::glm.control(maxit = maxit)
    )),
    error = function(e) NULL
  )
  if (is.null(fit) || is.null(fit$coefficients)) return(NULL)
  co <- fit$coefficients
  co[!is.finite(co)] <- 0
  list(coefficients = co, columns = colnames(X2), converged = isTRUE(fit$converged))
}

.tu_predict_binomial_matrix <- function(fit, X, fallback_prob) {
  if (is.null(fit)) return(rep(fallback_prob, nrow(X)))
  cols <- fit$columns
  missing <- setdiff(cols, colnames(X))
  if (length(missing) > 0L) {
    add <- matrix(0, nrow = nrow(X), ncol = length(missing))
    colnames(add) <- missing
    X <- cbind(X, add)
  }
  X <- X[, cols, drop = FALSE]
  eta <- drop(X %*% fit$coefficients)
  eta[!is.finite(eta)] <- stats::qlogis(fallback_prob)
  .tu_clip_probability(stats::plogis(eta))
}

.tu_add_time_basis <- function(ld, ns_df = 3L) {
  ns_df <- as.integer(ns_df)
  if (!is.finite(ns_df) || ns_df < 1L) ns_df <- 3L

  # For very small samples or degenerate tstart values, ns() can fail. In that
  # case we use simple polynomial time columns. Full runs should normally use ns.
  tb <- tryCatch(
    splines::ns(ld$tstart, df = ns_df),
    error = function(e) NULL
  )
  if (is.null(tb)) {
    t_scaled <- as.numeric(scale(ld$tstart))
    t_scaled[!is.finite(t_scaled)] <- 0
    tb <- sapply(seq_len(ns_df), function(k) t_scaled^k)
  }
  tb <- as.data.frame(tb)
  tb_names <- paste0("tu_t", seq_len(ncol(tb)))
  names(tb) <- tb_names
  data.table::setDT(tb)
  ld[, (tb_names) := tb]

  # Precompute interactions explicitly. This avoids formula parsing problems
  # with namespace-qualified functions inside interaction terms.
  for (nm in tb_names) {
    ld[, paste0("tu_z_", nm) := z_obs * get(nm)]
    ld[, paste0("tu_A_", nm) := A0 * get(nm)]
  }
  ld[, tu_A_z_obs := A0 * z_obs]

  list(ld = ld, tb_names = tb_names)
}

fit_ipcw_risk_timeupdated <- function(long_dt,
                                      pp_dt,
                                      t_cut,
                                      mis_spec = 0L,
                                      estimator = c("hajek", "ht"),
                                      w_floor = 0.05,
                                      ps_floor = 1e-3,
                                      dev_floor = 1e-6,
                                      ns_df = 3L) {
  estimator <- match.arg(estimator)
  stopifnot("id" %in% names(long_dt), "id" %in% names(pp_dt))

  d <- .tu_make_binary_horizon(pp_dt, t_cut = t_cut)
  W <- .tu_get_baseline_W(mis_spec)
  W <- W[W %in% names(d)]

  # Treatment model g(A|W), aligned with the primary standardized IPCW.
  f_g <- stats::as.formula(paste0("A0 ~ ", paste(W, collapse = " + ")))
  fit_g <- .tu_safe_glm_binomial(f_g, data = d)
  fallback_g <- mean(d$A0 == 1L, na.rm = TRUE)
  if (!is.finite(fallback_g) || fallback_g <= 0 || fallback_g >= 1) fallback_g <- 0.5
  g1 <- .tu_safe_predict_response(fit_g, newdata = d, fallback_prob = fallback_g)
  g1 <- pmin(pmax(g1, ps_floor), 1 - ps_floor)

  # Person-period deviation model.
  base_cols <- c("id", W)
  base_d <- unique(d[, base_cols, with = FALSE], by = "id")

  ld <- data.table::copy(long_dt)
  ld <- ld[tstart < t_cut]
  if (nrow(ld) == 0) stop("No person-period rows available before t_cut.")
  ld <- merge(ld, base_d, by = "id", all.x = TRUE)

  if (!("z_obs" %in% names(ld))) ld[, z_obs := 0]
  ld[!is.finite(z_obs), z_obs := 0]
  ld[!is.finite(tstart), tstart := 0]
  ld[, dev := as.integer(dev == 1L)]
  ld[, A0 := as.integer(A0 == 1L)]

  # Center/scale continuous covariates for numerical stability. Binary columns
  # are left unchanged. This stabilizes glm without changing the estimand.
  cont_cols <- intersect(c("age", "bmi", "egfr", "util", "z_obs"), names(ld))
  for (cc in cont_cols) {
    vv <- as.numeric(ld[[cc]])
    ss <- stats::sd(vv, na.rm = TRUE)
    mm <- mean(vv, na.rm = TRUE)
    if (is.finite(ss) && ss > 0) {
      ld[[cc]] <- (vv - mm) / ss
    } else {
      ld[[cc]] <- 0
    }
  }

  basis <- .tu_add_time_basis(ld, ns_df = ns_df)
  ld <- basis$ld
  tb_names <- basis$tb_names

  fallback_dev <- mean(ld$dev == 1L, na.rm = TRUE)
  if (!is.finite(fallback_dev)) fallback_dev <- 0.01
  fallback_dev <- .tu_clip_probability(fallback_dev, eps = dev_floor)

  if (length(unique(ld$dev[is.finite(ld$dev)])) < 2L) {
    p_dev <- rep(fallback_dev, nrow(ld))
    fit_dev <- NULL
    dev_model_status <- "constant_fallback_no_event_variation"
  } else {
    # Richer regime: baseline W + time-updated z + flexible time + z-by-time
    # and treatment-by-z terms. Reduced regime: reduced W + z_obs + flexible
    # time, without interactions.
    rhs <- c("A0", W, "z_obs", tb_names)
    if (as.integer(mis_spec) == 0L) {
      rhs <- c(rhs, paste0("tu_z_", tb_names), "tu_A_z_obs")
    }
    rhs <- unique(rhs[rhs %in% names(ld)])

    # v7 robustness patch:
    # Use an explicit numeric design matrix and glm.fit rather than formula-based
    # glm(). This avoids platform-specific formula/model.frame failures and keeps
    # the sensitivity estimator from silently falling back in every replicate.
    X_dev <- .tu_numeric_design_matrix(ld, rhs)
    fit_dev <- .tu_fit_binomial_matrix(ld$dev, X_dev)
    p_dev <- .tu_predict_binomial_matrix(fit_dev, X_dev, fallback_prob = fallback_dev)
    dev_model_status <- if (is.null(fit_dev)) {
      "glm_fallback"
    } else if (isFALSE(fit_dev$converged)) {
      if (as.integer(mis_spec) == 0L) "glm_timeupdated_richer_nonconverged" else "glm_timeupdated_reduced_nonconverged"
    } else {
      if (as.integer(mis_spec) == 0L) "glm_timeupdated_richer" else "glm_timeupdated_reduced"
    }
  }

  p_dev <- pmin(pmax(p_dev, dev_floor), 1 - dev_floor)
  ld[, p_no_dev := 1 - p_dev]

  c_by_id <- ld[, .(
    c_hat = prod(p_no_dev, na.rm = TRUE),
    n_intervals = .N,
    mean_p_dev = mean(1 - p_no_dev, na.rm = TRUE)
  ), by = id]
  d <- merge(d, c_by_id, by = "id", all.x = TRUE)
  d[is.na(c_hat) | !is.finite(c_hat), c_hat := 1]
  d[, c_hat := pmin(pmax(c_hat, w_floor), 1)]

  A <- d$A0
  Y <- d$Y
  U <- d$U

  H1 <- U * (A == 1L) / pmax(g1 * d$c_hat, ps_floor * w_floor)
  H0 <- U * (A == 0L) / pmax((1 - g1) * d$c_hat, ps_floor * w_floor)

  if (estimator == "ht") {
    risk1 <- mean(H1 * Y)
    risk0 <- mean(H0 * Y)
    ic1 <- H1 * Y - risk1
    ic0 <- H0 * Y - risk0
    den1 <- mean(H1)
    den0 <- mean(H0)
  } else {
    den1 <- mean(H1)
    den0 <- mean(H0)
    risk1 <- sum(H1 * Y) / pmax(sum(H1), 1e-12)
    risk0 <- sum(H0 * Y) / pmax(sum(H0), 1e-12)
    ic1 <- H1 * (Y - risk1) / pmax(den1, 1e-12)
    ic0 <- H0 * (Y - risk0) / pmax(den0, 1e-12)
  }

  rd <- risk1 - risk0
  ic_rd <- ic1 - ic0
  se_rd <- stats::sd(ic_rd, na.rm = TRUE) / sqrt(nrow(d))

  risk1_c <- pmax(risk1, 1e-9)
  risk0_c <- pmax(risk0, 1e-9)
  logrr <- log(risk1_c / risk0_c)
  ic_logrr <- (ic1 / risk1_c) - (ic0 / risk0_c)
  se_logrr <- stats::sd(ic_logrr, na.rm = TRUE) / sqrt(nrow(d))

  list(
    t = t_cut,
    method = paste0("ipcw_timeupdated_", estimator),
    dev_model_status = dev_model_status,
    mis_spec = as.integer(mis_spec),
    risk1 = as.numeric(risk1),
    risk0 = as.numeric(risk0),
    rd = as.numeric(rd),
    se_rd = as.numeric(se_rd),
    logrr = as.numeric(logrr),
    se_logrr = as.numeric(se_logrr),
    den1 = as.numeric(den1),
    den0 = as.numeric(den0),
    min_g1 = min(g1, na.rm = TRUE),
    max_g1 = max(g1, na.rm = TRUE),
    min_c_hat = min(d$c_hat, na.rm = TRUE),
    median_c_hat = stats::median(d$c_hat, na.rm = TRUE),
    median_intervals = stats::median(d$n_intervals, na.rm = TRUE),
    mean_interval_dev_prob = mean(d$mean_p_dev, na.rm = TRUE)
  )
}
