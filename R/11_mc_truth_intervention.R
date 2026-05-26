# R/11_mc_truth_intervention.R
# ------------------------------------------------------------
# Intervention-based Monte Carlo truth for fixed-horizon risks.
#
# Why this file exists:
# For the JCI submission, the true risk contrast should not be defined by
# applying one of the evaluated estimators to a very large simulated dataset.
# That strategy is common in applied simulation work, but a causal-inference
# reviewer can object that the truth is estimator-dependent.
#
# This module computes the target risks by direct intervention on treatment
# assignment under the known data-generating mechanism:
#   psi_a(t) = P(Y^a <= t) under sustained adherence/no artificial deviation.
#
# The deviation mechanism is deliberately disabled for the truth calculation.
# This matches the per-protocol counterfactual risk under sustained adherence.
# Deviation mechanisms still determine the observed-data problem and therefore
# the operating characteristics of IPCW and DR-AIPW in the finite samples.
# ------------------------------------------------------------

.simulate_intervention_risk_vectorized <- function(N,
                                                   A,
                                                   t_max = 5,
                                                   t_cut = NULL,
                                                   dt = 0.25,
                                                   beta_true = log(0.80),
                                                   lambda_event_base = 0.10,
                                                   theta_z = 0.50,
                                                   rho_z = 0.80,
                                                   seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  A <- as.integer(A)
  if (!A %in% c(0L, 1L)) stop("A must be 0 or 1.")
  if (is.null(t_cut)) t_cut <- t_max
  t_cut <- as.numeric(t_cut)

  # Same baseline covariate distribution as simulate_one_dataset().
  age <- pmin(pmax(stats::rnorm(N, 62, 10), 40), 90)
  sex <- stats::rbinom(N, 1, 0.45)
  bmi <- pmin(pmax(stats::rnorm(N, 31, 6), 18), 60)
  egfr <- pmin(pmax(stats::rnorm(N, 75, 18), 15), 120)
  util <- stats::rgamma(N, shape = 2, rate = 0.5)
  gall <- stats::rbinom(N, 1, stats::plogis(-2 + 0.02 * (stats::rnorm(N, 62, 10) - 60)))

  lp_prog0 <- 0.02 * (age - 60) + 0.10 * (bmi - 30) - 0.01 * (egfr - 75) + 0.25 * gall + 0.15 * util
  z <- as.numeric(scale(lp_prog0))

  t_grid <- seq(0, t_max, by = dt)
  if (tail(t_grid, 1) < t_max) t_grid <- c(t_grid, t_max)
  tstart_vec <- t_grid[-length(t_grid)]
  K <- length(tstart_vec)

  event_time <- rep(Inf, N)
  alive <- rep(TRUE, N)
  eps_dt <- 1e-6

  for (k in seq_len(K)) {
    idx <- which(alive)
    if (length(idx) == 0L) break

    haz_y <- lambda_event_base * exp(theta_z * z[idx] + beta_true * A)
    haz_y[!is.finite(haz_y) | haz_y <= 0] <- NA_real_
    ty <- stats::rexp(length(idx), rate = haz_y)
    ty <- pmax(ty, eps_dt)

    ev <- is.finite(ty) & ty < dt
    if (any(ev)) {
      event_time[idx[ev]] <- tstart_vec[k] + ty[ev]
      alive[idx[ev]] <- FALSE
    }

    if (k < K) {
      # Update latent prognosis for everyone still alive. This reproduces the
      # AR(1) latent prognosis process used in the observed-data DGP.
      idx2 <- which(alive)
      if (length(idx2) > 0L) {
        z[idx2] <- rho_z * z[idx2] + sqrt(1 - rho_z^2) * stats::rnorm(length(idx2))
      }
    }
  }

  list(
    risk = mean(event_time <= t_cut),
    event_time = event_time
  )
}

compute_intervention_truth_risk <- function(scenario_id,
                                            truth = c("non_null", "null"),
                                            N_truth = 200000L,
                                            B_truth = 5L,
                                            t_cut = 5,
                                            t_max = 5,
                                            dt = 0.25,
                                            beta_true_pp = log(0.80),
                                            seed_base = 2026L) {
  truth <- match.arg(as.character(truth), c("non_null", "null"))
  beta <- if (truth == "non_null") beta_true_pp else 0

  # Under the causal null, use exact zero contrasts. This avoids unnecessary
  # Monte Carlo noise in Type I error calculations.
  if (truth == "null" || isTRUE(beta == 0)) {
    return(data.table::data.table(
      scenario_id = scenario_id,
      truth = truth,
      t_cut = t_cut,
      beta_true_pp = beta_true_pp,
      N_truth = as.integer(N_truth),
      B_truth = as.integer(B_truth),
      truth_method = "intervention_exact_null",
      risk1_true = NA_real_,
      risk0_true = NA_real_,
      logrr_true = 0,
      rd_true = 0,
      logrr_true_sd = 0,
      rd_true_sd = 0,
      n_ok = as.integer(B_truth)
    ))
  }

  risk1_vec <- numeric(B_truth)
  risk0_vec <- numeric(B_truth)
  logrr_vec <- numeric(B_truth)
  rd_vec <- numeric(B_truth)

  for (b in seq_len(B_truth)) {
    seed_b <- as.integer(seed_base + 1009L * b)

    # Use common random numbers for the two treatment interventions. This
    # reduces Monte Carlo noise in the truth contrast and is acceptable because
    # the estimand is defined by the marginal risks, not by independent draws.
    r1 <- .simulate_intervention_risk_vectorized(
      N = N_truth, A = 1L, t_max = t_max, t_cut = t_cut, dt = dt,
      beta_true = beta, seed = seed_b
    )$risk
    r0 <- .simulate_intervention_risk_vectorized(
      N = N_truth, A = 0L, t_max = t_max, t_cut = t_cut, dt = dt,
      beta_true = beta, seed = seed_b
    )$risk

    risk1_vec[b] <- r1
    risk0_vec[b] <- r0
    logrr_vec[b] <- log(pmax(r1, 1e-12) / pmax(r0, 1e-12))
    rd_vec[b] <- r1 - r0
  }

  data.table::data.table(
    scenario_id = scenario_id,
    truth = truth,
    t_cut = t_cut,
    beta_true_pp = beta_true_pp,
    N_truth = as.integer(N_truth),
    B_truth = as.integer(B_truth),
    truth_method = "intervention_sustained_adherence",
    risk1_true = mean(risk1_vec),
    risk0_true = mean(risk0_vec),
    logrr_true = mean(logrr_vec),
    rd_true = mean(rd_vec),
    logrr_true_sd = if (length(logrr_vec) > 1L) stats::sd(logrr_vec) else NA_real_,
    rd_true_sd = if (length(rd_vec) > 1L) stats::sd(rd_vec) else NA_real_,
    n_ok = as.integer(length(logrr_vec))
  )
}
