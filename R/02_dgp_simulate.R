# R/02_dgp_simulate.R
# ------------------------------------------------------------
# Minimal but flexible DGP for:
#   - time-varying selection pressure via gamma(t) path (Axis A)
#   - nonlinear threshold-jump selection (Axis B)
#   - imperfect prognostic measurement (Axis C)
#
# Output:
#   - long_dt: counting-process style data with deviation (dev) and outcome (event)
#   - pp_dt:  one row per id for per-protocol analysis (censor at deviation)
# ------------------------------------------------------------

make_gamma_path <- function(axisA, tstart_vec, gamma_max = 1.0, t0 = NULL) {
  axisA <- match.arg(axisA, c("flat", "increasing", "late_surge"))

  if (axisA == "flat") {
    gamma <- rep(0, length(tstart_vec))
  } else if (axisA == "increasing") {
    # linear from 0 to gamma_max over follow-up
    gamma <- gamma_max * (tstart_vec / max(tstart_vec))
    gamma[!is.finite(gamma)] <- 0
  } else {
    # late surge after t0 (default: 70% of follow-up)
    if (is.null(t0)) t0 <- 0.7 * max(tstart_vec)
    gamma <- ifelse(tstart_vec < t0, 0, gamma_max)
  }
  gamma
}

simulate_one_dataset <- function(N = 1000,
                                 t_max = 5,
                                 dt = 0.25,
                                 beta_true = log(0.80),
                                 # hazards
                                 lambda_event_base = 0.10,
                                 lambda_dev_base   = 0.25,
                                 # outcome prognosis effect
                                 theta_z = 0.50,
                                 # Axis A
                                 axisA = c("flat","increasing","late_surge"),
                                 gamma_max = 1.0,
                                 t0 = NULL,
                                 # Axis B
                                 axisB = c("none","threshold_jump"),
                                 thresh_c = 1.0,
                                 thresh_mult = 3.0,
                                 # Axis C
                                 rho_meas = 1.0,
                                 # latent z dynamics
                                 rho_z = 0.80,
                                 seed = NULL,
                                 force_A = NULL,
                                 disable_deviation = FALSE) {
  if (!is.null(seed)) set.seed(seed)

  axisA <- match.arg(axisA, c("flat","increasing","late_surge"))
  axisB <- match.arg(axisB, c("none","threshold_jump"))

  # --- baseline covariates (kept similar to your prior simulation style) ---
  base <- data.table::data.table(
    id  = 1:N,
    age = pmin(pmax(stats::rnorm(N, 62, 10), 40), 90),
    sex = stats::rbinom(N, 1, 0.45),
    bmi = pmin(pmax(stats::rnorm(N, 31, 6), 18), 60),
    egfr= pmin(pmax(stats::rnorm(N, 75, 18), 15), 120),
    util= stats::rgamma(N, shape = 2, rate = 0.5),
    gall= stats::rbinom(N, 1, stats::plogis(-2 + 0.02 * (stats::rnorm(N, 62, 10) - 60)))
  )

  # treatment assignment (baseline)
  lin <- with(base,
              -1.0 + 0.02*(age-60) + 0.08*(bmi-30) - 0.015*(egfr-75) + 0.25*util + 0.7*gall)
  ps <- stats::plogis(lin)
  if (is.null(force_A)) {
    base[, A0 := stats::rbinom(.N, 1, ps)]
  } else {
    force_A <- as.integer(force_A)
    if (!force_A %in% c(0L, 1L)) stop("force_A must be NULL, 0, or 1.")
    base[, A0 := force_A]
  }

  # baseline prognosis driver
  lp_prog0 <- with(base,
                   0.02*(age-60) + 0.10*(bmi-30) - 0.01*(egfr-75) + 0.25*gall + 0.15*util)
  base[, lp_prog0 := lp_prog0]
  base[, z0_true := as.numeric(scale(lp_prog0))]

  # time grid
  t_grid <- seq(0, t_max, by = dt)
  if (tail(t_grid, 1) < t_max) t_grid <- c(t_grid, t_max)
  tstart_vec <- t_grid[-length(t_grid)]
  tstop_vec  <- t_grid[-1L]
  K <- length(tstart_vec)

  gamma_vec <- make_gamma_path(axisA, tstart_vec, gamma_max = gamma_max, t0 = t0)

  # allocate containers
  long_list <- vector("list", N)

  for (i in seq_len(N)) {
    id <- base$id[i]
    A0 <- base$A0[i]

    # latent z(t)
    z_true <- numeric(K)
    z_true[1] <- base$z0_true[i]
    if (K >= 2) {
      for (k in 2:K) {
        z_true[k] <- rho_z * z_true[k-1] + sqrt(1 - rho_z^2) * stats::rnorm(1)
      }
    }

    # observed z (measurement error)
    z_obs <- rho_meas * z_true + sqrt(1 - rho_meas^2) * stats::rnorm(K)

    # simulate competing times within each interval
    alive <- TRUE
    deviated <- FALSE
    rows <- list()

    current_t <- 0

    for (k in seq_len(K)) {
      if (!alive || deviated) break

      t0k <- tstart_vec[k]
      t1k <- tstop_vec[k]

      # hazards at interval start (piecewise constant within dt)
      haz_y <- lambda_event_base * exp(theta_z * z_true[k] + beta_true * A0)

      if (isTRUE(disable_deviation)) {
        haz_d <- 0
      } else {
        haz_d <- lambda_dev_base * exp(gamma_vec[k] * z_true[k])
        if (axisB == "threshold_jump" && z_true[k] > thresh_c) {
          haz_d <- haz_d * thresh_mult
        }
      }

      # draw times within interval
      # NOTE: rexp() can (very rarely) return values extremely close to 0.
      # survival::aeqSurv() treats near-zero intervals as having effective length 0.
      # We therefore enforce a tiny minimum time within an interval.
      eps_dt <- 1e-6
      ty <- max(stats::rexp(1, rate = haz_y), eps_dt)
      td <- if (isTRUE(disable_deviation) || !is.finite(haz_d) || haz_d <= 0) Inf else max(stats::rexp(1, rate = haz_d), eps_dt)

      t_event <- min(ty, td)

      if (t_event < dt) {
        # event happens inside interval
        t_end <- t0k + t_event
        if ((t_end - t0k) <= eps_dt) t_end <- t0k + eps_dt
        dev <- as.integer(td < ty)
        ev  <- as.integer(ty <= td)
        rows[[length(rows)+1L]] <- data.table::data.table(
          id = id,
          tstart = t0k,
          tstop  = t_end,
          interval = k,
          A0 = A0,
          z_true = z_true[k],
          z_obs  = z_obs[k],
          gamma_true = gamma_vec[k],
          dev   = dev,
          event = ev
        )
        if (dev == 1) deviated <- TRUE
        if (ev  == 1) alive    <- FALSE
        break
      } else {
        # no event in interval; carry forward
        rows[[length(rows)+1L]] <- data.table::data.table(
          id = id,
          tstart = t0k,
          tstop  = t1k,
          interval = k,
          A0 = A0,
          z_true = z_true[k],
          z_obs  = z_obs[k],
          gamma_true = gamma_vec[k],
          dev   = 0L,
          event = 0L
        )
      }
    }

    if (length(rows) == 0) {
      # degenerate safety
      rows[[1]] <- data.table::data.table(
        id = id, tstart = 0, tstop = dt, interval = 1L,
        A0 = A0, z_true = base$z0_true[i], z_obs = base$z0_true[i],
        gamma_true = gamma_vec[1], dev = 0L, event = 0L
      )
    }

    long_list[[i]] <- data.table::rbindlist(rows)
  }

  long_dt <- data.table::rbindlist(long_list)


  # per-protocol summary per id (censor at deviation)
  # dev_time: first deviation time; event_time: first outcome time
  dev_time <- long_dt[dev == 1, .(dev_time = min(tstop)), by = id]
  ev_time  <- long_dt[event == 1, .(event_time = min(tstop)), by = id]

  # build pp_dt from baseline (ensures all ids are present)
  pp_dt <- data.table::copy(base)[, .(id, A0, age, sex, bmi, egfr, util, gall)]

  pp_dt <- merge(pp_dt, dev_time, by = "id", all.x = TRUE)
  pp_dt <- merge(pp_dt, ev_time,  by = "id", all.x = TRUE)

  pp_dt[is.na(dev_time), dev_time := Inf]
  pp_dt[is.na(event_time), event_time := Inf]

  pp_dt[, time_pp := pmin(event_time, dev_time, t_max)]
  pp_dt[, delta_pp := as.integer(event_time <= dev_time & event_time <= t_max)]
  pp_dt[, dev_ind  := as.integer(dev_time < pmin(event_time, t_max))]

  # baseline observed z (for censoring models if desired)
  z0_obs <- long_dt[interval == 1, .(id, z0_obs = z_obs)]
  pp_dt <- merge(pp_dt, z0_obs, by = "id", all.x = TRUE)

  list(
    base   = base,
    long_dt= long_dt,
    pp_dt  = pp_dt,
    grid   = data.table::data.table(
      k = seq_len(K),
      tstart = tstart_vec,
      tstop  = tstop_vec,
      gamma_true = gamma_vec
    )
  )
}
