# R/05_tipping.R
# ------------------------------------------------------------
# Tipping diagnostics based on weight instability:
#   - rESS(t) = ESS(t) / N(t)
#   - tail share (top 1% weight mass)
# Returns first tipping time t*.
# ------------------------------------------------------------


# ---- helpers ----
ess_kish <- function(w) {
  (sum(w)^2) / sum(w^2)
}

tail_share <- function(w, q = 0.99) {
  w <- w[is.finite(w)]
  if (length(w) == 0) return(NA_real_)
  s <- sum(w)
  if (!is.finite(s) || s <= 0) return(NA_real_)
  cut <- stats::quantile(w, probs = q, names = FALSE, type = 7, na.rm = TRUE)
  sum(w[w >= cut]) / s
}

compute_weight_diagnostics <- function(long_dt,
                                       weights_dt,
                                       time_col = "tstart",
                                       id_col = "id",
                                       w_col = "w",
                                       q_tail = 0.99) {
  dt <- data.table::copy(long_dt)
  w  <- data.table::copy(weights_dt)

  # merge weights at time grid; expects weights_dt to have (id, time, w)
  data.table::setnames(w, c(id_col, "t", w_col), c(id_col, "t_merge", w_col))
  dt[, t_merge := get(time_col)]
  dt <- merge(dt, w, by = c(id_col, "t_merge"), all.x = TRUE)

  # keep one row per id per time (at risk at that time)
  dt <- dt[!is.na(get(w_col))]
  dt <- unique(dt, by = c(id_col, "t_merge"))

  diag <- dt[, {
    ww <- get(w_col)
    N  <- .N
    ESS <- ess_kish(ww)
    rESS <- ifelse(is.finite(ESS), ESS / N, NA_real_)
    tail <- tail_share(ww, q = q_tail)
    .(N = N, ESS = ESS, rESS = rESS, tail_share = tail)
  }, by = .(t = t_merge)]

  data.table::setorder(diag, t)
  diag[]
}

detect_tipping <- function(diag_dt, c_ess = 0.25, c_tail = 0.10) {
  d <- data.table::copy(diag_dt)

  d[, tip_ess  := as.integer(!is.na(rESS) & rESS < c_ess)]
  d[, tip_tail := as.integer(!is.na(tail_share) & tail_share > c_tail)]
  d[, tipped   := as.integer(tip_ess == 1 | tip_tail == 1)]

  if (any(d$tipped == 1)) {
    t_star <- d[tipped == 1, min(t)]
  } else {
    t_star <- NA_real_
  }

  list(
    t_star = t_star,
    diag   = d
  )
}
