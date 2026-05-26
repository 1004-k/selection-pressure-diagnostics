# R/06_performance.R
# ------------------------------------------------------------
# Performance summaries for simulation outputs
# ------------------------------------------------------------

summarize_perf <- function(beta_hat, se_hat, beta_true) {
  z  <- stats::qnorm(0.975)
  lo <- beta_hat - z * se_hat
  hi <- beta_hat + z * se_hat
  list(
    bias     = mean(beta_hat - beta_true, na.rm = TRUE),
    rmse     = sqrt(mean((beta_hat - beta_true)^2, na.rm = TRUE)),
    cover    = mean(lo <= beta_true & beta_true <= hi, na.rm = TRUE),
    sign_rev = mean(beta_hat < 0, na.rm = TRUE),
    # Under causal null (beta_true = 0), CI exclusion is type I error
    ci_excl0 = mean(!(lo <= 0 & 0 <= hi), na.rm = TRUE)
  )
}
