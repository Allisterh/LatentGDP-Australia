# Production fit: m13 only (m12's structural Sigma + SV on the measurement innovations +
# AR(1) measurement errors, plus stochastic volatility on the latent innovation).
# Runs in CI (weekly; picks up each quarterly ABS 5206.0 release).
# Fits the production model, saves the fit for downstream use, and then applies a
# strict convergence gate: if any check fails the script stop()s, so CI never
# publishes an estimate from a fit that has not converged.
#
# Why m13 (NOTES.md, "2026-10-18"): it replaced m12 under a rule fixed before m13's LFO
# results were known - pass this gate on two seeds, a reliable LFO score, beat m12 on
# PSIS-LFO-CV by z_NW >= 2 (also with each refit's first attempt only), and a local fit
# <= 30 min. Observed: both seeds pass; m13 - m12 = +59.65 nats (Newey-West SE 28.89,
# z 2.06; first attempt only +59.79, SE 29.00, z 2.06); local fit 23.1 min (seed 1).
# m13 was proposed after the re-specified m10's results were seen, so its gain is post hoc.
# (m12 was production from 2026-10-17 to 2026-10-18; m9 before that.)
#
# m13 is sampled with mu integrated out by a Kalman filter on quasi-differenced data,
# conditional on the latent log-variance path g and the measurement log-variances h (see
# the header of stan/m13_serial_latent_sv.stan); mu is drawn exactly by FFBS in generated
# quantities. g_t = log(tau^2) for all t gives m12 exactly.

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
})

dir.create("output", showWarnings = FALSE, recursive = TRUE)

source("R/production_config.R")   # PRODUCTION_START = 1980 Q1 (RDP 2014-12)
gdp_growth <- load_production_growth()

stan_data <- list(
  T = nrow(gdp_growth),
  K = 3,
  y = as.matrix(gdp_growth[, c("g_E", "g_I", "g_P")])
)

mod <- cmdstan_model("stan/m13_serial_latent_sv.stan")
fit <- mod$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 2000,
  seed = 1,
  refresh = 500,
  adapt_delta = 0.99,
  max_treedepth = 12
)

# Persist for the estimate writer; not committed (see .gitignore). Saved before
# the gate so a failing fit can still be inspected locally.
fit$save_object(file = "output/fit_production.rds")

cat("\n=== m13 production fit ===\n")
print(fit$summary(variables = c("mu_bar", "rho", "tau", "mu_g", "phi_g", "omega_g",
                                "rho_EI", "psi")))
cat(sprintf("Wall time: %.0f s\n", fit$time()$total))

# ---------------------------------------------------------------------------
# Convergence gate. Thresholds are fixed targets, not tuned to the result.
#   divergences          < 0.1% of post-warmup draws
#   min E-BFMI (chains)  >= 0.3   (0.2 is the hard floor at which Stan itself
#                                  warns; anything below 0.3 already fails)
#   max R-hat            <= 1.01 on hyperparameters, <= 1.05 on every mu[t], h[t,k], g[t]
#   min bulk & tail ESS  >= 400  on hyperparameters and on every mu[t]
# Chain and draw counts are read from the fit, never hard-coded.
# ---------------------------------------------------------------------------
GATE <- list(
  max_div_frac   = 0.001,
  min_ebfmi      = 0.3,
  ebfmi_floor    = 0.2,
  max_rhat_hyper = 1.01,
  max_rhat_state = 1.05,
  min_ess        = 400
)
# tau = exp(mu_g / 2) is a generated quantity in m13 (the latent innovation sd at the
# long-run log-variance); mu_g, phi_g, omega_g are the latent SV hyperparameters.
hyper_vars <- c("c", "rho", "tau", "mu_h", "phi", "omega", "rho_EI", "psi", "mu_bar",
                "mu_g", "phi_g", "omega_g")

n_chains <- fit$num_chains()
sdiag    <- fit$sampler_diagnostics(format = "draws_df")
n_draws  <- nrow(sdiag)                       # post-warmup draws, all chains
n_div    <- sum(sdiag$divergent__)
n_tdmax  <- sum(sdiag$treedepth__ >= fit$metadata()$max_treedepth)
ebfmi    <- fit$diagnostic_summary(diagnostics = "ebfmi", quiet = TRUE)$ebfmi

summ_conv <- function(vars) {
  posterior::summarise_draws(
    fit$draws(variables = vars),
    "rhat", "ess_bulk", "ess_tail"
  )
}
s_hyper <- summ_conv(hyper_vars)
s_mu    <- summ_conv("mu")
s_h     <- summ_conv("h")
s_g     <- summ_conv("g")

stopifnot(nrow(s_mu) == stan_data$T, nrow(s_h) == stan_data$T * stan_data$K,
          nrow(s_g) == stan_data$T)

# NA diagnostics (e.g. a constant draw) count as failures, never as passes.
worst_max <- function(x) if (anyNA(x)) NA_real_ else max(x)
worst_min <- function(x) if (anyNA(x)) NA_real_ else min(x)
which_at  <- function(s, col, fn) {
  v <- s[[col]]
  if (anyNA(v)) return(paste(s$variable[is.na(v)][1], "(NA)"))
  s$variable[fn(v)]
}

div_frac <- n_div / n_draws
min_eb   <- worst_min(ebfmi)

checks <- data.frame(
  metric = c(
    sprintf("divergences (%d of %d draws, %d chains)", n_div, n_draws, n_chains),
    "min E-BFMI over chains",
    "max R-hat, hyperparameters",
    "max R-hat, mu[t]",
    "max R-hat, h[t,k]",
    "max R-hat, g[t]",
    "min bulk ESS, hyperparameters",
    "min tail ESS, hyperparameters",
    "min bulk ESS, mu[t]",
    "min tail ESS, mu[t]"
  ),
  value = c(
    div_frac,
    min_eb,
    worst_max(s_hyper$rhat),
    worst_max(s_mu$rhat),
    worst_max(s_h$rhat),
    worst_max(s_g$rhat),
    worst_min(s_hyper$ess_bulk),
    worst_min(s_hyper$ess_tail),
    worst_min(s_mu$ess_bulk),
    worst_min(s_mu$ess_tail)
  ),
  op = c("<", ">=", "<=", "<=", "<=", "<=", ">=", ">=", ">=", ">="),
  threshold = c(
    GATE$max_div_frac, GATE$min_ebfmi,
    GATE$max_rhat_hyper, GATE$max_rhat_state, GATE$max_rhat_state, GATE$max_rhat_state,
    GATE$min_ess, GATE$min_ess, GATE$min_ess, GATE$min_ess
  ),
  worst = c(
    "",
    sprintf("chain %d", if (anyNA(ebfmi)) which(is.na(ebfmi))[1] else which.min(ebfmi)),
    which_at(s_hyper, "rhat", which.max),
    which_at(s_mu, "rhat", which.max),
    which_at(s_h, "rhat", which.max),
    which_at(s_g, "rhat", which.max),
    which_at(s_hyper, "ess_bulk", which.min),
    which_at(s_hyper, "ess_tail", which.min),
    which_at(s_mu, "ess_bulk", which.min),
    which_at(s_mu, "ess_tail", which.min)
  ),
  stringsAsFactors = FALSE
)
checks$pass <- mapply(function(v, op, th) {
  if (is.na(v)) return(FALSE)
  switch(op, "<" = v < th, "<=" = v <= th, ">=" = v >= th)
}, checks$value, checks$op, checks$threshold)

fmt_val <- function(v, i) {
  if (is.na(v)) return("NA")
  if (i == 1) return(sprintf("%.3f%%", 100 * v))
  if (checks$op[i] == "<=") return(sprintf("%.4f", v))
  if (grepl("ESS", checks$metric[i])) return(sprintf("%.0f", v))
  sprintf("%.3f", v)
}
fmt_th <- function(i) {
  th <- checks$threshold[i]
  if (i == 1) return(sprintf("%s %.1f%%", checks$op[i], 100 * th))
  sprintf("%s %s", checks$op[i], format(th))
}

cat("\n=== Convergence gate ===\n")
cat(sprintf("Per-chain E-BFMI: %s\n", paste(sprintf("%.3f", ebfmi), collapse = ", ")))
cat(sprintf("Max-treedepth hits: %d of %d draws (informational)\n", n_tdmax, n_draws))
cat(sprintf("%-46s %10s  %-10s %-6s %s\n", "metric", "value", "threshold", "status", "worst"))
for (i in seq_len(nrow(checks))) {
  cat(sprintf("%-46s %10s  %-10s %-6s %s\n",
              checks$metric[i], fmt_val(checks$value[i], i), fmt_th(i),
              if (checks$pass[i]) "PASS" else "FAIL", checks$worst[i]))
}

failures <- character(0)
for (i in which(!checks$pass)) {
  failures <- c(failures, sprintf(
    "%s = %s (required %s; worst: %s)",
    checks$metric[i], fmt_val(checks$value[i], i), fmt_th(i), checks$worst[i]
  ))
}
if (!is.na(min_eb) && min_eb < GATE$ebfmi_floor) {
  failures <- c(failures, sprintf(
    "min E-BFMI %.3f is below the hard floor %.1f", min_eb, GATE$ebfmi_floor))
}

if (length(failures) > 0) {
  stop(paste0(
    "Fit quality gate FAILED (", length(failures), " problem",
    if (length(failures) > 1) "s" else "", "):\n  - ",
    paste(failures, collapse = "\n  - ")
  ), call. = FALSE)
}
cat("\nFit quality gate passed: all checks within thresholds.\n")
