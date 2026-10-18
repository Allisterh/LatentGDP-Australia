# Writes machine-readable outputs from the production fit:
#   <docs>/estimates.csv             - full mu_t time series with 90% CI
#   <docs>/estimates_latest.json     - latest quarter headline numbers
#   <docs>/diagnostics.json          - run metadata, sampler / convergence
#                                      diagnostics and the input fingerprint
#   <docs>/model_summary.json        - posterior summaries quoted on the site (latent
#                                      volatility by decade, psi), annualised units
# These are the stable contract for downstream consumers of the estimate.
#
# Environment overrides (for smoke tests that must not touch docs/):
#   LATENTGDP_FIT   path to the saved CmdStanMCMC fit (default output/fit_production.rds)
#   LATENTGDP_DOCS  output directory                   (default docs)

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
  library(jsonlite)
})

source("R/input_hash.R")

fit_path <- Sys.getenv("LATENTGDP_FIT",  "output/fit_production.rds")
docs_dir <- Sys.getenv("LATENTGDP_DOCS", "docs")

fit        <- readRDS(fit_path)
gdp_growth <- readRDS("data/processed/gdp_growth.rds")

mu_draws <- fit$draws("mu", format = "draws_matrix")
stopifnot(ncol(mu_draws) == nrow(gdp_growth))

# Display convention: report quarter-on-quarter growth (100 x dlog), the
# standard Australian presentation. The model is fit on annualised growth
# (400 x dlog); QoQ is exactly that / 4, so rescaling the posterior draws and
# the observed series here is an exact transform of the same fit - no re-fit.
# Everything downstream (headline, CIs, latest JSON, variance ratio) derives
# from these, so it all carries the QoQ scale automatically.
SCALE <- 0.25
mu_draws       <- mu_draws * SCALE
gdp_growth$g_E <- gdp_growth$g_E * SCALE
gdp_growth$g_I <- gdp_growth$g_I * SCALE
gdp_growth$g_P <- gdp_growth$g_P * SCALE

estimates <- tibble(
  date          = gdp_growth$date,
  g_E           = gdp_growth$g_E,
  g_I           = gdp_growth$g_I,
  g_P           = gdp_growth$g_P,
  headline_mean = (gdp_growth$g_E + gdp_growth$g_I + gdp_growth$g_P) / 3,
  mu_median     = apply(mu_draws, 2, median),
  mu_q05        = apply(mu_draws, 2, quantile, 0.05),
  mu_q95        = apply(mu_draws, 2, quantile, 0.95),
  mu_q25        = apply(mu_draws, 2, quantile, 0.25),
  mu_q75        = apply(mu_draws, 2, quantile, 0.75)
)

dir.create(docs_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(estimates, file.path(docs_dir, "estimates.csv"), row.names = FALSE)

# Variance of the latent path relative to the headline (mean of E, I, P).
# Computed per posterior draw: var(mu^(s)) / var(headline). The variance of the
# posterior-median path would understate this, because taking the pointwise
# median shrinks the path towards its mean.
headline_var   <- var(estimates$headline_mean)
var_ratio_draw <- apply(unclass(mu_draws), 1, var) / headline_var
var_ratio_q    <- quantile(var_ratio_draw, c(0.05, 0.5, 0.95), names = FALSE)

# Latest-quarter headline JSON (the thing the site shows up top).
latest <- tail(estimates, 1)
# Build "YYYY Qn" by hand: format(.., "%q") is not supported by strftime.
latest_month <- as.integer(format(as.Date(latest$date), "%m"))
latest_year  <- as.integer(format(as.Date(latest$date), "%Y"))

latest_json <- list(
  quarter           = sprintf("%d Q%d", latest_year, ((latest_month - 1) %/% 3) + 1),
  date              = as.character(latest$date),
  mu_median         = round(latest$mu_median, 2),
  mu_q05            = round(latest$mu_q05, 2),
  mu_q95            = round(latest$mu_q95, 2),
  headline_mean     = round(latest$headline_mean, 2),
  observed = list(
    g_E = round(latest$g_E, 2),
    g_I = round(latest$g_I, 2),
    g_P = round(latest$g_P, 2)
  ),
  # Posterior median of the per-draw ratio (key kept for backward
  # compatibility) plus its 90% credible interval. Four decimals, so the
  # site's whole-percent display is not a double rounding (0.78529 -> 0.785
  # -> "78%" instead of "79%").
  variance_ratio_vs_headline     = round(var_ratio_q[2], 4),
  variance_ratio_vs_headline_q05 = round(var_ratio_q[1], 4),
  variance_ratio_vs_headline_q95 = round(var_ratio_q[3], 4),
  units = "quarter-on-quarter growth, % (100 x dlog of chain volume measures)",
  model = "m13: structural Sigma (Cov(eps_P, .) = 0) + stochastic volatility + AR(1) measurement errors + stochastic volatility on latent growth",
  data_source = "ABS 5206.0 Table 24, seasonally adjusted, chain volume measures"
)

write_json(latest_json, file.path(docs_dir, "estimates_latest.json"),
           pretty = TRUE, auto_unbox = TRUE, digits = 4)

# --- Model summary ----------------------------------------------------------
# Posterior summaries that the site's model description quotes, in the model's
# annualised units (no SCALE): the latent innovation sd tau_t = exp(g_t / 2), as the
# posterior median of each decade's average, and at 2020 Q2; and the median psi.
tau_t_draws <- unclass(posterior::as_draws_matrix(fit$draws("tau_t")))
psi_draws   <- unclass(posterior::as_draws_matrix(fit$draws("psi")))
stopifnot(ncol(tau_t_draws) == nrow(gdp_growth), ncol(psi_draws) == 3)
yr      <- as.integer(format(as.Date(gdp_growth$date), "%Y"))
decades <- sort(unique(yr[yr >= 1960] %/% 10 * 10))
decade_avg <- setNames(lapply(decades, function(d) {
  cols <- which(yr %/% 10 * 10 == d)
  round(median(rowMeans(tau_t_draws[, cols, drop = FALSE])), 4)
}), paste0(decades, "s"))
i_2020q2 <- which(as.Date(gdp_growth$date) == as.Date("2020-06-01"))
stopifnot(length(i_2020q2) == 1)
model_summary <- list(
  units = "annualised growth, % (400 x dlog)",
  latent_sd_decade_avg = decade_avg,
  latent_sd_2020q2 = round(median(tau_t_draws[, i_2020q2]), 4),
  psi = list(E = round(median(psi_draws[, 1]), 4), I = round(median(psi_draws[, 2]), 4),
             P = round(median(psi_draws[, 3]), 4))
)
write_json(model_summary, file.path(docs_dir, "model_summary.json"),
           pretty = TRUE, auto_unbox = TRUE, digits = NA)

# --- Diagnostics ------------------------------------------------------------
# Sampler settings come from the fit itself, not hard-coded constants.
meta            <- fit$metadata()
n_chains        <- fit$num_chains()
draws_per_chain <- posterior::niterations(fit$draws("lp__"))
total_draws     <- n_chains * draws_per_chain

diag <- fit$diagnostic_summary(quiet = TRUE)

# Convergence metrics over the hyperparameters and the latent path mu.
# Hyperparameter names are read defensively: candidates absent from the fit
# are skipped, so a renamed / reparameterised model still yields diagnostics for
# whatever exists (and "hyperparameters" in the JSON records what was used).
# psi (AR(1) measurement-error coefficients) is absent from an m9 fit; mu_g, phi_g and
# omega_g (latent-growth SV) are m13's and absent from an m12 fit. In m13, tau is the
# generated quantity exp(mu_g / 2). The set matches hyper_vars in R/fit_production.R (the
# gated quantities), including the derived mu_bar = c / (1 - rho).
HYPER_CANDIDATES <- c("c", "rho", "tau", "rho_EI", "mu_h", "phi", "omega", "psi", "mu_bar",
                      "mu_g", "phi_g", "omega_g")
stan_vars  <- fit$metadata()$stan_variables
hyper_used <- intersect(HYPER_CANDIDATES, stan_vars)
if (length(hyper_used) == 0) {
  warning("None of the expected production hyperparameters (",
          paste(HYPER_CANDIDATES, collapse = ", "),
          ") are in the fit; R-hat/ESS reported for mu only.", call. = FALSE)
}

conv_fun <- function(vars) {
  if (length(vars) == 0) return(NULL)
  posterior::summarise_draws(
    fit$draws(variables = vars),
    rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk,
    ess_tail = posterior::ess_tail
  )
}
conv_h  <- conv_fun(hyper_used)
conv_mu <- conv_fun("mu")
conv    <- rbind(conv_h, conv_mu)

# min / max that return NA (JSON null) rather than +-Inf for empty input.
smax <- function(x) { x <- x[is.finite(x)]; if (length(x) == 0) NA else max(x) }
smin <- function(x) { x <- x[is.finite(x)]; if (length(x) == 0) NA else min(x) }

diag_json <- list(
  run_time_utc          = format(Sys.time(), tz = "UTC"),
  input_hash            = compute_input_hash(),
  cmdstan_version       = cmdstan_version(),
  T                     = nrow(gdp_growth),
  sample_start          = as.character(min(estimates$date)),
  sample_end            = as.character(max(estimates$date)),
  chains                = n_chains,
  iter_warmup           = meta$iter_warmup,
  iter_sampling         = meta$iter_sampling,
  thin                  = meta$thin,
  draws_per_chain       = draws_per_chain,
  iter_per_chain        = draws_per_chain,   # legacy key: post-warmup draws per chain
  total_draws           = total_draws,
  adapt_delta           = meta$adapt_delta,
  max_treedepth         = meta$max_treedepth,
  num_divergent         = sum(diag$num_divergent),
  divergent_pct         = round(100 * sum(diag$num_divergent) / total_draws, 3),
  num_max_treedepth     = sum(diag$num_max_treedepth),
  ebfmi_min             = round(smin(diag$ebfmi), 3),
  hyperparameters       = I(hyper_used),
  rhat_max_hyper        = round(smax(conv_h$rhat), 4),
  rhat_max_mu           = round(smax(conv_mu$rhat), 4),
  ess_bulk_min          = round(smin(conv$ess_bulk), 0),
  ess_tail_min          = round(smin(conv$ess_tail), 0),
  ess_bulk_min_hyper    = round(smin(conv_h$ess_bulk), 0),
  ess_tail_min_hyper    = round(smin(conv_h$ess_tail), 0)
)
write_json(diag_json, file.path(docs_dir, "diagnostics.json"),
           pretty = TRUE, auto_unbox = TRUE, digits = NA)

cat(sprintf(
  "\nLatest estimate: %s mu = %.2f%% [%.2f, %.2f] (90%% CI)\n",
  latest_json$quarter, latest_json$mu_median,
  latest_json$mu_q05, latest_json$mu_q95
))
cat(sprintf("Headline (mean of E,I,P): %.2f%%\n", latest_json$headline_mean))
cat(sprintf("Variance ratio (mu vs headline): %.3f [%.3f, %.3f] (median, 90%% CI over draws)\n",
            var_ratio_q[2], var_ratio_q[1], var_ratio_q[3]))
cat(sprintf(
  "Diagnostics: %d chains x %d draws; %d divergent; min E-BFMI %.3f; max R-hat (hyper) %.4f; min ESS bulk/tail %s/%s\n",
  n_chains, draws_per_chain, diag_json$num_divergent, diag_json$ebfmi_min,
  diag_json$rhat_max_hyper, diag_json$ess_bulk_min, diag_json$ess_tail_min
))
cat(sprintf("Input hash: %s\n", diag_json$input_hash))
