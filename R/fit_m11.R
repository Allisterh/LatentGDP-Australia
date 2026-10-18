# Fit m11 (stan/m11_mixed_freq.stan: monthly latent growth, quarterly GDP(E, I, P) and
# monthly employment) and save it to <root>/fit_m11.rds.
#
#   Rscript R/fit_m11.R
#   LATENTGDP_START=1980-03-01 LATENTGDP_OUT=output/sample_1980 Rscript R/fit_m11.R
#
# Sample window and output root: R/sample_window.R. With LATENTGDP_START set, the monthly
# employment series starts in the first month of the start quarter (1980-01 for
# 1980-03-01) and the quarterly series in the start quarter, so the first quarter is
# complete; the alignment below keeps only quarters whose three months are all present.
# Unset = the full sample into output/, exactly as before. The diagnostics row for
# <root>/tables/fit_diagnostics.csv: Rscript R/02_fit_models.R --diag-only m11.

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
})

source("R/sample_window.R")
cat(sw_describe(), "\n")
emp_growth <- sw_trim_m(readRDS("data/processed/employment_growth_monthly.rds"))
gdp_growth <- sw_trim_q(readRDS("data/processed/gdp_growth.rds"))

# Align samples. Quarterly GDP date convention: first day of last month of quarter.
# Monthly employment: first day of month. We need each quarter's 3 months to be
# present in the monthly series, and we keep only quarters where that holds.

emp_growth <- emp_growth |>
  filter(date <= max(gdp_growth$date))   # trim post-quarterly months

# Quarter q_end_m is the row in emp_growth matching the quarterly date.
# We need q_end_m >= 3 (so months m-2, m-1, m exist).
m_dates <- emp_growth$date
candidate <- gdp_growth |>
  mutate(q_end_m = match(date, m_dates)) |>
  filter(!is.na(q_end_m), q_end_m >= 3)

gdp_growth_m11 <- candidate |> select(-q_end_m)
q_end_m <- candidate$q_end_m

stan_data_m11 <- list(
  T_m     = length(m_dates),
  T_q     = nrow(gdp_growth_m11),
  q_end_m = q_end_m,
  y_q     = as.matrix(gdp_growth_m11[, c("g_E", "g_I", "g_P")]),
  y_emp   = emp_growth$g_emp
)

cat(sprintf("m11 data: T_m = %d months (%s to %s), T_q = %d quarters\n",
            stan_data_m11$T_m, min(m_dates), max(m_dates), stan_data_m11$T_q))

mod <- cmdstan_model("stan/m11_mixed_freq.stan")
fit <- mod$sample(
  data = stan_data_m11, chains = 4, parallel_chains = 4,
  iter_warmup = 2000, iter_sampling = 2000, seed = 1, refresh = 1000,
  adapt_delta = 0.99, max_treedepth = 12
)
dir.create(sw_out_root(), showWarnings = FALSE, recursive = TRUE)
fit$save_object(sw_out("fit_m11.rds"))

print(fit$summary(variables = c(
  "mu_bar_m", "c", "rho", "tau", "sigma_q",
  "alpha_emp", "beta_emp", "sigma_emp"
)))
print(fit$diagnostic_summary())
