# Research suite: fit m1-m10 and m12 on the current data, one model at a time.
#   Rscript R/02_fit_models.R                  # all of m1..m10, m12
#   Rscript R/02_fit_models.R m9 m10           # a subset
#   Rscript R/02_fit_models.R --diag-only m12  # no fitting: tabulate the diagnostics of
#                                              # the saved output/fit_<m>.rds
# m11 (mixed frequency) is fitted by R/fit_m11.R. R/fit_m12.R fits m12 with the same
# settings as here and also applies the production gate (and has a warmup-stall harness).
#
# Writes output/fit_m{1..10,12}.rds (CmdStanFit objects via save_object) and
# output/tables/fit_diagnostics.csv (one row per model; rows for models not
# refitted in this run are kept).
#
# Samples:
#   m1-m7, m9, m12 : data/processed/gdp_growth.rds             (3 measures, full sample)
#   m8, m10        : data/processed/gdp_growth_with_hours.rds  (3 measures + hours,
#                    post-1978 subsample; written by R/01b_fetch_research_data.R)
#
# Sampler settings: 4 chains x (2000 warmup + 2000 draws), seed 1. m9 uses the
# production settings of R/fit_production.R (adapt_delta 0.99, max_treedepth 12);
# the other SV models (m4, m10, m12) use the same settings so that they are fitted
# like-for-like. Diagnostics are printed and tabulated, never used to stop the
# run: this is the research suite, and a model that samples badly is reported as
# such rather than silently dropped. Parallel chains: --parallel=K (default 4).
# m10 is the re-specified model with stochastic volatility on the latent innovation
# as well (tau is the generated quantity exp(mu_g / 2)); see the header of
# stan/m10_labour_sv.stan.

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
})

dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)

gdp_growth <- readRDS("data/processed/gdp_growth.rds")
hours_path <- "data/processed/gdp_growth_with_hours.rds"
if (!file.exists(hours_path)) {
  stop(hours_path, " not found - run Rscript R/01b_fetch_research_data.R first.",
       call. = FALSE)
}
gdp_growth_4 <- readRDS(hours_path)
stopifnot(!anyNA(gdp_growth[, c("g_E", "g_I", "g_P")]),
          !anyNA(gdp_growth_4[, c("g_E", "g_I", "g_P", "g_U")]))

# Standard 3-measure data block (m1-m7, m9)
stan_data <- list(
  T = nrow(gdp_growth),
  K = 3,
  y = as.matrix(gdp_growth[, c("g_E", "g_I", "g_P")])
)

# m7 needs the COVID quarter indices.
covid_dates <- as.Date(c("2020-06-01", "2020-09-01", "2020-12-01", "2021-03-01"))
covid_idx <- which(gdp_growth$date %in% covid_dates)
stopifnot(length(covid_idx) == length(covid_dates))
stan_data_m7 <- c(stan_data, list(n_covid = length(covid_idx), covid_idx = covid_idx))

# m8 / m10 use the 4-measure hours subsample.
stan_data_4 <- list(
  T = nrow(gdp_growth_4),
  y_gdp = as.matrix(gdp_growth_4[, c("g_E", "g_I", "g_P")]),
  y_U   = gdp_growth_4$g_U
)

sv_settings <- list(adapt_delta = 0.99, max_treedepth = 12)   # = R/fit_production.R

specs <- list(
  m1 = list(file = "stan/m1_baseline.stan",         data = stan_data,    extras = list(),
            params = c("mu_bar", "tau", "sigma")),
  m2 = list(file = "stan/m2_ar1.stan",              data = stan_data,    extras = list(),
            params = c("mu_bar", "c", "rho", "tau", "sigma")),
  m3 = list(file = "stan/m3_corr.stan",             data = stan_data,    extras = list(adapt_delta = 0.99),
            params = c("mu_bar", "c", "rho", "tau", "sigma", "Omega[2,1]", "Omega[3,1]", "Omega[3,2]")),
  m4 = list(file = "stan/m4_sv.stan",               data = stan_data,    extras = sv_settings,
            params = c("mu_bar", "c", "rho", "tau", "mu_h", "phi", "omega")),
  m5 = list(file = "stan/m5_serial.stan",           data = stan_data,    extras = list(adapt_delta = 0.99),
            params = c("mu_bar", "c", "rho", "tau", "sigma", "psi", "eps_sd")),
  m6 = list(file = "stan/m6_structural_corr.stan",  data = stan_data,    extras = list(adapt_delta = 0.99),
            params = c("mu_bar", "c", "rho", "tau", "sigma", "rho_EI")),
  m7 = list(file = "stan/m7_covid_scale.stan",      data = stan_data_m7, extras = list(adapt_delta = 0.99),
            params = c("mu_bar", "c", "rho", "tau", "sigma", "kappa")),
  m8 = list(file = "stan/m8_labour.stan",           data = stan_data_4,  extras = list(adapt_delta = 0.99),
            params = c("mu_bar", "c", "rho", "tau", "sigma_gdp", "sigma_U", "alpha", "beta")),
  m9 = list(file = "stan/m9_struct_sv.stan",        data = stan_data,    extras = sv_settings,
            params = c("mu_bar", "c", "rho", "tau", "mu_h", "phi", "omega", "rho_EI")),
  m10 = list(file = "stan/m10_labour_sv.stan",      data = stan_data_4,  extras = sv_settings,
             params = c("mu_bar", "c", "rho", "tau", "alpha", "beta", "mu_g", "phi_g", "omega_g",
                        "mu_h", "phi", "omega")),
  m12 = list(file = "stan/m12_struct_sv_serial.stan", data = stan_data,  extras = sv_settings,
             params = c("mu_bar", "c", "rho", "tau", "mu_h", "phi", "omega", "rho_EI", "psi"))
)

diag_row <- function(name, fit, spec, secs) {
  sd   <- fit$sampler_diagnostics(format = "draws_df")
  eb   <- fit$diagnostic_summary(diagnostics = "ebfmi", quiet = TRUE)$ebfmi
  s_hy <- summarise_draws(fit$draws(spec$params), "rhat", "ess_bulk", "ess_tail")
  s_mu <- summarise_draws(fit$draws("mu"), "rhat", "ess_bulk")
  data.frame(
    model           = name,
    T               = spec$data$T,
    chains          = fit$num_chains(),
    draws           = nrow(sd),
    adapt_delta     = fit$metadata()$adapt_delta,
    max_treedepth   = fit$metadata()$max_treedepth,
    divergences     = sum(sd$divergent__),
    treedepth_hits  = sum(sd$treedepth__ >= fit$metadata()$max_treedepth),
    min_ebfmi       = min(eb),
    max_rhat_hyper  = max(s_hy$rhat),
    worst_rhat_var  = s_hy$variable[which.max(s_hy$rhat)],
    min_ess_bulk_hyper = min(s_hy$ess_bulk),
    min_ess_tail_hyper = min(s_hy$ess_tail),
    max_rhat_mu     = max(s_mu$rhat),
    min_ess_bulk_mu = min(s_mu$ess_bulk),
    wall_secs       = round(secs),
    fitted_at       = format(Sys.time(), "%Y-%m-%d %H:%M"),
    stringsAsFactors = FALSE
  )
}

fit_one <- function(name, spec) {
  cat("\n===== ", name, " (T = ", spec$data$T, ") =====\n", sep = "")
  mod <- cmdstan_model(spec$file)
  args <- c(
    list(data = spec$data, chains = 4, parallel_chains = PARALLEL,
         iter_warmup = 2000, iter_sampling = 2000, seed = 1, refresh = 1000,
         show_exceptions = FALSE),
    spec$extras
  )
  t0 <- Sys.time()
  fit <- do.call(mod$sample, args)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  fit$save_object(file = sprintf("output/fit_%s.rds", name))
  print(fit$summary(variables = spec$params), n = Inf)
  dr <- diag_row(name, fit, spec, secs)
  print(dr)
  rm(fit); gc(verbose = FALSE)
  dr
}

diag_saved <- function(name, spec) {
  path <- sprintf("output/fit_%s.rds", name)
  if (!file.exists(path)) stop(path, " not found", call. = FALSE)
  fit <- readRDS(path)
  md <- fit$metadata()
  if (md$iter_sampling != 2000 || md$iter_warmup != 2000 || fit$num_chains() != 4 ||
      md$seed != 1)
    warning(name, ": saved fit does not use this script's settings (4 x 2000 + 2000, seed 1)",
            call. = FALSE)
  cat("\n===== ", name, " (saved fit ", path, ", ", format(file.mtime(path), "%Y-%m-%d %H:%M"),
      ") =====\n", sep = "")
  # Wall time of a saved fit: CmdStan's total sampling time recorded in the fit (it does
  # not include compilation; with parallel_chains < 4 it is longer than a 4-core run).
  dr <- diag_row(name, fit, spec, tryCatch(fit$time()$total, error = function(e) NA_real_))
  dr$fitted_at <- format(file.mtime(path), "%Y-%m-%d %H:%M")
  print(dr)
  dr
}

args_in <- commandArgs(trailingOnly = TRUE)
diag_only <- "--diag-only" %in% args_in
par_arg <- grep("^--parallel=", args_in, value = TRUE)
PARALLEL <- if (length(par_arg)) as.integer(sub("^--parallel=", "", par_arg)) else 4L
args_in <- setdiff(args_in, c("--diag-only", par_arg))
to_fit <- if (length(args_in)) args_in else names(specs)
bad <- setdiff(to_fit, names(specs))
if (length(bad)) stop("unknown model(s): ", paste(bad, collapse = ", "), call. = FALSE)

diag_path <- "output/tables/fit_diagnostics.csv"
for (nm in to_fit) {
  dr <- if (diag_only) diag_saved(nm, specs[[nm]]) else fit_one(nm, specs[[nm]])
  # Merge after every model so an interrupted run keeps what it finished.
  old <- if (file.exists(diag_path)) read.csv(diag_path, stringsAsFactors = FALSE) else NULL
  if (!is.null(old)) old <- old[old$model != nm, , drop = FALSE]
  all_d <- bind_rows(old, dr)
  all_d <- all_d[order(as.integer(sub("m", "", all_d$model))), , drop = FALSE]
  write.csv(all_d, diag_path, row.names = FALSE)
}
cat("\nDiagnostics written to ", diag_path, "\n", sep = "")
print(read.csv(diag_path))
