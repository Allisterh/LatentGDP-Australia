# Research suite: fit m1-m10, m12 and m13 on the current data, one model at a time.
#   Rscript R/02_fit_models.R                  # all of m1..m10, m12, m13
#   Rscript R/02_fit_models.R m9 m10           # a subset
#   Rscript R/02_fit_models.R --diag-only m12  # no fitting: tabulate the diagnostics of
#                                              # the saved output/fit_<m>.rds
# m11 (mixed frequency) is fitted by R/fit_m11.R. R/fit_m12.R fits m12 with the same
# settings as here and also applies the production gate (and has a warmup-stall harness).
# R/fit_m13.R does the same for m13 (m12 + SV on the latent innovation; tau is the
# generated quantity exp(mu_g / 2)); `--diag-only m13` tabulates its saved fit.
#
# Writes output/fit_m{1..10,12,13}.rds (CmdStanFit objects via save_object) and
# output/tables/fit_diagnostics.csv (one row per model; rows for models not
# refitted in this run are kept).
#
# Samples:
#   m1-m7, m9, m12, m13 : data/processed/gdp_growth.rds             (3 measures, full sample)
#   m8, m10        : data/processed/gdp_growth_with_hours.rds  (3 measures + hours,
#                    post-1978 subsample; written by R/01b_fetch_research_data.R)
#   Both are cut to start at LATENTGDP_START when it is set (see below); with 1980-03-01
#   all of them, hours included, run 1980 Q1 - latest (T = 186 at 2026 Q2).
#
# Sampler settings: 4 chains x (2000 warmup + 2000 draws), seed 1. m9 uses the
# production settings of R/fit_production.R (adapt_delta 0.99, max_treedepth 12);
# the other SV models (m4, m10, m12, m13) use the same settings so that they are fitted
# like-for-like. Diagnostics are printed and tabulated, never used to stop the
# run: this is the research suite, and a model that samples badly is reported as
# such rather than silently dropped. Parallel chains: --parallel=K (default 4).
# m10 is the re-specified model with stochastic volatility on the latent innovation
# as well (tau is the generated quantity exp(mu_g / 2)); see the header of
# stan/m10_labour_sv.stan.
#
# Sample window and output root (R/sample_window.R): LATENTGDP_START=YYYY-MM-DD sets the
# first growth observation of every sample (3-measure and hours alike); LATENTGDP_OUT sets
# the output root (default output/). Both unset = the full sample into output/, exactly as
# before. Example, the RBA RDP 2014-12 sample start:
#   LATENTGDP_START=1980-03-01 LATENTGDP_OUT=output/sample_1980 Rscript R/02_fit_models.R
# writes <root>/fit_m*.rds and <root>/tables/fit_diagnostics.csv.
#
# Extra modes (no fitting):
#   Rscript R/02_fit_models.R --diag-only m11   # m11's row, from <root>/fit_m11.rds
#                                               # (fitted by R/fit_m11.R)
#   Rscript R/02_fit_models.R --vardecomp       # per-draw variance decomposition of the
#       saved fits -> <root>/tables/variance_decomposition.csv (all quarters and excluding
#       2020 Q1 - 2021 Q1). Refuses to write output/tables/variance_decomposition.csv,
#       which R/03_compare.R owns.

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
})

source("R/sample_window.R")
cat(sw_describe(), "\n")
dir.create(sw_out("tables"), showWarnings = FALSE, recursive = TRUE)

gdp_growth <- sw_trim_q(readRDS("data/processed/gdp_growth.rds"))
hours_path <- "data/processed/gdp_growth_with_hours.rds"
if (!file.exists(hours_path)) {
  stop(hours_path, " not found - run Rscript R/01b_fetch_research_data.R first.",
       call. = FALSE)
}
gdp_growth_4 <- sw_trim_q(readRDS(hours_path))
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
             params = c("mu_bar", "c", "rho", "tau", "mu_h", "phi", "omega", "rho_EI", "psi")),
  m13 = list(file = "stan/m13_serial_latent_sv.stan", data = stan_data,  extras = sv_settings,
             params = c("mu_bar", "c", "rho", "tau", "mu_g", "phi_g", "omega_g", "mu_h", "phi",
                        "omega", "rho_EI", "psi"))
)

# m11 (fitted by R/fit_m11.R) can only be tabulated from its saved fit (--diag-only m11);
# it is not in `specs`, so it is never fitted here. Its latent state is the quarterly
# aggregate mu_q, and T is its number of quarters T_q.
diag_extra <- list(
  m11 = list(params = c("mu_bar_m", "c", "rho", "tau", "sigma_q", "alpha_emp", "beta_emp",
                        "sigma_emp"), mu_var = "mu_q")
)

diag_row <- function(name, fit, spec, secs) {
  sd   <- fit$sampler_diagnostics(format = "draws_df")
  eb   <- fit$diagnostic_summary(diagnostics = "ebfmi", quiet = TRUE)$ebfmi
  s_hy <- summarise_draws(fit$draws(spec$params), "rhat", "ess_bulk", "ess_tail")
  mu_var <- if (is.null(spec$mu_var)) "mu" else spec$mu_var
  s_mu <- summarise_draws(fit$draws(mu_var), "rhat", "ess_bulk")
  data.frame(
    model           = name,
    T               = if (is.null(spec$data)) nrow(s_mu) else spec$data$T,
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
  fit$save_object(file = sw_out(sprintf("fit_%s.rds", name)))
  print(fit$summary(variables = spec$params), n = Inf)
  dr <- diag_row(name, fit, spec, secs)
  print(dr)
  rm(fit); gc(verbose = FALSE)
  dr
}

diag_saved <- function(name, spec) {
  path <- sw_out(sprintf("fit_%s.rds", name))
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

# ---- variance decomposition of the saved fits (--vardecomp) ---------------------------
# ratio_s = var_t(mu_t^(s)) / var_t(headline_t), headline = mean of E, I, P, per posterior
# draw s; the table gives the median and the 5% / 95% quantiles (90% interval), on all
# quarters of each model's sample and on the same quarters without 2020 Q1 - 2021 Q1.
# Same columns as R/03_compare.R's table.
COVID_EXCL <- as.Date(c("2020-03-01", "2020-06-01", "2020-09-01", "2020-12-01", "2021-03-01"))
vardecomp <- function() {
  out_csv <- sw_out("tables", "variance_decomposition.csv")
  if (normalizePath(dirname(out_csv), mustWork = FALSE) ==
      normalizePath("output/tables", mustWork = FALSE))
    stop("--vardecomp does not write output/tables/variance_decomposition.csv (R/03_compare.R ",
         "owns it); set LATENTGDP_OUT to another root", call. = FALSE)
  yr <- function(d) format(d, "%Y")
  qlab <- function(d) sprintf("%sQ%d", yr(d), (as.integer(format(d, "%m")) + 2) %/% 3)
  lab <- function(kind, dates, excl) sprintf("%s %s-%s%s", kind, qlab(min(dates)), qlab(max(dates)),
                                             if (excl) " excl. 2020Q1-2021Q1" else "")
  vrow <- function(series, sample, T, v, hv)
    tibble(series = series, sample = sample, T = T, sd = sqrt(v), var = v,
           ratio_median = v / hv, ratio_q05 = v / hv, ratio_q95 = v / hv)
  draw_ratio <- function(path, vname, keep, hv, label, sample) {
    if (!file.exists(path)) { cat("  (", path, " not found: row skipped)\n", sep = ""); return(NULL) }
    d <- unclass(as.matrix(readRDS(path)$draws(vname, format = "draws_matrix")))
    if (ncol(d) != length(keep)) stop(path, ": ", ncol(d), " ", vname, " columns but the data give ",
                                      length(keep), " quarters", call. = FALSE)
    vs <- apply(d[, keep, drop = FALSE], 1, stats::var)
    r <- vs / hv
    tibble(series = label, sample = sample, T = sum(keep), sd = median(sqrt(vs)),
           var = median(vs), ratio_median = median(r),
           ratio_q05 = unname(quantile(r, 0.05)), ratio_q95 = unname(quantile(r, 0.95)))
  }
  block <- function(dates, y, kind, models, vname, path_fun, measures) {
    bind_rows(lapply(c(FALSE, TRUE), function(excl) {
      keep <- if (excl) !(dates %in% COVID_EXCL) else rep(TRUE, length(dates))
      hl <- rowMeans(y[keep, 1:3, drop = FALSE]); hv <- var(hl); smp <- lab(kind, dates, excl)
      bind_rows(
        vrow("Headline (mean of E, I, P)", smp, sum(keep), hv, hv),
        if (measures) bind_rows(lapply(1:3, function(k)
          vrow(c("GDP(E)", "GDP(I)", "GDP(P)")[k], smp, sum(keep), var(y[keep, k]), hv))),
        bind_rows(lapply(models, function(m)
          draw_ratio(path_fun(m), vname, keep, hv, sprintf("Latent %s (%s)", vname, m), smp))))
    }))
  }
  y3 <- as.matrix(gdp_growth[, c("g_E", "g_I", "g_P")])
  y4 <- as.matrix(gdp_growth_4[, c("g_E", "g_I", "g_P")])
  fp <- function(m) sw_out(sprintf("fit_%s.rds", m))
  vt <- bind_rows(
    block(gdp_growth$date, y3, "3-measure sample", c("m1", "m2", "m3", "m4", "m5", "m6", "m7",
                                                     "m9", "m12", "m13"), "mu", fp, TRUE),
    block(gdp_growth_4$date, y4, "hours sample", c("m8", "m10"), "mu", fp, FALSE))
  # m11: quarterly aggregate mu_q against headline on m11's quarters (alignment of R/fit_m11.R)
  em <- sw_trim_m(readRDS("data/processed/employment_growth_monthly.rds"))
  em <- em[em$date <= max(gdp_growth$date), ]
  qm <- match(gdp_growth$date, em$date)
  q_ok <- which(!is.na(qm) & qm >= 3)
  vt <- bind_rows(vt, block(gdp_growth$date[q_ok], y3[q_ok, , drop = FALSE], "mixed-frequency quarters",
                            "m11", "mu_q", fp, FALSE))
  cat("\n=== Variance of growth series relative to headline (per-draw ratios) ===\n")
  print(as.data.frame(vt), digits = 3)
  write.csv(vt, out_csv, row.names = FALSE)
  cat("Written to ", out_csv, "\n", sep = "")
}

args_in <- commandArgs(trailingOnly = TRUE)
if ("--vardecomp" %in% args_in) { vardecomp(); quit(save = "no") }
diag_only <- "--diag-only" %in% args_in
# --skip-existing: resume an interrupted run (models whose <root>/fit_<m>.rds exists and
# that have a row in <root>/tables/fit_diagnostics.csv are not refitted).
skip_existing <- "--skip-existing" %in% args_in
par_arg <- grep("^--parallel=", args_in, value = TRUE)
PARALLEL <- if (length(par_arg)) as.integer(sub("^--parallel=", "", par_arg)) else 4L
args_in <- setdiff(args_in, c("--diag-only", "--skip-existing", par_arg))
to_fit <- if (length(args_in)) args_in else names(specs)
known <- if (diag_only) c(names(specs), names(diag_extra)) else names(specs)
bad <- setdiff(to_fit, known)
if (length(bad)) stop("unknown model(s): ", paste(bad, collapse = ", "), call. = FALSE)

diag_path <- sw_out("tables", "fit_diagnostics.csv")
for (nm in to_fit) {
  if (skip_existing && !diag_only && file.exists(sw_out(sprintf("fit_%s.rds", nm))) &&
      file.exists(diag_path) && nm %in% read.csv(diag_path, stringsAsFactors = FALSE)$model) {
    cat("\n===== ", nm, ": fit and diagnostics row exist, skipped (--skip-existing) =====\n", sep = "")
    next
  }
  spec <- if (nm %in% names(specs)) specs[[nm]] else diag_extra[[nm]]
  dr <- if (diag_only) diag_saved(nm, spec) else fit_one(nm, spec)
  # Merge after every model so an interrupted run keeps what it finished.
  old <- if (file.exists(diag_path)) read.csv(diag_path, stringsAsFactors = FALSE) else NULL
  if (!is.null(old)) old <- old[old$model != nm, , drop = FALSE]
  all_d <- bind_rows(old, dr)
  all_d <- all_d[order(as.integer(sub("m", "", all_d$model))), , drop = FALSE]
  write.csv(all_d, diag_path, row.names = FALSE)
}
cat("\nDiagnostics written to ", diag_path, "\n", sep = "")
print(read.csv(diag_path))
