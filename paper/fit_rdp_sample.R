# paper/fit_rdp_sample.R
#
# Re-estimates selected specifications on the sample of RBA RDP 2014-12 (Rees, Lancaster &
# Finlay): first growth observation 1980 Q1, last 2013 Q2 (T = 134), on the current data
# vintage (data/processed/gdp_growth.rds, ABS 5206.0, June quarter 2026 release).
#
#   Rscript paper/fit_rdp_sample.R               # m2 m13 m5 m4 m9 m12, in that order
#   Rscript paper/fit_rdp_sample.R m2 m13        # a subset
#   Rscript paper/fit_rdp_sample.R --end=2019-12-01 --out=output/sample_rdp/end2019 m2 m13
#                                                # another sample end (1980 Q1 start kept)
#
# Settings are those of R/02_fit_models.R / R/fit_production.R: 4 chains x (2000 warmup +
# 2000 draws), seed 1, each model's own adapt_delta / max_treedepth. Fits are saved to
# output/sample_rdp/fit_<model>.rds (gitignored); one diagnostics row per model goes to
# output/sample_rdp/tables/fit_diagnostics.csv. Resumable: a model whose fit file exists is
# not refitted (its diagnostics row is recomputed from the saved fit).
#
# The truncation is done here, so R/sample_window.R and the research scripts are untouched.
# Nothing under output/ other than output/sample_rdp/ is written.

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
})
stopifnot(file.exists("R/production_config.R"))   # run from the repo root

args_all <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args_all, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[length(hit)]) else default
}
START <- as.Date("1980-03-01"); END <- as.Date(opt("end", "2013-06-01"))
OUT <- opt("out", "output/sample_rdp")
stopifnot(startsWith(normalizePath(OUT, mustWork = FALSE),
                     normalizePath("output/sample_rdp", mustWork = FALSE)))   # never elsewhere
dir.create(file.path(OUT, "tables"), showWarnings = FALSE, recursive = TRUE)

g <- readRDS("data/processed/gdp_growth.rds")
g <- g[g$date >= START & g$date <= END, ]
stopifnot(min(g$date) == START, max(g$date) == END, END != as.Date("2013-06-01") || nrow(g) == 134)
stan_data <- list(T = nrow(g), K = 3, y = as.matrix(g[, c("g_E", "g_I", "g_P")]))

sv <- list(adapt_delta = 0.99, max_treedepth = 12)
specs <- list(
  m2  = list(file = "stan/m2_ar1.stan", extras = list(),
             params = c("mu_bar", "c", "rho", "tau", "sigma")),
  m4  = list(file = "stan/m4_sv.stan", extras = sv,
             params = c("mu_bar", "c", "rho", "tau", "mu_h", "phi", "omega")),
  m5  = list(file = "stan/m5_serial.stan", extras = list(adapt_delta = 0.99),
             params = c("mu_bar", "c", "rho", "tau", "sigma", "psi", "eps_sd")),
  m9  = list(file = "stan/m9_struct_sv.stan", extras = sv,
             params = c("mu_bar", "c", "rho", "tau", "mu_h", "phi", "omega", "rho_EI")),
  m12 = list(file = "stan/m12_struct_sv_serial.stan", extras = sv,
             params = c("mu_bar", "c", "rho", "tau", "mu_h", "phi", "omega", "rho_EI", "psi")),
  m13 = list(file = "stan/m13_serial_latent_sv.stan", extras = sv,
             params = c("mu_bar", "c", "rho", "tau", "mu_g", "phi_g", "omega_g", "mu_h", "phi",
                        "omega", "rho_EI", "psi"))
)

args <- args_all[!startsWith(args_all, "--")]
todo <- if (length(args)) args else c("m2", "m13", "m5", "m4", "m9", "m12")
stopifnot(all(todo %in% names(specs)))

diag_row <- function(name, fit, spec, secs) {
  sd <- fit$sampler_diagnostics(format = "draws_df")
  eb <- fit$diagnostic_summary(diagnostics = "ebfmi", quiet = TRUE)$ebfmi
  s_hy <- summarise_draws(fit$draws(spec$params), "rhat", "ess_bulk", "ess_tail")
  s_mu <- summarise_draws(fit$draws("mu"), "rhat", "ess_bulk")
  data.frame(model = name, T = stan_data$T, first = format(START), last = format(END),
             chains = fit$num_chains(), draws = nrow(sd),
             adapt_delta = fit$metadata()$adapt_delta, max_treedepth = fit$metadata()$max_treedepth,
             divergences = sum(sd$divergent__),
             treedepth_hits = sum(sd$treedepth__ >= fit$metadata()$max_treedepth),
             min_ebfmi = min(eb), max_rhat_hyper = max(s_hy$rhat),
             worst_rhat_var = s_hy$variable[which.max(s_hy$rhat)],
             min_ess_bulk_hyper = min(s_hy$ess_bulk), min_ess_tail_hyper = min(s_hy$ess_tail),
             max_rhat_mu = max(s_mu$rhat), min_ess_bulk_mu = min(s_mu$ess_bulk),
             wall_secs = round(secs), fitted_at = format(Sys.time(), "%Y-%m-%d %H:%M"))
}

diag_csv <- file.path(OUT, "tables", "fit_diagnostics.csv")
for (nm in todo) {
  spec <- specs[[nm]]
  path <- file.path(OUT, sprintf("fit_%s.rds", nm))
  if (file.exists(path)) {
    cat(sprintf("[%s] %s exists, not refitted\n", format(Sys.time(), "%H:%M:%S"), path))
    fit <- readRDS(path); secs <- tryCatch(fit$time()$total, error = function(e) NA_real_)
  } else {
    cat(sprintf("[%s] start %s (T = %d)\n", format(Sys.time(), "%H:%M:%S"), nm, stan_data$T))
    mod <- cmdstan_model(spec$file)
    t0 <- Sys.time()
    fit <- do.call(mod$sample, c(list(data = stan_data, chains = 4, parallel_chains = 4,
                                      iter_warmup = 2000, iter_sampling = 2000, seed = 1,
                                      refresh = 1000, show_exceptions = FALSE), spec$extras))
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    fit$save_object(file = path)
  }
  dr <- diag_row(nm, fit, spec, secs)
  print(dr)
  old <- if (file.exists(diag_csv)) read.csv(diag_csv, stringsAsFactors = FALSE) else NULL
  if (!is.null(old)) old <- old[old$model != nm, ]
  write.csv(rbind(old, dr), diag_csv, row.names = FALSE)
  cat(sprintf("[%s] done %s\n", format(Sys.time(), "%H:%M:%S"), nm))
  rm(fit); invisible(gc(verbose = FALSE))
}
cat("ALL DONE\n")
