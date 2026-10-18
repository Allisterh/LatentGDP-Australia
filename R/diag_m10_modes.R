# Mode check for the re-specified m10 (stan/m10_labour_sv.stan, latent SV). Research only.
#
#   Rscript R/diag_m10_modes.R dispersed [--out=PATH]   # 16 chains, dispersed inits
#   Rscript R/diag_m10_modes.R modeA     [--out=PATH]   # 4 chains started in old mode A
#   Rscript R/diag_m10_modes.R summarise PATH [LABEL]   # per-chain table of a saved fit
#   options: --parallel=N (default 2)
#
# dispersed: 16 chains x (1000 warmup + 1000 draws), seed 777, adapt_delta 0.99,
#   max_treedepth 12. Inits: c ~ U(0, 3), rho ~ U(-0.5, 0.9), alpha ~ N(-1, 1),
#   beta ~ U(0.2, 1.5), mu_g and mu_h ~ U(-1, 4), phi_g and phi ~ U(0, 0.99), omega_g and
#   omega ~ U(0.05, 1.5), g_raw and h_raw ~ N(0, 0.5) (set.seed(4242)). This is the design of
#   the 2026-10-17 diagnosis (NOTES.md, "m10: two modes, diagnosis and re-specification").
#   The non-centred g_raw / h_raw start near 0, so no chain starts in the old mode A.
# modeA: 4 chains x (700 + 500), seed 2024, started at an old-mode-A configuration: all four
#   measurement sds about 10 in 2020 Q2 (about 8 in 2020 Q3), flat latent volatility
#   (omega_g 0.10-0.16), tau = exp(mu_g / 2) = 2.1, rho 0.39, phi_E 0.66, omega_E 0.67.
#   If mode A still existed under the new specification, these chains would stay there.
# summarise: one row per chain (medians of lp__, rho, tau, omega_g, phi[1], omega[1], mu,
#   sigma_t[, 1] and tau_t in 2020 Q2; divergences; E-BFMI) appended to
#   output/tables/m10_mode_check.csv under LABEL, plus R-hat across all chains of the run.
# The fit is saved to --out (default output/fit_m10_<mode>.rds, gitignored) and summarised.

suppressPackageStartupMessages({ library(cmdstanr); library(posterior) })
args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[length(hit)]) else default
}
pos <- args[!startsWith(args, "--")]
mode <- if (length(pos)) pos[1] else "dispersed"
stopifnot(mode %in% c("dispersed", "modeA", "summarise"))

g4 <- readRDS("data/processed/gdp_growth_with_hours.rds")
dat <- list(T = nrow(g4), y_gdp = as.matrix(g4[, c("g_E", "g_I", "g_P")]), y_U = g4$g_U)
T <- dat$T
iq <- which(g4$date == as.Date("2020-06-01"))
stan_file <- "stan/m10_labour_sv.stan"
tab_path <- "output/tables/m10_mode_check.csv"

summarise_fit <- function(fit, label) {
  v <- c("lp__", "rho", "tau", "omega_g", "phi[1]", "omega[1]", sprintf("mu[%d]", iq),
         sprintf("sigma_t[%d,1]", iq), sprintf("tau_t[%d]", iq))
  d <- as_draws_df(fit$draws(v))
  sdg <- fit$sampler_diagnostics(format = "draws_df")
  med <- aggregate(as.data.frame(d)[, v], by = list(chain = d$.chain), FUN = median)
  names(med) <- c("chain", "lp", "rho", "tau", "omega_g", "phi_E", "omega_E", "mu_2020Q2",
                  "sd_E_2020Q2", "tau_t_2020Q2")
  med$divergences <- as.vector(tapply(sdg$divergent__, sdg$.chain, sum))
  med$ebfmi <- fit$diagnostic_summary("ebfmi", quiet = TRUE)$ebfmi
  rh <- summarise_draws(fit$draws(c("rho", "mu_g", "phi_g", "omega_g", "beta", "mu_h", "phi", "omega")),
                        "rhat")
  med$run <- label
  med$stan_file <- stan_file
  med$n_chains <- fit$num_chains()
  med$max_rhat_all_chains <- max(rh$rhat)
  med$worst_rhat <- rh$variable[which.max(rh$rhat)]
  old <- if (file.exists(tab_path)) read.csv(tab_path, stringsAsFactors = FALSE) else NULL
  if (!is.null(old)) old <- old[old$run != label, , drop = FALSE]
  out <- rbind(old, med[, c("run", "stan_file", "n_chains", "chain", setdiff(names(med),
                            c("run", "stan_file", "n_chains", "chain")))])
  write.csv(out, tab_path, row.names = FALSE)
  print(format(med, digits = 3), row.names = FALSE)
  invisible(med)
}

if (mode == "summarise") {
  stopifnot(length(pos) >= 2)
  summarise_fit(readRDS(pos[2]), if (length(pos) >= 3) pos[3] else basename(pos[2]))
  quit(save = "no")
}

mod <- cmdstan_model(stan_file)
par <- as.integer(opt("parallel", "2"))
out <- opt("out", sprintf("output/fit_m10_%s.rds", mode))
if (mode == "dispersed") {
  set.seed(4242)
  inits <- lapply(1:16, function(i) list(
    c = runif(1, 0, 3), rho = runif(1, -0.5, 0.9), alpha = rnorm(1, -1, 1), beta = runif(1, 0.2, 1.5),
    mu_g = runif(1, -1, 4), phi_g = runif(1, 0, 0.99), omega_g = runif(1, 0.05, 1.5),
    g_raw = rnorm(T, 0, 0.5), mu_h = runif(4, -1, 4), phi = runif(4, 0, 0.99),
    omega = runif(4, 0.05, 1.5), h_raw = matrix(rnorm(T * 4, 0, 0.5), T, 4)))
  fit <- mod$sample(data = dat, chains = 16, parallel_chains = par, iter_warmup = 1000,
                    iter_sampling = 1000, seed = 777, init = inits, refresh = 0,
                    adapt_delta = 0.99, max_treedepth = 12, show_exceptions = FALSE)
} else {
  set.seed(99)
  inits <- lapply(1:4, function(j) {
    hr <- matrix(rnorm(T * 4, 0, 0.3), T, 4)
    phi <- c(0.66, 0.9, 0.9, 0.6); om <- c(0.67, 0.6, 0.6, 0.6); muh <- c(2, 1.8, 1.4, 1.8)
    hr[iq, ] <- (log(100) - muh) / om                                 # sd ~ 10 in 2020 Q2
    hr[iq + 1, ] <- (log(64) - muh - phi * (log(100) - muh)) / om     # sd ~ 8 in 2020 Q3
    list(c = 1.3 + 0.1 * j, rho = 0.39, alpha = -1, beta = 0.9, mu_g = 2 * log(2.1),
         phi_g = 0.5, omega_g = 0.08 + 0.02 * j, g_raw = rnorm(T, 0, 0.1),
         mu_h = muh, phi = phi, omega = om, h_raw = hr)
  })
  fit <- mod$sample(data = dat, chains = 4, parallel_chains = par, iter_warmup = 700,
                    iter_sampling = 500, seed = 2024, init = inits, refresh = 0,
                    adapt_delta = 0.99, max_treedepth = 12, show_exceptions = FALSE)
}
fit$save_object(out)
summarise_fit(fit, mode)
