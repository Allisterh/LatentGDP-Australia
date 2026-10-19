# paper/ppc_m13.R
#
# Prior and posterior predictive checks of the production model m13 (1980 Q1 - 2026 Q2).
#
#   Rscript paper/ppc_m13.R
#   Rscript paper/ppc_m13.R --variant=tight   # the tight-priors research variant
#                                             # (stan/m13_tight_priors.stan): prior draws from
#                                             # its priors, posterior from
#                                             # output/sample_1980/fit_m13_tight.rds; writes
#                                             # output/ppc/ppc_tight_draws.rds and
#                                             # paper/data/ppc_tight_{prior,posterior,coverage}.csv
#
# Reads data/processed/gdp_growth.rds and output/fit_production.rds. Writes
#   output/ppc/ppc_draws.rds           replicated test statistics, per draw (gitignored)
#   paper/data/ppc_prior.csv           prior predictive summaries
#   paper/data/ppc_posterior.csv       posterior predictive checks (observed statistic,
#                                      replicated quantiles, posterior predictive p-value)
#   paper/data/ppc_coverage.csv        coverage of the observed measures by 50% / 90%
#                                      posterior predictive intervals
# No Stan model is run: replicated series are simulated in R from the model equations
# (stan/m13_serial_latent_sv.stan), vectorised over draws.
#
# Two kinds of replication:
#   "prior"       hyperparameters from the priors (priors.md, m13 section), then g, h, mu,
#                 eps and y simulated forward from their stationary starts;
#   "full"        hyperparameters from the posterior, then g, h, mu, eps and y simulated
#                 forward (a fresh latent history: checks the whole generative model);
#   "conditional" posterior draws of the latent path mu and the measurement log-variances
#                 h kept, only the measurement errors eps simulated (checks the error model
#                 given the latent states).
# The fit's own y_rep is not used: it conditions on the observed y_{t-1} (one-step
# replication), which makes autocorrelation checks uninformative.
# Units: annualised per cent (400 x dlog), as everywhere in the project.

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
})
stopifnot(file.exists("R/production_config.R"))
source("R/production_config.R")
set.seed(20261019)
dir.create("output/ppc", showWarnings = FALSE, recursive = TRUE)
args <- commandArgs(trailingOnly = TRUE)
variant <- if (length(hit <- grep("^--variant=", args, value = TRUE))) sub("^--variant=", "", hit) else "m13"
stopifnot(variant %in% c("m13", "tight"))
tight <- variant == "tight"
fit_path <- if (tight) "output/sample_1980/fit_m13_tight.rds" else "output/fit_production.rds"
sfx <- if (tight) "ppc_tight_" else "ppc_"

g <- load_production_growth()
Y <- as.matrix(g[, c("g_E", "g_I", "g_P")]); T <- nrow(Y)
dates <- as.Date(g$date)

# ---- test statistics ---------------------------------------------------------------------
ac1 <- function(x) cor(x[-1], x[-length(x)])
pre <- dates < as.Date("2020-01-01")
stats_of <- function(y) {           # y: T x 3 matrix
  hl <- rowMeans(y)
  dEI <- y[, 1] - y[, 2]; dEP <- y[, 1] - y[, 3]; dIP <- y[, 2] - y[, 3]
  c(sd_E = sd(y[, 1]), sd_I = sd(y[, 2]), sd_P = sd(y[, 3]), sd_head = sd(hl),
    sd_head_pre2020 = sd(hl[pre]),
    ac1_E = ac1(y[, 1]), ac1_I = ac1(y[, 2]), ac1_P = ac1(y[, 3]), ac1_head = ac1(hl),
    sd_EI = sd(dEI), sd_EP = sd(dEP), sd_IP = sd(dIP),
    ac1_EI = ac1(dEI), ac1_EP = ac1(dEP), ac1_IP = ac1(dIP),
    cor_EI = cor(y[, 1], y[, 2]),
    min_head = min(hl), max_head = max(hl),
    maxabs_EI = max(abs(dEI)),
    n_head_below_m4 = sum(hl < -4),                      # quarters with headline below -4% ann.
    sd_EI_2010s = sd(dEI[dates >= as.Date("2010-03-01") & dates <= as.Date("2019-12-01")]),
    sd_EI_1980s = sd(dEI[dates < as.Date("1990-01-01")]))
}
obs <- stats_of(Y)

# ---- simulator (vectorised over S draws) -------------------------------------------------
# par: list of vectors of length S (c, rho, mu_g, phi_g, omega_g, rho_EI) and S x 3 matrices
# (mu_h, phi, omega, psi). If mu / h are given (S x T, S x T x 3) they are used instead of
# being simulated ("conditional").
simulate <- function(par, S, mu = NULL, h = NULL) {
  out <- array(NA_real_, c(S, T, 3))
  r <- par$rho_EI
  if (is.null(h)) {
    h <- array(NA_real_, c(S, T, 3))
    for (k in 1:3) {
      h[, 1, k] <- par$mu_h[, k] + par$omega[, k] / sqrt(1 - par$phi[, k]^2) * rnorm(S)
      for (t in 2:T) h[, t, k] <- par$mu_h[, k] + par$phi[, k] * (h[, t - 1, k] - par$mu_h[, k]) +
          par$omega[, k] * rnorm(S)
    }
  }
  if (is.null(mu)) {
    gg <- matrix(NA_real_, S, T); mu <- matrix(NA_real_, S, T)
    gg[, 1] <- par$mu_g + par$omega_g / sqrt(1 - par$phi_g^2) * rnorm(S)
    for (t in 2:T) gg[, t] <- par$mu_g + par$phi_g * (gg[, t - 1] - par$mu_g) + par$omega_g * rnorm(S)
    mu[, 1] <- par$c / (1 - par$rho) + exp(gg[, 1] / 2) / sqrt(1 - par$rho^2) * rnorm(S)
    for (t in 2:T) mu[, t] <- par$c + par$rho * mu[, t - 1] + exp(gg[, t] / 2) * rnorm(S)
  }
  psi <- par$psi
  # period 1: eps_1 ~ N(0, V_1), V_1[j,k] = R_1[j,k] / (1 - psi_j psi_k)
  s1 <- exp(h[, 1, ] / 2)
  vE <- s1[, 1]^2 / (1 - psi[, 1]^2); vI <- s1[, 2]^2 / (1 - psi[, 2]^2); vP <- s1[, 3]^2 / (1 - psi[, 3]^2)
  cEI <- r * s1[, 1] * s1[, 2] / (1 - psi[, 1] * psi[, 2])
  z1 <- rnorm(S); z2 <- rnorm(S)
  eE <- sqrt(vE) * z1
  eI <- cEI / sqrt(vE) * z1 + sqrt(pmax(vI - cEI^2 / vE, 0)) * z2
  eP <- sqrt(vP) * rnorm(S)
  out[, 1, ] <- mu[, 1] + cbind(eE, eI, eP)
  for (t in 2:T) {
    st <- exp(h[, t, ] / 2)
    z1 <- rnorm(S); z2 <- rnorm(S); z3 <- rnorm(S)
    uE <- st[, 1] * z1; uI <- st[, 2] * (r * z1 + sqrt(1 - r^2) * z2); uP <- st[, 3] * z3
    eE <- psi[, 1] * eE + uE; eI <- psi[, 2] * eI + uI; eP <- psi[, 3] * eP + uP
    out[, t, ] <- mu[, t] + cbind(eE, eI, eP)
  }
  list(y = out, mu = mu)
}
stat_draws <- function(yarr) t(apply(yarr, 1, function(y) stats_of(y)))

# ---- prior predictive --------------------------------------------------------------------
S0 <- 4000
rtn <- function(n, m, s) { x <- rnorm(n, m, s); while (any(b <- abs(x) >= 1)) x[b] <- rnorm(sum(b), m, s); x }
pri <- list(c = rnorm(S0, 1.5, 1), rho = rtn(S0, 0.5, 0.3), mu_g = rnorm(S0, 2, 1),
            phi_g = runif(S0, -1, 1), omega_g = abs(rnorm(S0, 0, 0.5)),
            mu_h = matrix(rnorm(3 * S0, 2, 1), S0), phi = matrix(runif(3 * S0, -1, 1), S0),
            omega = matrix(abs(rnorm(3 * S0, 0, 0.5)), S0), rho_EI = rtn(S0, 0, 0.5),
            psi = matrix(rtn(3 * S0, 0, 0.5), S0))
if (tight) {
  # stan/m13_tight_priors.stan: rho ~ Beta(2.625, 2.625), mu_bar ~ N(3, 2), c = mu_bar (1 - rho);
  # stationary sds s_g, s_h ~ half-N(0, 1), phi ~ U(-1, 1), omega = s sqrt(1 - phi^2).
  rho_t <- rbeta(S0, 2.625, 2.625)
  pri$rho <- rho_t; pri$c <- rnorm(S0, 3, 2) * (1 - rho_t)
  pri$omega_g <- abs(rnorm(S0, 0, 1)) * sqrt(1 - pri$phi_g^2)
  pri$omega <- matrix(abs(rnorm(3 * S0, 0, 1)), S0) * sqrt(1 - pri$phi^2)
}
sp <- simulate(pri, S0)
st_prior <- stat_draws(sp$y)
lr_sd_mu <- apply(sp$mu, 1, sd)
tau_prior <- exp(pri$mu_g / 2)
eps_bar_prior <- exp(pri$mu_h / 2) / sqrt(1 - pri$psi^2)
mean_head_prior <- apply(sp$y, 1, function(y) mean(rowMeans(y)))
qq <- function(x) quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95, 0.99), names = FALSE, na.rm = TRUE)
prior_rows <- rbind(
  data.frame(quantity = "long-run mean of latent growth c/(1-rho)", t(qq(pri$c / (1 - pri$rho)))),
  data.frame(quantity = "mean headline growth over the sample", t(qq(mean_head_prior))),
  data.frame(quantity = "tau = exp(mu_g/2)", t(qq(tau_prior))),
  data.frame(quantity = "sd of the latent path", t(qq(lr_sd_mu))),
  data.frame(quantity = "long-run error sd, E", t(qq(eps_bar_prior[, 1]))),
  data.frame(quantity = "sd of headline growth", t(qq(st_prior[, "sd_head"]))),
  data.frame(quantity = "sd of E - I difference", t(qq(st_prior[, "sd_EI"]))),
  data.frame(quantity = "max |headline| (min and max)", t(qq(pmax(abs(st_prior[, "min_head"]), abs(st_prior[, "max_head"]))))),
  data.frame(quantity = "first-order autocorrelation of headline", t(qq(st_prior[, "ac1_head"]))))
names(prior_rows)[-1] <- c("q05", "q25", "q50", "q75", "q95", "q99")
prior_rows$observed <- c(NA, mean(rowMeans(Y)), NA, NA, NA, obs["sd_head"], obs["sd_EI"],
                         max(abs(obs[c("min_head", "max_head")])), obs["ac1_head"])
prior_rows$p_obs_below <- c(NA,
  mean(mean_head_prior <= mean(rowMeans(Y))), NA, NA, NA,
  mean(st_prior[, "sd_head"] <= obs["sd_head"]), mean(st_prior[, "sd_EI"] <= obs["sd_EI"]),
  mean(pmax(abs(st_prior[, "min_head"]), abs(st_prior[, "max_head"])) <= max(abs(obs[c("min_head", "max_head")]))),
  mean(st_prior[, "ac1_head"] <= obs["ac1_head"]))
prior_rows$share_sd_head_10_or_more <- c(rep(NA, 5), mean(st_prior[, "sd_head"] > 10), rep(NA, 3))
write.csv(prior_rows, paste0("paper/data/", sfx, "prior.csv"), row.names = FALSE)

# ---- posterior predictive ----------------------------------------------------------------
fit <- readRDS(fit_path)
S <- 1000
dm <- as_draws_matrix(fit$draws(c("c", "rho", "mu_g", "phi_g", "omega_g", "mu_h", "phi", "omega",
                                  "rho_EI", "psi")))
idx <- round(seq(1, nrow(dm), length.out = S))
dm <- unclass(dm)[idx, ]
post <- list(c = dm[, "c"], rho = dm[, "rho"], mu_g = dm[, "mu_g"], phi_g = dm[, "phi_g"],
             omega_g = dm[, "omega_g"], rho_EI = dm[, "rho_EI"],
             mu_h = dm[, paste0("mu_h[", 1:3, "]")], phi = dm[, paste0("phi[", 1:3, "]")],
             omega = dm[, paste0("omega[", 1:3, "]")], psi = dm[, paste0("psi[", 1:3, "]")])
mu_post <- unclass(fit$draws("mu", format = "draws_matrix"))[idx, ]
hmat <- unclass(fit$draws("h", format = "draws_matrix"))[idx, ]
h_post <- array(NA_real_, c(S, T, 3))
for (k in 1:3) h_post[, , k] <- hmat[, paste0("h[", 1:T, ",", k, "]")]
rm(fit, hmat); invisible(gc(verbose = FALSE))

sf <- simulate(post, S)                                  # full replication
sc <- simulate(post, S, mu = mu_post, h = h_post)        # conditional on mu, h
st_full <- stat_draws(sf$y); st_cond <- stat_draws(sc$y)
saveRDS(list(obs = obs, prior = st_prior, full = st_full, conditional = st_cond, idx = idx),
        paste0("output/ppc/", sfx, "draws.rds"))

ppc_rows <- do.call(rbind, lapply(c("full", "conditional"), function(kind) {
  st <- if (kind == "full") st_full else st_cond
  do.call(rbind, lapply(names(obs), function(nm) {
    x <- st[, nm]
    data.frame(kind = kind, statistic = nm, observed = obs[[nm]],
               q05 = quantile(x, 0.05, names = FALSE, na.rm = TRUE),
               q50 = quantile(x, 0.5, names = FALSE, na.rm = TRUE),
               q95 = quantile(x, 0.95, names = FALSE, na.rm = TRUE),
               p_rep_ge_obs = mean(x >= obs[[nm]], na.rm = TRUE))
  }))
}))
write.csv(ppc_rows, paste0("paper/data/", sfx, "posterior.csv"), row.names = FALSE)

# coverage of each observed measure by the conditional replication's pointwise intervals
cov_rows <- do.call(rbind, lapply(1:3, function(k) {
  lo90 <- apply(sc$y[, , k], 2, quantile, 0.05); hi90 <- apply(sc$y[, , k], 2, quantile, 0.95)
  lo50 <- apply(sc$y[, , k], 2, quantile, 0.25); hi50 <- apply(sc$y[, , k], 2, quantile, 0.75)
  inside90 <- Y[, k] >= lo90 & Y[, k] <= hi90; inside50 <- Y[, k] >= lo50 & Y[, k] <= hi50
  data.frame(measure = c("E", "I", "P")[k], cover50 = mean(inside50), cover90 = mean(inside90),
             cover90_pre2020 = mean(inside90[pre]), cover90_post2020 = mean(inside90[!pre]))
}))
write.csv(cov_rows, paste0("paper/data/", sfx, "coverage.csv"), row.names = FALSE)
cat("done\n"); print(prior_rows); print(ppc_rows); print(cov_rows)
