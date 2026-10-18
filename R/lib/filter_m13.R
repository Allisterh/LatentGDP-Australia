# Marginal-likelihood filters for m13 (stan/m13_serial_latent_sv.stan: m12 + stochastic
# volatility on the latent innovation), with the latent states integrated out.
# Rcpp code: R/lib/filter_m13.cpp.
#
#   source("R/lib/filters.R")        # dcol(), as_dm(), model_defs, filter_loglik()
#   source("R/lib/filter_m13.R")
#   ll <- filter_loglik_m13(draws, y, N = 2000L, seed = 1)   # list(full = S x T, sub, ess_min, ess_frac)
#
# filter_loglik_m13() has the signature and return value of filter_loglik() in
# R/lib/filters.R, minus the `model` argument: `full` is the S x T matrix of
# log p-hat(y_t | y_{1:t-1}, theta_s) from a Rao-Blackwellised particle filter (particles
# over (h, g), exact Kalman over (mu_t, mu_{t-1}) given them), `sub` is NULL (all three
# measures are the target), `ess_min` the minimum ESS of the particle weights per draw,
# `ess_frac` the S x T matrix ESS_t / N_t. N may be a length-T schedule (R/lib/lfo.R,
# pf_schedule), exactly as for m12 and m10.
#
# Proposal for the latent log-variance g_t (see filter_m13.cpp): M13_PROP = 1 is the locally
# adapted proposal of m10's filter (t_4 around the conditional mode, defensive share
# M13_EPS_DEF from the transition); 0 is the bootstrap proposal. The weights carry the exact
# correction, so the estimate is unbiased either way.
#
# Registering m13 with the LFO code (R/lib/lfo.R, R/03_compare.R): call register_m13()
# after sourcing R/lib/filters.R and this file. It adds model_defs$m13 (type "sv", so
# filter_par() runs it in forked workers and psis_lfo() builds a particle schedule for it;
# pars = M13_PARS plus z_EI, as 03_compare.R sets for m12) and wraps filter_loglik() so
# that filter_loglik("m13", ...) dispatches here; every other model goes to the previous
# function unchanged (it composes with register_m12() and register_m10_filter() in any
# order). Call load_filter_m13() in the parent before forking (03_compare.R does this for
# the other filters). The 03_compare.R spec entry is
#   m13 = list(group = 1, file = "stan/m13_serial_latent_sv.stan", data = data_g1, extras = sv_set)
# with "m13" in group_models$`1`. On 12 posterior draws of the full-sample fit, a
# pilot-based schedule raised N in 23 of 267 quarters (to 100000 only at 2020 Q2; cost
# x1.31 of N = 2000), and 4 scheduled passes over the 12 draws took a few seconds on 4
# cores, so the default SV_DRAWS should be affordable.
#
# Initial condition: exactly the Stan program's (eps_1 ~ N(0, V_1), V_1[j,k] = R_1[j,k] /
# (1 - psi_j psi_k), g_1 and h_1 stationary, mu_1 ~ N(c/(1-rho), exp(g_1)/(1-rho^2))).

suppressPackageStartupMessages(library(Rcpp))

M13_PARS    <- c("c", "rho", "mu_g", "phi_g", "omega_g", "mu_h", "phi", "omega", "rho_EI", "psi")
M13_PROP    <- 1L
M13_EPS_DEF <- 0.1

.m13_env <- new.env()
load_filter_m13 <- function() {
  if (!exists("rbpf_m13_cpp", envir = .m13_env)) {
    src <- "R/lib/filter_m13.cpp"
    if (!file.exists(src)) stop("run from the repository root (R/lib/filter_m13.cpp not found)")
    Rcpp::sourceCpp(src, env = .m13_env, verbose = FALSE)
  }
  invisible(.m13_env)
}

# 1 + r and 1 - r: from z_EI = atanh(r) when available (exact, as in the Stan program),
# otherwise from r itself (as m12_epem in R/lib/filter_m12.R).
m13_epem <- function(p) {
  if (!is.null(p$z_EI)) {
    z <- p$z_EI
    c(ep = 2 * plogis(2 * z), em = 2 * plogis(-2 * z))
  } else {
    r <- p$rho_EI
    c(ep = 1 + r, em = 1 - r)
  }
}

# Per-draw parameter lists from a draws matrix (Stan column names); z_EI is used if present.
m13_draw_list <- function(draws) {
  dm <- as_dm(draws)
  cols <- lapply(setNames(M13_PARS, M13_PARS), function(v) dcol(dm, v))
  if ("z_EI" %in% colnames(dm)) cols$z_EI <- dm[, "z_EI", drop = FALSE]
  lapply(seq_len(nrow(dm)), function(s) lapply(cols, function(m) unname(m[s, ])))
}

filter_loglik_m13 <- function(draws, y, ctx = list(), N = 1000L, ess_thr = 0.5,
                              seed = NULL, h1_sd_max = 3, prop = M13_PROP,
                              eps_def = M13_EPS_DEF) {
  env <- load_filter_m13()
  y <- unname(as.matrix(y)); storage.mode(y) <- "double"
  stopifnot(ncol(y) == 3, !anyNA(y))
  if (!is.null(seed)) set.seed(seed)
  pl <- m13_draw_list(draws)
  S <- length(pl); T <- nrow(y)
  full <- matrix(NA_real_, S, T)
  ess_min <- numeric(S); ess_frac <- matrix(NA_real_, S, T)
  for (s in seq_len(S)) {
    p <- pl[[s]]
    e <- m13_epem(p)
    o <- env$rbpf_m13_cpp(y, p$c, p$rho, e[["ep"]], e[["em"]], p$psi,
                          p$mu_g, p$phi_g, p$omega_g, p$mu_h, p$phi, p$omega,
                          as.integer(N), ess_thr, h1_sd_max, as.integer(prop), eps_def)
    full[s, ] <- o[, 1]
    ess_min[s] <- min(o[, 3])
    ess_frac[s, ] <- o[, 3] / o[, 4]
  }
  list(full = full, sub = NULL, ess_min = ess_min, ess_frac = ess_frac)
}

# Conditional-on-(g, h) Kalman filter for one draw (p: named list as from m13_draw_list),
# one h path (T x 3) and one g path (length T). Returns the T x 6 matrix of Stan's kf_m13():
# log p(y_t | y_{1:t-1}, g, h, theta), m_t, P_t and the backward-sampling coefficients.
kf_m13_given <- function(p, y, h, g, tau_t = exp(0.5 * g)) {
  env <- load_filter_m13()
  e <- m13_epem(p)
  y <- unname(as.matrix(y)); storage.mode(y) <- "double"
  h <- unname(as.matrix(h)); storage.mode(h) <- "double"
  env$kf_m13_cpp(y, h, as.double(tau_t), p$c, p$rho, e[["ep"]], e[["em"]], p$psi)
}

# Install m13 into the LFO machinery (see the header). Idempotent.
register_m13 <- function(envir = globalenv()) {
  if (!exists("filter_loglik", envir = envir) || !exists("model_defs", envir = envir))
    stop("source R/lib/filters.R before register_m13()")
  md <- get("model_defs", envir = envir)
  md$m13 <- list(type = "sv", pars = c(M13_PARS, "z_EI"))
  assign("model_defs", md, envir = envir)
  base <- get("filter_loglik", envir = envir)
  if (isTRUE(attr(base, "m13_wrapped"))) return(invisible(TRUE))
  wrapped <- function(model, draws, y, ctx = list(), N = 1000L, ess_thr = 0.5,
                      seed = NULL, h1_sd_max = 3) {
    if (identical(model, "m13"))
      return(filter_loglik_m13(draws, y, ctx, N = N, ess_thr = ess_thr, seed = seed,
                               h1_sd_max = h1_sd_max))
    base(model, draws, y, ctx, N = N, ess_thr = ess_thr, seed = seed, h1_sd_max = h1_sd_max)
  }
  # keep the other models' wrap markers so their register functions stay idempotent
  for (a in c("m12_wrapped", "m10_wrapped"))
    if (isTRUE(attr(base, a))) attr(wrapped, a) <- TRUE
  attr(wrapped, "m13_wrapped") <- TRUE
  assign("filter_loglik", wrapped, envir = envir)
  invisible(TRUE)
}
