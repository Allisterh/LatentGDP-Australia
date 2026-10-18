# LFO filter for the re-specified m10 (stan/m10_labour_sv.stan: 4 signals, SV on every
# measurement error and on the latent innovation). Rcpp code: R/lib/filter_m10.cpp.
#
#   source("R/lib/filters.R"); source("R/lib/filter_m10.R")
#   register_m10_filter()      # route filter_loglik("m10", ...) to this filter
#
# register_m10_filter() replaces model_defs$m10 (new parameter list) and wraps
# filter_loglik so that model "m10" uses m10_filter_loglik(); every other model is
# untouched. R/lib/lfo.R and R/03_compare.R look both up by name at call time, so
# calling it once after sourcing filters.R is enough.
#
# Proposal for the latent log-variance g_t (see filter_m10.cpp): M10_PROP = 1 is the
# locally adapted proposal (Laplace mode of prior x conditional likelihood, t_4 with a
# defensive share M10_EPS_DEF from the transition); 0 is the bootstrap proposal. The
# weights carry the exact correction, so the estimate is unbiased either way.

M10_PROP    <- 1L
M10_EPS_DEF <- 0.1

.m10_env <- new.env()
load_filter_m10 <- function() {
  if (!exists("m10_rbpf_cpp", envir = .m10_env)) {
    if (!file.exists("R/lib/filter_m10.cpp"))
      stop("run from the repository root (R/lib/filter_m10.cpp not found)")
    Rcpp::sourceCpp("R/lib/filter_m10.cpp", env = .m10_env, verbose = FALSE)
  }
  invisible(.m10_env)
}

M10_PARS <- c("c", "rho", "alpha", "beta", "mu_g", "phi_g", "omega_g", "mu_h", "phi", "omega")

m10_draw_list <- function(draws) {
  dm <- as_dm(draws)
  cols <- lapply(setNames(M10_PARS, M10_PARS), function(v) dcol(dm, v))
  lapply(seq_len(nrow(dm)), function(s) lapply(cols, function(m) unname(m[s, ])))
}

# Same return shape as filter_loglik(): list(full = S x T, sub = S x T (GDP-only
# marginal, hours in the past only), ess_min = S vector, ess_frac = S x T ESS_t / N_t).
# N may be a length-T particle-count schedule (see R/lib/lfo.R).
m10_filter_loglik <- function(draws, y, ctx = list(), N = 1000L, ess_thr = 0.5,
                              seed = NULL, h1_sd_max = 3,
                              prop = M10_PROP, eps_def = M10_EPS_DEF) {
  env <- load_filter_m10()
  if (!is.null(seed)) set.seed(seed)
  pl <- m10_draw_list(draws)
  S <- length(pl); T <- nrow(y)
  full <- sub <- ess_frac <- matrix(NA_real_, S, T); ess_min <- numeric(S)
  for (s in seq_len(S)) {
    p <- pl[[s]]
    o <- env$m10_rbpf_cpp(y, c(1, 1, 1, p$beta), c(0, 0, 0, p$alpha), 3L, p$c, p$rho,
                          p$mu_g, p$phi_g, p$omega_g, p$mu_h, p$phi, p$omega,
                          as.integer(N), ess_thr, h1_sd_max, as.integer(prop), eps_def)
    full[s, ] <- o[, 1]; sub[s, ] <- o[, 2]; ess_min[s] <- min(o[, 3])
    ess_frac[s, ] <- o[, 3] / o[, 4]
  }
  list(full = full, sub = sub, ess_min = ess_min, ess_frac = ess_frac)
}

# Conditional (given g and h paths) Kalman filter: [log p(y_t | ...), GDP-only column].
m10_kf_given <- function(p, y, h, g) {
  env <- load_filter_m10()
  env$m10_kf_given_cpp(y, h, g, c(1, 1, 1, p$beta), c(0, 0, 0, p$alpha), 3L, p$c, p$rho)
}

register_m10_filter <- function() {
  if (!exists("model_defs") || !exists("filter_loglik"))
    stop("source R/lib/filters.R before R/lib/filter_m10.R")
  md <- get("model_defs", envir = globalenv())
  md$m10 <- list(type = "sv", sv = NULL, pars = M10_PARS)
  assign("model_defs", md, envir = globalenv())
  if (!isTRUE(attr(get("filter_loglik", envir = globalenv()), "m10_wrapped"))) {
    base_fl <- get("filter_loglik", envir = globalenv())
    wrapped <- function(model, draws, y, ctx = list(), N = 1000L, ess_thr = 0.5,
                        seed = NULL, h1_sd_max = 3) {
      if (identical(model, "m10"))
        return(m10_filter_loglik(draws, y, ctx, N = N, ess_thr = ess_thr, seed = seed,
                                 h1_sd_max = h1_sd_max))
      base_fl(model, draws, y, ctx, N = N, ess_thr = ess_thr, seed = seed,
              h1_sd_max = h1_sd_max)
    }
    attr(wrapped, "m10_wrapped") <- TRUE
    assign("filter_loglik", wrapped, envir = globalenv())
  }
  invisible(TRUE)
}
