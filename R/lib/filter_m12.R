# Marginal-likelihood filters for m12 (stan/m12_struct_sv_serial.stan: m9 + AR(1)
# measurement errors), with the latent states integrated out.
#
#   source("R/lib/filters.R")        # dcol(), as_dm(), model_defs, filter_loglik()
#   source("R/lib/filter_m12.R")
#   ll <- filter_loglik_m12(draws, y, N = 2000L, seed = 1)   # list(full = S x T, sub, ess_min)
#
# filter_loglik_m12() has the signature and return value of filter_loglik() in
# R/lib/filters.R, minus the `model` argument: `full` is the S x T matrix of
# log p-hat(y_t | y_{1:t-1}, theta_s) from a Rao-Blackwellised particle filter (particles
# over h, exact Kalman over (mu_t, mu_{t-1}) given h), `sub` is NULL (all three measures
# are the target), `ess_min` the minimum ESS of the particle weights per draw, `ess_frac`
# the S x T matrix ESS_t / N_t. N may be a length-T schedule (see R/lib/lfo.R).
#
# Registering m12 with the LFO code (R/lib/lfo.R, R/03_compare.R): call register_m12()
# after sourcing both files. It adds model_defs$m12 (type "sv", so filter_par() runs it in
# parallel workers) and wraps filter_loglik() so that filter_loglik("m12", ...) dispatches
# here; every other model goes to the original function unchanged. Equivalently, add to
# R/lib/filters.R
#   model_defs$m12 <- list(type = "sv", pars = M12_PARS)
# and, as the first line of filter_loglik(),
#   if (model == "m12") return(filter_loglik_m12(draws, y, ctx, N, ess_thr, seed, h1_sd_max))
# The 03_compare.R spec entry is
#   m12 = list(group = 1, file = "stan/m12_struct_sv_serial.stan", data = data_g1, extras = sv_set)
#
# Initial condition: exactly the Stan program's (eps_1 ~ N(0, V_1),
# V_1[j,k] = R_1[j,k] / (1 - psi_j psi_k), mu_1 stationary, h_1 stationary).

suppressPackageStartupMessages(library(Rcpp))

M12_PARS <- c("c", "rho", "tau", "mu_h", "phi", "omega", "rho_EI", "psi")

.m12_env <- new.env()
load_filter_m12 <- function() {
  if (!exists("rbpf_m12_cpp", envir = .m12_env)) {
    src <- "R/lib/filter_m12.cpp"
    if (!file.exists(src)) stop("run from the repository root (R/lib/filter_m12.cpp not found)")
    Rcpp::sourceCpp(src, env = .m12_env, verbose = FALSE)
  }
  invisible(.m12_env)
}

# 1 + r and 1 - r. From z_EI = atanh(r) when available (exact, as in the Stan program),
# otherwise from r itself (fine away from |r| = 1; posterior rho_EI is about 0.1-0.2).
m12_epem <- function(p) {
  if (!is.null(p$z_EI)) {
    z <- p$z_EI
    c(ep = 2 * plogis(2 * z), em = 2 * plogis(-2 * z))
  } else {
    r <- p$rho_EI
    c(ep = 1 + r, em = 1 - r)
  }
}

# Per-draw parameter lists from a draws matrix (Stan column names); z_EI is used if present.
m12_draw_list <- function(draws) {
  dm <- as_dm(draws)
  cols <- lapply(setNames(M12_PARS, M12_PARS), function(v) dcol(dm, v))
  if ("z_EI" %in% colnames(dm)) cols$z_EI <- dm[, "z_EI", drop = FALSE]
  lapply(seq_len(nrow(dm)), function(s) lapply(cols, function(m) unname(m[s, ])))
}

filter_loglik_m12 <- function(draws, y, ctx = list(), N = 1000L, ess_thr = 0.5,
                              seed = NULL, h1_sd_max = 3) {
  env <- load_filter_m12()
  y <- unname(as.matrix(y)); storage.mode(y) <- "double"
  stopifnot(ncol(y) == 3, !anyNA(y))
  if (!is.null(seed)) set.seed(seed)
  pl <- m12_draw_list(draws)
  S <- length(pl); T <- nrow(y)
  full <- matrix(NA_real_, S, T)
  ess_min <- numeric(S); ess_frac <- matrix(NA_real_, S, T)
  for (s in seq_len(S)) {
    p <- pl[[s]]
    e <- m12_epem(p)
    o <- env$rbpf_m12_cpp(y, p$c, p$rho, p$tau, e[["ep"]], e[["em"]], p$psi,
                          p$mu_h, p$phi, p$omega, as.integer(N), ess_thr, h1_sd_max)
    full[s, ] <- o[, 1]
    ess_min[s] <- min(o[, 3])
    ess_frac[s, ] <- o[, 3] / o[, 4]
  }
  list(full = full, sub = NULL, ess_min = ess_min, ess_frac = ess_frac)
}

# Conditional-on-h Kalman filter for one draw (p: named list as from m12_draw_list) and one
# h path (T x 3). Returns the T x 6 matrix of Stan's kf_m12(): log p(y_t | y_{1:t-1}, h,
# theta), m_t, P_t and the backward-sampling coefficients.
kf_m12_given_h <- function(p, y, h) {
  env <- load_filter_m12()
  e <- m12_epem(p)
  y <- unname(as.matrix(y)); storage.mode(y) <- "double"
  h <- unname(as.matrix(h)); storage.mode(h) <- "double"
  env$kf_m12_cpp(y, h, p$c, p$rho, p$tau, e[["ep"]], e[["em"]], p$psi)
}

# Install m12 into the LFO machinery (see the header). Idempotent.
register_m12 <- function(envir = globalenv()) {
  if (!exists("filter_loglik", envir = envir) || !exists("model_defs", envir = envir))
    stop("source R/lib/filters.R before register_m12()")
  md <- get("model_defs", envir = envir)
  md$m12 <- list(type = "sv", pars = M12_PARS)
  assign("model_defs", md, envir = envir)
  base <- get("filter_loglik", envir = envir)
  if (isTRUE(attr(base, "m12_wrapped"))) return(invisible(TRUE))
  wrapped <- function(model, draws, y, ctx = list(), N = 1000L, ess_thr = 0.5,
                      seed = NULL, h1_sd_max = 3) {
    if (identical(model, "m12"))
      return(filter_loglik_m12(draws, y, ctx, N = N, ess_thr = ess_thr, seed = seed,
                               h1_sd_max = h1_sd_max))
    base(model, draws, y, ctx, N = N, ess_thr = ess_thr, seed = seed, h1_sd_max = h1_sd_max)
  }
  attr(wrapped, "m12_wrapped") <- TRUE
  assign("filter_loglik", wrapped, envir = envir)
  invisible(TRUE)
}
