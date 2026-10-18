# Marginal-likelihood filters for the LatentGDP specs (latent states integrated out).
#
#   source("R/lib/filters.R")
#   ll <- filter_loglik("m6", draws, y)          # S x T matrix of log p(y_t | y_{1:t-1}, theta_s)
#
# Linear-Gaussian specs (m1, m2, m3, m5, m6, m7, m8, m11) use the exact Kalman filter
# (kf_loglik_cpp); the SV specs (m4, m9) use a Rao-Blackwellised particle filter
# (rbpf_cpp): particles over the log-variance paths h, Kalman over mu given h.
# m12 (R/lib/filter_m12.R) and the re-specified m10 (R/lib/filter_m10.R) have their own
# RBPFs; R/03_compare.R installs them with register_m12() / register_m10_filter().
# model_defs$m10 / sv_m10 below is the OLD constant-tau m10, kept only because
# tests/test_filters.R (ii) and tests/test_m3_m10.R (h) use it as a reference filter.
# Each ssm_* function below transcribes one Stan file: same parameter names, same
# observation equation, same distribution for the initial state.
#
# draws: a posterior::draws_matrix (or anything coercible to a matrix with Stan's
# column names, e.g. "sigma[2]", "Omega[3,1]") holding the hyperparameters.

suppressPackageStartupMessages({
  library(Rcpp)
  library(posterior)
})

.filters_env <- new.env()
load_filters <- function() {
  if (!exists("kf_loglik_cpp", envir = .filters_env)) {
    here <- if (file.exists("R/lib/filters.cpp")) "R/lib/filters.cpp" else
      stop("run from the repository root (R/lib/filters.cpp not found)")
    Rcpp::sourceCpp(here, env = .filters_env, verbose = FALSE)
  }
  invisible(.filters_env)
}

# ---- helpers ---------------------------------------------------------------
# Column(s) of a draws matrix for a Stan variable, in Stan's (column-major) order.
dcol <- function(dm, name) {
  cn <- colnames(dm)
  if (name %in% cn) return(dm[, name, drop = FALSE])
  hit <- startsWith(cn, paste0(name, "["))
  if (!any(hit)) stop("variable '", name, "' not in draws")
  dm[, hit, drop = FALSE]
}
as_dm <- function(draws) {
  m <- unclass(as.matrix(draws))
  storage.mode(m) <- "double"
  m
}
lse <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }

# ---- linear-Gaussian state-space forms (one draw -> SSM) --------------------
# Every function takes p = a named list of numeric vectors for one draw, y, and ctx
# (context: e.g. covid indices), and returns list(y, Z, d, H, Tm, cst, Q, a1, P1, blk).
# Z: k x m x (1 or T) cube; H: k x k x (1 or T) cube; d: (1 or T) x k matrix.

cube1 <- function(M) array(M, dim = c(nrow(M), ncol(M), 1))

ar1_scalar <- function(c, rho, tau) {
  list(Tm = matrix(rho), cst = c, Q = matrix(tau^2),
       a1 = c / (1 - rho), P1 = matrix(tau^2 / (1 - rho^2)))
}

ssm_m1 <- function(p, y, ctx) {
  # mu_t ~ N(mu_bar, tau^2) iid;  y_t^k = mu_t + N(0, sigma_k^2)
  K <- ncol(y)
  c(list(y = y, Z = cube1(matrix(1, K, 1)), d = matrix(0, 1, K),
         H = cube1(diag(p$sigma^2, K)), blk = rep(1L, K)),
    list(Tm = matrix(0), cst = p$mu_bar, Q = matrix(p$tau^2),
         a1 = p$mu_bar, P1 = matrix(p$tau^2)))
}
ssm_m2 <- function(p, y, ctx) {
  K <- ncol(y)
  c(list(y = y, Z = cube1(matrix(1, K, 1)), d = matrix(0, 1, K),
         H = cube1(diag(p$sigma^2, K)), blk = rep(1L, K)),
    ar1_scalar(p$c, p$rho, p$tau))
}
ssm_m3 <- function(p, y, ctx) {
  K <- ncol(y)
  Om <- matrix(p$Omega, K, K)                  # Stan writes Omega column-major
  H <- diag(p$sigma, K) %*% Om %*% diag(p$sigma, K)
  c(list(y = y, Z = cube1(matrix(1, K, 1)), d = matrix(0, 1, K),
         H = cube1(H), blk = rep(1L, K)),
    ar1_scalar(p$c, p$rho, p$tau))
}
ssm_m5 <- function(p, y, ctx) {
  # Serially correlated errors eps_t^k = psi_k eps_{t-1}^k + u_t^k, stationary start.
  # Stan: y_1^k ~ N(mu_1, sigma_k / sqrt(1 - psi_k^2)),
  #       y_t^k ~ N(mu_t + psi_k (y_{t-1}^k - mu_{t-1}), sigma_k).
  # Quasi-differenced observation ytil_t = y_t - psi * y_{t-1} = mu_t - psi mu_{t-1} + u_t
  # (unit Jacobian given y_{t-1}), state alpha_t = (mu_t, mu_{t-1}).
  T <- nrow(y); K <- ncol(y)
  psi <- p$psi
  yt <- y
  yt[-1, ] <- y[-1, , drop = FALSE] - sweep(y[-T, , drop = FALSE], 2, psi, `*`)
  Z <- array(0, dim = c(K, 2, T))
  Z[, 1, ] <- 1
  Z[, 2, -1] <- -psi
  H <- array(0, dim = c(K, K, T))
  for (t in seq_len(T)) H[, , t] <- diag(p$sigma^2, K)
  H[, , 1] <- diag(p$sigma^2 / (1 - psi^2), K)
  m0 <- p$c / (1 - p$rho); v0 <- p$tau^2 / (1 - p$rho^2)
  list(y = yt, Z = Z, d = matrix(0, 1, K), H = H, blk = rep(1L, K),
       Tm = matrix(c(p$rho, 1, 0, 0), 2, 2), cst = c(p$c, 0),
       Q = diag(c(p$tau^2, 0)),
       # mu_0 (second slot at t = 1) never enters: Z_1[, 2] = 0 and it is
       # overwritten by mu_1 in the first transition.
       a1 = c(m0, m0), P1 = diag(c(v0, v0)))
}
ssm_m6 <- function(p, y, ctx) {
  s <- p$sigma; r <- p$rho_EI
  H <- matrix(c(s[1]^2, r * s[1] * s[2], 0,
                r * s[1] * s[2], s[2]^2, 0,
                0, 0, s[3]^2), 3, 3)
  c(list(y = y, Z = cube1(matrix(1, 3, 1)), d = matrix(0, 1, 3),
         H = cube1(H), blk = rep(1L, 3)),
    ar1_scalar(p$c, p$rho, p$tau))
}
ssm_m7 <- function(p, y, ctx) {
  T <- nrow(y); K <- ncol(y)
  kap <- rep(1, T)
  if (length(ctx$covid_idx)) kap[ctx$covid_idx] <- p$kappa_full
  H <- array(0, dim = c(K, K, T))
  for (t in seq_len(T)) H[, , t] <- diag((p$sigma * kap[t])^2, K)
  c(list(y = y, Z = cube1(matrix(1, K, 1)), d = matrix(0, 1, K), H = H, blk = rep(1L, K)),
    ar1_scalar(p$c, p$rho, p$tau))
}
ssm_m8 <- function(p, y, ctx) {
  # y = (E, I, P, hours); hours = alpha + beta mu + N(0, sigma_U^2).
  # Block 1 = GDP measures, block 2 = hours: column 1 of the filter output is the
  # GDP-only predictive log p(y^GDP_t | y_{1:t-1}).
  c(list(y = y, Z = cube1(matrix(c(1, 1, 1, p$beta), 4, 1)),
         d = matrix(c(0, 0, 0, p$alpha), 1, 4),
         H = cube1(diag(c(p$sigma_gdp^2, p$sigma_U^2), 4)), blk = c(1L, 1L, 1L, 2L)),
    ar1_scalar(p$c, p$rho, p$tau))
}
ssm_m11 <- function(p, y, ctx) {
  # Monthly state alpha_m = (mu_m, mu_{m-1}, mu_{m-2}), mu_m AR(1) with stationary start.
  # Observation at month m: y[m, 1] = employment growth = alpha_emp + beta_emp mu_m + e,
  # y[m, 2:4] = GDP(E, I, P) of the quarter ending in month m (NA in other months)
  # = (mu_m + mu_{m-1} + mu_{m-2}) / 3 + e_k. Block 1 = employment, block 2 = GDP,
  # so column 2 is log p(y^GDP_q | all GDP before q, employment up to the quarter's last month).
  # Initial lags: the Stan program's mu_m[1] is stationary; the slots for months 0
  # and -1 never load on anything because q_end_m >= 3, so their prior is irrelevant.
  m0 <- p$c / (1 - p$rho); v0 <- p$tau^2 / (1 - p$rho^2)
  Zq <- matrix(c(p$beta_emp, 0, 0,
                 1/3, 1/3, 1/3,
                 1/3, 1/3, 1/3,
                 1/3, 1/3, 1/3), 4, 3, byrow = TRUE)
  list(y = y, Z = cube1(Zq), d = matrix(c(p$alpha_emp, 0, 0, 0), 1, 4),
       H = cube1(diag(c(p$sigma_emp^2, p$sigma_q^2), 4)), blk = c(1L, 2L, 2L, 2L),
       Tm = rbind(c(p$rho, 0, 0), c(1, 0, 0), c(0, 1, 0)), cst = c(p$c, 0, 0),
       Q = diag(c(p$tau^2, 0, 0)), a1 = rep(m0, 3), P1 = diag(v0, 3))
}

# ---- SV forms for the RBPF ---------------------------------------------------
sv_m4 <- function(p, y, ctx) {
  K <- ncol(y)
  list(y = y, Z = rep(1, K), d = rep(0, K), C = diag(K), c = p$c, rho = p$rho, tau = p$tau,
       mu_h = p$mu_h, phi = p$phi, omega = p$omega, nsub = K)
}
sv_m9 <- function(p, y, ctx) {
  r <- p$rho_EI
  C <- matrix(c(1, r, 0, r, 1, 0, 0, 0, 1), 3, 3)
  list(y = y, Z = rep(1, 3), d = rep(0, 3), C = C, c = p$c, rho = p$rho, tau = p$tau,
       mu_h = p$mu_h, phi = p$phi, omega = p$omega, nsub = 3)
}
sv_m10 <- function(p, y, ctx) {
  list(y = y, Z = c(1, 1, 1, p$beta), d = c(0, 0, 0, p$alpha), C = diag(4),
       c = p$c, rho = p$rho, tau = p$tau, mu_h = p$mu_h, phi = p$phi, omega = p$omega,
       nsub = 3)
}

model_defs <- list(
  m1  = list(type = "kf", ssm = ssm_m1, pars = c("mu_bar", "tau", "sigma")),
  m2  = list(type = "kf", ssm = ssm_m2, pars = c("c", "rho", "tau", "sigma")),
  m3  = list(type = "kf", ssm = ssm_m3, pars = c("c", "rho", "tau", "sigma", "Omega")),
  m5  = list(type = "kf", ssm = ssm_m5, pars = c("c", "rho", "tau", "sigma", "psi")),
  m6  = list(type = "kf", ssm = ssm_m6, pars = c("c", "rho", "tau", "sigma", "rho_EI")),
  m7  = list(type = "kf", ssm = ssm_m7, pars = c("c", "rho", "tau", "sigma", "kappa")),
  m8  = list(type = "kf", ssm = ssm_m8, pars = c("c", "rho", "tau", "sigma_gdp", "sigma_U", "alpha", "beta")),
  m11 = list(type = "kf", ssm = ssm_m11, pars = c("c", "rho", "tau", "sigma_q", "alpha_emp", "beta_emp", "sigma_emp")),
  m4  = list(type = "sv", sv = sv_m4, pars = c("c", "rho", "tau", "mu_h", "phi", "omega")),
  m9  = list(type = "sv", sv = sv_m9, pars = c("c", "rho", "tau", "mu_h", "phi", "omega", "rho_EI")),
  m10 = list(type = "sv", sv = sv_m10, pars = c("c", "rho", "tau", "alpha", "beta", "mu_h", "phi", "omega"))
)

# Split a draws matrix into per-draw parameter lists (only the variables we need).
# m7: kappa for COVID quarters outside the training sample (index > ctx$n_train) is
# not in the fit; it is drawn from its prior lognormal(0, 1), one draw per posterior
# draw, which is the model's own predictive for those quarters.
draw_list <- function(model, draws, ctx = list()) {
  dm <- as_dm(draws)
  pars <- model_defs[[model]]$pars
  cols <- list()
  for (v in pars) {
    if (model == "m7" && v == "kappa") next
    cols[[v]] <- dcol(dm, v)
  }
  S <- nrow(dm)
  if (model == "m7") {
    n_cov <- length(ctx$covid_idx)
    kf <- matrix(NA_real_, S, n_cov)
    n_in <- sum(ctx$covid_idx <= ctx$n_train)
    if (n_in > 0) kf[, seq_len(n_in)] <- dcol(dm, "kappa")[, seq_len(n_in)]
    if (n_in < n_cov) {
      kf[, (n_in + 1):n_cov] <- matrix(rlnorm(S * (n_cov - n_in), 0, 1), S)
    }
    cols$kappa_full <- kf
  }
  lapply(seq_len(S), function(s) lapply(cols, function(m) unname(m[s, ])))
}

# Per-draw, per-period log predictive densities.
# Returns list(full = S x T, sub = S x T or NULL, ess_min = S vector or NULL,
#              ess_frac = S x T or NULL).
#   full: log p(y_t | y_{1:t-1}, theta_s)
#   sub : for m8/m10 the GDP-only marginal log p(y^GDP_t | y_{1:t-1}, theta_s);
#         for m11 log p(y^GDP_q | y_{<q}, employment through the quarter, theta_s).
#   ess_frac (SV models): ESS of the particle weights at t divided by N_t.
# N (SV models): a particle count, or a length-T schedule N_t fixed in advance (see
# pf_schedule() in R/lib/lfo.R); the estimate is unbiased for any fixed schedule.
filter_loglik <- function(model, draws, y, ctx = list(), N = 1000L, ess_thr = 0.5,
                          seed = NULL, h1_sd_max = 3) {
  env <- load_filters()
  def <- model_defs[[model]]
  if (is.null(ctx$n_train)) ctx$n_train <- nrow(y)
  if (!is.null(seed)) set.seed(seed)
  pl <- draw_list(model, draws, ctx)
  S <- length(pl); T <- nrow(y)
  full <- matrix(NA_real_, S, T)
  sub <- NULL; ess_min <- NULL; ess_frac <- NULL
  if (def$type == "kf") {
    nb <- NULL
    for (s in seq_len(S)) {
      m <- def$ssm(pl[[s]], y, ctx)
      if (s == 1) {
        nb <- max(m$blk)
        if (nb > 1) check_block_diag(m$H, m$blk)
        if (nb > 1) sub <- matrix(NA_real_, S, T)
      }
      o <- env$kf_loglik_cpp(m$y, m$Z, m$d, m$H, m$Tm, m$cst, m$Q, m$a1, m$P1,
                             as.integer(m$blk), nb)
      full[s, ] <- rowSums(o)
      if (nb > 1) sub[s, ] <- if (model == "m11") o[, 2] else o[, 1]
    }
  } else {
    m1 <- def$sv(pl[[1]], y, ctx)
    if (m1$nsub < ncol(y)) sub <- matrix(NA_real_, S, T)
    ess_min <- numeric(S); ess_frac <- matrix(NA_real_, S, T)
    for (s in seq_len(S)) {
      m <- sv_prep(def$sv(pl[[s]], y, ctx))
      o <- env$rbpf_cpp(m$y, m$Z, m$d, m$Linv, m$logdetC, m$logdetC_sub,
                        as.integer(m$nsub), m$c, m$rho, m$tau, m$mu_h, m$phi, m$omega,
                        as.integer(N), ess_thr, h1_sd_max)
      full[s, ] <- o[, 1]
      if (!is.null(sub)) sub[s, ] <- o[, 2]
      ess_min[s] <- min(o[, 3])
      ess_frac[s, ] <- o[, 3] / o[, 4]
    }
  }
  list(full = full, sub = sub, ess_min = ess_min, ess_frac = ess_frac)
}

# Inverse Cholesky factor of C and log|C| (and of C's leading nsub x nsub block, whose
# Cholesky factor is the leading block of chol(C)).
sv_prep <- function(m) {
  ns <- m$nsub
  Lc <- t(chol(m$C))
  m$Linv <- forwardsolve(Lc, diag(nrow(Lc)))
  m$logdetC <- 2 * sum(log(diag(Lc)))
  m$logdetC_sub <- 2 * sum(log(diag(Lc)[seq_len(ns)]))
  m
}

# Conditional-on-h Kalman (the RBPF's inner step) for one draw and one h path.
kf_given_h <- function(model, p, y, h, ctx = list()) {
  env <- load_filters()
  m <- sv_prep(model_defs[[model]]$sv(p, y, ctx))
  env$kf_given_h_cpp(y, h, m$Z, m$d, m$Linv, m$logdetC, m$logdetC_sub,
                     as.integer(m$nsub), m$c, m$rho, m$tau)
}

check_block_diag <- function(H, blk) {
  for (i in seq_along(blk)) for (j in seq_along(blk)) {
    if (blk[i] != blk[j] && any(abs(H[i, j, ]) > 0))
      stop("H is not block-diagonal across observation blocks")
  }
}

# Single-draw convenience: total log-likelihood of y under the Kalman filter.
kf_total <- function(model, p, y, ctx = list()) {
  env <- load_filters()
  m <- model_defs[[model]]$ssm(p, y, ctx)
  sum(env$kf_loglik_cpp(m$y, m$Z, m$d, m$H, m$Tm, m$cst, m$Q, m$a1, m$P1,
                        as.integer(m$blk), max(m$blk)))
}
