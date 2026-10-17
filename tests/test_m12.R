# Tests for m12 (stan/m12_struct_sv_serial.stan) and its filters (R/lib/filter_m12.R,
# R/lib/filter_m12.cpp). Run from the repository root:
#   Rscript tests/test_m12.R
#
# The Stan program's own functions (kf_m12, ffbs_m12_rng) are compiled from a temporary
# copy of the .stan file (so the model executable in stan/ is not touched) and called from R.
#
# (i)   Stan kf_m12 log-likelihood == dense multivariate-normal log-density of the stacked
#       y (simulated from the model, given h) for several psi, including psi = 0,
#       psi_E = psi_I, large |psi|, rho_EI near -1 (z_EI = -5) and tiny measurement
#       variances (h = -15). The dense covariance is written from the model definition
#       (eps a VAR(1) with diag(psi), time-varying innovation covariance R_t, eps_1 ~ N(0, V_1)),
#       independently of the quasi-differenced filter.
# (ii)  psi = 0: m12's filter equals m9's (Stan kf_m9 and the RBPF's conditional step in
#       R/lib/filters.R) given the same h: per-period log-likelihood and filtered moments.
# (iii) FFBS (Stan ffbs_m12_rng) matches the dense Gaussian smoother p(mu | y, h, theta):
#       means, variances and lag-1 covariances within Monte Carlo error.
# (iv)  The C++ conditional filter (used inside the RBPF) equals Stan's kf_m12 (all 6 columns).
# (v)   RBPF with omega -> 0 reproduces the exact Kalman result; at psi = 0 it reproduces the
#       m9 RBPF of R/lib/filters.R for the same seed; with omega > 0 its estimate of
#       p(y_{1:T} | theta) agrees with a brute-force Monte Carlo average over h paths.
# (vi)  Numerical robustness: finite filter output with rho_EI rounding to -1 (z_EI = -20).
# (vii) The C++ conditional filter equals the fitted Stan program's log_lik_pred on
#       posterior draws from output/fit_m12.rds (skipped if that fit is absent).

suppressPackageStartupMessages({
  library(testthat)
  library(mvtnorm)
  library(cmdstanr)
})
options(cmdstanr_no_ver_check = TRUE)
source("R/lib/filters.R")
source("R/lib/filter_m12.R")
load_filters(); load_filter_m12()
set.seed(20261016)

# ---- Stan functions -----------------------------------------------------------
expose_copy <- function(path) {
  td <- tempfile("stanfn"); dir.create(td)
  f <- file.path(td, basename(path))
  file.copy(path, f)
  mod <- cmdstan_model(f, compile_standalone = TRUE, quiet = TRUE)
  mod$functions
}
cat("Compiling Stan functions (temporary copies)...\n")
fn12 <- expose_copy("stan/m12_struct_sv_serial.stan")
fn9  <- expose_copy("stan/m9_struct_sv.stan")

# ---- dense references -----------------------------------------------------------
ar1_cov <- function(T, rho, tau) {
  v0 <- tau^2 / (1 - rho^2)
  v0 * rho^abs(outer(seq_len(T), seq_len(T), "-"))
}
stack <- function(y) as.vector(t(y))
Cmat <- function(r) matrix(c(1, r, 0, r, 1, 0, 0, 0, 1), 3, 3)
R_t <- function(h_t, r) { D <- diag(exp(h_t / 2)); D %*% Cmat(r) %*% D }

# Joint law of the stacked y = (y_1', ..., y_T')' given h and theta, from the model definition.
dense_m12_moments <- function(p, h) {
  T <- nrow(h); r <- tanh(p$z_EI); Psi <- diag(p$psi)
  V <- vector("list", T)
  V[[1]] <- R_t(h[1, ], r) / (1 - outer(p$psi, p$psi))        # stationary analogue
  for (t in seq_len(T)[-1]) V[[t]] <- Psi %*% V[[t - 1]] %*% Psi + R_t(h[t, ], r)
  Se <- matrix(0, 3 * T, 3 * T)
  for (t in seq_len(T)) for (s in seq_len(t)) {
    B <- diag(p$psi^(t - s), 3) %*% V[[s]]                    # Cov(eps_t, eps_s)
    it <- (t - 1) * 3 + 1:3; is <- (s - 1) * 3 + 1:3
    Se[it, is] <- B; Se[is, it] <- t(B)
  }
  Smu <- ar1_cov(T, p$rho, p$tau)
  A <- kronecker(diag(T), matrix(1, 3, 1))
  list(m = rep(p$c / (1 - p$rho), 3 * T), S = A %*% Smu %*% t(A) + Se, A = A, Smu = Smu,
       mmu = rep(p$c / (1 - p$rho), T))
}
dense_m12 <- function(y, p, h) {
  d <- dense_m12_moments(p, h)
  dmvnorm(stack(y), d$m, d$S, log = TRUE)
}
dense_smoother <- function(y, p, h) {
  d <- dense_m12_moments(p, h)
  K <- d$Smu %*% t(d$A) %*% solve(d$S)
  list(mean = as.vector(d$mmu + K %*% (stack(y) - d$m)),
       cov = d$Smu - K %*% d$A %*% d$Smu)
}
# Simulate y from the model given h.
sim_m12 <- function(T, p, h) {
  r <- tanh(p$z_EI)
  mu <- numeric(T)
  mu[1] <- rnorm(1, p$c / (1 - p$rho), p$tau / sqrt(1 - p$rho^2))
  for (t in 2:T) mu[t] <- p$c + p$rho * mu[t - 1] + p$tau * rnorm(1)
  eps <- matrix(0, T, 3)
  V1 <- R_t(h[1, ], r) / (1 - outer(p$psi, p$psi))
  eps[1, ] <- as.vector(t(chol(V1)) %*% rnorm(3))
  for (t in 2:T) eps[t, ] <- p$psi * eps[t - 1, ] + as.vector(t(chol(R_t(h[t, ], r))) %*% rnorm(3))
  mu + eps
}
sim_h <- function(T, mu_h = c(2.2, 1.8, 1.4), sd = 0.8) {
  sweep(matrix(rnorm(T * 3, 0, sd), T, 3), 2, mu_h, `+`)
}
stan_kf <- function(y, h, p) fn12$kf_m12(y, h, p$c, p$rho, p$tau, p$z_EI, p$psi)
cpp_kf  <- function(y, h, p) kf_m12_given_h(p, y, h)

base_p <- list(c = 1.3, rho = 0.45, tau = 3.2, z_EI = atanh(0.3))
psi_cases <- list(
  zero      = c(0, 0, 0),
  m5_like   = c(-0.5, -0.4, -0.2),
  equal_EI  = c(0.35, 0.35, -0.6),
  large     = c(0.95, -0.9, 0.5),
  mixed     = c(-0.8, 0.6, 0.1)
)

test_that("(i) Stan kf_m12 == dense MVN log-density", {
  T <- 9
  for (nm in names(psi_cases)) {
    p <- c(base_p, list(psi = psi_cases[[nm]]))
    h <- sim_h(T)
    y <- sim_m12(T, p, h)
    expect_equal(sum(stan_kf(y, h, p)[, 1]), dense_m12(y, p, h), tolerance = 1e-9,
                 info = nm)
    # per-period decomposition: the first period alone is the marginal of y_1
    d <- dense_m12_moments(p, h)
    expect_equal(stan_kf(y, h, p)[1, 1],
                 dmvnorm(y[1, ], d$m[1:3], d$S[1:3, 1:3], log = TRUE), tolerance = 1e-10,
                 info = nm)
    # last period: conditional of y_T given y_{1:T-1}
    past <- 1:(3 * (T - 1)); cur <- 3 * (T - 1) + 1:3; x <- stack(y)
    Scp <- d$S[cur, past] %*% solve(d$S[past, past])
    ref_T <- dmvnorm(x[cur], as.vector(d$m[cur] + Scp %*% (x[past] - d$m[past])),
                     d$S[cur, cur] - Scp %*% d$S[past, cur], log = TRUE)
    expect_equal(stan_kf(y, h, p)[T, 1], ref_T, tolerance = 1e-8, info = nm)
  }
  # rho_EI near -1 (z_EI = -5: 1 + rho_EI = 9.1e-5)
  for (nm in c("zero", "m5_like", "large")) {
    p <- list(c = 1.3, rho = 0.45, tau = 3.2, z_EI = -5, psi = psi_cases[[nm]])
    h <- sim_h(T)
    y <- sim_m12(T, p, h)
    expect_equal(sum(stan_kf(y, h, p)[, 1]), dense_m12(y, p, h), tolerance = 1e-7,
                 info = paste("rho_EI near -1,", nm))
  }
  # tiny measurement variances in some periods (h = -15), as in test_filters.R (vi)
  p <- c(base_p, list(psi = psi_cases$m5_like))
  h <- matrix(c(-15, -14, -16, 1, 2, 0.5), 6, 3, byrow = TRUE)[rep(1:2, 3), ]
  y <- sim_m12(6, p, h)
  expect_equal(sum(stan_kf(y, h, p)[, 1]), dense_m12(y, p, h), tolerance = 1e-7)
})

test_that("(ii) psi = 0: m12 filter == m9 filter given the same h", {
  T <- 40
  p <- c(base_p, list(psi = c(0, 0, 0)))
  h <- sim_h(T)
  y <- sim_m12(T, p, h)
  f12 <- stan_kf(y, h, p)
  f9 <- fn9$kf_m9(y, h, p$c, p$rho, p$tau, p$z_EI)
  expect_equal(f12[, 1], f9[, 1], tolerance = 1e-10)
  expect_equal(f12[, 2], f9[, 2], tolerance = 1e-10)
  expect_equal(f12[, 3], f9[, 3], tolerance = 1e-10)
  # and the conditional step of the m9 RBPF in R/lib/filters.R
  k9 <- kf_given_h("m9", list(c = p$c, rho = p$rho, tau = p$tau, rho_EI = tanh(p$z_EI)), y, h)
  expect_equal(f12[, 1], k9[, 1], tolerance = 1e-10)
})

test_that("(iii) FFBS matches the dense smoother", {
  T <- 8
  nsim <- 20000
  for (nm in c("m5_like", "large")) {
    p <- c(base_p, list(psi = psi_cases[[nm]]))
    h <- sim_h(T)
    y <- sim_m12(T, p, h)
    f <- stan_kf(y, h, p)
    draws <- t(vapply(seq_len(nsim), function(i) fn12$ffbs_m12_rng(f, seed = i), numeric(T)))
    ref <- dense_smoother(y, p, h)
    # means: |difference| < 4.5 Monte Carlo se
    se_m <- sqrt(diag(ref$cov) / nsim)
    expect_true(all(abs(colMeans(draws) - ref$mean) < 4.5 * se_m), info = nm)
    # covariances (incl. variances and lag-1): se of a sample covariance ~
    # sqrt((S_ii S_jj + S_ij^2) / n)
    Sh <- cov(draws)
    se_S <- sqrt((outer(diag(ref$cov), diag(ref$cov)) + ref$cov^2) / nsim)
    lag01 <- abs(outer(1:T, 1:T, "-")) <= 1
    expect_true(all(abs(Sh - ref$cov)[lag01] < 4.5 * se_S[lag01]), info = nm)
    expect_true(max(abs(Sh - ref$cov) / se_S) < 5, info = nm)
  }
})

test_that("(iv) C++ conditional filter == Stan kf_m12", {
  T <- 30
  for (nm in names(psi_cases)) {
    p <- c(base_p, list(psi = psi_cases[[nm]]))
    h <- sim_h(T)
    y <- sim_m12(T, p, h)
    expect_equal(cpp_kf(y, h, p), stan_kf(y, h, p), tolerance = 1e-11, info = nm)
    # also when only rho_EI (not z_EI) is supplied, as in the LFO draws
    p2 <- p; p2$rho_EI <- tanh(p$z_EI); p2$z_EI <- NULL
    expect_equal(cpp_kf(y, h, p2)[, 1], stan_kf(y, h, p)[, 1], tolerance = 1e-10, info = nm)
  }
})

mk_dm <- function(p, with_z = FALSE) {
  v <- c(c = p$c, rho = p$rho, tau = p$tau,
         setNames(p$mu_h, sprintf("mu_h[%d]", 1:3)), setNames(p$phi, sprintf("phi[%d]", 1:3)),
         setNames(p$omega, sprintf("omega[%d]", 1:3)), rho_EI = tanh(p$z_EI),
         setNames(p$psi, sprintf("psi[%d]", 1:3)))
  if (with_z) v <- c(v, z_EI = p$z_EI)
  matrix(v, 1, dimnames = list(NULL, names(v)))
}

test_that("(v) RBPF: omega -> 0 gives the Kalman result; psi = 0 gives the m9 RBPF", {
  T <- 60
  mu_h <- c(2.4, 1.9, 1.5); phi <- c(0.9, 0.5, -0.3)
  for (nm in c("m5_like", "large", "zero")) {
    p <- c(base_p, list(psi = psi_cases[[nm]], mu_h = mu_h, phi = phi, omega = rep(1e-10, 3)))
    hc <- matrix(mu_h, T, 3, byrow = TRUE)
    y <- sim_m12(T, p, hc)
    o <- filter_loglik_m12(mk_dm(p), y, N = 50L, seed = 1)
    ex <- stan_kf(y, hc, p)[, 1]
    expect_equal(as.vector(o$full), ex, tolerance = 1e-6, info = nm)
    expect_null(o$sub)
  }
  # psi = 0, omega > 0: identical particle scheme and RNG stream to the m9 RBPF
  p <- c(base_p, list(psi = c(0, 0, 0), mu_h = mu_h, phi = c(0.95, 0.8, 0.9),
                      omega = c(0.3, 0.5, 0.4)))
  h <- sim_h(T)
  y <- sim_m12(T, p, h)
  dm <- mk_dm(p)
  o12 <- filter_loglik_m12(dm, y, N = 500L, seed = 7)
  dm9 <- dm[, !grepl("^psi", colnames(dm)), drop = FALSE]
  o9 <- filter_loglik("m9", dm9, y, N = 500L, seed = 7)
  expect_equal(o12$full, o9$full, tolerance = 1e-8)
  expect_equal(o12$ess_min, o9$ess_min, tolerance = 1e-8)
})

test_that("(v) RBPF estimate of p(y | theta) agrees with brute-force MC over h", {
  T <- 5
  p <- c(base_p, list(psi = c(-0.5, 0.4, -0.2), mu_h = c(2.2, 1.8, 1.4),
                      phi = c(0.9, 0.6, 0.8), omega = c(0.4, 0.5, 0.3)))
  hh <- sim_h(T)
  y <- sim_m12(T, p, hh)
  # brute force: log mean_j p(y | h^(j), theta) over h^(j) from the stationary SV prior
  M <- 100000
  sim_path <- function() {
    h <- matrix(0, T, 3)
    h[1, ] <- p$mu_h + p$omega / sqrt(1 - p$phi^2) * rnorm(3)
    for (t in 2:T) h[t, ] <- p$mu_h + p$phi * (h[t - 1, ] - p$mu_h) + p$omega * rnorm(3)
    h
  }
  set.seed(11)
  lv <- vapply(seq_len(M), function(j) sum(cpp_kf(y, sim_path(), p)[, 1]), 0)
  bf <- lse(lv) - log(M)
  w <- exp(lv - max(lv)); bf_se <- sd(w) / mean(w) / sqrt(M)      # delta-method se of log mean
  # RBPF: average of 20 independent runs at N = 5000 (log of mean p-hat is unbiased-ish)
  rb <- vapply(1:20, function(s) sum(filter_loglik_m12(mk_dm(p), y, N = 5000L, seed = s)$full), 0)
  rb_mean <- lse(rb) - log(length(rb))
  cat(sprintf("\n  brute force log p(y|theta) = %.4f (se %.4f); RBPF = %.4f (sd over runs %.4f)\n",
              bf, bf_se, rb_mean, sd(rb)))
  expect_lt(abs(rb_mean - bf), 4 * sqrt(bf_se^2 + var(rb) / length(rb)) + 0.005)
})

test_that("(vi) numerical robustness at rho_EI = -1 in double precision", {
  T <- 12
  p <- list(c = 1.3, rho = 0.45, tau = 3.2, z_EI = -20, psi = c(-0.5, -0.4, -0.2))
  expect_equal(tanh(p$z_EI), -1)
  h <- sim_h(T)
  y <- sim_m12(T, list(c = 1.3, rho = 0.45, tau = 3.2, z_EI = -5, psi = p$psi), h)
  f <- stan_kf(y, h, p)
  expect_true(all(is.finite(f)))
  expect_true(all(f[, 3] > 0) && all(f[-1, 6] > 0))
  expect_equal(cpp_kf(y, h, p), f, tolerance = 1e-10)
})

test_that("(vii) C++ conditional filter == Stan log_lik_pred on fitted draws", {
  fp <- "output/fit_m12.rds"
  if (!file.exists(fp)) skip("output/fit_m12.rds not available")
  fit <- readRDS(fp)
  gdp <- readRDS("data/processed/gdp_growth.rds")
  y <- as.matrix(gdp[, c("g_E", "g_I", "g_P")])
  T <- nrow(y)
  dm <- as_dm(posterior::subset_draws(
    fit$draws(c(M12_PARS, "z_EI", "h", "log_lik_pred"), format = "draws_matrix"),
    draw = c(1, 1500, 3000, 4500, 7999)))
  for (s in seq_len(nrow(dm))) {
    p <- m12_draw_list(dm[s, , drop = FALSE])[[1]]
    h <- matrix(dcol(dm, "h")[s, ], T, 3)
    expect_equal(kf_m12_given_h(p, y, h)[, 1], unname(dcol(dm, "log_lik_pred")[s, ]),
                 tolerance = 1e-6, info = paste("draw", s))
  }
})

cat("\nAll m12 tests finished.\n")
