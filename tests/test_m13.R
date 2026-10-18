# Tests for m13 (stan/m13_serial_latent_sv.stan: m12 + SV on the latent innovation) and its
# filters (R/lib/filter_m13.R, R/lib/filter_m13.cpp). Run from the repository root:
#   Rscript tests/test_m13.R
#
# The Stan programs' functions (kf_m13, ffbs_m13_rng; m12's kf_m12) are compiled from
# temporary copies of the .stan files (the model executables in stan/ are not touched).
#
# (i)    Stan kf_m13 log-likelihood == dense multivariate-normal log-density of the stacked y
#        (simulated from the model given non-constant g and h) for several psi, including
#        psi = 0, psi_E = psi_I and large |psi|, a COVID-like jump in g, rho_EI near -1 and
#        tiny measurement variances. The dense covariance is written from the model
#        definition (mu an AR(1) with time-varying innovation variance exp(g_t) and
#        Var(mu_1) = exp(g_1)/(1 - rho^2); eps a VAR(1) with diag(psi), innovation covariance
#        R_t, eps_1 ~ N(0, V_1)), independently of the quasi-differenced filter.
# (ii)   g constant at log(tau^2): kf_m13 reproduces m12's kf_m12 exactly (bit for bit with
#        tau_t = tau supplied directly; to rounding through exp(g / 2)), in Stan and in C++.
# (iii)  FFBS (Stan ffbs_m13_rng) matches the dense Gaussian smoother p(mu | y, g, h, theta)
#        with non-constant g: means, variances and lag-1 covariances within Monte Carlo error.
# (iv)   The C++ conditional filter (used inside the RBPF) == Stan kf_m13 (all 6 columns).
# (v)    RBPF: omega, omega_g -> 0 reproduces the exact Kalman result (bootstrap proposal,
#        scalar N and a length-T particle schedule; adapted proposal with omega_g = 0, and
#        within Monte Carlo error at omega_g = 1e-10); with omega_g = 0 and
#        mu_g = log(tau^2) it reproduces m12's RBPF for the same seed (both proposals).
# (vi)   RBPF estimate of p(y_{1:T} | theta) with omega, omega_g > 0 agrees with a
#        brute-force Monte Carlo average over (g, h) paths, for both g proposals.
# (vii)  Registration with the LFO machinery (R/lib/lfo.R): pf_schedule() and filter_par()
#        run m13 through filter_loglik("m13", ...); m12 still dispatches to its own filter.
# (viii) Numerical robustness: finite output with rho_EI rounding to -1 (z_EI = -20).
# (ix)   The C++ conditional filter == the fitted Stan program's log_lik_pred on posterior
#        draws from output/fit_m13.rds (skipped if that fit is absent).

suppressPackageStartupMessages({
  library(testthat)
  library(mvtnorm)
  library(cmdstanr)
})
options(cmdstanr_no_ver_check = TRUE)
source("R/lib/filters.R")
source("R/lib/lfo.R")
source("R/lib/filter_m12.R")
source("R/lib/filter_m13.R")
load_filters(); load_filter_m12(); load_filter_m13()
set.seed(20261017)

# ---- Stan functions -----------------------------------------------------------
expose_copy <- function(path) {
  td <- tempfile("stanfn"); dir.create(td)
  f <- file.path(td, basename(path))
  file.copy(path, f)
  mod <- cmdstan_model(f, compile_standalone = TRUE, quiet = TRUE)
  mod$functions
}
cat("Compiling Stan functions (temporary copies)...\n")
fn13 <- expose_copy("stan/m13_serial_latent_sv.stan")
fn12 <- expose_copy("stan/m12_struct_sv_serial.stan")

# ---- dense references -----------------------------------------------------------
stack <- function(y) as.vector(t(y))
Cmat <- function(r) matrix(c(1, r, 0, r, 1, 0, 0, 0, 1), 3, 3)
R_t <- function(h_t, r) { D <- diag(exp(h_t / 2)); D %*% Cmat(r) %*% D }
# Cov of mu_{1:T}: Var(mu_1) = exp(g_1)/(1-rho^2), Var(mu_t) = rho^2 Var(mu_{t-1}) + exp(g_t),
# Cov(mu_t, mu_s) = rho^(t-s) Var(mu_s) for s <= t.
tv_ar1_cov <- function(rho, g) {
  T <- length(g)
  v <- numeric(T)
  v[1] <- exp(g[1]) / (1 - rho^2)
  for (t in seq_len(T)[-1]) v[t] <- rho^2 * v[t - 1] + exp(g[t])
  S <- matrix(0, T, T)
  for (t in seq_len(T)) for (s in seq_len(t)) S[t, s] <- S[s, t] <- rho^(t - s) * v[s]
  S
}
dense_m13_moments <- function(p, h, g) {
  T <- nrow(h); r <- tanh(p$z_EI); Psi <- diag(p$psi)
  V <- vector("list", T)
  V[[1]] <- R_t(h[1, ], r) / (1 - outer(p$psi, p$psi))
  for (t in seq_len(T)[-1]) V[[t]] <- Psi %*% V[[t - 1]] %*% Psi + R_t(h[t, ], r)
  Se <- matrix(0, 3 * T, 3 * T)
  for (t in seq_len(T)) for (s in seq_len(t)) {
    B <- diag(p$psi^(t - s), 3) %*% V[[s]]                    # Cov(eps_t, eps_s)
    it <- (t - 1) * 3 + 1:3; is <- (s - 1) * 3 + 1:3
    Se[it, is] <- B; Se[is, it] <- t(B)
  }
  Smu <- tv_ar1_cov(p$rho, g)
  A <- kronecker(diag(T), matrix(1, 3, 1))
  list(m = rep(p$c / (1 - p$rho), 3 * T), S = A %*% Smu %*% t(A) + Se, A = A, Smu = Smu,
       mmu = rep(p$c / (1 - p$rho), T))
}
dense_m13 <- function(y, p, h, g) {
  d <- dense_m13_moments(p, h, g)
  dmvnorm(stack(y), d$m, d$S, log = TRUE)
}
dense_smoother <- function(y, p, h, g) {
  d <- dense_m13_moments(p, h, g)
  K <- d$Smu %*% t(d$A) %*% solve(d$S)
  list(mean = as.vector(d$mmu + K %*% (stack(y) - d$m)),
       cov = d$Smu - K %*% d$A %*% d$Smu)
}
sim_m13 <- function(T, p, h, g) {
  r <- tanh(p$z_EI)
  tt <- exp(g / 2)
  mu <- numeric(T)
  mu[1] <- rnorm(1, p$c / (1 - p$rho), tt[1] / sqrt(1 - p$rho^2))
  for (t in 2:T) mu[t] <- p$c + p$rho * mu[t - 1] + tt[t] * rnorm(1)
  eps <- matrix(0, T, 3)
  V1 <- R_t(h[1, ], r) / (1 - outer(p$psi, p$psi))
  eps[1, ] <- as.vector(t(chol(V1)) %*% rnorm(3))
  for (t in 2:T) eps[t, ] <- p$psi * eps[t - 1, ] + as.vector(t(chol(R_t(h[t, ], r))) %*% rnorm(3))
  mu + eps
}
sim_h <- function(T, mu_h = c(2.2, 1.8, 1.4), sd = 0.8) {
  sweep(matrix(rnorm(T * 3, 0, sd), T, 3), 2, mu_h, `+`)
}
sim_g <- function(T, mu_g = 1.5, sd = 1) mu_g + rnorm(T, 0, sd)
stan_kf <- function(y, h, g, p) fn13$kf_m13(y, h, exp(0.5 * g), p$c, p$rho, p$z_EI, p$psi)
cpp_kf  <- function(y, h, g, p) kf_m13_given(p, y, h, g)

base_p <- list(c = 1.3, rho = 0.45, z_EI = atanh(0.3))
psi_cases <- list(
  zero      = c(0, 0, 0),
  m5_like   = c(-0.5, -0.4, -0.2),
  equal_EI  = c(0.35, 0.35, -0.6),
  large     = c(0.95, -0.9, 0.5),
  mixed     = c(-0.8, 0.6, 0.1)
)

test_that("(i) Stan kf_m13 == dense MVN log-density (non-constant g and h)", {
  T <- 9
  for (nm in names(psi_cases)) {
    p <- c(base_p, list(psi = psi_cases[[nm]]))
    h <- sim_h(T); g <- sim_g(T)
    y <- sim_m13(T, p, h, g)
    expect_equal(sum(stan_kf(y, h, g, p)[, 1]), dense_m13(y, p, h, g), tolerance = 1e-9,
                 info = nm)
    d <- dense_m13_moments(p, h, g)
    expect_equal(stan_kf(y, h, g, p)[1, 1],
                 dmvnorm(y[1, ], d$m[1:3], d$S[1:3, 1:3], log = TRUE), tolerance = 1e-10,
                 info = nm)
    # every period: conditional of y_t given y_{1:t-1}
    x <- stack(y)
    for (t in 2:T) {
      past <- 1:(3 * (t - 1)); cur <- 3 * (t - 1) + 1:3
      Scp <- d$S[cur, past] %*% solve(d$S[past, past])
      ref_t <- dmvnorm(x[cur], as.vector(d$m[cur] + Scp %*% (x[past] - d$m[past])),
                       d$S[cur, cur] - Scp %*% d$S[past, cur], log = TRUE)
      expect_equal(stan_kf(y, h, g, p)[t, 1], ref_t, tolerance = 1e-8,
                   info = paste(nm, "t =", t))
    }
  }
  # COVID-like latent jump: g jumps by 5 in one period (latent sd x 12)
  p <- c(base_p, list(psi = psi_cases$m5_like))
  h <- sim_h(10); g <- rep(1, 10); g[6] <- 6; g[7] <- 4
  y <- sim_m13(10, p, h, g)
  expect_equal(sum(stan_kf(y, h, g, p)[, 1]), dense_m13(y, p, h, g), tolerance = 1e-9)
  # very small latent variance in some periods
  g2 <- g; g2[c(2, 3, 8)] <- -12
  expect_equal(sum(stan_kf(y, h, g2, p)[, 1]), dense_m13(y, p, h, g2), tolerance = 1e-8)
  # rho_EI near -1 (z_EI = -5)
  for (nm in c("zero", "m5_like", "large")) {
    p <- list(c = 1.3, rho = 0.45, z_EI = -5, psi = psi_cases[[nm]])
    h <- sim_h(T); g <- sim_g(T)
    y <- sim_m13(T, p, h, g)
    expect_equal(sum(stan_kf(y, h, g, p)[, 1]), dense_m13(y, p, h, g), tolerance = 1e-7,
                 info = paste("rho_EI near -1,", nm))
  }
  # tiny measurement variances in some periods (h = -15)
  p <- c(base_p, list(psi = psi_cases$m5_like))
  h <- matrix(c(-15, -14, -16, 1, 2, 0.5), 6, 3, byrow = TRUE)[rep(1:2, 3), ]
  g <- sim_g(6)
  y <- sim_m13(6, p, h, g)
  expect_equal(sum(stan_kf(y, h, g, p)[, 1]), dense_m13(y, p, h, g), tolerance = 1e-7)
})

test_that("(ii) g constant at log(tau^2): kf_m13 reproduces m12's kf_m12 exactly", {
  T <- 40
  tau <- 3.2
  for (nm in names(psi_cases)) {
    p <- c(base_p, list(psi = psi_cases[[nm]]))
    h <- sim_h(T)
    y <- sim_m13(T, p, h, rep(log(tau^2), T))
    f12 <- fn12$kf_m12(y, h, p$c, p$rho, tau, p$z_EI, p$psi)
    # bit for bit when tau_t = tau is passed directly
    expect_identical(fn13$kf_m13(y, h, rep(tau, T), p$c, p$rho, p$z_EI, p$psi), f12, info = nm)
    # through g = log(tau^2), exp(g / 2) = tau to rounding
    expect_equal(stan_kf(y, h, rep(log(tau^2), T), p), f12, tolerance = 1e-13, info = nm)
    # C++: m13 conditional filter == m12 conditional filter
    p12 <- c(p, list(tau = tau))
    expect_identical(kf_m13_given(p, y, h, tau_t = rep(tau, T)), kf_m12_given_h(p12, y, h),
                     info = nm)
  }
})

test_that("(iii) FFBS matches the dense smoother (non-constant g)", {
  T <- 8
  nsim <- 20000
  for (nm in c("m5_like", "large")) {
    p <- c(base_p, list(psi = psi_cases[[nm]]))
    h <- sim_h(T); g <- sim_g(T, sd = 1.2)
    y <- sim_m13(T, p, h, g)
    f <- stan_kf(y, h, g, p)
    draws <- t(vapply(seq_len(nsim), function(i) fn13$ffbs_m13_rng(f, seed = i), numeric(T)))
    ref <- dense_smoother(y, p, h, g)
    se_m <- sqrt(diag(ref$cov) / nsim)
    expect_true(all(abs(colMeans(draws) - ref$mean) < 4.5 * se_m), info = nm)
    Sh <- cov(draws)
    se_S <- sqrt((outer(diag(ref$cov), diag(ref$cov)) + ref$cov^2) / nsim)
    lag01 <- abs(outer(1:T, 1:T, "-")) <= 1
    expect_true(all(abs(Sh - ref$cov)[lag01] < 4.5 * se_S[lag01]), info = nm)
    expect_true(max(abs(Sh - ref$cov) / se_S) < 5, info = nm)
    cat(sprintf("\n  FFBS %s: worst mean error %.2f MC se, worst cov error %.2f MC se\n", nm,
                max(abs(colMeans(draws) - ref$mean) / se_m), max(abs(Sh - ref$cov) / se_S)))
  }
})

test_that("(iv) C++ conditional filter == Stan kf_m13", {
  T <- 30
  for (nm in names(psi_cases)) {
    p <- c(base_p, list(psi = psi_cases[[nm]]))
    h <- sim_h(T); g <- sim_g(T)
    y <- sim_m13(T, p, h, g)
    expect_equal(cpp_kf(y, h, g, p), stan_kf(y, h, g, p), tolerance = 1e-11, info = nm)
    p2 <- p; p2$rho_EI <- tanh(p$z_EI); p2$z_EI <- NULL
    expect_equal(cpp_kf(y, h, g, p2)[, 1], stan_kf(y, h, g, p)[, 1], tolerance = 1e-10, info = nm)
  }
})

mk_dm <- function(p, with_z = FALSE) {
  v <- c(c = p$c, rho = p$rho, mu_g = p$mu_g, phi_g = p$phi_g, omega_g = p$omega_g,
         setNames(p$mu_h, sprintf("mu_h[%d]", 1:3)), setNames(p$phi, sprintf("phi[%d]", 1:3)),
         setNames(p$omega, sprintf("omega[%d]", 1:3)), rho_EI = tanh(p$z_EI),
         setNames(p$psi, sprintf("psi[%d]", 1:3)))
  if (with_z) v <- c(v, z_EI = p$z_EI)
  matrix(v, 1, dimnames = list(NULL, names(v)))
}
mk_dm12 <- function(p, tau) {
  v <- c(c = p$c, rho = p$rho, tau = tau,
         setNames(p$mu_h, sprintf("mu_h[%d]", 1:3)), setNames(p$phi, sprintf("phi[%d]", 1:3)),
         setNames(p$omega, sprintf("omega[%d]", 1:3)), rho_EI = tanh(p$z_EI),
         setNames(p$psi, sprintf("psi[%d]", 1:3)))
  matrix(v, 1, dimnames = list(NULL, names(v)))
}

test_that("(v) RBPF: omega, omega_g -> 0 gives the Kalman result; omega_g = 0 gives m12's RBPF", {
  T <- 60
  mu_h <- c(2.4, 1.9, 1.5); phi <- c(0.9, 0.5, -0.3)
  for (nm in c("m5_like", "large", "zero")) {
    p <- c(base_p, list(psi = psi_cases[[nm]], mu_h = mu_h, phi = phi, omega = rep(1e-10, 3),
                        mu_g = 2.1, phi_g = 0.7, omega_g = 1e-10))
    hc <- matrix(mu_h, T, 3, byrow = TRUE); gc <- rep(p$mu_g, T)
    y <- sim_m13(T, p, hc, gc)
    ex <- stan_kf(y, hc, gc, p)[, 1]
    # bootstrap g proposal, scalar N and a length-T schedule (as from pf_schedule)
    o <- filter_loglik_m13(mk_dm(p), y, N = 50L, seed = 1, prop = 0L)
    expect_equal(as.vector(o$full), ex, tolerance = 1e-6, info = nm)
    expect_null(o$sub)
    sched <- as.integer(rep(c(50L, 120L, 80L), length.out = T))
    o <- filter_loglik_m13(mk_dm(p), y, N = sched, seed = 2, prop = 0L)
    expect_equal(as.vector(o$full), ex, tolerance = 1e-6, info = paste(nm, "schedule"))
    expect_equal(as.vector(o$ess_frac), rep(1, T), tolerance = 1e-6)
    # adapted proposal with omega_g = 0 exactly: no g draws, exact
    p0 <- p; p0$omega_g <- 0
    o <- filter_loglik_m13(mk_dm(p0), y, N = sched, seed = 3, prop = 1L)
    expect_equal(as.vector(o$full), ex, tolerance = 1e-6, info = paste(nm, "adapted, omega_g = 0"))
    # adapted proposal at omega_g = 1e-10: unbiased but not exact (the importance ratio
    # N / mixture-proposal is not constant), so within Monte Carlo error
    o <- filter_loglik_m13(mk_dm(p), y, N = 20000L, seed = 4, prop = 1L)
    dif <- as.vector(o$full) - ex
    cat(sprintf("\n  adapted proposal, omega_g = 1e-10, %s: max |per-period error| %.4f, total %.4f\n",
                nm, max(abs(dif)), sum(dif)))
    expect_lt(max(abs(dif)), 0.02)
    expect_lt(abs(sum(dif)), 0.1)
  }
  # omega_g = 0, mu_g = log(tau^2), omega > 0: same particle scheme and RNG stream as m12's RBPF
  tau <- 3.2
  p <- c(base_p, list(psi = psi_cases$m5_like, mu_h = mu_h, phi = c(0.95, 0.8, 0.9),
                      omega = c(0.3, 0.5, 0.4), mu_g = log(tau^2), phi_g = 0.6, omega_g = 0))
  h <- sim_h(T)
  y <- sim_m13(T, p, h, rep(log(tau^2), T))
  o12 <- filter_loglik_m12(mk_dm12(p, tau), y, N = 500L, seed = 7)
  for (pr in 0:1) {
    o13 <- filter_loglik_m13(mk_dm(p), y, N = 500L, seed = 7, prop = pr)
    expect_equal(o13$full, o12$full, tolerance = 1e-10, info = paste("prop", pr))
    expect_equal(o13$ess_min, o12$ess_min, tolerance = 1e-8, info = paste("prop", pr))
  }
  sched <- as.integer(rep(c(500L, 1500L, 800L), length.out = T))
  o12 <- filter_loglik_m12(mk_dm12(p, tau), y, N = sched, seed = 8)
  o13 <- filter_loglik_m13(mk_dm(p), y, N = sched, seed = 8)
  expect_equal(o13$full, o12$full, tolerance = 1e-10, info = "schedule")
})

test_that("(vi) RBPF estimate of p(y | theta) agrees with brute-force MC over (g, h)", {
  T <- 5
  p <- c(base_p, list(psi = c(-0.5, 0.4, -0.2), mu_h = c(2.2, 1.8, 1.4),
                      phi = c(0.9, 0.6, 0.8), omega = c(0.4, 0.5, 0.3),
                      mu_g = 1.6, phi_g = 0.8, omega_g = 0.9))
  hh <- sim_h(T); gg <- c(1.5, 1.2, 4.5, 2.5, 1.8)          # a latent-volatility jump at t = 3
  y <- sim_m13(T, p, hh, gg)
  M <- 200000
  sim_paths <- function() {
    h <- matrix(0, T, 3); g <- numeric(T)
    h[1, ] <- p$mu_h + p$omega / sqrt(1 - p$phi^2) * rnorm(3)
    g[1] <- p$mu_g + p$omega_g / sqrt(1 - p$phi_g^2) * rnorm(1)
    for (t in 2:T) {
      h[t, ] <- p$mu_h + p$phi * (h[t - 1, ] - p$mu_h) + p$omega * rnorm(3)
      g[t] <- p$mu_g + p$phi_g * (g[t - 1] - p$mu_g) + p$omega_g * rnorm(1)
    }
    list(h = h, g = g)
  }
  set.seed(11)
  lv <- vapply(seq_len(M), function(j) {
    s <- sim_paths(); sum(cpp_kf(y, s$h, s$g, p)[, 1])
  }, 0)
  bf <- lse(lv) - log(M)
  w <- exp(lv - max(lv)); bf_se <- sd(w) / mean(w) / sqrt(M)
  for (pr in 0:1) {
    rb <- vapply(1:20, function(s) sum(filter_loglik_m13(mk_dm(p), y, N = 5000L, seed = s,
                                                         prop = pr)$full), 0)
    rb_mean <- lse(rb) - log(length(rb))
    cat(sprintf("\n  prop %d: brute force log p(y|theta) = %.4f (se %.4f); RBPF = %.4f (sd over runs %.4f)\n",
                pr, bf, bf_se, rb_mean, sd(rb)))
    expect_lt(abs(rb_mean - bf), 4 * sqrt(bf_se^2 + var(rb) / length(rb)) + 0.005)
  }
})

test_that("(vii) registration: pf_schedule / filter_par run m13; m12 still dispatches", {
  register_m12(); register_m13()
  register_m13()                                     # idempotent
  expect_identical(model_defs$m13$type, "sv")
  expect_true(all(c(M13_PARS, "z_EI") %in% model_defs$m13$pars))
  T <- 40
  p <- c(base_p, list(psi = psi_cases$m5_like, mu_h = c(2.2, 1.8, 1.4), phi = c(0.9, 0.6, 0.8),
                      omega = c(0.3, 0.4, 0.3), mu_g = 1.6, phi_g = 0.8, omega_g = 0.6))
  y <- sim_m13(T, p, sim_h(T, sd = 0.4), sim_g(T, 1.6, 0.6))
  dm <- do.call(rbind, lapply(1:6, function(i) { q <- p; q$c <- p$c + 0.02 * i; mk_dm(q, TRUE) }))
  pf <- list(n_base = 200L, ess_target = 100, n_max = 2000L, pilot_draws = 4L, pilot_n = 200L, q = 0.1)
  Nt <- pf_schedule("m13", dm, y, list(), pf, seed = 5L, cores = 2L)
  expect_length(Nt, T)
  expect_true(all(Nt >= pf$n_base & Nt <= pf$n_max))
  r <- filter_par("m13", dm, y, list(), Nt, seed = 6L, cores = 2L)
  expect_equal(dim(r$full), c(nrow(dm), T))
  expect_true(all(is.finite(r$full)))
  expect_equal(dim(r$ess_frac), c(nrow(dm), T))
  # same result as calling the m13 filter directly on one chunk
  d1 <- filter_loglik("m13", dm[1:2, , drop = FALSE], y, N = 300L, seed = 9L)
  d2 <- filter_loglik_m13(dm[1:2, , drop = FALSE], y, N = 300L, seed = 9L)
  expect_identical(d1$full, d2$full)
  # m12 is unaffected by the m13 wrapper
  p12 <- p; p12$tau <- 3
  dm12 <- mk_dm12(p, 3)
  expect_identical(filter_loglik("m12", dm12, y, N = 300L, seed = 3L)$full,
                   filter_loglik_m12(dm12, y, N = 300L, seed = 3L)$full)
})

test_that("(viii) numerical robustness at rho_EI = -1 in double precision", {
  T <- 12
  p <- list(c = 1.3, rho = 0.45, z_EI = -20, psi = c(-0.5, -0.4, -0.2))
  expect_equal(tanh(p$z_EI), -1)
  h <- sim_h(T); g <- sim_g(T)
  y <- sim_m13(T, list(c = 1.3, rho = 0.45, z_EI = -5, psi = p$psi), h, g)
  f <- stan_kf(y, h, g, p)
  expect_true(all(is.finite(f)))
  expect_true(all(f[, 3] > 0) && all(f[-1, 6] > 0))
  expect_equal(cpp_kf(y, h, g, p), f, tolerance = 1e-10)
})

test_that("(ix) C++ conditional filter == Stan log_lik_pred on fitted draws", {
  fp <- "output/fit_m13.rds"
  if (!file.exists(fp)) skip("output/fit_m13.rds not available")
  fit <- readRDS(fp)
  gdp <- readRDS("data/processed/gdp_growth.rds")
  y <- as.matrix(gdp[, c("g_E", "g_I", "g_P")])
  T <- nrow(y)
  dm <- as_dm(posterior::subset_draws(
    fit$draws(c(M13_PARS, "z_EI", "h", "g", "log_lik_pred"), format = "draws_matrix"),
    draw = c(1, 1500, 3000, 4500, 7999)))
  for (s in seq_len(nrow(dm))) {
    p <- m13_draw_list(dm[s, , drop = FALSE])[[1]]
    h <- matrix(dcol(dm, "h")[s, ], T, 3)
    g <- as.vector(dcol(dm, "g")[s, ])
    expect_equal(kf_m13_given(p, y, h, g)[, 1], unname(dcol(dm, "log_lik_pred")[s, ]),
                 tolerance = 1e-6, info = paste("draw", s))
  }
})

cat("\nAll m13 tests finished.\n")
