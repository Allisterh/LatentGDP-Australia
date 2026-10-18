# Tests for R/lib/filters.R.  Run from the repository root:
#   Rscript tests/test_filters.R
#
# (i)   Kalman log-likelihood == dense multivariate-normal log-density of the stacked
#       observations (small T, simulated data) for every linear-Gaussian spec
#       (m1, m2, m3, m5, m6, m7, m8, m11). The dense references below are written from
#       each model's definition, independently of the state-space code.
# (ii)  RBPF with omega -> 0 reproduces the Kalman result (m4 -> m2, m9 -> m6,
#       m10 -> m8), including the GDP-only marginal column for m10.
# (iii) RBPF Monte Carlo variability across seeds is small relative to the model
#       differences the comparison reports.
# (iv)  The RBPF's conditional-on-h Kalman step equals the Stan programs'
#       log_lik_pred (m4, m9) on posterior draws from output/fit_m*.rds.
# (v)   Per-period decomposition: the sum of per-period predictive densities is the
#       joint density, and the block (GDP-only) column is the correct marginal.
# (vi)  Numerical robustness of the RBPF step: exact at extreme log-variances
#       (dense MVN at h = -15), and finite / sensible on the m9 posterior draw with
#       phi_P = 0.9998 that broke the first implementation (log p-hat = 1e74 from
#       cancellation in Q - P q^2 / g with particles at h ~ -50).
# (vii) Per-period particle-count schedules N_t (used by R/03_compare.R), for all three
#       RBPFs (rbpf_cpp for m4/m9, R/lib/filter_m12.cpp, R/lib/filter_m10.cpp):
#       a constant schedule reproduces the scalar-N output exactly (same seed); with
#       omega -> 0 a varying schedule reproduces the dense multivariate-normal
#       log-density of simulated data (m9, m12 with AR(1) errors, the re-specified m10
#       with 4 signals); with omega > 0 the estimate of p(y_{1:T} | theta) under a varying
#       schedule (N_t from 5 to 400) is unbiased: the mean of p-hat over many runs agrees
#       with a brute-force Monte Carlo average over the volatility paths; and
#       pf_schedule() raises N only where the pilot ESS collapses.
# (iv) covers m4 and m9 only: the re-specified m10 has its own filter, checked against
#       its log_lik_pred in tests/test_m3_m10.R (i), and m12 in tests/test_m12.R (vii).

suppressPackageStartupMessages({
  library(testthat)
  library(mvtnorm)
})
source("R/lib/filters.R")
source("R/lib/lfo.R")
source("R/lib/filter_m12.R")     # filter_loglik_m12(), kf_m12_given_h()
source("R/lib/filter_m10.R")     # m10_filter_loglik(), m10_kf_given() (not registered here:
                                 # (ii) uses the old constant-tau m10 of filters.R)
env <- load_filters()
set.seed(20261016)

# ---- dense references --------------------------------------------------------
ar1_cov <- function(T, rho, tau) {
  v0 <- tau^2 / (1 - rho^2)
  v0 * rho^abs(outer(seq_len(T), seq_len(T), "-"))
}
# Stacked vec(t(y)): y[1,], y[2,], ...
stack <- function(y) as.vector(t(y))
dense_scalar <- function(y, mu_mean, Smu, Z, d, Hlist) {
  T <- nrow(y); k <- ncol(y)
  A <- kronecker(diag(T), matrix(Z, k, 1))                     # (T k) x T
  m <- as.vector(A %*% rep(mu_mean, length.out = T)) + rep(d, T)
  S <- A %*% Smu %*% t(A)
  for (t in seq_len(T)) {
    ii <- (t - 1) * k + seq_len(k)
    S[ii, ii] <- S[ii, ii] + Hlist[[t]]
  }
  dmvnorm(stack(y), m, S, log = TRUE)
}

gdp <- readRDS("data/processed/gdp_growth.rds")
gdp4 <- readRDS("data/processed/gdp_growth_with_hours.rds")
y_full <- as.matrix(gdp[, c("g_E", "g_I", "g_P")])
y_full4 <- as.matrix(gdp4[, c("g_E", "g_I", "g_P", "g_U")])

T0 <- 9
sim_y <- function(T, k, scale = 4) matrix(rnorm(T * k, 3, scale), T, k)
y3 <- sim_y(T0, 3); y4 <- sim_y(T0, 4)

test_that("(i) Kalman == dense MVN for every linear-Gaussian spec", {
  # m1: iid latent
  p <- list(mu_bar = 2.7, tau = 3.1, sigma = c(4, 3, 2.5))
  ref <- dense_scalar(y3, p$mu_bar, diag(p$tau^2, T0), rep(1, 3), rep(0, 3),
                      rep(list(diag(p$sigma^2)), T0))
  expect_equal(kf_total("m1", p, y3), ref, tolerance = 1e-9)

  # m2: AR(1), stationary start
  p <- list(c = 1.4, rho = 0.45, tau = 3.3, sigma = c(4, 3, 2.5))
  ref <- dense_scalar(y3, p$c / (1 - p$rho), ar1_cov(T0, p$rho, p$tau), rep(1, 3),
                      rep(0, 3), rep(list(diag(p$sigma^2)), T0))
  expect_equal(kf_total("m2", p, y3), ref, tolerance = 1e-9)

  # m3: full error correlation (Omega passed column-major as in Stan output)
  Om <- matrix(c(1, 0.3, -0.2, 0.3, 1, 0.25, -0.2, 0.25, 1), 3)
  p <- list(c = 1.4, rho = -0.3, tau = 3.3, sigma = c(4, 3, 2.5), Omega = as.vector(Om))
  H <- diag(p$sigma) %*% Om %*% diag(p$sigma)
  ref <- dense_scalar(y3, p$c / (1 - p$rho), ar1_cov(T0, p$rho, p$tau), rep(1, 3),
                      rep(0, 3), rep(list(H), T0))
  expect_equal(kf_total("m3", p, y3), ref, tolerance = 1e-9)

  # m5: AR(1) measurement errors, stationary start; dense covariance written directly
  p <- list(c = 1.4, rho = 0.6, tau = 3.3, sigma = c(4, 3, 2.5), psi = c(-0.5, 0.4, -0.2))
  T <- T0; k <- 3
  Smu <- ar1_cov(T, p$rho, p$tau)
  S <- kronecker(Smu, matrix(1, k, k))
  for (j in 1:k) {
    Se <- ar1_cov(T, p$psi[j], p$sigma[j])
    ii <- seq(j, T * k, by = k)
    S[ii, ii] <- S[ii, ii] + Se
  }
  ref <- dmvnorm(stack(y3), rep(p$c / (1 - p$rho), T * k), S, log = TRUE)
  expect_equal(kf_total("m5", p, y3), ref, tolerance = 1e-9)
  # and against the Stan program's own factorisation, given mu (sanity check of the
  # quasi-differencing transcription): integrate mu by brute-force Monte Carlo is
  # unnecessary; the dense form above is the model's joint density.

  # m6: structural Sigma
  p <- list(c = 1.4, rho = 0.2, tau = 3.3, sigma = c(4, 3, 2.5), rho_EI = 0.35)
  s <- p$sigma
  H <- matrix(c(s[1]^2, p$rho_EI * s[1] * s[2], 0, p$rho_EI * s[1] * s[2], s[2]^2, 0,
                0, 0, s[3]^2), 3)
  ref <- dense_scalar(y3, p$c / (1 - p$rho), ar1_cov(T0, p$rho, p$tau), rep(1, 3),
                      rep(0, 3), rep(list(H), T0))
  expect_equal(kf_total("m6", p, y3), ref, tolerance = 1e-9)

  # m7: COVID scale on quarters 4 and 6
  ctx <- list(covid_idx = c(4L, 6L))
  p <- list(c = 1.4, rho = 0.2, tau = 3.3, sigma = c(4, 3, 2.5), kappa_full = c(5, 0.6))
  kap <- rep(1, T0); kap[ctx$covid_idx] <- p$kappa_full
  ref <- dense_scalar(y3, p$c / (1 - p$rho), ar1_cov(T0, p$rho, p$tau), rep(1, 3),
                      rep(0, 3), lapply(kap, function(kk) diag((p$sigma * kk)^2)))
  expect_equal(kf_total("m7", p, y3, ctx), ref, tolerance = 1e-9)

  # m8: 4 signals with hours loading
  p <- list(c = 1.4, rho = 0.2, tau = 3.3, sigma_gdp = c(4, 3, 2.5), sigma_U = 2,
            alpha = -0.8, beta = 0.8)
  ref <- dense_scalar(y4, p$c / (1 - p$rho), ar1_cov(T0, p$rho, p$tau),
                      c(1, 1, 1, p$beta), c(0, 0, 0, p$alpha),
                      rep(list(diag(c(p$sigma_gdp^2, p$sigma_U^2))), T0))
  expect_equal(kf_total("m8", p, y4), ref, tolerance = 1e-9)

  # m11: monthly latent, monthly employment, quarterly GDP = 3-month average
  Tm <- 14
  q_end <- c(3L, 6L, 9L, 12L)
  p <- list(c = 0.3, rho = 0.7, tau = 4, sigma_q = c(3, 2.5, 2), alpha_emp = -1.2,
            beta_emp = 0.6, sigma_emp = 5)
  emp <- rnorm(Tm, 1, 5); yq <- sim_y(length(q_end), 3)
  ym <- matrix(NA_real_, Tm, 4); ym[, 1] <- emp; ym[q_end, 2:4] <- yq
  # dense: observations = A mu + b + noise
  rows <- list(); b <- c(); nv <- c()
  for (m in 1:Tm) {
    a <- rep(0, Tm); a[m] <- p$beta_emp
    rows[[length(rows) + 1]] <- a; b <- c(b, p$alpha_emp); nv <- c(nv, p$sigma_emp^2)
    if (m %in% q_end) for (k in 1:3) {
      a <- rep(0, Tm); a[(m - 2):m] <- 1/3
      rows[[length(rows) + 1]] <- a; b <- c(b, 0); nv <- c(nv, p$sigma_q[k]^2)
    }
  }
  A <- do.call(rbind, rows)
  obs <- as.vector(t(ym)); obs <- obs[!is.na(obs)]           # same month-major order
  mu_mean <- rep(p$c / (1 - p$rho), Tm)
  ref <- dmvnorm(obs, as.vector(A %*% mu_mean) + b,
                 A %*% ar1_cov(Tm, p$rho, p$tau) %*% t(A) + diag(nv), log = TRUE)
  expect_equal(kf_total("m11", p, ym), ref, tolerance = 1e-9)
})

test_that("(v) per-period and block decomposition", {
  p <- list(c = 1.4, rho = 0.2, tau = 3.3, sigma_gdp = c(4, 3, 2.5), sigma_U = 2,
            alpha = -0.8, beta = 0.8)
  dm <- matrix(unlist(p), 1, dimnames = list(NULL, c("c", "rho", "tau",
         "sigma_gdp[1]", "sigma_gdp[2]", "sigma_gdp[3]", "sigma_U", "alpha", "beta")))
  o <- filter_loglik("m8", dm, y4)
  expect_equal(sum(o$full), kf_total("m8", p, y4), tolerance = 1e-10)
  # GDP-only column = Kalman on the 3 GDP measures with hours *in the past* only:
  # check period 1 (no past) against the 3-dim marginal N(mean, Z Z' v0 + H).
  v0 <- p$tau^2 / (1 - p$rho^2)
  ref1 <- dmvnorm(y4[1, 1:3], rep(p$c / (1 - p$rho), 3),
                  matrix(v0, 3, 3) + diag(p$sigma_gdp^2), log = TRUE)
  expect_equal(o$sub[1, 1], ref1, tolerance = 1e-10)
  # and the last period's GDP marginal against the dense conditional
  T <- nrow(y4)
  A <- kronecker(diag(T), matrix(c(1, 1, 1, p$beta), 4, 1))
  m <- as.vector(A %*% rep(p$c / (1 - p$rho), T)) + rep(c(0, 0, 0, p$alpha), T)
  S <- A %*% ar1_cov(T, p$rho, p$tau) %*% t(A) +
    kronecker(diag(T), diag(c(p$sigma_gdp^2, p$sigma_U^2)))
  past <- 1:(4 * (T - 1)); cur <- 4 * (T - 1) + 1:3
  x <- stack(y4)
  Scp <- S[cur, past] %*% solve(S[past, past])
  cm <- m[cur] + Scp %*% (x[past] - m[past])
  cS <- S[cur, cur] - Scp %*% S[past, cur]
  expect_equal(o$sub[1, T], dmvnorm(x[cur], as.vector(cm), cS, log = TRUE), tolerance = 1e-9)
})

test_that("(ii) RBPF with omega -> 0 reproduces the Kalman filter", {
  Tt <- 40
  y3l <- sim_y(Tt, 3); y4l <- sim_y(Tt, 4)
  mk <- function(v) { m <- matrix(unlist(v), 1); colnames(m) <- names(unlist(v)); m }
  # m4 -> m2
  mu_h <- c(2.5, 2.0, 1.6); phi <- c(0.9, 0.5, -0.3); om <- rep(1e-10, 3)
  p4 <- c(c = 1.2, rho = 0.4, tau = 3.2, `mu_h[1]` = mu_h[1], `mu_h[2]` = mu_h[2],
          `mu_h[3]` = mu_h[3], `phi[1]` = phi[1], `phi[2]` = phi[2], `phi[3]` = phi[3],
          `omega[1]` = om[1], `omega[2]` = om[2], `omega[3]` = om[3])
  dm4 <- matrix(p4, 1, dimnames = list(NULL, names(p4)))
  o4 <- filter_loglik("m4", dm4, y3l, N = 50L, seed = 1)
  k2 <- kf_total("m2", list(c = 1.2, rho = 0.4, tau = 3.2, sigma = exp(mu_h / 2)), y3l)
  expect_equal(sum(o4$full), k2, tolerance = 1e-6)
  # m9 -> m6
  dm9 <- cbind(dm4, rho_EI = -0.4)
  o9 <- filter_loglik("m9", dm9, y3l, N = 50L, seed = 1)
  k6 <- kf_total("m6", list(c = 1.2, rho = 0.4, tau = 3.2, sigma = exp(mu_h / 2),
                            rho_EI = -0.4), y3l)
  expect_equal(sum(o9$full), k6, tolerance = 1e-6)
  # m10 -> m8 (joint and GDP-only marginal)
  mh4 <- c(mu_h, 1.1)
  p10 <- c(c = 1.2, rho = 0.4, tau = 3.2, alpha = -0.7, beta = 0.9,
           setNames(mh4, sprintf("mu_h[%d]", 1:4)), setNames(c(phi, 0.2), sprintf("phi[%d]", 1:4)),
           setNames(rep(1e-10, 4), sprintf("omega[%d]", 1:4)))
  dm10 <- matrix(p10, 1, dimnames = list(NULL, names(p10)))
  o10 <- filter_loglik("m10", dm10, y4l, N = 50L, seed = 1)
  p8 <- list(c = 1.2, rho = 0.4, tau = 3.2, sigma_gdp = exp(mu_h / 2), sigma_U = exp(1.1 / 2),
             alpha = -0.7, beta = 0.9)
  dm8 <- matrix(unlist(p8), 1, dimnames = list(NULL, c("c", "rho", "tau",
          "sigma_gdp[1]", "sigma_gdp[2]", "sigma_gdp[3]", "sigma_U", "alpha", "beta")))
  o8 <- filter_loglik("m8", dm8, y4l)
  expect_equal(o10$full, o8$full, tolerance = 1e-6)
  expect_equal(o10$sub, o8$sub, tolerance = 1e-6)
})

test_that("(vi) RBPF step is stable at extreme log-variances", {
  # Conditional-on-h step vs dense MVN with tiny measurement variances (m9 structure).
  Tt <- 6
  yy <- sim_y(Tt, 3)
  p <- list(c = 1.2, rho = 0.4, tau = 3.2, rho_EI = 0.6)
  h <- matrix(c(-15, -14, -16, 1, 2, 0.5), Tt, 3, byrow = TRUE)[rep(1:2, 3), ]
  o <- kf_given_h("m9", p, yy, h)
  R <- lapply(seq_len(Tt), function(t) {
    s <- exp(h[t, ] / 2)
    matrix(c(s[1]^2, p$rho_EI * s[1] * s[2], 0, p$rho_EI * s[1] * s[2], s[2]^2, 0,
             0, 0, s[3]^2), 3)
  })
  ref <- dense_scalar(yy, p$c / (1 - p$rho), ar1_cov(Tt, p$rho, p$tau), rep(1, 3),
                      rep(0, 3), R)
  expect_equal(sum(o[, 1]), ref, tolerance = 1e-8)
  # The draw that produced log p-hat = 1e74 before the stable rewrite.
  bad <- c(c = 2.681707, rho = 0.2773076, tau = 3.552991, `mu_h[1]` = 2.125618,
           `mu_h[2]` = 2.291553, `mu_h[3]` = 2.288216, `phi[1]` = 0.9928036,
           `phi[2]` = 0.495468, `phi[3]` = 0.9997955, `omega[1]` = 0.1240946,
           `omega[2]` = 0.687032, `omega[3]` = 0.3701453, rho_EI = 0.2885642)
  dmb <- matrix(bad, 1, dimnames = list(NULL, names(bad)))
  tot <- sapply(1:5, function(seed) {
    ob <- filter_loglik("m9", dmb, y_full, N = 2000L, seed = seed)$full
    expect_true(all(is.finite(ob)))
    expect_lt(max(ob), 0)        # a density of a 3-vector in annualised % is < 1 here
    sum(ob)
  })
  expect_lt(sd(tot), 1)
})

test_that("(vii) particle-count schedules: exact at omega -> 0, unbiased otherwise", {
  # -- dense references with constant variances (omega -> 0) --
  Tt <- 12
  sched <- as.integer(c(40, 7, 300, 25, 25, 90, 3, 150, 60, 10, 200, 33))
  mh <- c(2.2, 1.7, 1.3); ph <- c(0.9, 0.5, -0.3)
  th <- list(c = 1.1, rho = 0.45, tau = 3.0)
  # m9 structure (E-I correlated, P independent)
  r <- 0.35
  y <- sim_y(Tt, 3)
  v9 <- c(c = th$c, rho = th$rho, tau = th$tau, setNames(mh, sprintf("mu_h[%d]", 1:3)),
          setNames(ph, sprintf("phi[%d]", 1:3)), setNames(rep(1e-10, 3), sprintf("omega[%d]", 1:3)),
          rho_EI = r)
  dm9 <- matrix(v9, 1, dimnames = list(NULL, names(v9)))
  s_ <- exp(mh / 2)
  H9 <- matrix(c(s_[1]^2, r * s_[1] * s_[2], 0, r * s_[1] * s_[2], s_[2]^2, 0, 0, 0, s_[3]^2), 3)
  ref9 <- dense_scalar(y, th$c / (1 - th$rho), ar1_cov(Tt, th$rho, th$tau), rep(1, 3), rep(0, 3),
                       rep(list(H9), Tt))
  o9 <- filter_loglik("m9", dm9, y, N = sched, seed = 3)
  expect_equal(sum(o9$full), ref9, tolerance = 1e-6)
  expect_equal(dim(o9$ess_frac), c(1L, Tt))
  # constant schedule == scalar N, bit for bit (same seed, omega > 0)
  dm9b <- dm9; dm9b[, grepl("^omega", colnames(dm9b))] <- c(0.3, 0.5, 0.4)
  a <- filter_loglik("m9", dm9b, y, N = 300L, seed = 5)
  b <- filter_loglik("m9", dm9b, y, N = rep(300L, Tt), seed = 5)
  expect_identical(a$full, b$full)
  expect_error(filter_loglik("m9", dm9b, y, N = c(10L, 20L), seed = 5), "length")

  # m12: AR(1) measurement errors; eps_1 ~ N(0, V1), V1 = R / (1 - psi psi'), stationary
  psi <- c(-0.5, 0.4, -0.2)
  V1 <- H9 / (1 - outer(psi, psi))
  S <- kronecker(ar1_cov(Tt, th$rho, th$tau), matrix(1, 3, 3))
  for (t in 1:Tt) for (u in 1:Tt) {
    lag <- abs(t - u)
    blk <- if (t >= u) diag(psi^lag, 3) %*% V1 else V1 %*% diag(psi^lag, 3)
    S[(t - 1) * 3 + 1:3, (u - 1) * 3 + 1:3] <- S[(t - 1) * 3 + 1:3, (u - 1) * 3 + 1:3] + blk
  }
  ref12 <- dmvnorm(stack(y), rep(th$c / (1 - th$rho), 3 * Tt), S, log = TRUE)
  dm12 <- cbind(dm9, matrix(psi, 1, dimnames = list(NULL, sprintf("psi[%d]", 1:3))))
  o12 <- filter_loglik_m12(dm12, y, N = sched, seed = 3)
  expect_equal(sum(o12$full), ref12, tolerance = 1e-6)
  a <- filter_loglik_m12(cbind(dm9b, dm12[, 14:16, drop = FALSE]), y, N = 250L, seed = 8)
  b <- filter_loglik_m12(cbind(dm9b, dm12[, 14:16, drop = FALSE]), y, N = rep(250L, Tt), seed = 8)
  expect_identical(a$full, b$full)

  # m10 (re-specified): 4 signals, latent variance exp(g_t), all omegas -> 0
  y4s <- sim_y(Tt, 4)
  mh4 <- c(mh, 1.0); mu_g <- 2.1
  v10 <- c(c = th$c, rho = th$rho, alpha = -0.7, beta = 0.9, mu_g = mu_g, phi_g = 0.6,
           omega_g = 1e-10, setNames(mh4, sprintf("mu_h[%d]", 1:4)),
           setNames(c(ph, 0.2), sprintf("phi[%d]", 1:4)), setNames(rep(1e-10, 4), sprintf("omega[%d]", 1:4)))
  dm10 <- matrix(v10, 1, dimnames = list(NULL, names(v10)))
  ref10 <- dense_scalar(y4s, th$c / (1 - th$rho), ar1_cov(Tt, th$rho, exp(mu_g / 2)),
                        c(1, 1, 1, 0.9), c(0, 0, 0, -0.7), rep(list(diag(exp(mh4))), Tt))
  o10 <- m10_filter_loglik(dm10, y4s, N = sched, seed = 3, prop = 0L)
  expect_equal(sum(o10$full), ref10, tolerance = 1e-6)
  dm10b <- dm10; dm10b[, "omega_g"] <- 0.7; dm10b[, grepl("^omega\\[", colnames(dm10b))] <- 0.3
  a <- m10_filter_loglik(dm10b, y4s, N = 200L, seed = 4)
  b <- m10_filter_loglik(dm10b, y4s, N = rep(200L, Tt), seed = 4)
  expect_identical(a$full, b$full)

  # -- unbiasedness under a varying schedule (omega > 0), against brute force --
  T5 <- 5
  sch5 <- as.integer(c(30, 400, 5, 120, 40))
  sim_path <- function(M, mu, phi, om) {
    K <- length(mu); H <- array(0, c(M, T5, K))
    for (k in 1:K) {
      H[, 1, k] <- mu[k] + om[k] / sqrt(1 - phi[k]^2) * rnorm(M)
      for (t in 2:T5) H[, t, k] <- mu[k] + phi[k] * (H[, t - 1, k] - mu[k]) + om[k] * rnorm(M)
    }
    H
  }
  check_unbiased <- function(label, run, lik_given, M = 100000, R = 3000, paths) {
    lv <- vapply(seq_len(M), function(j) lik_given(paths, j), 0)
    mx <- max(lv); w <- exp(lv - mx); ref <- mx + log(mean(w)); ref_se <- sd(w) / mean(w) / sqrt(M)
    tot <- vapply(seq_len(R), run, 0)
    mt <- max(tot); wt <- exp(tot - mt)
    est <- mt + log(mean(wt)); est_se <- sd(wt) / mean(wt) / sqrt(R)
    cat(sprintf("\n  (vii) %s: brute force %.4f (se %.4f); schedule RBPF log mean p-hat %.4f (se %.4f); mean log p-hat %.4f\n",
                label, ref, ref_se, est, est_se, mean(tot)))
    expect_lt(abs(est - ref), 4 * sqrt(ref_se^2 + est_se^2) + 0.005, label = label)
  }
  set.seed(77)
  y5 <- sim_y(T5, 3, 3); y5[3, ] <- y5[3, ] + c(5, -4, 1)      # measurement disagreement
  p9 <- list(c = 1.1, rho = 0.45, tau = 3.0, rho_EI = 0.35, mu_h = c(2.0, 1.6, 1.2),
             phi = c(0.9, 0.6, 0.8), omega = c(0.5, 0.6, 0.4))
  dmp <- matrix(c(p9$c, p9$rho, p9$tau, p9$mu_h, p9$phi, p9$omega, p9$rho_EI), 1,
                dimnames = list(NULL, colnames(dm9)))
  P9 <- sim_path(100000, p9$mu_h, p9$phi, p9$omega)
  check_unbiased("m9", function(s) sum(filter_loglik("m9", dmp, y5, N = sch5, seed = s)$full),
                 function(H, j) sum(kf_given_h("m9", p9, y5, H[j, , ])[, 1]), paths = P9)
  p12 <- c(p9, list(psi = c(-0.5, 0.4, -0.2)))
  dmp12 <- cbind(dmp, matrix(p12$psi, 1, dimnames = list(NULL, sprintf("psi[%d]", 1:3))))
  check_unbiased("m12", function(s) sum(filter_loglik_m12(dmp12, y5, N = sch5, seed = s)$full),
                 function(H, j) sum(kf_m12_given_h(p12, y5, H[j, , ])[, 1]), paths = P9)
  y54 <- sim_y(T5, 4, 3); y54[3, ] <- y54[3, ] - 12                  # common shock (latent SV)
  p10 <- list(c = 1.1, rho = 0.45, alpha = -0.7, beta = 0.9, mu_g = 1.5, phi_g = 0.8,
              omega_g = 0.6, mu_h = c(2.0, 1.6, 1.2, 1.0), phi = c(0.9, 0.6, 0.8, 0.5),
              omega = c(0.5, 0.6, 0.4, 0.3))
  dmp10 <- matrix(unlist(p10), 1, dimnames = list(NULL, colnames(dm10)))
  P10 <- sim_path(100000, c(p10$mu_h, p10$mu_g), c(p10$phi, p10$phi_g), c(p10$omega, p10$omega_g))
  check_unbiased("m10 (adapted proposal)",
                 function(s) sum(m10_filter_loglik(dmp10, y54, N = sch5, seed = s, prop = 1L)$full),
                 function(H, j) sum(m10_kf_given(p10, y54, H[j, , 1:4], H[j, , 5])[, 1]), paths = P10)

  # -- pf_schedule: raises N only where the pilot ESS collapses --
  T40 <- 40
  y40 <- sim_y(T40, 3, 2); y40[25, ] <- y40[25, ] + c(25, -20, 0)   # one wild disagreement
  dm40 <- dm9b[rep(1, 20), , drop = FALSE]
  pf <- list(n_base = 200L, ess_target = 100, n_max = 5000L, pilot_draws = 20L, pilot_n = 200L, q = 0.1)
  Nt <- pf_schedule("m9", dm40, y40, list(), pf, seed = 1, cores = 1)
  expect_length(Nt, T40)
  expect_true(all(Nt >= pf$n_base & Nt <= pf$n_max))
  expect_equal(which.max(Nt), 25L)
  expect_gt(Nt[25], pf$n_base)
})

# ---- tests that use the fitted models ---------------------------------------
fit_path <- function(nm) sprintf("output/fit_%s.rds", nm)
has_pred <- function(nm) {
  if (!file.exists(fit_path(nm))) return(FALSE)
  f <- readRDS(fit_path(nm))
  "log_lik_pred" %in% f$metadata()$stan_variables
}

test_that("(iv) conditional Kalman inside the RBPF == Stan log_lik_pred", {
  for (nm in c("m4", "m9")) {
    if (!has_pred(nm)) { skip(paste(nm, "fit with log_lik_pred not available")) }
    f <- readRDS(fit_path(nm))
    y <- if (nm == "m10") y_full4 else y_full
    pars <- model_defs[[nm]]$pars
    dm <- as_dm(posterior::subset_draws(f$draws(c(pars, "h", "log_lik_pred"),
                                                format = "draws_matrix"),
                                        draw = c(1, 1500, 3000, 4500, 7999)))
    T <- nrow(y); K <- ncol(y)
    for (s in seq_len(nrow(dm))) {
      p <- draw_list(nm, dm[s, , drop = FALSE])[[1]]
      h <- matrix(dcol(dm, "h")[s, ], T, K)               # Stan column-major
      o <- kf_given_h(nm, p, y, h)
      expect_equal(o[, 1], unname(dcol(dm, "log_lik_pred")[s, ]), tolerance = 1e-6,
                   info = paste(nm, "draw", s))
    }
  }
})

test_that("(iii) RBPF Monte Carlo error is small relative to model differences", {
  if (!all(file.exists(fit_path(c("m4", "m9", "m6"))))) skip("fits not available")
  # Five posterior draws per SV model, 10 independent filter runs each (N = 2000,
  # the particle count used by R/03_compare.R).
  sds <- c(); means <- list()
  for (nm in c("m4", "m9")) {
    f <- readRDS(fit_path(nm))
    dm <- as_dm(posterior::subset_draws(f$draws(model_defs[[nm]]$pars, format = "draws_matrix"),
                                        draw = c(10, 2000, 4000, 6000, 8000)))
    tot <- sapply(1:10, function(seed)
      rowSums(filter_loglik(nm, dm, y_full, N = 2000L, seed = seed)$full))
    sds <- c(sds, apply(tot, 1, sd))
    means[[nm]] <- rowMeans(tot)
  }
  f6 <- readRDS(fit_path("m6"))
  dm6 <- as_dm(posterior::subset_draws(f6$draws(model_defs$m6$pars, format = "draws_matrix"),
                                       draw = c(10, 2000, 4000, 6000, 8000)))
  k6 <- rowSums(filter_loglik("m6", dm6, y_full)$full)
  cat(sprintf("\n  RBPF sd of log p-hat(y_1:T | theta) across seeds, per draw (T = %d): %s\n",
              nrow(y_full), paste(sprintf("%.2f", sds), collapse = " ")))
  cat(sprintf("  log p(y | theta) means: m9 %s | m4 %s | m6 (Kalman) %s\n",
              paste(round(means$m9), collapse = " "), paste(round(means$m4), collapse = " "),
              paste(round(k6), collapse = " ")))
  # Full-sample, single-draw MC sd must be small in absolute terms (< 1 nat) ...
  expect_lt(max(sds), 1)
  # ... and small relative to the SV vs constant-variance difference (m9 vs m6).
  expect_lt(max(sds), 0.05 * abs(mean(means$m9) - mean(k6)))
})

# (x) m13 as R/03_compare.R wires it: register_m12(), register_m10_filter(), register_m13()
# in that order. filter_loglik("m13", ...) must equal the dense MVN likelihood of m12 when
# the latent and measurement volatilities are constant (omega, omega_g -> 0, mu_g = log
# tau^2; under a varying particle schedule, through filter_par's forked workers too), equal
# filter_loglik_m13() bit for bit, and leave the m12, m10 and m9 dispatch unchanged.
# The m13 filter's own tests (dense MVN with time-varying g and h, brute force, FFBS) are
# in tests/test_m13.R.
test_that("(x) m13 registered as in R/03_compare.R: dense MVN at constant variances, dispatch", {
  source("R/lib/filter_m13.R")
  Tt <- 12
  sched <- as.integer(c(40, 7, 300, 25, 25, 90, 3, 150, 60, 10, 200, 33))
  mh <- c(2.2, 1.7, 1.3); ph <- c(0.9, 0.5, -0.3); psi <- c(-0.5, 0.4, -0.2)
  th <- list(c = 1.1, rho = 0.45, tau = 3.0); r <- 0.35
  y <- sim_y(Tt, 3)
  s_ <- exp(mh / 2)
  H <- matrix(c(s_[1]^2, r * s_[1] * s_[2], 0, r * s_[1] * s_[2], s_[2]^2, 0, 0, 0, s_[3]^2), 3)
  V1 <- H / (1 - outer(psi, psi))
  S <- kronecker(ar1_cov(Tt, th$rho, th$tau), matrix(1, 3, 3))
  for (t in 1:Tt) for (u in 1:Tt) {
    lag <- abs(t - u)
    blk <- if (t >= u) diag(psi^lag, 3) %*% V1 else V1 %*% diag(psi^lag, 3)
    S[(t - 1) * 3 + 1:3, (u - 1) * 3 + 1:3] <- S[(t - 1) * 3 + 1:3, (u - 1) * 3 + 1:3] + blk
  }
  ref <- dmvnorm(stack(y), rep(th$c / (1 - th$rho), 3 * Tt), S, log = TRUE)
  v12 <- c(c = th$c, rho = th$rho, tau = th$tau, setNames(mh, sprintf("mu_h[%d]", 1:3)),
           setNames(ph, sprintf("phi[%d]", 1:3)), setNames(rep(1e-10, 3), sprintf("omega[%d]", 1:3)),
           rho_EI = r, setNames(psi, sprintf("psi[%d]", 1:3)), z_EI = atanh(r))
  # omega_g = 0 exactly: the adapted g proposal (the LFO default) is then exact; at
  # omega_g = 1e-10 it is only unbiased (see tests/test_m13.R (v)), the bootstrap one exact.
  v13 <- c(v12[c("c", "rho")], mu_g = log(th$tau^2), phi_g = 0.6, omega_g = 0,
           v12[!names(v12) %in% c("c", "rho", "tau")])
  dm12 <- matrix(v12, 1, dimnames = list(NULL, names(v12)))
  dm13 <- matrix(v13, 1, dimnames = list(NULL, names(v13)))
  base_fl <- filter_loglik
  base_md <- model_defs
  register_m12(); register_m10_filter(); register_m13()
  on.exit({ assign("filter_loglik", base_fl, envir = globalenv())
            assign("model_defs", base_md, envir = globalenv()) }, add = TRUE)
  expect_identical(model_defs$m13$type, "sv")
  expect_true(all(c("mu_g", "phi_g", "omega_g", "psi", "rho_EI", "z_EI") %in% model_defs$m13$pars))
  o <- filter_loglik("m13", dm13, y, N = sched, seed = 3)            # adapted, omega_g = 0
  expect_equal(sum(o$full), ref, tolerance = 1e-6)
  expect_equal(dim(o$ess_frac), c(1L, Tt))
  dm13e <- dm13; dm13e[, "omega_g"] <- 1e-10                          # bootstrap, omega_g > 0
  o <- filter_loglik_m13(dm13e, y, N = sched, seed = 3, prop = 0L)
  expect_equal(sum(o$full), ref, tolerance = 1e-6)
  # forked workers (filter_par) on 3 copies of the draw
  op <- filter_par("m13", dm13[c(1, 1, 1), , drop = FALSE], y, list(), sched, 9L, 2L)
  expect_equal(rowSums(op$full), rep(ref, 3), tolerance = 1e-6)
  # dispatch: m13 == direct call; m12 unchanged by registering m13 (omega, omega_g > 0)
  dm13b <- dm13; dm13b[, "omega_g"] <- 0.6; dm13b[, grepl("^omega\\[", colnames(dm13b))] <- 0.3
  expect_identical(filter_loglik("m13", dm13b, y, N = 300L, seed = 5)$full,
                   filter_loglik_m13(dm13b, y, N = 300L, seed = 5)$full)
  dm12b <- dm12; dm12b[, grepl("^omega\\[", colnames(dm12b))] <- 0.3
  expect_identical(filter_loglik("m12", dm12b, y, N = 300L, seed = 5)$full,
                   filter_loglik_m12(dm12b, y, N = 300L, seed = 5)$full)
  dm9 <- dm12b[, !grepl("^psi|z_EI", colnames(dm12b)), drop = FALSE]
  expect_identical(filter_loglik("m9", dm9, y, N = 300L, seed = 5)$full,
                   base_fl("m9", dm9, y, N = 300L, seed = 5)$full)
})

cat("\nAll filter tests finished.\n")
