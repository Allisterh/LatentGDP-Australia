# Tests for the marginalised (Kalman / FFBS) Stan programs of m3 and m10.
# Run from the repository root:
#   Rscript tests/test_m3_m10.R
#
# The Stan functions are compiled and exposed with cmdstanr (compile_standalone),
# so the code under test is exactly the code the sampler runs.
#
# m3 (stan/m3_corr.stan):
#   (a) kf_m3: sum of log p(y_t | y_{1:t-1}, theta) == dense multivariate-normal
#       log-density of the stacked observations (simulated data, several Omega),
#       including an Omega close to singular.
#   (b) ffbs_rng: draws of mu match the dense Gaussian smoother p(mu | y, theta)
#       in mean and covariance (Monte Carlo tolerance).
#   (c) The LFO Kalman filter in R/lib/filters.R (ssm_m3) gives the same
#       per-period densities as the Stan program's log_lik_pred on posterior
#       draws of output/fit_m3.rds (the LFO comparison scores the model that was
#       fitted).
# m10 tests are further down.

suppressPackageStartupMessages({
  library(testthat)
  library(mvtnorm)
  library(cmdstanr)
  library(posterior)
})
source("R/lib/filters.R")
set.seed(20261016)

# ---- dense references (written from the model definitions) ------------------
ar1_cov <- function(T, rho, tau) {
  v0 <- tau^2 / (1 - rho^2)
  v0 * rho^abs(outer(seq_len(T), seq_len(T), "-"))
}
stack <- function(y) as.vector(t(y))                      # y[1, ], y[2, ], ...
# y_t = d + Z mu_t + e_t, e_t ~ N(0, H_t); mu ~ N(mu_mean 1, Smu).
dense_parts <- function(T, mu_mean, Smu, Z, d, Hlist) {
  k <- length(Z)
  A <- kronecker(diag(T), matrix(Z, k, 1))                 # (T k) x T
  m <- as.vector(A %*% rep(mu_mean, length.out = T)) + rep(d, T)
  S <- A %*% Smu %*% t(A)
  for (t in seq_len(T)) {
    ii <- (t - 1) * k + seq_len(k)
    S[ii, ii] <- S[ii, ii] + Hlist[[t]]
  }
  list(A = A, m = m, S = S)
}
dense_ll <- function(y, ...) {
  dp <- dense_parts(nrow(y), ...)
  dmvnorm(stack(y), dp$m, dp$S, log = TRUE)
}
# Posterior of mu | y: mean and covariance from the joint Gaussian.
dense_smoother <- function(y, mu_mean, Smu, Z, d, Hlist) {
  dp <- dense_parts(nrow(y), mu_mean, Smu, Z, d, Hlist)
  Cxy <- Smu %*% t(dp$A)
  G <- Cxy %*% solve(dp$S)
  list(mean = as.vector(rep(mu_mean, length.out = nrow(y)) + G %*% (stack(y) - dp$m)),
       cov = Smu - G %*% t(Cxy))
}
corr_from_chol <- function(L) { O <- L %*% t(L); O }
ffbs_check <- function(draws, sm, label) {
  # Mean: max |z| over t; covariance: max abs error relative to the diagonal scale.
  n <- nrow(draws)
  se <- sqrt(diag(sm$cov) / n)
  z <- (colMeans(draws) - sm$mean) / se
  expect_lt(max(abs(z)), 4.5, label = paste(label, "FFBS mean z"))
  ec <- cov(draws)
  sc <- sqrt(outer(diag(sm$cov), diag(sm$cov)))
  expect_lt(max(abs(ec - sm$cov) / sc), 0.06, label = paste(label, "FFBS cov error"))
}

# ---- m3 ---------------------------------------------------------------------
# Expose the Stan program's own functions (kf_*, ffbs_rng) to R.
stan_functions <- function(file) {
  # Compiled into a temporary directory so the project's executables are untouched.
  mod <- cmdstan_model(file, dir = tempdir(), compile_standalone = TRUE,
                       force_recompile = TRUE)
  mod$functions
}
f3 <- stan_functions("stan/m3_corr.stan")

test_that("m3 (a): kf_m3 == dense MVN", {
  T0 <- 10
  y <- matrix(rnorm(T0 * 3, 3, 4), T0, 3)
  cases <- list(
    list(c = 1.4, rho = -0.3, tau = 3.3, sigma = c(4, 3, 2.5),
         Om = matrix(c(1, 0.3, -0.2, 0.3, 1, 0.25, -0.2, 0.25, 1), 3)),
    list(c = 0.2, rho = 0.85, tau = 1.1, sigma = c(1.5, 2.5, 0.7),
         Om = matrix(c(1, -0.6, 0.4, -0.6, 1, -0.1, 0.4, -0.1, 1), 3)),
    # nearly singular Omega and a tiny measurement variance
    list(c = 2.0, rho = 0.05, tau = 4.0, sigma = c(4.4, 3.1, 0.05),
         Om = matrix(c(1, 0.98, 0.97, 0.98, 1, 0.985, 0.97, 0.985, 1), 3))
  )
  for (i in seq_along(cases)) {
    p <- cases[[i]]
    L <- diag(p$sigma) %*% t(chol(p$Om))
    H <- L %*% t(L)
    ref <- dense_ll(y, p$c / (1 - p$rho), ar1_cov(T0, p$rho, p$tau), rep(1, 3),
                    rep(0, 3), rep(list(H), T0))
    o <- f3$kf_m3(y, L, p$c, p$rho, p$tau)
    expect_equal(sum(o[, 1]), ref, tolerance = 1e-8, label = paste("m3 case", i))
    # per-period decomposition agrees with filters.R (the LFO filter)
    lfo <- kf_total("m3", list(c = p$c, rho = p$rho, tau = p$tau, sigma = p$sigma,
                               Omega = as.vector(p$Om)), y)
    expect_equal(sum(o[, 1]), lfo, tolerance = 1e-8)
  }
})

test_that("m3 (b): ffbs_rng matches the dense smoother", {
  T0 <- 8
  y <- matrix(rnorm(T0 * 3, 3, 4), T0, 3)
  p <- list(c = 1.4, rho = 0.6, tau = 2.3, sigma = c(4, 3, 2.5),
            Om = matrix(c(1, 0.3, -0.2, 0.3, 1, 0.25, -0.2, 0.25, 1), 3))
  L <- diag(p$sigma) %*% t(chol(p$Om))
  H <- L %*% t(L)
  sm <- dense_smoother(y, p$c / (1 - p$rho), ar1_cov(T0, p$rho, p$tau), rep(1, 3),
                       rep(0, 3), rep(list(H), T0))
  fo <- f3$kf_m3(y, L, p$c, p$rho, p$tau)
  n <- 20000
  dr <- t(vapply(seq_len(n), function(i) as.vector(f3$ffbs_rng(fo, p$c, p$rho, p$tau, seed = i)),
                 numeric(T0)))
  ffbs_check(dr, sm, "m3")
  # filtered moments at T equal the smoother's at T
  expect_equal(fo[T0, 2], sm$mean[T0], tolerance = 1e-8)
  expect_equal(fo[T0, 3], sm$cov[T0, T0], tolerance = 1e-8)
})

test_that("m3 (c): LFO filter (filters.R) == Stan log_lik_pred on fitted draws", {
  if (!file.exists("output/fit_m3.rds")) skip("output/fit_m3.rds not available")
  fit <- readRDS("output/fit_m3.rds")
  if (!"log_lik_pred" %in% fit$metadata()$stan_variables) skip("fit predates log_lik_pred")
  gdp <- readRDS("data/processed/gdp_growth.rds")
  y <- as.matrix(gdp[, c("g_E", "g_I", "g_P")])
  dm <- as_dm(subset_draws(fit$draws(c(model_defs$m3$pars, "log_lik_pred"),
                                     format = "draws_matrix"),
                           draw = c(1, 2100, 4200, 6300, 7999)))
  o <- filter_loglik("m3", dm, y)
  expect_equal(o$full, unname(dcol(dm, "log_lik_pred")), tolerance = 1e-6)
})

# ---- m10 --------------------------------------------------------------------
# stan/m10_labour_sv.stan: 4 signals (E, I, P, hours), SV on every measurement error
# and on the latent innovation (log-variance g_t). LFO filter: R/lib/filter_m10.R/.cpp.
#   (d) Stan kf_tv == dense MVN given the (g, h) paths (incl. extreme h), and == the
#       Rcpp conditional filter m10_kf_given_cpp (the RBPF's inner step).
#   (e) Stan ffbs_tv_rng matches the dense smoother p(mu | y, g, h, theta).
#   (f) RBPF with every omega -> 0 reproduces the m8 Kalman filter (joint and GDP-only
#       marginal): exactly with the bootstrap proposal, within Monte Carlo error with
#       the adapted proposal for g.
#   (g) RBPF with latent SV active == brute-force Monte Carlo over g paths of the exact
#       conditional likelihood (checks the proposal correction), for both proposals.
#   (h) Old m10 as a special case (omega_g = 0, mu_g = log tau^2) agrees with the
#       R/lib/filters.R RBPF within Monte Carlo error.
#   (i) Stan log_lik_pred on posterior draws of output/fit_m10.rds == m10_kf_given.
source("R/lib/filter_m10.R")
register_m10_filter()
f10 <- stan_functions("stan/m10_labour_sv.stan")
env10 <- load_filter_m10()

tv_cov <- function(rho, qv) {
  T <- length(qv); V <- numeric(T); V[1] <- qv[1] / (1 - rho^2)
  for (t in 2:T) V[t] <- rho^2 * V[t - 1] + qv[t]
  S <- matrix(0, T, T)
  for (s in 1:T) for (t in s:T) { S[s, t] <- S[t, s] <- rho^(t - s) * V[s] }
  S
}
m10_dense <- function(y, p, g, h) {
  T <- nrow(y)
  dense_ll(y, p$c / (1 - p$rho), tv_cov(p$rho, exp(g)), c(1, 1, 1, p$beta),
           c(0, 0, 0, p$alpha), lapply(seq_len(T), function(t) diag(exp(h[t, ]))))
}
p10 <- list(c = 1.2, rho = 0.45, alpha = -0.7, beta = 0.85)

test_that("m10 (d): kf_tv == dense MVN and == m10_kf_given_cpp", {
  T0 <- 9
  y <- matrix(rnorm(T0 * 4, 3, 4), T0, 4)
  g <- rnorm(T0, 1.5, 1); g[5] <- 6
  h <- matrix(rnorm(T0 * 4, 1.5, 1), T0, 4)
  hx <- h; hx[3, ] <- c(-15, -14, -16, 1); hx[7, 2] <- -20       # tiny variances
  for (hh in list(h, hx)) {
    ref <- m10_dense(y, p10, g, hh)
    o <- f10$kf_tv(y, hh, c(1, 1, 1, p10$beta), c(0, 0, 0, p10$alpha), p10$c, p10$rho, exp(g))
    expect_equal(sum(o[, 1]), ref, tolerance = 1e-8)
    oc <- m10_kf_given(p10, y, hh, g)
    expect_equal(oc[, 1], o[, 1], tolerance = 1e-10)
  }
  # GDP-only column of the Rcpp filter at the last period == dense conditional
  ref_dp <- dense_parts(T0, p10$c / (1 - p10$rho), tv_cov(p10$rho, exp(g)),
                        c(1, 1, 1, p10$beta), c(0, 0, 0, p10$alpha),
                        lapply(seq_len(T0), function(t) diag(exp(h[t, ]))))
  past <- 1:(4 * (T0 - 1)); cur <- 4 * (T0 - 1) + 1:3
  x <- stack(y); Sm <- ref_dp$S
  Scp <- Sm[cur, past] %*% solve(Sm[past, past])
  cm <- ref_dp$m[cur] + Scp %*% (x[past] - ref_dp$m[past])
  cS <- Sm[cur, cur] - Scp %*% Sm[past, cur]
  expect_equal(m10_kf_given(p10, y, h, g)[T0, 2],
               dmvnorm(x[cur], as.vector(cm), cS, log = TRUE), tolerance = 1e-9)
})

test_that("m10 (e): ffbs_tv_rng matches the dense smoother", {
  T0 <- 8
  y <- matrix(rnorm(T0 * 4, 3, 4), T0, 4)
  g <- rnorm(T0, 1.5, 1); h <- matrix(rnorm(T0 * 4, 1.5, 1), T0, 4)
  sm <- dense_smoother(y, p10$c / (1 - p10$rho), tv_cov(p10$rho, exp(g)),
                       c(1, 1, 1, p10$beta), c(0, 0, 0, p10$alpha),
                       lapply(seq_len(T0), function(t) diag(exp(h[t, ]))))
  fo <- f10$kf_tv(y, h, c(1, 1, 1, p10$beta), c(0, 0, 0, p10$alpha), p10$c, p10$rho, exp(g))
  n <- 20000
  dr <- t(vapply(seq_len(n), function(i)
    as.vector(f10$ffbs_tv_rng(fo, p10$c, p10$rho, exp(g), seed = i)), numeric(T0)))
  ffbs_check(dr, sm, "m10")
})

mk_dm <- function(p, mu_g, phi_g, omega_g, mu_h, phi, omega) {
  v <- c(c = p$c, rho = p$rho, alpha = p$alpha, beta = p$beta, mu_g = mu_g,
         phi_g = phi_g, omega_g = omega_g, setNames(mu_h, sprintf("mu_h[%d]", 1:4)),
         setNames(phi, sprintf("phi[%d]", 1:4)), setNames(omega, sprintf("omega[%d]", 1:4)))
  matrix(v, 1, dimnames = list(NULL, names(v)))
}

test_that("m10 (f): RBPF with omega -> 0 reproduces the m8 Kalman filter", {
  Tt <- 40
  y <- matrix(rnorm(Tt * 4, 3, 4), Tt, 4)
  mh <- c(2.5, 2.0, 1.6, 1.1); tau <- 3.2
  dm <- mk_dm(p10, log(tau^2), 0.5, 1e-10, mh, c(0.9, 0.5, -0.3, 0.2), rep(1e-10, 4))
  p8 <- c(c = p10$c, rho = p10$rho, tau = tau, `sigma_gdp[1]` = exp(mh[1] / 2),
          `sigma_gdp[2]` = exp(mh[2] / 2), `sigma_gdp[3]` = exp(mh[3] / 2),
          sigma_U = exp(mh[4] / 2), alpha = p10$alpha, beta = p10$beta)
  o8 <- filter_loglik("m8", matrix(p8, 1, dimnames = list(NULL, names(p8))), y)
  o0 <- m10_filter_loglik(dm, y, N = 50L, seed = 1, prop = 0L)
  expect_equal(o0$full, o8$full, tolerance = 1e-6)
  expect_equal(o0$sub, o8$sub, tolerance = 1e-6)
  # adapted proposal (the default): the correction weights add Monte Carlo noise only
  tot <- sapply(1:10, function(s) sum(filter_loglik("m10", dm, y, N = 2000L, seed = s)$full))
  expect_lt(abs(mean(tot) - sum(o8$full)), 4 * sd(tot) / sqrt(10) + 0.01)
  expect_lt(sd(tot), 0.2)
})

test_that("m10 (g): RBPF == brute-force Monte Carlo over latent log-variance paths", {
  T0 <- 6
  y <- matrix(rnorm(T0 * 4, 3, 3), T0, 4)
  y[4, ] <- y[4, ] - 10                                  # a large common shock
  mh <- c(2.0, 1.8, 1.2, 1.0)
  mu_g <- 1.5; phi_g <- 0.8; om_g <- 0.6
  dm <- mk_dm(p10, mu_g, phi_g, om_g, mh, rep(0, 4), rep(1e-10, 4))
  h <- matrix(mh, T0, 4, byrow = TRUE)
  # brute force: g paths from their prior, exact conditional likelihood per path
  M <- 200000
  G <- matrix(0, M, T0)
  G[, 1] <- mu_g + om_g / sqrt(1 - phi_g^2) * rnorm(M)
  for (t in 2:T0) G[, t] <- mu_g + phi_g * (G[, t - 1] - mu_g) + om_g * rnorm(M)
  ll <- vapply(seq_len(M), function(i) sum(m10_kf_given(p10, y, h, G[i, ])[, 1]), 0)
  mx <- max(ll); w <- exp(ll - mx)
  ref <- mx + log(mean(w))
  ref_se <- sd(w) / sqrt(M) / mean(w)                    # delta-method se of log mean
  for (pr in 0:1) {
    tot <- sapply(1:20, function(s) sum(m10_filter_loglik(dm, y, N = 5000L, seed = s,
                                                          prop = pr)$full))
    se <- sqrt(var(tot) / 20 + ref_se^2)
    expect_lt(abs(mean(tot) - ref), 4 * se + 0.02, label = paste("proposal", pr))
  }
})

test_that("m10 (h): old m10 (constant tau) is the omega_g = 0 special case", {
  y4 <- as.matrix(readRDS("data/processed/gdp_growth_with_hours.rds")[1:60,
                  c("g_E", "g_I", "g_P", "g_U")])
  mh <- c(2.1, 2.1, 1.0, 1.5); ph <- c(0.9, 0.8, 0.95, 0.7); om <- c(0.3, 0.4, 0.3, 0.5)
  tau <- 2.5
  old <- c(c = p10$c, rho = p10$rho, tau = tau, alpha = p10$alpha, beta = p10$beta,
           setNames(mh, sprintf("mu_h[%d]", 1:4)), setNames(ph, sprintf("phi[%d]", 1:4)),
           setNames(om, sprintf("omega[%d]", 1:4)))
  dm_old <- matrix(old, 1, dimnames = list(NULL, names(old)))
  env_f <- load_filters()
  m_old <- sv_prep(list(y = y4, Z = c(1, 1, 1, p10$beta), d = c(0, 0, 0, p10$alpha),
                        C = diag(4), nsub = 3))
  set.seed(5)
  a <- replicate(20, sum(env_f$rbpf_cpp(y4, m_old$Z, m_old$d, m_old$Linv, m_old$logdetC,
                                         m_old$logdetC_sub, 3L, p10$c, p10$rho, tau, mh, ph,
                                         om, 2000L, 0.5, 3)[, 1]))
  dm <- mk_dm(p10, log(tau^2), 0, 0, mh, ph, om)
  b <- sapply(1:20, function(s) sum(m10_filter_loglik(dm, y4, N = 2000L, seed = s)$full))
  expect_lt(abs(mean(a) - mean(b)), 4 * sqrt(var(a) / 20 + var(b) / 20) + 0.02)
})

test_that("m10 (i): conditional filter == Stan log_lik_pred on fitted draws", {
  if (!file.exists("output/fit_m10.rds")) skip("output/fit_m10.rds not available")
  fit <- readRDS("output/fit_m10.rds")
  if (!all(c("g", "log_lik_pred") %in% fit$metadata()$stan_variables))
    skip("output/fit_m10.rds is not the latent-SV specification")
  y4 <- as.matrix(readRDS("data/processed/gdp_growth_with_hours.rds")[,
                  c("g_E", "g_I", "g_P", "g_U")])
  T <- nrow(y4)
  dm <- as_dm(subset_draws(fit$draws(c(M10_PARS, "g", "h", "log_lik_pred"),
                                     format = "draws_matrix"),
                           draw = c(1, 2100, 4200, 6300, 7999)))
  pl <- m10_draw_list(dm)
  for (s in seq_len(nrow(dm))) {
    h <- matrix(dcol(dm, "h")[s, ], T, 4)
    g <- unname(dcol(dm, "g")[s, ])
    o <- m10_kf_given(pl[[s]], y4, h, g)
    expect_equal(o[, 1], unname(dcol(dm, "log_lik_pred")[s, ]), tolerance = 1e-6)
  }
})

cat("\nAll m3 / m10 tests finished.\n")
