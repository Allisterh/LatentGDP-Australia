# Tests for the tight-priors research variant of m13 (stan/m13_tight_priors.stan). Run from
# the repository root:
#   Rscript tests/test_m13_tight.R
#
# The variant changes only priors and the sampled coordinates (mu_bar instead of c; the
# stationary sds s_g, s_h instead of omega_g, omega); its functions block is a verbatim copy
# of m13's. The checks:
#   (i)   the functions and data blocks are byte-identical to stan/m13_serial_latent_sv.stan;
#   (ii)  kf_m13 compiled from the variant returns exactly m13's kf_m13 output (all six
#         columns) for the same (y, h, g, c, rho, z_EI, psi), for several psi and rho, and
#         FFBS draws agree for the same RNG seed;
#   (iii) the full log density: at the same point, mapped between the two coordinate systems
#         (mu_bar = c / (1 - rho), s = omega / sqrt(1 - phi^2), same g_raw, h_raw), the
#         variant's log density minus m13's (both without Jacobians, on the 1980-2026 data)
#         equals the difference of the two prior densities computed in R. So g, h and the
#         likelihood are identical and only the priors differ.
# The model executables in stan/ are not touched (temporary copies are compiled).

suppressPackageStartupMessages({
  library(testthat)
  library(cmdstanr)
})
options(cmdstanr_no_ver_check = TRUE)
set.seed(20261019)
F13 <- "stan/m13_serial_latent_sv.stan"
FT  <- "stan/m13_tight_priors.stan"

block <- function(f) { x <- readLines(f); x[grep("^functions \\{", x):(grep("^parameters \\{", x) - 1)] }
test_that("(i) functions and data blocks are byte-identical to m13's", {
  expect_identical(block(FT), block(F13))
})

copy_to_tmp <- function(path) {
  td <- tempfile("stantight"); dir.create(td)
  f <- file.path(td, basename(path)); file.copy(path, f); f
}
cat("Compiling Stan functions (temporary copies)...\n")
fn13 <- cmdstan_model(copy_to_tmp(F13), compile_standalone = TRUE, quiet = TRUE)$functions
fnT  <- cmdstan_model(copy_to_tmp(FT),  compile_standalone = TRUE, quiet = TRUE)$functions

test_that("(ii) kf_m13 and FFBS identical to m13's", {
  T <- 40
  for (psi in list(c(0, 0, 0), c(-0.45, -0.35, -0.2), c(0.95, -0.9, 0.5))) for (rho in c(0.05, 0.4, 0.97)) {
    y <- matrix(rnorm(T * 3, 3, 4), T, 3)
    h <- sweep(matrix(rnorm(T * 3, 0, 0.8), T, 3), 2, c(2.2, 1.8, 1.4), `+`)
    g <- 1.5 + rnorm(T); g[30] <- g[30] + 5
    c0 <- 2.6 * (1 - rho); z <- atanh(0.3)
    a <- fn13$kf_m13(y, h, exp(0.5 * g), c0, rho, z, psi)
    b <- fnT$kf_m13(y, h, exp(0.5 * g), c0, rho, z, psi)
    expect_identical(a, b)
    expect_identical(fn13$ffbs_m13_rng(a, seed = 7), fnT$ffbs_m13_rng(b, seed = 7))
  }
})

test_that("(iii) log density differs from m13's only by the prior", {
  source("R/production_config.R")
  gr <- load_production_growth()
  sd <- list(T = nrow(gr), K = 3, y = as.matrix(gr[, c("g_E", "g_I", "g_P")]))
  m13 <- cmdstan_model(copy_to_tmp(F13), compile_model_methods = TRUE, quiet = TRUE)
  mT  <- cmdstan_model(copy_to_tmp(FT),  compile_model_methods = TRUE, quiet = TRUE)
  # one-iteration runs only to get objects with the model methods (E-BFMI warnings expected)
  f13 <- suppressWarnings(m13$sample(data = sd, chains = 1, iter_warmup = 1, iter_sampling = 1,
                    refresh = 0, seed = 1, show_messages = FALSE, show_exceptions = FALSE))
  fT  <- suppressWarnings(mT$sample(data = sd, chains = 1, iter_warmup = 1, iter_sampling = 1,
                   refresh = 0, seed = 1, show_messages = FALSE, show_exceptions = FALSE))
  f13$init_model_methods(); fT$init_model_methods()
  T <- sd$T
  for (i in 1:5) {
    rho <- runif(1, 0.05, 0.95); c0 <- rnorm(1, 1.5, 0.5); phi_g <- runif(1, -0.95, 0.99)
    om_g <- runif(1, 0.1, 1.2); phi <- runif(3, -0.9, 0.99); om <- runif(3, 0.05, 0.8)
    mu_g <- rnorm(1, 1.5, 0.5); mu_h <- rnorm(3, 1.8, 0.5); z <- rnorm(1, 0, 0.3)
    psi <- runif(3, -0.6, 0.3); g_raw <- rnorm(T); h_raw <- matrix(rnorm(T * 3), T, 3)
    p13 <- list(c = c0, rho = rho, mu_g = mu_g, phi_g = phi_g, omega_g = om_g, mu_h = mu_h,
                phi = phi, omega = om, z_EI = z, psi = psi, g_raw = g_raw, h_raw = h_raw)
    s_g <- om_g / sqrt(1 - phi_g^2); s_h <- om / sqrt(1 - phi^2); mu_bar <- c0 / (1 - rho)
    pT <- list(mu_bar = mu_bar, rho = rho, mu_g = mu_g, phi_g = phi_g, s_g = s_g, mu_h = mu_h,
               phi = phi, s_h = s_h, z_EI = z, psi = psi, g_raw = g_raw, h_raw = h_raw)
    lp13 <- f13$log_prob(f13$unconstrain_variables(p13), jacobian = FALSE)
    lpT  <- fT$log_prob(fT$unconstrain_variables(pT), jacobian = FALSE)
    # Stan drops constants in ~ statements; compare the parts that differ, with the same
    # convention (kernels up to constants that are the same at every point): the difference
    # of the two log densities must be constant across points.
    prior13 <- dnorm(c0, 1.5, 1, log = TRUE) + dnorm(rho, 0.5, 0.3, log = TRUE) +
      dnorm(om_g, 0, 0.5, log = TRUE) + sum(dnorm(om, 0, 0.5, log = TRUE))
    priorT <- dnorm(mu_bar, 3, 2, log = TRUE) + dbeta(rho, 2.625, 2.625, log = TRUE) +
      dnorm(s_g, 0, 1, log = TRUE) + sum(dnorm(s_h, 0, 1, log = TRUE))
    dd <- (lpT - lp13) - (priorT - prior13)
    if (i == 1) d0 <- dd else expect_equal(dd, d0, tolerance = 1e-8)
  }
  cat(sprintf("(iii) log-density difference minus prior difference: constant %.6f over 5 points\n", d0))
})
cat("test_m13_tight.R done\n")
