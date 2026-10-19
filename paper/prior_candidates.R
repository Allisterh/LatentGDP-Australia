# paper/prior_candidates.R
#
# Prior predictive choice of tighter priors for m13 (NOTES.md, "Tighter priors for m13
# (prior predictive check and refit), 2026-10-19"). Run from the repository root:
#   Rscript paper/prior_candidates.R          (about a minute)
#
# Simulates the model forward from its priors (paper/ppc_m13.R's simulator: stationary
# starts, T = 186 quarters, 1980 Q1 - 2026 Q2, annualised units), 20,000 draws per prior set
# (seed 1, in chunks of 2,500), for m13's priors and three candidate sets:
#   C1  rho ~ Beta(2.625, 2.625) on (0, 1) (RDP 2014-12: mean 0.5, sd 0.2); mu_bar ~ N(3, 2),
#       c = mu_bar (1 - rho); SV blocks as m13
#   C2  C1 + (phi + 1)/2 ~ Beta(20, 1.5) on phi_g and phi_k (Kim, Shephard & Chib 1998)
#   C3  C1 + stationary sd of each log-variance s = omega / sqrt(1 - phi^2) ~ half-N(0, 1),
#       phi ~ U(-1, 1) (chosen; stan/m13_tight_priors.stan)
# and evaluates them against a plausibility target fixed before any candidate was simulated:
#   T1  90% interval of the sd of headline growth within [1, 15]
#   T2  99th percentile of max_t |headline growth| <= 100
#   T3  90% interval of mean headline growth within [-3, 9]
#   T4  each observed statistic inside the central 95% of its prior predictive
# Selection rule (fixed in advance): the first of C1, C2, C3 meeting T1-T4; if none does, the
# candidate meeting T4 with the smallest T2 value. Writes
#   paper/data/prior_candidates.csv      prior predictive quantiles by candidate
#   paper/data/prior_target.csv          target checks by candidate
#   paper/data/prior_tail_diagnosis.csv  m13 prior: what drives the tails (conditional
#                                        prior predictive and tail-draw attribution)
# No model is fitted and no posterior result is used.

stopifnot(file.exists("R/production_config.R"))

# ---- simulator and statistics (as paper/ppc_m13.R) ---------------------------------------
source("R/production_config.R")
g0 <- load_production_growth()
Y <- as.matrix(g0[, c("g_E", "g_I", "g_P")]); T <- nrow(Y)
dates <- as.Date(g0$date)
ac1 <- function(x) cor(x[-1], x[-length(x)])
pre <- dates < as.Date("2020-01-01")
stats_small <- function(y) {
  hl <- rowMeans(y)
  c(sd_head = sd(hl), mean_head = mean(hl), maxabs_head = max(abs(hl)),
    ac1_head = ac1(hl), sd_EI = sd(y[, 1] - y[, 2]), sd_head_pre2020 = sd(hl[pre]),
    maxabs_any = max(abs(y)))
}
obs <- stats_small(Y)
simulate <- function(par, S) {
  out <- array(NA_real_, c(S, T, 3)); r <- par$rho_EI
  h <- array(NA_real_, c(S, T, 3))
  for (k in 1:3) {
    h[, 1, k] <- par$mu_h[, k] + par$omega[, k] / sqrt(1 - par$phi[, k]^2) * rnorm(S)
    for (t in 2:T) h[, t, k] <- par$mu_h[, k] + par$phi[, k] * (h[, t - 1, k] - par$mu_h[, k]) +
        par$omega[, k] * rnorm(S)
  }
  gg <- matrix(NA_real_, S, T); mu <- matrix(NA_real_, S, T)
  gg[, 1] <- par$mu_g + par$omega_g / sqrt(1 - par$phi_g^2) * rnorm(S)
  for (t in 2:T) gg[, t] <- par$mu_g + par$phi_g * (gg[, t - 1] - par$mu_g) + par$omega_g * rnorm(S)
  mu[, 1] <- par$c / (1 - par$rho) + exp(gg[, 1] / 2) / sqrt(1 - par$rho^2) * rnorm(S)
  for (t in 2:T) mu[, t] <- par$c + par$rho * mu[, t - 1] + exp(gg[, t] / 2) * rnorm(S)
  psi <- par$psi
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
  # summaries of the latent and log-variance paths per draw
  list(y = out, sd_mu = apply(mu, 1, sd), max_g = apply(gg, 1, max),
       max_h = apply(h, 1, max), maxabs_mu = apply(abs(mu), 1, max))
}
rtn <- function(n, m, s) { x <- rnorm(n, m, s); while (any(b <- abs(x) >= 1)) x[b] <- rnorm(sum(b), m, s); x }

# Prior specification: a list of functions; each returns the hyperparameter draws.
#   latent = "m13": c ~ N(1.5, 1), rho ~ N(0.5, 0.3) on (-1, 1)
#   latent = "rdp": rho ~ Beta(2.625, 2.625) (mean 0.5, sd 0.2) on (0, 1), mu_bar ~ N(3, 2), c = mu_bar (1 - rho)
#   sv = "unif": phi ~ U(-1, 1), omega ~ half-N(0, 0.5)   (m13)
#   sv = "beta": (phi + 1)/2 ~ Beta(a, b), omega ~ half-N(0, 0.5)
#   sv = "stat": phi ~ U(-1,1) or Beta, stationary sd s = omega / sqrt(1 - phi^2) ~ half-N(0, s_sd)
draw_prior <- function(S, latent = "m13", sv = "unif", ab = c(20, 1.5), s_sd = 1) {
  if (latent == "m13") { c <- rnorm(S, 1.5, 1); rho <- rtn(S, 0.5, 0.3) }
  else { rho <- rbeta(S, 2.625, 2.625); c <- rnorm(S, 3, 2) * (1 - rho) }
  dphi <- function(n) if (sv == "unif") runif(n, -1, 1) else if (sv == "stat_unif") runif(n, -1, 1) else 2 * rbeta(n, ab[1], ab[2]) - 1
  phi_g <- dphi(S); phi <- matrix(dphi(3 * S), S)
  if (sv %in% c("stat", "stat_unif")) {
    omega_g <- abs(rnorm(S, 0, s_sd)) * sqrt(1 - phi_g^2)
    omega <- matrix(abs(rnorm(3 * S, 0, s_sd)), S) * sqrt(1 - phi^2)
  } else {
    omega_g <- abs(rnorm(S, 0, 0.5)); omega <- matrix(abs(rnorm(3 * S, 0, 0.5)), S)
  }
  list(c = c, rho = rho, mu_g = rnorm(S, 2, 1), phi_g = phi_g, omega_g = omega_g,
       mu_h = matrix(rnorm(3 * S, 2, 1), S), phi = phi, omega = omega,
       rho_EI = rtn(S, 0, 0.5), psi = matrix(rtn(3 * S, 0, 0.5), S))
}
# Run in chunks (memory); returns a data frame: one row per draw, parameters + statistics
run_pp <- function(S = 20000, chunk = 2500, seed = 1, ...) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(S / chunk), function(i) {
    p <- draw_prior(chunk, ...)
    sm <- simulate(p, chunk)
    st <- t(apply(sm$y, 1, stats_small))
    data.frame(c = p$c, rho = p$rho, mu_bar = p$c / (1 - p$rho), mu_g = p$mu_g, phi_g = p$phi_g,
               omega_g = p$omega_g, s_g = p$omega_g / sqrt(1 - p$phi_g^2),
               mu_h = p$mu_h, phi = p$phi, omega = p$omega,
               s_h = p$omega / sqrt(1 - p$phi^2), max_s_h = apply(p$omega / sqrt(1 - p$phi^2), 1, max),
               psi = p$psi, rho_EI = p$rho_EI,
               sd_mu = sm$sd_mu, max_g = sm$max_g, max_h = sm$max_h, maxabs_mu = sm$maxabs_mu, st)
  }))
}
qv <- function(x, p = c(0.05, 0.5, 0.95, 0.99)) quantile(x, p, names = FALSE, na.rm = TRUE)
summ_pp <- function(d, label) {
  f <- function(v, nm) data.frame(candidate = label, quantity = nm, q05 = qv(v)[1], q50 = qv(v)[2],
                                  q95 = qv(v)[3], q99 = qv(v)[4])
  rbind(f(d$sd_head, "sd of headline growth"), f(d$mean_head, "mean headline growth"),
        f(d$maxabs_head, "max |headline|"), f(d$ac1_head, "ac1 of headline"),
        f(d$sd_EI, "sd of E - I"), f(d$sd_mu, "sd of latent path"),
        f(d$sd_head_pre2020, "sd of headline, pre-2020"))
}

# ---- candidates ----------------------------------------------------------------------------
cands <- list(
  m13 = list(latent = "m13", sv = "unif"),
  C1 = list(latent = "rdp", sv = "unif"),
  C2 = list(latent = "rdp", sv = "beta", ab = c(20, 1.5)),
  C3 = list(latent = "rdp", sv = "stat_unif", s_sd = 1))
summ <- list(); tgt <- list(); draws <- list()
for (nm in names(cands)) {
  d <- do.call(run_pp, c(list(S = 20000, seed = 1), cands[[nm]]))
  draws[[nm]] <- d
  summ[[nm]] <- summ_pp(d, nm)
  pct <- sapply(c("sd_head", "mean_head", "maxabs_head", "ac1_head", "sd_EI", "sd_head_pre2020"),
                function(v) 100 * mean(d[[v]] <= obs[[v]]))
  T1 <- unname(c(qv(d$sd_head)[1], qv(d$sd_head)[3])); T2 <- qv(d$maxabs_head)[4]
  T3 <- unname(c(qv(d$mean_head)[1], qv(d$mean_head)[3]))
  tgt[[nm]] <- data.frame(candidate = nm,
    T1_sd_q05 = T1[1], T1_sd_q95 = T1[2], T1 = T1[1] >= 1 & T1[2] <= 15,
    T2_maxabs_q99 = T2, T2 = T2 <= 100,
    T3_mean_q05 = T3[1], T3_mean_q95 = T3[2], T3 = T3[1] >= -3 & T3[2] <= 9,
    T4_min_pct = min(pct), T4_max_pct = max(pct), T4 = all(pct >= 2.5 & pct <= 97.5),
    share_sd_gt10 = mean(d$sd_head > 10), share_maxabs_gt100 = mean(d$maxabs_head > 100),
    t(round(pct, 2)))
  tgt[[nm]]$all_pass <- with(tgt[[nm]], T1 & T2 & T3 & T4)
}
tgt <- do.call(rbind, tgt)
cand_order <- c("C1", "C2", "C3")
ok <- cand_order[tgt$all_pass[match(cand_order, tgt$candidate)]]
chosen <- if (length(ok)) ok[1] else {
  x <- tgt[tgt$candidate %in% cand_order & tgt$T4, ]; x$candidate[which.min(x$T2_maxabs_q99)] }
tgt$chosen <- tgt$candidate == chosen
tgt$rule <- ifelse(tgt$chosen, if (length(ok)) "first candidate meeting T1-T4" else
  "fallback: no candidate meets T1-T4; meets T4 with the smallest T2", "")
cat("chosen:", chosen, "\n")
write.csv(do.call(rbind, summ), "paper/data/prior_candidates.csv", row.names = FALSE)
write.csv(tgt, "paper/data/prior_target.csv", row.names = FALSE)

# ---- diagnosis of the m13 prior's tails -------------------------------------------------------
d <- draws$m13
d$max_s <- pmax(d$s_g, d$max_s_h)
d$absphi_max <- pmax(abs(d$phi_g), abs(d$phi.1), abs(d$phi.2), abs(d$phi.3))
tail <- d$maxabs_head > 100
attr_row <- function(sel, lab) {
  x <- d[sel, ]
  data.frame(table = "attribution", set = lab, share = mean(sel),
             rho_gt_0.9 = mean(x$rho > 0.9), absmubar_gt_10 = mean(abs(x$mu_bar) > 10),
             any_absphi_gt_0.95 = mean(x$absphi_max > 0.95), s_g_gt_2 = mean(x$s_g > 2),
             max_s_h_gt_2 = mean(x$max_s_h > 2), any_s_gt_2 = mean(x$max_s > 2),
             latent_maxabs_gt_100 = mean(x$maxabs_mu > 100),
             sd_head_q95 = NA, maxabs_q99 = NA, mean_q95 = NA)
}
cond_row <- function(sel, lab) {
  data.frame(table = "conditional", set = lab, share = mean(sel),
             rho_gt_0.9 = NA, absmubar_gt_10 = NA, any_absphi_gt_0.95 = NA, s_g_gt_2 = NA,
             max_s_h_gt_2 = NA, any_s_gt_2 = NA, latent_maxabs_gt_100 = NA,
             sd_head_q95 = quantile(d$sd_head[sel], 0.95, names = FALSE),
             maxabs_q99 = quantile(d$maxabs_head[sel], 0.99, names = FALSE),
             mean_q95 = quantile(d$mean_head[sel], 0.95, names = FALSE))
}
diag <- rbind(
  attr_row(rep(TRUE, nrow(d)), "all draws"), attr_row(!tail, "max |headline| <= 100"),
  attr_row(tail, "max |headline| > 100"),
  cond_row(rep(TRUE, nrow(d)), "none"), cond_row(d$rho <= 0.9, "rho <= 0.9"),
  cond_row(d$rho > 0.9, "rho > 0.9"), cond_row(d$max_s <= 2, "all s <= 2"),
  cond_row(d$max_s > 2, "some s > 2"), cond_row(d$absphi_max <= 0.95, "all |phi| <= 0.95"),
  cond_row(d$rho <= 0.9 & d$max_s <= 2, "rho <= 0.9 and all s <= 2"))
# the chosen candidate's remaining tail (post-selection information, not used for the choice)
d <- draws[[chosen]]
d$max_s <- pmax(d$s_g, d$max_s_h)
d$absphi_max <- pmax(abs(d$phi_g), abs(d$phi.1), abs(d$phi.2), abs(d$phi.3))
tail <- d$maxabs_head > 100
diag <- rbind(diag, transform(attr_row(tail, "max |headline| > 100"), table = paste("attribution", chosen)),
              transform(attr_row(rep(TRUE, nrow(d)), "all draws"), table = paste("attribution", chosen)))
write.csv(diag, "paper/data/prior_tail_diagnosis.csv", row.names = FALSE)
print(tgt); print(diag, digits = 3)
cat("done\n")
