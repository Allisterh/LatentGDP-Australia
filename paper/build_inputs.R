# paper/build_inputs.R
#
# Extracts everything the working paper (paper/latentgdp_wp.qmd) needs from the local
# fit objects into small CSV files in paper/data/. The .qmd reads only paper/data/*.csv,
# output/tables/*.csv, output/sample_1980/tables/*.csv and docs/*.json|csv, so the paper
# can be rendered from committed files without the (gitignored, ~200 MB each) fits.
#
# Run from the repository root:
#   Rscript paper/build_inputs.R
#
# Inputs (local, gitignored):
#   data/processed/gdp_growth.rds                    ABS 5206.0 Table 24, 2026 Q2 vintage
#   output/fit_production.rds                        production m13, 1980 Q1 - 2026 Q2, seed 1
#   output/sample_1980/fit_m{1..13}.rds              1980 Q1 research fits
#   output/sample_1980/fit_m13_seed2.rds             1980 Q1 m13, seed 2
#   output/fit_m{1..13}.rds                          full-sample (1959 Q4) research fits
#   output/fit_m13_seed2.rds, output/fit_m13_omegag1.rds
#   output/sample_1980/fit_m13_omegag1.rds           1980 Q1 m13, omega_g ~ half-N(0, 1)
#   output/sample_1980/fit_m13_tight{,_seed2}.rds    1980 Q1 tight-priors variant of m13
#                                                    (stan/m13_tight_priors.stan; seeds 1, 2)
#   output/sample_rdp/fit_m13_tight.rds              the same variant on 1980 Q1 - 2013 Q2
#   output/sample_rdp/fit_m{2,4,5,9,12,13}.rds       RDP 2014-12 sample, 1980 Q1 - 2013 Q2
#                                                    (paper/fit_rdp_sample.R)
#   output/sample_rdp/end2019/fit_m{2,13}.rds        1980 Q1 - 2019 Q4 (paper/fit_rdp_sample.R --end)
#   output/sample_1980/logs/test_*.log               unit-test logs (counts of passed checks)
# Nothing outside paper/data/ is written. No model is fitted.
#
# Units: every model works in annualised growth, 400 x dlog (percent). All CSVs below keep
# that unit; the paper converts where it says so.

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
})

stopifnot(file.exists("R/production_config.R"))   # run from the repo root
out_dir <- "paper/data"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
wcsv <- function(x, name) {
  path <- file.path(out_dir, name)
  num <- vapply(x, is.numeric, TRUE)
  x[num] <- lapply(x[num], function(v) signif(v, 7))
  write.csv(x, path, row.names = FALSE)
  cat(sprintf("wrote %-34s %6d rows %8d bytes\n", name, nrow(x), file.size(path)))
}
qs <- function(x, p = c(0.05, 0.25, 0.5, 0.75, 0.95)) quantile(x, p, names = FALSE, na.rm = TRUE)

# ---------------------------------------------------------------------------------------
# 1. Data
# ---------------------------------------------------------------------------------------
g <- readRDS("data/processed/gdp_growth.rds")
g <- as.data.frame(g)
g$date <- as.Date(g$date)
g$headline <- rowMeans(g[, c("g_E", "g_I", "g_P")])
wcsv(data.frame(date = format(g$date), g_E = g$g_E, g_I = g$g_I, g_P = g$g_P,
                headline = g$headline), "growth.csv")
dates_full <- g$date
last_date  <- max(dates_full)
# Every fit except the RDP-sample fits ends in 2026 Q2; the dates of a fit with T quarters
# are the last T quarters up to its end date.
dates_for <- function(T, last = last_date) tail(dates_full[dates_full <= last], T)

# ---------------------------------------------------------------------------------------
# 2. Fits to read
# ---------------------------------------------------------------------------------------
models <- paste0("m", 1:13)
fit_list <- rbind(
  data.frame(key = paste0(models, "_full"), model = models, sample = "full",
             path = file.path("output", paste0("fit_", models, ".rds"))),
  data.frame(key = paste0(models, "_1980"), model = models, sample = "1980",
             path = file.path("output/sample_1980", paste0("fit_", models, ".rds"))),
  data.frame(key = c("m13_prod", "m13_1980_seed2", "m13_full_seed2", "m13_full_omegag1"),
             model = "m13", sample = c("1980", "1980", "full", "full"),
             path = c("output/fit_production.rds", "output/sample_1980/fit_m13_seed2.rds",
                      "output/fit_m13_seed2.rds", "output/fit_m13_omegag1.rds")),
  data.frame(key = "m13_1980_omegag1", model = "m13", sample = "1980",
             path = "output/sample_1980/fit_m13_omegag1.rds"),
  data.frame(key = paste0(c("m2", "m4", "m5", "m9", "m12", "m13"), "_rdp"),
             model = c("m2", "m4", "m5", "m9", "m12", "m13"), sample = "rdp",
             path = file.path("output/sample_rdp", paste0("fit_", c("m2", "m4", "m5", "m9", "m12", "m13"), ".rds")))
)
fit_list <- rbind(fit_list, data.frame(key = c("m2_e2019", "m13_e2019"), model = c("m2", "m13"),
  sample = "e2019", path = file.path("output/sample_rdp/end2019", c("fit_m2.rds", "fit_m13.rds"))))
# Tight-priors variant of m13 (stan/m13_tight_priors.stan; keys m13t_*). model = "m13" so the
# m13 detail (paths of mu, tau_t, sigma_t; decade means) is extracted for it too.
fit_list <- rbind(fit_list, data.frame(key = c("m13t_1980", "m13t_1980_seed2", "m13t_rdp"),
  model = "m13", sample = c("1980", "1980", "rdp"),
  path = c("output/sample_1980/fit_m13_tight.rds", "output/sample_1980/fit_m13_tight_seed2.rds",
           "output/sample_rdp/fit_m13_tight.rds")))
fit_list$last <- ifelse(fit_list$sample == "rdp", "2013-06-01",
                 ifelse(fit_list$sample == "e2019", "2019-12-01", format(last_date)))
fit_list <- fit_list[file.exists(fit_list$path), ]
cat("fits found:", nrow(fit_list), "\n")

hyper_patterns <- c("^c$", "^rho$", "^tau$", "^mu_bar$", "^mu_bar_m$", "^sigma\\[", "^sigma_q\\[",
                    "^sigma_gdp\\[", "^sigma_U$", "^sigma_emp$", "^alpha$", "^beta$",
                    "^alpha_emp$", "^beta_emp$", "^mu_h\\[", "^phi\\[", "^omega\\[",
                    "^rho_EI$", "^psi\\[", "^eps_sd_bar\\[", "^eps_sd\\[", "^mu_g$",
                    "^phi_g$", "^omega_g$", "^kappa\\[", "^Omega\\[", "^s_g$", "^s_h\\[")

params_rows <- list(); covid_rows <- list(); path_rows <- list(); vr_rows <- list()
m13_paths <- list(); decade_rows <- list(); stats_rows <- list(); comb_rows <- list()
med_paths <- list()

win <- function(d, from, to, excl_covid = FALSE) {
  k <- d >= as.Date(from) & d <= as.Date(to)
  if (excl_covid) k <- k & !(d >= as.Date("2020-03-01") & d <= as.Date("2021-03-01"))
  k
}
windows <- list(
  list(name = "1959Q4-2026Q2", from = "1959-12-01", to = "2026-06-01", ex = FALSE),
  list(name = "1959Q4-2026Q2 excl. 2020Q1-2021Q1", from = "1959-12-01", to = "2026-06-01", ex = TRUE),
  list(name = "1978Q4-2026Q2", from = "1978-12-01", to = "2026-06-01", ex = FALSE),
  list(name = "1978Q4-2026Q2 excl. 2020Q1-2021Q1", from = "1978-12-01", to = "2026-06-01", ex = TRUE),
  list(name = "1980Q1-2026Q2", from = "1980-03-01", to = "2026-06-01", ex = FALSE),
  list(name = "1980Q1-2026Q2 excl. 2020Q1-2021Q1", from = "1980-03-01", to = "2026-06-01", ex = TRUE),
  list(name = "1980Q1-2013Q2", from = "1980-03-01", to = "2013-06-01", ex = FALSE),
  list(name = "1980Q1-2019Q4", from = "1980-03-01", to = "2019-12-01", ex = FALSE),
  list(name = "1980Q1-2007Q4", from = "1980-03-01", to = "2007-12-01", ex = FALSE)
)

for (r in seq_len(nrow(fit_list))) {
  fl <- fit_list[r, ]
  t0 <- Sys.time()
  fit <- readRDS(fl$path)
  vars <- fit$metadata()$stan_variables
  muvar <- if ("mu_q" %in% vars) "mu_q" else "mu"
  mu <- unclass(fit$draws(muvar, format = "draws_matrix"))
  T <- ncol(mu); d <- dates_for(T, as.Date(fl$last))
  hl <- g$headline[match(d, g$date)]

  # 2a. hyperparameter summaries
  allv <- fit$metadata()$variables
  hv <- allv[Reduce(`|`, lapply(hyper_patterns, grepl, x = allv))]
  hv <- hv[!grepl("^sigma\\[[0-9]+,", hv)]
  if (length(hv)) {
    s <- summarise_draws(fit$draws(hv), mean, sd,
                         ~quantile(.x, c(0.05, 0.5, 0.95), names = FALSE),
                         rhat, ess_bulk, ess_tail)
    names(s)[4:6] <- c("q05", "q50", "q95")
    params_rows[[fl$key]] <- data.frame(key = fl$key, model = fl$model, sample = fl$sample,
                                        as.data.frame(s))
  }

  # 2b. latent growth in 2020 Q2, and the median path (all fits)
  j <- which(d == as.Date("2020-06-01"))
  if (length(j)) covid_rows[[fl$key]] <- data.frame(key = fl$key, model = fl$model, sample = fl$sample,
                                     T = T, first = format(d[1]),
                                     mu_mean = mean(mu[, j]), mu_q05 = qs(mu[, j])[1],
                                     mu_q50 = qs(mu[, j])[3], mu_q95 = qs(mu[, j])[5])
  if (fl$key %in% c(paste0(models, "_full"), paste0(models, "_1980"))) {
    qq <- apply(mu, 2, qs, p = c(0.05, 0.5, 0.95))
    path_rows[[fl$key]] <- data.frame(model = fl$model, sample = fl$sample, date = format(d),
                                      q05 = round(qq[1, ], 3), q50 = round(qq[2, ], 3),
                                      q95 = round(qq[3, ], 3))
  }

  # 2c. per-draw variance ratios var_t(mu^(s)) / var_t(headline) on several windows,
  #     and the (misleading) ratio of the posterior-mean and posterior-median paths
  for (w in windows) {
    k <- win(d, w$from, w$to, w$ex)
    if (d[1] > as.Date(w$from) || max(d) < as.Date(w$to)) next   # window not covered
    vh <- var(hl[k])
    ratio <- apply(mu[, k, drop = FALSE], 1, var) / vh
    vr_rows[[paste(fl$key, w$name)]] <- data.frame(
      key = fl$key, model = fl$model, sample = fl$sample, window = w$name, n = sum(k),
      ratio_q05 = qs(ratio)[1], ratio_q50 = qs(ratio)[3], ratio_q95 = qs(ratio)[5],
      ratio_meanpath = var(colMeans(mu[, k, drop = FALSE])) / vh,
      ratio_medianpath = var(apply(mu[, k, drop = FALSE], 2, median)) / vh)
  }

  # 2d. m13 detail: paths of mu, tau_t, sigma_t; decade means of tau_t
  if (fl$model == "m13") {
    tt <- unclass(fit$draws("tau_t", format = "draws_matrix"))
    st <- unclass(fit$draws("sigma_t", format = "draws_matrix"))
    qmu <- apply(mu, 2, qs)
    qtt <- apply(tt, 2, qs, p = c(0.05, 0.5, 0.95))
    df <- data.frame(key = fl$key, date = format(d), mu_mean = colMeans(mu),
                     mu_q05 = qmu[1, ], mu_q25 = qmu[2, ], mu_q50 = qmu[3, ],
                     mu_q75 = qmu[4, ], mu_q95 = qmu[5, ],
                     tau_q05 = qtt[1, ], tau_q50 = qtt[2, ], tau_q95 = qtt[3, ])
    for (kk in 1:3) {
      cols <- paste0("sigma_t[", seq_len(T), ",", kk, "]")
      qk <- apply(st[, cols, drop = FALSE], 2, qs, p = c(0.05, 0.5, 0.95))
      nm <- c("E", "I", "P")[kk]
      df[[paste0("sig", nm, "_q05")]] <- qk[1, ]
      df[[paste0("sig", nm, "_q50")]] <- qk[2, ]
      df[[paste0("sig", nm, "_q95")]] <- qk[3, ]
    }
    m13_paths[[fl$key]] <- df
    dec <- c("1960s" = 1960, "1970s" = 1970, "1980s" = 1980, "1990s" = 1990,
             "2000s" = 2000, "2010s" = 2010, "2020s" = 2020)
    yr <- as.integer(format(d, "%Y"))
    for (nm in names(dec)) {
      k <- yr >= dec[[nm]] & yr < dec[[nm]] + 10
      if (!any(k)) next
      dm <- rowMeans(tt[, k, drop = FALSE])
      decade_rows[[paste(fl$key, nm)]] <- data.frame(
        key = fl$key, decade = nm, n = sum(k), q05 = qs(dm)[1], q50 = qs(dm)[3],
        q95 = qs(dm)[5],
        mean_of_median_path = mean(apply(tt[, k, drop = FALSE], 2, median)))
    }
  }

  # 2e. RDP Table 5 / Table 4 analogues on 1980Q1-2013Q2 for the posterior-median path
  #     (quarterly units, /4), for the 1980-sample m2 (closest to RDP Model 1) and m13
  if (fl$sample %in% c("rdp", "e2019") ||
      fl$key %in% c("m13_prod", paste0(c("m2", "m4", "m5", "m9", "m12", "m13"), "_1980"))) {
    k <- win(d, "1980-03-01", "2013-06-01")
    med <- apply(mu, 2, median)
    med_paths[[fl$key]] <- data.frame(key = fl$key, date = format(d), median = med)
    x <- med[k] / 4
    stats_rows[[fl$key]] <- data.frame(series = fl$key, mean = mean(x), sd = sd(x),
                                       rho1 = cor(x[-1], x[-length(x)]))
    E <- g$g_E[match(d[k], g$date)] / 4; I <- g$g_I[match(d[k], g$date)] / 4
    P <- g$g_P[match(d[k], g$date)] / 4
    # min sum (a1 E + a2 I + (1 - a1 - a2) P - M)^2, as in RDP 2014-12 section 4.5.2
    yy <- x - P; X <- cbind(E - P, I - P)
    a <- solve(crossprod(X), crossprod(X, yy))
    comb_rows[[fl$key]] <- data.frame(series = fl$key, w_E = a[1], w_I = a[2],
                                      w_P = 1 - sum(a))
  }

  cat(sprintf("%-18s T=%3d  %5.1f s\n", fl$key, T, as.numeric(Sys.time() - t0, units = "secs")))
  rm(fit, mu); invisible(gc(verbose = FALSE))
}

# Same Table 5 statistics for the published series over 1980Q1-2013Q2 (quarterly units)
k <- g$date >= as.Date("1980-03-01") & g$date <= as.Date("2013-06-01")
for (nm in c("headline", "g_E", "g_I", "g_P")) {
  x <- g[[nm]][k] / 4
  stats_rows[[nm]] <- data.frame(series = nm, mean = mean(x), sd = sd(x),
                                 rho1 = cor(x[-1], x[-length(x)]))
}

params <- do.call(rbind, params_rows); rownames(params) <- NULL
wcsv(params, "params_all.csv")
wcsv(do.call(rbind, covid_rows), "mu_2020q2_by_fit.csv")
wcsv(do.call(rbind, path_rows), "mu_paths_all.csv")
wcsv(do.call(rbind, vr_rows), "variance_ratios.csv")
wcsv(do.call(rbind, m13_paths), "m13_paths.csv")
wcsv(do.call(rbind, decade_rows), "m13_tau_decades.csv")
wcsv(do.call(rbind, stats_rows), "rdp_table5_analogue.csv")
wcsv(do.call(rbind, comb_rows), "rdp_table4_analogue.csv")

# ---------------------------------------------------------------------------------------
# 3. Robustness of m13: path and parameter differences between fits
# ---------------------------------------------------------------------------------------
mp <- do.call(rbind, m13_paths)
cmp <- function(a, b, label) {
  A <- mp[mp$key == a, ]; B <- mp[mp$key == b, ]
  common <- intersect(A$date, B$date)
  A <- A[match(common, A$date), ]; B <- B[match(common, B$date), ]
  dmu <- abs(A$mu_q50 - B$mu_q50); dmean <- abs(A$mu_mean - B$mu_mean)
  rt <- B$tau_q50 / A$tau_q50 - 1
  data.frame(comparison = label, a = a, b = b, n_common = length(common),
             mu_median_mad = mean(dmu), mu_median_max = max(dmu),
             mu_median_max_date = common[which.max(dmu)],
             mu_mean_mad = mean(dmean), mu_mean_max = max(dmean),
             mu_mean_cor = cor(A$mu_mean, B$mu_mean),
             tau_median_rel_mean = mean(abs(rt)), tau_median_rel_max = max(abs(rt)),
             tau_median_rel_max_date = common[which.max(abs(rt))])
}
rob <- rbind(
  cmp("m13_full", "m13_full_omegag1", "prior: omega_g ~ half-N(0, 1) vs half-N(0, 0.5), full sample"),
  cmp("m13_full", "m13_full_seed2", "seed 2 vs seed 1, full sample"),
  cmp("m13_prod", "m13_1980_seed2", "seed 2 vs seed 1, 1980 sample"),
  cmp("m13_full", "m13_prod", "1980 sample vs full sample (common quarters)"),
  cmp("m13_prod", "m13_1980", "production fit vs 1980 research fit"),
  if ("m13_1980_omegag1" %in% names(m13_paths))
    cmp("m13_prod", "m13_1980_omegag1", "prior: omega_g ~ half-N(0, 1) vs half-N(0, 0.5), 1980 sample"),
  if ("m13_rdp" %in% names(m13_paths))
    cmp("m13_prod", "m13_rdp", "RDP sample (to 2013 Q2) vs 1980 sample (to 2026 Q2), common quarters"),
  if ("m13t_1980" %in% names(m13_paths))
    cmp("m13_prod", "m13t_1980", "tight priors vs m13 priors, 1980 sample (seed 1)"),
  if (all(c("m13t_1980", "m13t_1980_seed2") %in% names(m13_paths)))
    cmp("m13t_1980", "m13t_1980_seed2", "tight priors: seed 2 vs seed 1, 1980 sample"),
  if (all(c("m13_rdp", "m13t_rdp") %in% names(m13_paths)))
    cmp("m13_rdp", "m13t_rdp", "tight priors vs m13 priors, RDP sample (to 2013 Q2)")
)
wcsv(rob, "m13_robustness_paths.csv")

# Gate record of the 1980-sample prior-sensitivity fit (written by R/fit_m13.R)
gp <- "output/sample_1980/fit_m13_omegag1_gate.rds"
if (file.exists(gp)) {
  gr <- readRDS(gp)
  wcsv(data.frame(gr$checks, ebfmi_by_chain = paste(sprintf("%.3f", gr$ebfmi), collapse = " "),
                  sampling_secs_total = gr$time$total, wall_secs = gr$wall_secs,
                  omega_g_sd = gr$omega_g_sd), "m13_1980_omegag1_gate.csv")
}

# Gate records of the tight-priors variant (written by R/fit_m13.R --variant=tight), with
# the production fit's seed-1 record for comparison (output/sample_1980/fit_m13_gate.rds
# is the research fit that the production fit reproduces draw for draw)
tg <- list(m13_prod = "output/sample_1980/fit_m13_gate.rds",
           m13_1980_seed2 = "output/sample_1980/fit_m13_seed2_gate.rds",
           m13t_1980 = "output/sample_1980/fit_m13_tight_gate.rds",
           m13t_1980_seed2 = "output/sample_1980/fit_m13_tight_seed2_gate.rds")
tg <- tg[file.exists(unlist(tg))]
if (length(tg)) {
  wcsv(do.call(rbind, lapply(names(tg), function(k) {
    gr <- readRDS(tg[[k]])
    data.frame(key = k, seed = gr$seed, gr$checks, gate_pass = gr$gate_pass,
               ebfmi_by_chain = paste(sprintf("%.3f", gr$ebfmi), collapse = " "),
               treedepth_hits = gr$n_tdmax, sampling_secs_total = gr$time$total,
               wall_secs = gr$wall_secs)
  })), "m13_tight_gate.csv")
}

# ---------------------------------------------------------------------------------------
# 4. Unit-test logs: number of passed expectations per test file
# ---------------------------------------------------------------------------------------
tl <- list.files("output/sample_1980/logs", pattern = "^test_.*\\.log$", full.names = TRUE)
tests <- do.call(rbind, lapply(tl, function(f) {
  x <- readLines(f, warn = FALSE)
  n <- as.integer(sub(".*Test passed with ([0-9]+) success.*", "\\1",
                      grep("Test passed with [0-9]+ success", x, value = TRUE)))
  nf <- length(grep("Failure|Test failed", x))
  data.frame(file = sub("\\.log$", ".R", basename(f)), log = f, passed = sum(n),
             blocks = length(n), failed_blocks = nf)
}))
wcsv(tests, "test_counts.csv")

# ---------------------------------------------------------------------------------------
# 5. Transcribed values (not recomputable from committed tables)
# ---------------------------------------------------------------------------------------
# 5a. RBA RDP 2014-12 (Rees, Lancaster & Finlay), as published. Quarterly units
#     (per cent per quarter; variances in squared quarterly per cent). Source: the paper's
#     Tables 1, 2, 3, 4, 5 and 6 (sample 1980:Q1-2013:Q2).
rdp <- read.csv(text = "
table,model,parameter,stat,value
1,Model 1,mu,median,0.78
1,Model 1,rho,median,0.40
1,Model 1,sigma2_G,median,0.39
1,Model 1,sigma2_E,median,0.44
1,Model 1,sigma2_I,median,0.71
1,Model 1,sigma2_P,median,0.31
1,Model 1,sigma2_G,q05,0.29
1,Model 1,sigma2_G,q95,0.52
1,Model 1,sigma2_E,q05,0.34
1,Model 1,sigma2_E,q95,0.57
1,Model 1,sigma2_I,q05,0.56
1,Model 1,sigma2_I,q95,0.89
1,Model 1,sigma2_P,q05,0.24
1,Model 1,sigma2_P,q95,0.42
1,Model 1,rho,q05,0.24
1,Model 1,rho,q95,0.56
2,Model 2,rho,median,0.46
2,Model 2,sigma2_G,median,0.43
2,Model 2,sigma2_E,median,0.91
2,Model 2,sigma2_I,median,0.98
2,Model 2,sigma2_P,median,0.52
2,Model 2,sigma_EI,median,0.39
2,Model 2,sigma_EP,median,0.36
2,Model 2,sigma_IP,median,0.20
3,Model 3,rho,median,0.62
3,Model 3,sigma2_G,median,0.28
3,Model 3,sigma2_E,median,0.70
3,Model 3,sigma2_I,median,0.87
3,Model 3,sigma2_P,median,0.47
3,Model 3,sigma_EI,median,0.23
3,Model 3,sigma_EP,median,0.22
3,Model 3,sigma_IP,median,0.11
3,Model 3,lambda,median,-0.64
4,Model 1,w_E,value,0.30
4,Model 1,w_I,value,0.21
4,Model 1,w_P,value,0.49
4,Model 2,w_E,value,0.09
4,Model 2,w_I,value,0.26
4,Model 2,w_P,value,0.65
4,Model 3,w_E,value,0.22
4,Model 3,w_I,value,0.24
4,Model 3,w_P,value,0.54
5,GDP(A),sd,value,0.76
5,GDP(E),sd,value,0.93
5,GDP(I),sd,value,1.03
5,GDP(P),sd,value,0.84
5,Model 1,sd,value,0.57
5,Model 2,sd,value,0.56
5,Model 3,sd,value,0.49
5,GDP(A),rho1,value,0.21
5,GDP(E),rho1,value,-0.03
5,GDP(I),rho1,value,-0.19
5,GDP(P),rho1,value,0.31
5,Model 1,rho1,value,0.47
5,Model 2,rho1,value,0.79
5,Model 3,rho1,value,0.68
5,GDP(A),mean,value,0.79
5,GDP(E),mean,value,0.80
5,GDP(I),mean,value,0.78
5,GDP(P),mean,value,0.80
6,GDP(A),mean_abs_revision,value,0.29
6,Model 1,mean_abs_revision,value,0.20
6,Model 2,mean_abs_revision,value,0.19
6,Model 3,mean_abs_revision,value,0.18
", strip.white = TRUE)
wcsv(rdp, "rdp2014_reported.csv")

# 5b. Facts recorded only in the research log (NOTES.md) or Stan file headers, with the
#     section they come from. They describe runs whose outputs are not kept as tables.
facts <- read.csv(text = "
id,value,source
m13_stall_init2_chains,100,NOTES.md 'm13: specification' (warmup-stall harness)
m13_stall_init2_stalled,0,NOTES.md 'm13: specification'
m13_stall_init05_chains,40,NOTES.md 'm13: specification'
m13_stall_init05_stalled,0,NOTES.md 'm13: specification'
m9_stall_bound_default_stalled,3,stan/m9_struct_sv.stan header (of 100 chains; Stan's built-in bound on rho_EI)
m9_stall_bound_init05_stalled,2,stan/m9_struct_sv.stan header (of 100 chains)
m9_stall_tanh_chains,200,stan/m9_struct_sv.stan header
m9_stall_tanh_stalled,0,stan/m9_struct_sv.stan header
m9_old_div_local,26,stan/m9_struct_sv.stan header (mu sampled; T = 267; of 8000)
m9_old_ebfmi_local,0.11,stan/m9_struct_sv.stan header
m9_old_div_ci,44,NOTES.md 'm9 reparameterisation' (CI fit published 2026-10-14; of 8000)
m9_old_ebfmi_ci,0.087,NOTES.md 'm9 reparameterisation'
m3_old_ebfmi,0.03,NOTES.md 'm3: rewrite and identification' (mu sampled)
m3_old_rhat,1.029,NOTES.md 'm3: rewrite and identification'
m3_old_ess,132,NOTES.md 'm3: rewrite and identification'
m4_rw_rhat,1.6,NOTES.md 'Identification issue in m4' (random-walk SV)
m4_ar_musampled_rhat_phi3,1.30,NOTES.md 'Identification issue in m4'
m10_old_chains,16,NOTES.md 'm10: two modes'
m10_old_modeA_chains,9,NOTES.md 'm10: two modes'
m10_old_modeB_chains,3,NOTES.md 'm10: two modes'
m10_old_stranded_chains,4,NOTES.md 'm10: two modes'
m10_old_modeA_rho,0.39,NOTES.md 'm10: two modes'
m10_old_modeA_tau,2.1,NOTES.md 'm10: two modes'
m10_old_modeA_mu2020q2,-1.3,NOTES.md 'm10: two modes'
m10_old_modeA_ratio,0.36,NOTES.md 'm10: two modes'
m10_old_modeA_lp,-2287,NOTES.md 'm10: two modes'
m10_old_modeB_rho,0.12,NOTES.md 'm10: two modes'
m10_old_modeB_tau,3.6,NOTES.md 'm10: two modes'
m10_old_modeB_mu2020q2,-26.9,NOTES.md 'm10: two modes'
m10_old_modeB_ratio,0.86,NOTES.md 'm10: two modes'
m10_old_modeB_lp,-2330,NOTES.md 'm10: two modes'
loo_ploo_min,205,NOTES.md 'Why the earlier PSIS-LOO comparison was invalid' (T = 267 models)
loo_ploo_max,390,NOTES.md 'Why the earlier PSIS-LOO comparison was invalid'
loo_khigh_min,187,NOTES.md 'Why the earlier PSIS-LOO comparison was invalid' (of 267)
loo_khigh_max,242,NOTES.md 'Why the earlier PSIS-LOO comparison was invalid' (of 267)
loo_old_m9_m4,6.3,NOTES.md 'Why the earlier PSIS-LOO comparison was invalid'
loo_old_sv_gain,240,NOTES.md 'Why the earlier PSIS-LOO comparison was invalid'
meanpath_ratio_old,0.65,NOTES.md 'Headline replication result' (posterior-mean path; most models)
pf_brute_force_ref,-37.496,NOTES.md 'm13: specification' (200000 (g; h) paths; se 0.005)
pf_brute_force_boot,-37.488,NOTES.md 'm13: specification'
pf_brute_force_adapt,-37.497,NOTES.md 'm13: specification'
pf_adapt_covid_runsd,0.65,NOTES.md 'm13: LFO' (N = 2000; 12 draws x 4 runs)
pf_boot_covid_runsd,1.38,NOTES.md 'm13: LFO'
pf_adapt_covid_bias,-0.27,NOTES.md 'm13: LFO'
pf_boot_covid_bias,-0.96,NOTES.md 'm13: LFO'
pf_sched_raised_quarters,23,NOTES.md 'm13: LFO' (of 267)
pf_sched_cost,1.31,NOTES.md 'm13: LFO'
lfo_full_m13_hours,3.9,NOTES.md 'm13: LFO' (about 3 h 55 min incl. particle checks)
lfo_1980_hours,7.1,NOTES.md '1980Q1 sample' (about 7 h 05 min)
fits_1980_minutes,42,NOTES.md '1980Q1 sample' (13 fits; 55 with m13 seed 2)
mode_check_1980_first_iter_min,8,NOTES.md 'Posterior: what changes' (m4; m12 warmup)
mode_check_1980_first_iter_max,52,NOTES.md 'Posterior: what changes'
m12_research_fit_parallel_chains,2,NOTES.md 'm12: fits and posterior' (seed-1 research fit)
m12_prod_fit_secs,1260,NOTES.md 'm12: fits and posterior' (production fit; 4 parallel chains; CmdStan sampling total)
m12_stall_chains,100,NOTES.md 'm12: specification' (warmup-stall harness; default inits)
m12_stall_stalled,0,NOTES.md 'm12: specification'
ffbs_wrong_mcse,51.6,NOTES.md 'm12: specification' (m9's backward step reused for m12; worst error in MC SEs)
m10_covid_latent_sigmas,13,NOTES.md 'm10: two modes' (2020 Q2 as a multiple of the constant latent sd; 'roughly 13-sigma')
ci_job_minutes,27,README.md 'Sampler settings and diagnostics' (whole CI job for the 1980-sample production fit on GitHub's runners)
m9_prod_fit_minutes,9,NOTES.md 'm9 reparameterisation' (about 8-9 minutes on 4 local cores)
", strip.white = TRUE)
wcsv(facts, "record_facts.csv")

cat("done\n")
