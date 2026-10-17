# Model comparison on the marginal likelihood (latent states integrated out):
# PSIS-LFO-CV of the one-step-ahead predictive density p(y_t | y_{1:t-1}).
#
#   Rscript R/03_compare.R               # run / resume every model, then summarise
#   Rscript R/03_compare.R m4 m9         # run / resume only these models' LFO
#   Rscript R/03_compare.R summary       # only rebuild tables from checkpoints
#
# Why not PSIS-LOO on log_lik: every spec has a per-quarter latent mu_t (and the SV
# specs a per-quarter h_t), so the saved pointwise log_lik[t] conditions on states
# fitted to the left-out observation. The old LOO objects had p_loo of 208-336 for
# T = 265 and 194-259 of 265 Pareto k > 0.7 for every model: the m9-over-m4 and
# SV-over-non-SV elpd margins they reported were not evidence. Here the latent
# states are integrated out by R/lib/filters.R (Kalman filter for the linear-Gaussian
# specs, Rao-Blackwellised particle filter for the SV specs), so the predictive
# density of each quarter conditions only on the past.
#
# Groups (different data, so compared separately):
#   group 1: 3 measures (E, I, P), full sample  -> m1, m2, m3, m4, m5, m6, m7, m9
#   group 2: 3 measures + hours, post-1978      -> m8, m10
# Group 2 also gets the GDP-only predictive log p(y^GDP_t | y_{1:t-1}) (hours in the
# past, not in the target), which is on the same scale as group 1's density over
# the same evaluation window, and is reported as a supplementary comparison.
# group 3: m11 (mixed frequency: monthly latent, monthly employment). Its LFO periods
# are quarters; employment inside the target quarter is conditioned on (it is
# released before GDP), so its GDP-only score has a larger information set.
#
# Bounds (fixed in advance; everything bounded is logged to output/lfo/lfo_log.txt):
#   evaluation window : 1990Q1 - last quarter (same calendar window for every group)
#   Pareto-k threshold: 0.7 (refit when exceeded)
#   refit cap         : 12 refits per model after the initial fit; beyond that,
#                       high-k steps are used as they are and counted in the table
#   refit sampler     : 4 chains x (1000 warmup + 1000 draws), each model's own
#                       adapt_delta / max_treedepth from R/02_fit_models.R
#   refit accept rule (tightened 2026-10-16, in the spirit of the production gate but
#                       deliberately looser, because refits have 4000 draws, not 8000):
#                         divergences            <= 1% of post-warmup draws
#                         max R-hat (hyperparams) <= 1.05
#                         min E-BFMI (chains)    >= 0.2
#                         min bulk ESS (hyperparams) >= 100
#                       "hyperparameters" = model_defs[[model]]$pars (what the filter
#                       uses). A failing attempt is retried with a new seed up to 3 times
#                       (4 attempts, seeds 10000 + i + (a - 1) * 7919); if none passes, the
#                       best attempt (fewest failed criteria, then lowest R-hat, then
#                       fewest divergences) is kept, marked below_rule = TRUE and counted
#                       in model_comparison.csv (n_fits_below_rule).
#                       Retrying until a seed passes is seed selection: a posterior that
#                       is sampled unreliably at that data end point can pass by chance.
#                       Every attempt, failed or kept, is therefore recorded
#                       (output/tables/lfo_refit_attempts.csv; n_attempts_failed and
#                       n_fits_retried in model_comparison.csv), and `retrycheck` mode
#                       (below) rescores the targets of every retried fit with the first
#                       attempt's draws and with all attempts' draws pooled.
#   refit reuse       : a refit with the same model, Stan file, data end point and seed as
#                       one stored in an earlier checkpoint (REFIT_CACHE) gives the same
#                       draws, so it is taken from that checkpoint instead of re-run; it is
#                       still checked against the rule above (logged as "cached").
#   RBPF (SV models)  : per-fit particle schedule N_t from a pilot pass (R/lib/lfo.R,
#                       pf_schedule): base PF$n_base = 2000 particles, raised to
#                       ceiling(PF$ess_target / f_t) where the pilot's 10% quantile of the
#                       ESS fraction f_t is low (the COVID quarters, the 1970s), capped at
#                       PF$n_max = 100000; adaptive systematic resampling at ESS < N/2.
#                       One filter pass per posterior draw on SV_DRAWS (2000) of the 4000
#                       draws (every 2nd draw of each chain); m10, whose filter is 2.5x
#                       slower, uses 1000 (every 4th draw).
#   Monte Carlo error : the main pass plus PF_RERUNS = 3 independent re-runs of every
#                       filter pass (same draws, same refit and particle schedules, new
#                       seeds). The reported elpd of an SV model is the mean over the 4
#                       passes and its MCSE is sd(pass totals) / sqrt(4); the pointwise
#                       values used for differences are the pass means. The finite-N
#                       downward bias is measured separately (`pfcheck` mode: N = 2000 vs
#                       the schedule vs a high-N reference on the COVID targets, with 3
#                       reference passes for 2020Q2 (`pfref` mode);
#                       output/tables/lfo_pf_covid_check.csv; for group 2 also the
#                       GDP-only target, columns elpd_sub_* / bias_sub_*).
#   SEs of differences: sqrt(n var(d)) treats the quarters as independent. The
#                       pointwise differences are autocorrelated - positively for most
#                       SV-versus-constant-variance pairs, negatively for some (m12 - m9:
#                       lag-1 autocorrelation about -0.4, so there the Newey-West SE is
#                       the smaller) - so a Newey-West SE (Bartlett kernel, NW_LAG = 8
#                       quarters) is reported alongside (se_*_nw, z_nw; lag1_acf in
#                       model_comparison_pairwise.csv).
#
#   Rscript R/03_compare.R m12 --cores=4   # filter workers / parallel chains (default 4)
#   Rscript R/03_compare.R pfcheck m9 m12  # particle-count bias check on the COVID targets
#   Rscript R/03_compare.R pfref m9 m12    #   + 2 more high-N reference passes for 2020Q2
#   Rscript R/03_compare.R pfcheck m9 m12 --precovid  # the same check on 2015Q1-2019Q4
#   Rscript R/03_compare.R retrycheck m9 m12  # rescore the targets of retried refits (seed selection)
#   Rscript R/03_compare.R m9 --no-mcse    # main pass only; then, in parallel processes,
#   Rscript R/03_compare.R rerun m9 1      #   re-runs 1..PF_RERUNS (each stored in order)

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(loo)
  library(dplyr)
})
source("R/lib/filters.R")
source("R/lib/lfo.R")
source("R/lib/filter_m12.R"); register_m12()        # m12: m9 + AR(1) measurement errors
source("R/lib/filter_m10.R"); register_m10_filter() # re-specified m10 (latent SV)
model_defs$m12$pars <- c(M12_PARS, "z_EI")          # z_EI: exact 1 +- rho_EI in the filter
# Compile every filter in the parent: forked filter workers that compile on first use
# race on Rcpp's cache files ("load(file = index_file): empty input file").
load_filters(); load_filter_m12(); load_filter_m10()

WINDOW_START   <- as.Date("1990-03-01")
K_THRES        <- 0.7
MAX_REFITS     <- 12L
N_PARTICLES    <- 2000L     # base particle count (= PF$n_base)
PF <- list(n_base = N_PARTICLES, ess_target = 500, n_max = 100000L,
           pilot_draws = 200L, pilot_n = 2000L, q = 0.1)
PF_RERUNS      <- 3L        # independent re-runs of every filter pass (MCSE)
SV_DRAWS       <- 2000L     # SV models: filter every 2nd draw of each chain (cost)
SV_DRAWS_M     <- c(m10 = 1000L)
REFIT_WARMUP   <- 1000L
REFIT_SAMPLING <- 1000L
NW_LAG         <- 8L        # Newey-West lag (quarters) for autocorrelation-robust SEs
args_all <- commandArgs(trailingOnly = TRUE)
cores_arg <- grep("^--cores=", args_all, value = TRUE)
CORES <- if (length(cores_arg)) as.integer(sub("^--cores=", "", cores_arg)) else 4L
NO_MCSE <- "--no-mcse" %in% args_all    # main LFO pass only; re-runs via `rerun <model> <r>`
args_in <- setdiff(args_all, c(cores_arg, "--no-mcse"))
# Earlier checkpoints whose refits may be reused (same Stan file; see the header).
REFIT_CACHE <- c(m4 = "output/lfo/superseded/lfo_m4_N2000_oldrule.rds",
                 m9 = "output/lfo/superseded/lfo_m9_N2000_oldrule.rds")

dir.create("output/lfo", showWarnings = FALSE, recursive = TRUE)
dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)
log_path <- "output/lfo/lfo_log.txt"
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", sprintf(...))
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = log_path, append = TRUE)
}
log_line <- function(x) log_msg("%s", x)

# ---- data ----------------------------------------------------------------------
gdp   <- readRDS("data/processed/gdp_growth.rds")
gdp4  <- readRDS("data/processed/gdp_growth_with_hours.rds")
y3    <- as.matrix(gdp[, c("g_E", "g_I", "g_P")])
y4    <- as.matrix(gdp4[, c("g_E", "g_I", "g_P", "g_U")])
covid_dates <- as.Date(c("2020-06-01", "2020-09-01", "2020-12-01", "2021-03-01"))
covid_idx <- which(gdp$date %in% covid_dates)
stopifnot(length(covid_idx) == 4)

L1 <- which(gdp$date == WINDOW_START) - 1L     # train on y[1:L], first target = 1990Q1
L2 <- which(gdp4$date == WINDOW_START) - 1L
stopifnot(length(L1) == 1, length(L2) == 1)

data_g1 <- function(i) list(T = i, K = 3, y = y3[seq_len(i), , drop = FALSE])
data_g2 <- function(i) list(T = i, y_gdp = y4[seq_len(i), 1:3, drop = FALSE],
                            y_U = y4[seq_len(i), 4])
data_m7 <- function(i) {
  ci <- covid_idx[covid_idx <= i]
  c(data_g1(i), list(n_covid = length(ci), covid_idx = as.array(as.integer(ci))))
}
# m11: monthly employment + quarterly GDP (same alignment as R/fit_m11.R). LFO periods
# are quarters; the filter runs monthly and ll_post_m11 aggregates months to quarters.
emp <- readRDS("data/processed/employment_growth_monthly.rds")
emp <- emp[emp$date <= max(gdp$date), ]
q_cand <- match(gdp$date, emp$date)
q_keep <- which(!is.na(q_cand) & q_cand >= 3)
q_end <- q_cand[q_keep]
q_dates <- gdp$date[q_keep]
y11 <- matrix(NA_real_, nrow(emp), 4)
y11[, 1] <- emp$g_emp
y11[q_end, 2:4] <- y3[q_keep, ]
L3 <- which(q_dates == WINDOW_START) - 1L
stopifnot(length(L3) == 1, all(diff(q_end) == 3))
data_m11 <- function(i) list(T_m = q_end[i], T_q = i, q_end_m = as.array(q_end[seq_len(i)]),
                             y_q = y3[q_keep[seq_len(i)], , drop = FALSE],
                             y_emp = emp$g_emp[seq_len(q_end[i])])
# Monthly filter output -> quarters. full = log p(employment in the quarter's months,
# GDP of the quarter | past); sub = log p(GDP_q | past, employment through the quarter);
# pre = log p(employment in the quarter's months | past).
ll_post_m11 <- function(ll) {
  emp_ll <- ll$full - ll$sub                     # sub is 0 in months without GDP
  starts <- c(1L, head(q_end, -1) + 1L)
  pre <- sapply(seq_along(q_end), function(j) rowSums(emp_ll[, starts[j]:q_end[j], drop = FALSE]))
  if (!is.matrix(pre)) pre <- matrix(pre, nrow = nrow(ll$full))
  tgt <- ll$sub[, q_end, drop = FALSE]
  list(full = pre + tgt, sub = tgt, pre = pre)
}

sv_set <- list(adapt_delta = 0.99, max_treedepth = 12)
ad99   <- list(adapt_delta = 0.99)

specs <- list(
  m1  = list(group = 1, file = "stan/m1_baseline.stan",        data = data_g1, extras = list()),
  m2  = list(group = 1, file = "stan/m2_ar1.stan",             data = data_g1, extras = list()),
  m3  = list(group = 1, file = "stan/m3_corr.stan",            data = data_g1, extras = ad99),
  m4  = list(group = 1, file = "stan/m4_sv.stan",              data = data_g1, extras = sv_set),
  m5  = list(group = 1, file = "stan/m5_serial.stan",          data = data_g1, extras = ad99),
  m6  = list(group = 1, file = "stan/m6_structural_corr.stan", data = data_g1, extras = ad99),
  m7  = list(group = 1, file = "stan/m7_covid_scale.stan",     data = data_m7, extras = ad99),
  m9  = list(group = 1, file = "stan/m9_struct_sv.stan",       data = data_g1, extras = sv_set),
  m8  = list(group = 2, file = "stan/m8_labour.stan",          data = data_g2, extras = ad99),
  m10 = list(group = 2, file = "stan/m10_labour_sv.stan",      data = data_g2, extras = sv_set),
  m11 = list(group = 3, file = "stan/m11_mixed_freq.stan",     data = data_m11, extras = sv_set),
  m12 = list(group = 1, file = "stan/m12_struct_sv_serial.stan", data = data_g1, extras = sv_set)
)
group_models <- list(`1` = c("m1", "m2", "m3", "m4", "m5", "m6", "m7", "m9", "m12"),
                     `2` = c("m8", "m10"), `3` = "m11")

ctx_fun_for <- function(model) {
  if (model == "m7") function(i) list(covid_idx = covid_idx, n_train = i)
  else function(i) list(n_train = i)
}
pars_for <- function(model, i) {
  p <- model_defs[[model]]$pars
  if (model == "m7" && !any(covid_idx <= i)) p <- setdiff(p, "kappa")
  p
}

csv_dir <- file.path(tempdir(), "lfo_csv")
dir.create(csv_dir, showWarnings = FALSE)

# Refit accept rule (see the header). rule_fails() returns the failed criteria.
RULE <- list(max_div_frac = 0.01, max_rhat = 1.05, min_ebfmi = 0.2, min_ess_bulk = 100)
MAX_ATTEMPTS <- 4L                   # first attempt + up to 3 retries with new seeds
rule_fails <- function(d, n_draws) {
  f <- c(divergences = !(d$divergences <= RULE$max_div_frac * n_draws),
         rhat = !(d$max_rhat <= RULE$max_rhat),
         ebfmi = !(d$min_ebfmi >= RULE$min_ebfmi),
         ess_bulk = !(d$min_ess_bulk >= RULE$min_ess_bulk))
  f[is.na(f)] <- TRUE
  names(f)[f]
}
attempt_seed <- function(i, attempt) 10000L + i + (attempt - 1L) * 7919L
isTRUE_vec <- function(x) !is.na(x) & x

load_cache <- function(model) {
  path <- REFIT_CACHE[model]
  if (is.na(path) || !file.exists(path)) return(NULL)
  old <- readRDS(path)
  # Same Stan program: the source must be older than the cached run.
  if (file.mtime(specs[[model]]$file) >= file.mtime(path)) {
    log_msg("[%s] refit cache %s ignored: %s changed after it was written", model, path,
            specs[[model]]$file)
    return(NULL)
  }
  # old$fits[[k]]$seed is the filter seed; the Stan seed of the kept attempt is in fit_diag.
  lapply(old$fits, function(f) {
    row <- old$fit_diag[old$fit_diag$i_fit == f$i, , drop = FALSE]
    if (nrow(row) != 1) return(NULL)
    list(i = f$i, seed = as.integer(row$seed), dm = f$dm, chain_id = f$chain_id,
         diag = row[, c("divergences", "min_ebfmi", "max_rhat", "min_ess_bulk", "worst_rhat")])
  })
}

# One sampling attempt of a refit on y_{1:i} with the given seed: list(dm, chain_id, diag)
# (draws thinned to the filter's draw count for the SV models), or the error message.
sample_attempt <- function(model, mod, i, seed, attempt) {
  sp <- specs[[model]]
  n_draws <- 4L * REFIT_SAMPLING
  sv_draws <- if (model %in% names(SV_DRAWS_M)) SV_DRAWS_M[[model]] else SV_DRAWS
  fit <- tryCatch(do.call(mod$sample, c(list(
    data = sp$data(i), chains = 4, parallel_chains = min(4L, CORES),
    iter_warmup = REFIT_WARMUP, iter_sampling = REFIT_SAMPLING, seed = seed,
    refresh = 0, show_messages = FALSE, show_exceptions = FALSE,
    output_dir = csv_dir), sp$extras)), error = function(e) e)
  ok <- !inherits(fit, "error") && length(fit$return_codes()) == 4 &&
    all(fit$return_codes() == 0)
  if (!ok) {
    msg <- if (inherits(fit, "error")) conditionMessage(fit) else
      paste("return codes", paste(fit$return_codes(), collapse = ","))
    if (!inherits(fit, "error")) unlink(fit$output_files())
    return(msg)
  }
  dr <- fit$draws(pars_for(model, i), format = "draws_df")
  dm_all <- as_dm(as_draws_matrix(dr))
  keep <- apply(dm_all, 2, sd) > 0                    # e.g. Omega diagonal
  sm <- summarise_draws(subset_draws(as_draws_matrix(dr), variable = colnames(dm_all)[keep]),
                        "rhat", "ess_bulk")
  sd_ <- fit$sampler_diagnostics(format = "draws_df")
  stopifnot(nrow(sd_) == n_draws)
  diag <- data.frame(
    divergences = sum(sd_$divergent__),
    min_ebfmi = min(fit$diagnostic_summary("ebfmi", quiet = TRUE)$ebfmi),
    max_rhat = max(sm$rhat), min_ess_bulk = min(sm$ess_bulk),
    worst_rhat = sm$variable[which.max(sm$rhat)], seed = seed, attempt = attempt,
    cached = FALSE)
  unlink(fit$output_files())
  rm(fit); gc(verbose = FALSE)
  if (model_defs[[model]]$type == "sv") {
    dr <- thin_draws(dr, thin = as.integer(n_draws / sv_draws))
  }
  list(dm = as_dm(as_draws_matrix(dr)), chain_id = dr$.chain, diag = diag)
}

make_refit <- function(model) {
  sp <- specs[[model]]
  mod <- cmdstan_model(sp$file)
  cache <- Filter(Negate(is.null), load_cache(model) %||% list())
  n_draws <- 4L * REFIT_SAMPLING
  function(i) {
    best <- NULL; best_key <- NULL
    attempts <- list()                       # every attempt's diagnostics, failed or kept
    for (attempt in seq_len(MAX_ATTEMPTS)) {
      seed <- attempt_seed(i, attempt)
      hit <- Filter(function(f) f$i == i && f$seed == seed, cache)
      if (length(hit)) {
        h <- hit[[1]]
        diag <- cbind(h$diag, data.frame(seed = seed, attempt = attempt, cached = TRUE))
        cand <- list(dm = h$dm, chain_id = h$chain_id, diag = diag)
        log_msg("[%s] fit on 1:%d attempt %d (seed %d) taken from the refit cache", model, i,
                attempt, seed)
      } else {
        cand <- sample_attempt(model, mod, i, seed, attempt)
        if (is.character(cand)) {
          log_msg("[%s] fit on 1:%d attempt %d failed (%s); retrying with a new seed",
                  model, i, attempt, cand)
          attempts[[length(attempts) + 1]] <- data.frame(seed = seed, attempt = attempt,
                                                         cached = FALSE, fails = "sampler error")
          next
        }
      }
      fails <- rule_fails(cand$diag, n_draws)
      cand$diag$fails <- if (length(fails)) paste(fails, collapse = "+") else ""
      key <- c(length(fails), cand$diag$max_rhat, cand$diag$divergences)
      if (is.null(best) || key[1] < best_key[1] ||
          (key[1] == best_key[1] && (key[2] < best_key[2] ||
                                     (key[2] == best_key[2] && key[3] < best_key[3])))) {
        best <- cand; best_key <- key
      }
      attempts[[length(attempts) + 1]] <- cand$diag
      if (!length(fails)) {
        cand$diag$below_rule <- FALSE
        cand$attempts <- dplyr::bind_rows(attempts)
        cand$attempts$kept <- cand$attempts$attempt == attempt
        return(cand)
      }
      log_msg("[%s] fit on 1:%d attempt %d (seed %d) fails the refit rule [%s]: divergences %d, min E-BFMI %.3f, max R-hat %.3f (%s), min bulk ESS %.0f%s",
              model, i, attempt, seed, paste(fails, collapse = ", "), cand$diag$divergences,
              cand$diag$min_ebfmi, cand$diag$max_rhat, cand$diag$worst_rhat,
              cand$diag$min_ess_bulk, if (attempt < MAX_ATTEMPTS) "; retrying with a new seed" else "")
    }
    if (is.null(best)) stop(sprintf("[%s] fit on 1:%d failed %d times", model, i, MAX_ATTEMPTS))
    log_msg("[%s] fit on 1:%d: no attempt met the refit rule; KEEPING attempt %d (seed %d; fails %s) and counting it as below_rule",
            model, i, best$diag$attempt, best$diag$seed, best$diag$fails)
    best$diag$below_rule <- TRUE
    best$attempts <- dplyr::bind_rows(attempts)
    best$attempts$kept <- best$attempts$attempt == best$diag$attempt
    best
  }
}
`%||%` <- function(a, b) if (is.null(a)) b else a

run_lfo <- function(model) {
  sp <- specs[[model]]
  g <- sp$group
  y <- switch(g, y3, y4, y11)
  dates <- switch(g, gdp$date, gdp4$date, q_dates)
  L <- switch(g, L1, L2, L3)
  log_msg("[%s] PSIS-LFO start: group %d, train on periods 1:%d (to %s), %d one-step predictions %s-%s",
          model, g, L, format(dates[L]), length(dates) - L, format(dates[L + 1]),
          format(dates[length(dates)]))
  psis_lfo(model, y, dates, L, refit_fun = make_refit(model), ctx_fun = ctx_fun_for(model),
           k_thres = K_THRES, max_refits = MAX_REFITS, N = N_PARTICLES, cores = CORES,
           ckpt = sprintf("output/lfo/lfo_%s.rds", model), log_fun = log_line,
           ll_post = if (model == "m11") ll_post_m11 else identity,
           pf = if (model_defs[[model]]$type == "sv") PF else NULL)
}

# Monte Carlo error of the RBPF-based elpd: PF_RERUNS independent re-runs of every filter
# pass with new seeds (same draws, refit schedule and particle schedules). Resumable:
# output/lfo/lfo_<model>_rerun.rds holds the list of finished re-run pointwise tables.
run_mcse <- function(model) {
  path <- sprintf("output/lfo/lfo_%s_rerun.rds", model)
  st <- readRDS(sprintf("output/lfo/lfo_%s.rds", model))
  reruns <- if (file.exists(path)) readRDS(path) else list()
  if (is.data.frame(reruns)) reruns <- list()       # old single-rerun format: start over
  y <- switch(specs[[model]]$group, y3, y4, y11)
  while (length(reruns) < PF_RERUNS) {
    r <- length(reruns) + 1L
    t0 <- Sys.time()
    pw2 <- lfo_rerun_filter(st, y, ctx_fun_for(model), seed_offset = 500000L * r, cores = CORES,
                            ll_post = if (model == "m11") ll_post_m11 else identity)
    log_msg("[%s] filter re-run %d/%d with new seeds: elpd %.2f vs main pass %.2f (%.0f s)", model,
            r, PF_RERUNS, sum(pw2$elpd), sum(st$pointwise$elpd),
            as.numeric(difftime(Sys.time(), t0, units = "secs")))
    reruns[[r]] <- pw2
    saveRDS(reruns, path)
  }
  reruns
}

# Particle-count bias check on the COVID targets (`Rscript R/03_compare.R pfcheck m9 m12`).
# For every fit whose prediction window contains a COVID quarter, the pointwise elpd of
# those targets is recomputed with (a) a constant N = 2000 (2 seeds, the old setting),
# (b) the stored schedule (the main pass and the re-runs, already in the checkpoints),
# (c) a high-N reference: schedule with n_base 8000 (4x), ESS target 4000 (8x), cap 400000.
PF_REF <- list(n_base = 8000L, ess_target = 4000, n_max = 400000L,
               pilot_draws = 200L, pilot_n = 2000L, q = 0.1)
# targets = "precovid": the same check on the 20 quarters 2015Q1-2019Q4 (a calm stretch, to
# see whether the bias is confined to the COVID quarters); output/lfo/pfcheck_<model>_pre.rds,
# output/tables/lfo_pf_precovid_check.csv (`Rscript R/03_compare.R pfcheck m9 m12 --precovid`).
run_pfcheck <- function(model, targets = c("covid", "precovid")) {
  targets <- match.arg(targets)
  st <- readRDS(sprintf("output/lfo/lfo_%s.rds", model))
  y <- switch(specs[[model]]$group, y3, y4, y11)
  reruns <- if (file.exists(sprintf("output/lfo/lfo_%s_rerun.rds", model)))
    readRDS(sprintf("output/lfo/lfo_%s_rerun.rds", model)) else list()
  pw <- st$pointwise
  # by date, so that group 2 (whose period index differs from group 1's) gets the same quarters
  tg <- if (targets == "covid") which(as.Date(pw$date) %in% c(as.Date("2020-03-01"), covid_dates)) # 2020Q1-2021Q1
        else which(as.Date(pw$date) >= as.Date("2015-01-01") & as.Date(pw$date) <= as.Date("2019-12-31"))
  # Resumable: each finished fit is stored in output/lfo/pfcheck_<model>[_pre].rds.
  ck_path <- sprintf("output/lfo/pfcheck_%s%s.rds", model, if (targets == "covid") "" else "_pre")
  out <- if (file.exists(ck_path)) readRDS(ck_path) else list()
  for (f in st$fits) {
    rows <- tg[pw$i_fit[tg] == f$i]
    if (!length(rows)) next
    if (any(vapply(out, function(o) o$i_fit[1] == f$i, TRUE))) {
      log_msg("[%s] pfcheck fit 1:%d taken from %s", model, f$i, ck_path); next
    }
    t_last <- max(pw$t[rows])
    yy <- y[seq_len(t_last), , drop = FALSE]
    ctx <- ctx_fun_for(model)(f$i)
    # rows x 2: the target's elpd and, for group 2, the GDP-only elpd_sub (NA otherwise)
    step_both <- function(ll) matrix(t(vapply(rows, function(r) {
      s <- lfo_step(ll, f$i, pw$i_train[r], f$chain_id); c(s$elpd, s$elpd_sub)
    }, numeric(2))), ncol = 2)
    t0 <- Sys.time()
    b2000 <- lapply(1:2, function(k) step_both(filter_par(model, f$dm, yy, ctx, 2000L, f$seed + 610000L + k, CORES)))
    e2000 <- cbind(b2000[[1]][, 1], b2000[[2]][, 1]); s2000 <- cbind(b2000[[1]][, 2], b2000[[2]][, 2])
    Nref <- pf_schedule(model, f$dm, yy, ctx, PF_REF, f$seed + 620000L, CORES)
    bref <- step_both(filter_par(model, f$dm, yy, ctx, Nref, f$seed + 630000L, CORES))
    eref <- bref[, 1]; sref <- bref[, 2]
    esch <- matrix(c(pw$elpd[rows], unlist(lapply(reruns, function(r2) r2$elpd[rows]))), nrow = length(rows))
    ssch <- matrix(c(pw$elpd_sub[rows], unlist(lapply(reruns, function(r2) r2$elpd_sub[rows]))), nrow = length(rows))
    log_msg("[%s] pfcheck fit 1:%d (%d %s targets, %.0f s): reference %s", model, f$i, length(rows), targets,
            as.numeric(difftime(Sys.time(), t0, units = "secs")), sched_summary(Nref))
    out[[length(out) + 1]] <- data.frame(
      model = model, i_fit = f$i, date = pw$date[rows], t = pw$t[rows],
      N_sched_at_t = if (length(f$Nt) > 1) f$Nt[pw$t[rows]] else f$Nt,
      N_ref_at_t = Nref[pw$t[rows]],
      elpd_N2000_a = e2000[, 1], elpd_N2000_b = e2000[, 2],
      elpd_sched_mean = rowMeans(esch), elpd_sched_sd = if (ncol(esch) > 1) apply(esch, 1, sd) else NA_real_,
      n_sched_passes = ncol(esch), elpd_ref = eref,
      bias_N2000 = rowMeans(e2000) - eref, bias_sched = rowMeans(esch) - eref,
      elpd_sub_N2000_mean = rowMeans(s2000), elpd_sub_sched_mean = rowMeans(ssch), elpd_sub_ref = sref,
      bias_sub_N2000 = rowMeans(s2000) - sref, bias_sub_sched = rowMeans(ssch) - sref)
    saveRDS(out, ck_path)
  }
  do.call(rbind, out)
}

# Extra reference passes for the 2020Q2 target (`Rscript R/03_compare.R pfref m12 m9`), run
# after pfcheck. 2020Q2 carries almost all of the particle-filter spread (per-pass sd
# 0.4-0.7 nats against < 0.03 elsewhere), so one reference pass is as noisy as the bias it
# is meant to measure. PF_REF_EXTRA more reference passes (new pilot and filter seeds)
# are added for that target only; the table then compares N = 2000 and the schedule with
# the mean of all reference passes (elpd_ref_mean, elpd_ref_se = sd / sqrt(passes)).
PF_REF_EXTRA <- 2L
run_pfref <- function(model) {
  st <- readRDS(sprintf("output/lfo/lfo_%s.rds", model))
  y <- switch(specs[[model]]$group, y3, y4, y11)
  pw <- st$pointwise
  r <- which(as.Date(pw$date) == covid_dates[1])         # 2020Q2 (by date: any group)
  f <- st$fits[[which(vapply(st$fits, function(x) as.integer(x$i), 0L) == pw$i_fit[r])]]
  yy <- y[seq_len(pw$t[r]), , drop = FALSE]
  ctx <- ctx_fun_for(model)(f$i)
  ex_path <- sprintf("output/lfo/pfref_%s.rds", model)
  ex <- if (file.exists(ex_path)) readRDS(ex_path) else list()
  while (length(ex) < PF_REF_EXTRA) {
    k <- length(ex) + 1L
    t0 <- Sys.time()
    Nref <- pf_schedule(model, f$dm, yy, ctx, PF_REF, f$seed + 620000L + 1000L * k, CORES)
    ll <- filter_par(model, f$dm, yy, ctx, Nref, f$seed + 630000L + 1000L * k, CORES)
    s_ <- lfo_step(ll, f$i, pw$i_train[r], f$chain_id)
    ex[[k]] <- data.frame(model = model, i_fit = f$i, t = pw$t[r], pass = k + 1L,
                          N_ref_at_t = Nref[pw$t[r]], elpd_ref = s_$elpd, elpd_sub_ref = s_$elpd_sub)
    saveRDS(ex, ex_path)
    log_msg("[%s] pfref 2020Q2 (fit 1:%d) reference pass %d: elpd %.3f at N_t %d (%.0f s)", model, f$i,
            k + 1L, ex[[k]]$elpd_ref, ex[[k]]$N_ref_at_t, as.numeric(difftime(Sys.time(), t0, units = "secs")))
  }
  do.call(rbind, ex)
}
# Adds the reference-mean columns to output/tables/lfo_pf_covid_check.csv.
pf_table_ref_mean <- function(tab, extra) {
  tab$n_ref_passes <- 1L; tab$elpd_ref_mean <- tab$elpd_ref; tab$elpd_ref_se <- NA_real_
  for (m in unique(extra$model)) {
    e <- extra[extra$model == m, ]
    j <- which(tab$model == m & tab$t == e$t[1])
    if (length(j) != 1) next
    v <- c(tab$elpd_ref[j], e$elpd_ref)
    tab$n_ref_passes[j] <- length(v); tab$elpd_ref_mean[j] <- mean(v)
    tab$elpd_ref_se[j] <- sd(v) / sqrt(length(v))
  }
  tab$bias_N2000_vs_refmean <- (tab$elpd_N2000_a + tab$elpd_N2000_b) / 2 - tab$elpd_ref_mean
  tab$bias_sched_vs_refmean <- tab$elpd_sched_mean - tab$elpd_ref_mean
  # GDP-only target (group 2), where the checks recorded it
  if ("elpd_sub_ref" %in% names(tab)) {
    tab$elpd_sub_ref_mean <- tab$elpd_sub_ref
    if ("elpd_sub_ref" %in% names(extra)) for (m in unique(extra$model)) {
      e <- extra[extra$model == m, ]
      j <- which(tab$model == m & tab$t == e$t[1])
      if (length(j) == 1 && !anyNA(e$elpd_sub_ref)) tab$elpd_sub_ref_mean[j] <- mean(c(tab$elpd_sub_ref[j], e$elpd_sub_ref))
    }
    tab$bias_sub_sched_vs_refmean <- tab$elpd_sub_sched_mean - tab$elpd_sub_ref_mean
  }
  tab
}

# Seed-selection check (`Rscript R/03_compare.R retrycheck m9 m12`). For every stored fit
# that was kept only after earlier attempts failed the refit rule, the failed attempts are
# re-created (same seeds, so the same draws; taken from the refit cache where it has them),
# and the targets that fit serves are rescored three ways, each with the fit's stored
# particle schedule and one common filter seed:
#   kept   - the kept attempt's draws (the main pass, re-filtered with that seed);
#   first  - attempt 1's draws only: what the score would be with no retry at all;
#   pooled - all attempts' chains pooled, the same number of draws from each attempt (no
#            selection among seeds).
# The refit schedule (which fit serves which target) is held fixed. Resumable:
# output/lfo/retrycheck_<model>.rds; table output/tables/lfo_retry_sensitivity.csv. The
# effect on the pairwise differences is summarised by the `summary` step
# (output/tables/lfo_retry_sensitivity_summary.csv).
RETRY_FILTER_SEED <- 640000L
run_retrycheck <- function(model) {
  st <- readRDS(sprintf("output/lfo/lfo_%s.rds", model))
  y <- switch(specs[[model]]$group, y3, y4, y11)
  fd <- st$fit_diag
  retried <- fd$i_fit[isTRUE_vec(fd$attempt > 1)]
  ck_path <- sprintf("output/lfo/retrycheck_%s.rds", model)
  out <- if (file.exists(ck_path)) readRDS(ck_path) else list()
  if (!length(retried)) { log_msg("[%s] retrycheck: no retried fits", model); return(NULL) }
  mod <- cmdstan_model(specs[[model]]$file)
  cache <- Filter(Negate(is.null), load_cache(model) %||% list())
  pw <- st$pointwise
  for (i in retried) {
    if (any(vapply(out, function(o) o$rows$i_fit[1] == i, TRUE))) {
      log_msg("[%s] retrycheck fit 1:%d taken from %s", model, i, ck_path); next
    }
    t0 <- Sys.time()
    f <- st$fits[[which(vapply(st$fits, function(x) as.integer(x$i), 0L) == i)]]
    ka <- fd$attempt[fd$i_fit == i]
    att <- vector("list", ka)
    for (a in seq_len(ka - 1L)) {
      seed <- attempt_seed(i, a)
      hit <- Filter(function(h) h$i == i && h$seed == seed, cache)
      att[[a]] <- if (length(hit))
        list(dm = hit[[1]]$dm, chain_id = hit[[1]]$chain_id,
             diag = cbind(hit[[1]]$diag, data.frame(seed = seed, attempt = a, cached = TRUE)))
      else sample_attempt(model, mod, i, seed, a)
      if (is.character(att[[a]])) stop(sprintf("[%s] retrycheck 1:%d attempt %d: %s", model, i, a, att[[a]]))
      log_msg("[%s] retrycheck fit 1:%d attempt %d (seed %d) re-created%s: divergences %d, min E-BFMI %.3f, max R-hat %.3f (%s), min bulk ESS %.0f",
              model, i, a, seed, if (isTRUE(att[[a]]$diag$cached)) " from the refit cache" else "",
              att[[a]]$diag$divergences, att[[a]]$diag$min_ebfmi, att[[a]]$diag$max_rhat,
              att[[a]]$diag$worst_rhat, att[[a]]$diag$min_ess_bulk)
    }
    att[[ka]] <- list(dm = f$dm, chain_id = f$chain_id, diag = fd[fd$i_fit == i, , drop = FALSE])
    # Pooled: the same number of draws from every chain of every attempt (evenly spaced
    # within each chain), chains renumbered 1..4 * ka.
    per_chain <- min(table(f$chain_id)) %/% ka
    pool <- lapply(seq_len(ka), function(a) {
      ci <- att[[a]]$chain_id
      rows <- unlist(lapply(sort(unique(ci)), function(c) {
        r <- which(ci == c); r[unique(round(seq(1, length(r), length.out = per_chain)))]
      }))
      list(dm = att[[a]]$dm[rows, , drop = FALSE], chain_id = (a - 1L) * 4L + ci[rows])
    })
    pooled <- list(dm = do.call(rbind, lapply(pool, `[[`, "dm")),
                   chain_id = unlist(lapply(pool, `[[`, "chain_id")))
    # Convergence across all pooled chains (the attempts' thinned draws).
    arr <- function(a) {
      ci <- att[[a]]$chain_id; n <- min(table(ci))
      lapply(sort(unique(ci)), function(c) att[[a]]$dm[which(ci == c)[seq_len(n)], , drop = FALSE])
    }
    chains <- unlist(lapply(seq_len(ka), arr), recursive = FALSE)
    A <- array(NA_real_, c(nrow(chains[[1]]), length(chains), ncol(f$dm)),
               dimnames = list(NULL, NULL, colnames(f$dm)))
    for (c in seq_along(chains)) A[, c, ] <- chains[[c]]
    A <- A[, , apply(f$dm, 2, sd) > 0, drop = FALSE]
    smp <- summarise_draws(as_draws_array(A), "rhat", "ess_bulk")
    rows <- which(pw$i_fit == i)
    ctx <- ctx_fun_for(model)(i)
    score <- function(dm, chain_id) {
      ll <- filter_par(model, dm, y, ctx, f$Nt, f$seed + RETRY_FILTER_SEED, CORES)
      t(vapply(rows, function(r) {
        s <- lfo_step(ll, i, pw$i_train[r], chain_id); c(s$elpd, s$elpd_sub, s$k)
      }, numeric(3)))
    }
    e_kept <- score(f$dm, f$chain_id)
    e_first <- score(att[[1]]$dm, att[[1]]$chain_id)
    e_pool <- score(pooled$dm, pooled$chain_id)
    res <- data.frame(model = model, i_fit = i, date_fit = fd$date_fit[fd$i_fit == i],
                      kept_attempt = ka, t = pw$t[rows], date = pw$date[rows],
                      elpd_main = pw$elpd[rows], elpd_kept = e_kept[, 1], elpd_first = e_first[, 1],
                      elpd_pooled = e_pool[, 1], elpd_sub_kept = e_kept[, 2],
                      elpd_sub_first = e_first[, 2], elpd_sub_pooled = e_pool[, 2],
                      k_kept = e_kept[, 3], k_first = e_first[, 3], k_pooled = e_pool[, 3])
    fitrow <- data.frame(model = model, i_fit = i, kept_attempt = ka, n_draws_pooled = nrow(pooled$dm),
                         pooled_chains = length(chains), pooled_max_rhat = max(smp$rhat),
                         pooled_worst_rhat = smp$variable[which.max(smp$rhat)],
                         pooled_min_ess_bulk = min(smp$ess_bulk))
    out[[length(out) + 1]] <- list(rows = res, fit = fitrow,
                                   attempts = bind_rows(lapply(att, `[[`, "diag")))
    saveRDS(out, ck_path)
    log_msg("[%s] retrycheck fit 1:%d (%d targets, %.0f s): pooled %d chains max R-hat %.3f (%s), min bulk ESS %.0f; sum elpd kept %.3f, first attempt %.3f, pooled %.3f (main pass %.3f)",
            model, i, length(rows), as.numeric(difftime(Sys.time(), t0, units = "secs")),
            length(chains), fitrow$pooled_max_rhat, fitrow$pooled_worst_rhat, fitrow$pooled_min_ess_bulk,
            sum(res$elpd_kept), sum(res$elpd_first), sum(res$elpd_pooled), sum(res$elpd_main))
  }
  out
}

# ---- run ----------------------------------------------------------------------------
# `Rscript R/03_compare.R rerun m9 3 --cores=3`: compute filter re-run number 3 of a
# finished LFO run in a separate process (same seeds as run_mcse would use), so that
# re-runs can proceed in parallel. It is stored at position 3 of lfo_<model>_rerun.rds once
# re-runs 1-2 are there (it waits for them).
if (length(args_in) && args_in[1] == "rerun") {
  m <- args_in[2]; r <- as.integer(args_in[3])
  stopifnot(m %in% names(specs), !is.na(r), r >= 1)
  st <- readRDS(sprintf("output/lfo/lfo_%s.rds", m))
  stopifnot(isTRUE(st$done))
  y <- switch(specs[[m]]$group, y3, y4, y11)
  t0 <- Sys.time()
  pw2 <- lfo_rerun_filter(st, y, ctx_fun_for(m), seed_offset = 500000L * r, cores = CORES,
                          ll_post = if (m == "m11") ll_post_m11 else identity)
  path <- sprintf("output/lfo/lfo_%s_rerun.rds", m)
  repeat {
    rr <- if (file.exists(path)) tryCatch(readRDS(path), error = function(e) NULL) else list()
    if (!is.null(rr) && length(rr) >= r - 1) break
    Sys.sleep(20)
  }
  if (length(rr) >= r) stop(sprintf("[%s] re-run %d is already stored", m, r))
  rr[[r]] <- pw2
  tmp <- paste0(path, ".tmp", Sys.getpid())
  saveRDS(rr, tmp); file.rename(tmp, path)          # atomic: other re-run processes poll it
  log_msg("[%s] filter re-run %d/%d with new seeds (separate process): elpd %.2f vs main pass %.2f (%.0f s)",
          m, r, PF_RERUNS, sum(pw2$elpd), sum(st$pointwise$elpd),
          as.numeric(difftime(Sys.time(), t0, units = "secs")))
  quit(save = "no")
}
if (length(args_in) && args_in[1] == "pfref") {
  ms <- args_in[-1]
  stopifnot(length(ms) > 0, all(ms %in% names(specs)))
  log_msg("=== extra reference passes for the 2020Q2 target: %s (%d cores)", paste(ms, collapse = ", "), CORES)
  extra <- do.call(rbind, lapply(ms, run_pfref))
  path <- "output/tables/lfo_pf_covid_check.csv"
  if (!file.exists(path)) { log_msg("pfref: %s not found (run pfcheck first); extra passes kept in output/lfo/pfref_*.rds", path); quit(save = "no") }
  tab <- read.csv(path, stringsAsFactors = FALSE)
  tab <- tab[, setdiff(names(tab), c("n_ref_passes", "elpd_ref_mean", "elpd_ref_se",
                                     "bias_N2000_vs_refmean", "bias_sched_vs_refmean",
                                     "elpd_sub_ref_mean", "bias_sub_sched_vs_refmean"))]
  all_extra <- dplyr::bind_rows(lapply(unique(tab$model), function(m) {
    pth <- sprintf("output/lfo/pfref_%s.rds", m)
    if (file.exists(pth)) do.call(rbind, readRDS(pth)) else NULL
  }))
  tab <- if (is.null(all_extra)) tab else pf_table_ref_mean(tab, all_extra)
  write.csv(tab, path, row.names = FALSE)
  for (m in ms) {
    tm <- tab[tab$model == m, ]
    q2 <- as.Date(tm$date) == covid_dates[1]
    log_msg("[%s] pfref: 2020Q2 reference mean %.3f (se %.3f, %d passes); schedule mean %.3f, N=2000 mean %.3f; sum over COVID targets vs reference mean: bias N=2000 %.3f, bias schedule %.3f",
            m, tm$elpd_ref_mean[q2], tm$elpd_ref_se[q2],
            tm$n_ref_passes[q2], tm$elpd_sched_mean[q2],
            mean(c(tm$elpd_N2000_a[q2], tm$elpd_N2000_b[q2])),
            sum(tm$bias_N2000_vs_refmean), sum(tm$bias_sched_vs_refmean))
  }
  quit(save = "no")
}
if (length(args_in) && args_in[1] == "retrycheck") {
  ms <- args_in[-1]
  stopifnot(length(ms) > 0, all(ms %in% names(specs)))
  log_msg("=== seed-selection check of retried refits: %s (%d cores)", paste(ms, collapse = ", "), CORES)
  for (m in ms) run_retrycheck(m)
  tab <- bind_rows(lapply(names(specs), function(m) {
    pth <- sprintf("output/lfo/retrycheck_%s.rds", m)
    if (file.exists(pth)) bind_rows(lapply(readRDS(pth), function(o) cbind(o$rows, o$fit[, -(1:3)])))
  }))
  write.csv(tab, "output/tables/lfo_retry_sensitivity.csv", row.names = FALSE)
  quit(save = "no")
}
if (length(args_in) && args_in[1] == "pfcheck") {
  pre <- "--precovid" %in% args_in
  ms <- setdiff(args_in[-1], "--precovid")
  stopifnot(length(ms) > 0, all(ms %in% names(specs)))
  log_msg("=== particle-count check on the %s targets: %s (%d cores)", if (pre) "2015-2019" else "COVID",
          paste(ms, collapse = ", "), CORES)
  path <- if (pre) "output/tables/lfo_pf_precovid_check.csv" else "output/tables/lfo_pf_covid_check.csv"
  for (m in ms) {
    res <- run_pfcheck(m, if (pre) "precovid" else "covid")
    old <- if (file.exists(path)) read.csv(path, stringsAsFactors = FALSE) else NULL
    if (!is.null(old)) old <- old[old$model != m, , drop = FALSE]
    write.csv(dplyr::bind_rows(old, res), path, row.names = FALSE)
    log_msg("[%s] pfcheck: sum over %s targets: bias N=2000 %.3f, bias schedule %.3f (vs reference)%s",
            m, if (pre) "2015-2019" else "COVID", sum(res$bias_N2000), sum(res$bias_sched),
            if (!is.null(res$bias_sub_sched) && !anyNA(res$bias_sub_sched))
              sprintf("; GDP-only target: bias N=2000 %.3f, bias schedule %.3f",
                      sum(res$bias_sub_N2000), sum(res$bias_sub_sched)) else "")
  }
  quit(save = "no")
}
summary_only <- identical(args_in, "summary")
to_run <- if (summary_only) character(0) else if (length(args_in)) args_in else names(specs)
bad <- setdiff(to_run, names(specs))
if (length(bad)) stop("unknown model(s): ", paste(bad, collapse = ", "))
if (length(to_run)) {
  log_msg("=== PSIS-LFO run (%s; %d cores): window from %s, k threshold %.1f, refit cap %d, refits %d+%d x 4 chains, rule div<=%.0f%% R-hat<=%.2f E-BFMI>=%.1f bulk ESS>=%d (%d attempts), RBPF N_t %d-%d (ESS target %.0f), %d re-runs",
          paste(to_run, collapse = ", "), CORES, format(WINDOW_START), K_THRES, MAX_REFITS,
          REFIT_WARMUP, REFIT_SAMPLING, 100 * RULE$max_div_frac, RULE$max_rhat, RULE$min_ebfmi,
          RULE$min_ess_bulk, MAX_ATTEMPTS, PF$n_base, PF$n_max, PF$ess_target, PF_RERUNS)
}
for (m in to_run) {
  run_lfo(m)
  if (model_defs[[m]]$type == "sv" && !NO_MCSE) run_mcse(m)
}

# ---- summarise ------------------------------------------------------------------------
all_models <- names(specs)
have <- all_models[file.exists(sprintf("output/lfo/lfo_%s.rds", all_models))]
st_all <- lapply(setNames(have, have), function(m) readRDS(sprintf("output/lfo/lfo_%s.rds", m)))
done <- names(Filter(function(s) isTRUE(s$done), st_all))
if (!length(done)) { log_msg("no completed LFO runs yet"); quit(save = "no") }

# SV models: the pointwise elpd is the mean over the main filter pass and its re-runs
# (elpd_main keeps the main pass, which made the refit decisions); MCSE of the total =
# sd(pass totals) / sqrt(n_passes).
n_passes <- setNames(rep(1L, length(done)), done)
mcse_rows <- list()
for (m in done) {
  p0 <- st_all[[m]]$pointwise
  p0$elpd_main <- p0$elpd; p0$elpd_sub_main <- p0$elpd_sub
  path <- sprintf("output/lfo/lfo_%s_rerun.rds", m)
  if (model_defs[[m]]$type == "sv" && file.exists(path)) {
    rr <- readRDS(path)
    if (is.data.frame(rr)) rr <- list(rr)
    E <- cbind(p0$elpd, sapply(rr, `[[`, "elpd"))
    Es <- cbind(p0$elpd_sub, sapply(rr, `[[`, "elpd_sub"))
    p0$elpd <- rowMeans(E); p0$elpd_sub <- rowMeans(Es)
    p0$elpd_mc_sd <- apply(E, 1, sd)
    n_passes[m] <- ncol(E)
    tot <- colSums(E); tot_s <- colSums(Es)
    mcse_rows[[m]] <- tibble(model = m, n_pf_passes = ncol(E),
      pass_totals = paste(sprintf("%.2f", tot), collapse = " "),
      mcse_single_pass = sd(tot), mcse_elpd = sd(tot) / sqrt(ncol(E)),
      mcse_elpd_gdp = if (anyNA(tot_s)) NA_real_ else sd(tot_s) / sqrt(ncol(E)))
  } else {
    p0$elpd_mc_sd <- if (model_defs[[m]]$type == "sv") NA_real_ else 0
  }
  st_all[[m]]$pointwise <- p0
}
pw <- bind_rows(lapply(st_all[done], `[[`, "pointwise"))
mcse_tab <- bind_rows(mcse_rows)
if (!nrow(mcse_tab)) mcse_tab <- tibble(model = character(), n_pf_passes = integer(),
                                         pass_totals = character(), mcse_single_pass = numeric(),
                                         mcse_elpd = numeric(), mcse_elpd_gdp = numeric())
mcse_of <- function(m, gdp = FALSE) {
  if (model_defs[[m]]$type != "sv") return(0)
  v <- if (gdp && specs[[m]]$group != 1) mcse_tab$mcse_elpd_gdp[mcse_tab$model == m] else
    mcse_tab$mcse_elpd[mcse_tab$model == m]
  if (length(v)) v else NA_real_
}

# Newey-West (Bartlett kernel) SE of a sum of n pointwise values: sqrt(n * long-run
# variance). Reported alongside the iid sqrt(n var(x)), which ignores autocorrelation.
se_nw <- function(x, L = NW_LAG) {
  n <- length(x); e <- x - mean(x)
  lrv <- sum(e * e) / n
  for (l in seq_len(min(L, n - 1)))
    lrv <- lrv + 2 * (1 - l / (L + 1)) * sum(e[-(1:l)] * e[1:(n - l)]) / n
  sqrt(n * lrv)
}
# Every stored fit, checked against the current refit rule (also the fits of runs made
# under the earlier rule, which tested only divergences and R-hat).
fit_diag_all <- bind_rows(lapply(st_all[done], `[[`, "fit_diag"))
if (!"cached" %in% names(fit_diag_all)) fit_diag_all$cached <- NA
if (!"attempt" %in% names(fit_diag_all)) fit_diag_all$attempt <- NA_integer_
fit_diag_all$fails_rule <- vapply(seq_len(nrow(fit_diag_all)), function(r)
  paste(rule_fails(fit_diag_all[r, ], 4L * REFIT_SAMPLING), collapse = "+"), "")
fit_diag_all <- fit_diag_all |>
  mutate(meets_rule = fails_rule == "",
         weak_mixing = min_ebfmi < 0.3 | min_ess_bulk < 400)
# Every refit attempt, failed ones included (output/tables/lfo_refit_attempts.csv).
# Checkpoints written since attempts are stored carry st$attempts. Older checkpoints (all
# of the 2026-10-17 runs) store only the kept attempt; their failed attempts are recovered
# from the "fails the refit rule" lines of output/lfo/lfo_log.txt, keeping only lines
# logged after the model's Stan file was last changed, for a stored fit, an attempt before
# the kept one and that attempt's deterministic seed. n_unrecorded counts failed attempts
# implied by the kept attempt number but found in neither place.
log_attempts <- function(m) {
  if (!file.exists(log_path)) return(NULL)
  L <- readLines(log_path)
  pat <- sprintf(paste0("^(\\S+ \\S+)  \\[%s\\] fit on 1:([0-9]+) attempt ([0-9]+) \\(seed ([0-9]+)\\) ",
                        "fails the refit rule \\[([^]]*)\\]: divergences ([0-9]+), min E-BFMI ([0-9.]+), ",
                        "max R-hat ([0-9.]+) \\((.*)\\), min bulk ESS ([0-9.]+)"), m)
  hit <- Filter(length, regmatches(L, regexec(pat, L)))
  if (!length(hit)) return(NULL)
  d <- do.call(rbind, lapply(hit, function(h) data.frame(
    logged_at = h[2], i_fit = as.integer(h[3]), attempt = as.integer(h[4]), seed = as.integer(h[5]),
    fails = gsub(", ", "+", h[6]), divergences = as.numeric(h[7]), min_ebfmi = as.numeric(h[8]),
    max_rhat = as.numeric(h[9]), worst_rhat = h[10], min_ess_bulk = as.numeric(h[11]))))
  # Keep only attempts whose seed matches the deterministic retry scheme, and the
  # latest log line per (refit, attempt). No time cutoff: the caller keeps only
  # attempts below the checkpoint's kept attempt, and the current run always
  # logged those itself, so its lines (the latest) replace any stale ones from an
  # older model version. (A file.mtime() cutoff used to drop every line once a
  # .stan file was touched, and then crashed on a 0-row frame.)
  d <- d[d$seed == attempt_seed(d$i_fit, d$attempt), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  d[!duplicated(d[, c("i_fit", "attempt")], fromLast = TRUE), , drop = FALSE]
}
attempt_cols <- c("model", "i_fit", "date_fit", "attempt", "seed", "kept", "cached", "divergences",
                  "min_ebfmi", "max_rhat", "worst_rhat", "min_ess_bulk", "fails", "source")
attempts_all <- bind_rows(lapply(done, function(m) {
  fd <- fit_diag_all[fit_diag_all$model == m, , drop = FALSE]
  a <- st_all[[m]]$attempts
  if (!is.null(a) && nrow(a)) {
    a$source <- "checkpoint"
  } else {
    kept <- fd |> mutate(kept = TRUE, fails = fails_rule, source = "checkpoint (kept attempt)")
    lg <- log_attempts(m)
    if (!is.null(lg)) {
      lg <- lg[lg$i_fit %in% fd$i_fit, , drop = FALSE]
      ka <- fd$attempt[match(lg$i_fit, fd$i_fit)]
      lg <- lg[!is.na(ka) & lg$attempt < ka, , drop = FALSE]
      if (!nrow(lg)) lg <- NULL
    }
    if (!is.null(lg)) {
      lg$model <- m; lg$kept <- FALSE; lg$cached <- NA; lg$source <- "lfo_log.txt"
      lg$date_fit <- fd$date_fit[match(lg$i_fit, fd$i_fit)]
    }
    a <- bind_rows(kept, lg)
  }
  a <- a[, intersect(attempt_cols, names(a)), drop = FALSE]
  a[order(a$i_fit, a$attempt), , drop = FALSE]
}))
n_unrecorded <- function(m) {
  fd <- fit_diag_all[fit_diag_all$model == m & isTRUE_vec(fit_diag_all$attempt > 1), , drop = FALSE]
  a <- attempts_all[attempts_all$model == m & !attempts_all$kept, , drop = FALSE]
  sum(vapply(seq_len(nrow(fd)), function(r) (fd$attempt[r] - 1L) - sum(a$i_fit == fd$i_fit[r]), 0))
}
write.csv(attempts_all, "output/tables/lfo_refit_attempts.csv", row.names = FALSE)

for (m in done) {
  bad <- fit_diag_all[fit_diag_all$model == m & !fit_diag_all$meets_rule, ]
  if (nrow(bad)) log_msg("summary: %s has %d fit(s) below the refit rule: %s", m, nrow(bad),
                         paste(sprintf("1:%d (%s) [%s]", bad$i_fit, bad$date_fit, bad$fails_rule), collapse = "; "))
}

summ_group <- function(g) {
  ms <- intersect(group_models[[g]], done)
  if (!length(ms)) return(NULL)
  P <- sapply(ms, function(m) st_all[[m]]$pointwise$elpd)
  if (!is.matrix(P)) P <- matrix(P, ncol = length(ms), dimnames = list(NULL, ms))
  ts <- lapply(ms, function(m) st_all[[m]]$pointwise$t)
  stopifnot(all(vapply(ts, identical, TRUE, ts[[1]])))
  n <- nrow(P)
  tot <- colSums(P)
  best <- names(which.max(tot))
  tibble(
    group = as.integer(g), model = ms,
    sample = switch(g, `1` = "E, I, P; full sample", `2` = "E, I, P, hours; post-1978",
                    `3` = "E, I, P quarterly + monthly employment; post-1978"),
    window = sprintf("%s to %s", st_all[[ms[1]]]$pointwise$date[1], st_all[[ms[1]]]$pointwise$date[n]),
    n_pred = n,
    elpd_lfo = tot,
    se_elpd = apply(P, 2, function(x) sqrt(n * var(x))),
    elpd_diff_vs_best = tot - tot[best],
    se_diff_vs_best = apply(P, 2, function(x) sqrt(n * var(x - P[, best]))),
    se_diff_vs_best_nw = apply(P, 2, function(x) se_nw(x - P[, best])),
    best_in_group = best,
    n_fits = vapply(ms, function(m) length(st_all[[m]]$fits), 0L),
    n_refits = vapply(ms, function(m) st_all[[m]]$n_refits, 0L),
    refit_cap = vapply(ms, function(m) st_all[[m]]$max_refits, 0L),
    n_k_high_unrefit = vapply(ms, function(m) sum(st_all[[m]]$pointwise$k_high_unrefit), 0L),
    n_fits_below_rule = vapply(ms, function(m) sum(!fit_diag_all$meets_rule[fit_diag_all$model == m]), 0L),
    n_fits_retried = vapply(ms, function(m) sum(isTRUE_vec(fit_diag_all$attempt[fit_diag_all$model == m] > 1)), 0L),
    n_attempts_failed = vapply(ms, function(m) sum(!attempts_all$kept[attempts_all$model == m]), 0L),
    n_attempts_unrecorded = vapply(ms, function(m) as.integer(n_unrecorded(m)), 0L),
    n_fits_cached = vapply(ms, function(m) sum(isTRUE_vec(fit_diag_all$cached[fit_diag_all$model == m])), 0L),
    n_fits_weak_mixing = vapply(ms, function(m) sum(fit_diag_all$weak_mixing[fit_diag_all$model == m]), 0L),
    max_k_used = vapply(ms, function(m) {
      k <- st_all[[m]]$pointwise$pareto_k[!st_all[[m]]$pointwise$refit]
      if (all(is.na(k))) NA_real_ else max(k, na.rm = TRUE)
    }, 0),
    filter = vapply(ms, function(m) {
      if (model_defs[[m]]$type != "sv") return("Kalman (exact)")
      Ns <- unlist(lapply(st_all[[m]]$fits, `[[`, "Nt"))
      sprintf("RBPF N_t %d-%d (mean %.0f), %d draws, %d passes", min(Ns), max(Ns), mean(Ns),
              nrow(st_all[[m]]$fits[[1]]$dm), n_passes[[m]])
    }, ""),
    elpd_lfo_gdp_only = if (g == "1") tot else
      vapply(ms, function(m) sum(st_all[[m]]$pointwise$elpd_sub), 0)
  ) |>
    left_join(mcse_tab, by = "model") |>
    arrange(desc(elpd_lfo))
}
cmp <- bind_rows(summ_group("1"), summ_group("2"), summ_group("3"))

# Pairwise differences (row model minus column model) with SE, within each group, and
# for the GDP-only target across all models on the common window.
pairwise <- function(P, label, gdp = FALSE) {
  n <- nrow(P); out <- list()
  for (a in colnames(P)) for (b in colnames(P)) if (a != b) {
    d <- P[, a] - P[, b]
    out[[length(out) + 1]] <- tibble(comparison = label, model_a = a, model_b = b,
      n_pred = n, elpd_diff = sum(d), se_diff = sqrt(n * var(d)),
      z = sum(d) / sqrt(n * var(d)), se_diff_nw = se_nw(d), z_nw = sum(d) / se_nw(d),
      lag1_acf = if (var(d) > 0) cor(d[-1], d[-n]) else NA_real_,
      # particle-filter Monte Carlo error of the difference (independent filter runs)
      mcse_diff = sqrt(mcse_of(a, gdp)^2 + mcse_of(b, gdp)^2))
  }
  bind_rows(out)
}
pw_tabs <- list()
for (g in names(group_models)) {
  ms <- intersect(group_models[[g]], done)
  if (length(ms) > 1) pw_tabs[[g]] <- pairwise(
    sapply(ms, function(m) st_all[[m]]$pointwise$elpd), paste0("group ", g, " (joint)"))
}
# GDP-only target on the common window: group 1 joint density (= GDP), group 2 and
# m11 GDP-only column. Different information sets (hours / employment in the past,
# and for m11 employment within the target quarter), same target and window.
g1 <- intersect(group_models$`1`, done)
g23 <- intersect(c(group_models$`2`, group_models$`3`), done)
if (length(g1) && length(g23)) {
  Pg <- cbind(sapply(g1, function(m) st_all[[m]]$pointwise$elpd),
              sapply(g23, function(m) st_all[[m]]$pointwise$elpd_sub))
  colnames(Pg) <- c(g1, g23)
  dd <- lapply(c(g1, g23), function(m) st_all[[m]]$pointwise$date)
  stopifnot(all(vapply(dd, identical, TRUE, dd[[1]])))
  pw_tabs$gdp <- pairwise(Pg, "GDP-only target, common window (all groups)", gdp = TRUE)
}
pairs_tab <- bind_rows(pw_tabs)

# Sub-period sums (group 1) and the key differences by sub-period, with iid and
# Newey-West SEs (NW lag min(8, n - 1); with n = 4 in 2020 the NW SE is not meaningful).
periods <- list(`1990-99` = c("1990-01-01", "1999-12-31"), `2000-07` = c("2000-01-01", "2007-12-31"),
                `2008-19` = c("2008-01-01", "2019-12-31"), `2020` = c("2020-01-01", "2020-12-31"),
                `2021-22` = c("2021-01-01", "2022-12-31"), `2023-26` = c("2023-01-01", "2026-12-31"),
                `pre-2020` = c("1990-01-01", "2019-12-31"), `2020 on` = c("2020-01-01", "2026-12-31"),
                `all excl. COVID` = c("1990-01-01", "2026-12-31"))
sub_rows <- list()
key_pairs <- list(c("m12", "m9"), c("m12", "m4"), c("m12", "m5"), c("m9", "m5"), c("m9", "m4"),
                  c("m5", "m2"), c("m9", "m2"), c("m12", "m2"))
if (length(g1)) {
  P1 <- sapply(g1, function(m) st_all[[m]]$pointwise$elpd)
  dts <- as.Date(st_all[[g1[1]]]$pointwise$date)
  for (pn in names(periods)) {
    sel <- dts >= as.Date(periods[[pn]][1]) & dts <= as.Date(periods[[pn]][2])
    if (pn == "all excl. COVID") sel <- sel & !(st_all[[g1[1]]]$pointwise$t %in% covid_idx)
    for (m in g1) sub_rows[[length(sub_rows) + 1]] <- tibble(period = pn, n_pred = sum(sel),
      model_a = m, model_b = NA_character_, elpd = sum(P1[sel, m]),
      se_iid = NA_real_, se_nw = NA_real_)
    for (kp in key_pairs) if (all(kp %in% g1)) {
      d <- P1[sel, kp[1]] - P1[sel, kp[2]]; n <- length(d)
      sub_rows[[length(sub_rows) + 1]] <- tibble(period = pn, n_pred = n, model_a = kp[1],
        model_b = kp[2], elpd = sum(d), se_iid = sqrt(n * var(d)),
        se_nw = se_nw(d, L = min(NW_LAG, n - 1L)))
    }
  }
}
sub_tab <- bind_rows(sub_rows)
write.csv(sub_tab, "output/tables/model_comparison_subperiod.csv", row.names = FALSE)

# Seed-selection sensitivity (from `retrycheck`): the key group-1 differences with the
# targets of every retried fit rescored from attempt 1 alone ("first") or from all
# attempts pooled ("pooled"), each as a correction (variant - kept, same filter seed) to
# the pass-mean pointwise values; and with those targets dropped ("drop").
rc_path <- "output/tables/lfo_retry_sensitivity.csv"
if (file.exists(rc_path) && length(g1)) {
  rc <- read.csv(rc_path, stringsAsFactors = FALSE)
  rc <- rc[rc$model %in% g1, , drop = FALSE]
  tt <- st_all[[g1[1]]]$pointwise$t
  adj <- function(v) {
    P <- P1
    for (m in unique(rc$model)) {
      r <- rc[rc$model == m, , drop = FALSE]; j <- match(r$t, tt)
      P[j, m] <- P[j, m] + r[[paste0("elpd_", v)]] - r$elpd_kept
    }
    P
  }
  affected <- tt %in% rc$t
  variants <- list(main = P1, first = adj("first"), pooled = adj("pooled"))
  sens_rows <- list()
  for (kp in key_pairs) if (all(kp %in% g1)) {
    for (v in c(names(variants), "drop")) {
      d <- if (v == "drop") (P1[, kp[1]] - P1[, kp[2]])[!affected] else
        variants[[v]][, kp[1]] - variants[[v]][, kp[2]]
      n <- length(d)
      sens_rows[[length(sens_rows) + 1]] <- tibble(model_a = kp[1], model_b = kp[2], variant = v,
        n_pred = n, n_targets_rescored = if (v %in% c("first", "pooled")) sum(affected) else 0L,
        elpd_diff = sum(d), se_iid = sqrt(n * var(d)), se_nw = se_nw(d), z_nw = sum(d) / se_nw(d))
    }
  }
  sens_tab <- bind_rows(sens_rows)
  write.csv(sens_tab, "output/tables/lfo_retry_sensitivity_summary.csv", row.names = FALSE)
  cat("\n=== Seed-selection sensitivity of retried refits (a - b) ===\n")
  print(as.data.frame(sens_tab), digits = 4)
}

write.csv(cmp, "output/tables/model_comparison.csv", row.names = FALSE)
write.csv(pairs_tab, "output/tables/model_comparison_pairwise.csv", row.names = FALSE)
write.csv(pw, "output/tables/lfo_pointwise.csv", row.names = FALSE)
write.csv(fit_diag_all, "output/tables/lfo_fit_diagnostics.csv", row.names = FALSE)

cat("\n=== PSIS-LFO-CV, one-step-ahead, marginal likelihood (states integrated out) ===\n")
print(as.data.frame(cmp |> select(group, model, n_pred, elpd_lfo, se_elpd, elpd_diff_vs_best,
                                  se_diff_vs_best, se_diff_vs_best_nw, n_refits, n_k_high_unrefit, max_k_used,
                                  n_fits_below_rule, n_fits_retried, mcse_elpd, elpd_lfo_gdp_only)), digits = 4)
if (nrow(sub_tab)) {
  cat("\n=== Key differences by sub-period (a - b) ===\n")
  print(as.data.frame(sub_tab |> filter(!is.na(model_b))), digits = 3)
}
if (nrow(pairs_tab)) {
  cat("\n=== Pairwise differences (a - b), selected ===\n")
  print(as.data.frame(pairs_tab |> filter(elpd_diff >= 0) |> arrange(comparison, desc(abs(z)))),
        digits = 3)
}
# What was bounded, per model
cat("\n=== Bounds actually hit ===\n")
for (m in done) {
  s <- st_all[[m]]
  hk <- s$pointwise[s$pointwise$k_high_unrefit, ]
  line <- sprintf("%s: %d fits (initial + %d refits, cap %d); %d high-k steps used without refit%s",
                  m, length(s$fits), s$n_refits, s$max_refits, nrow(hk),
                  if (nrow(hk)) sprintf(" (max k %.2f at %s)", max(hk$pareto_k),
                                        hk$date[which.max(hk$pareto_k)]) else "")
  cat(line, "\n")
}
if (!"m11" %in% done) cat("m11 (mixed frequency): LFO not run/completed - not compared\n")

# ---- variance decomposition, per posterior draw -------------------------------------
# ratio_s = var_t(mu_t^(s)) / var_t(headline_t), headline = mean of E, I, P.
# Uses the full-sample fits from R/02_fit_models.R and R/fit_m11.R.
headline_full <- rowMeans(y3)
headline_sub  <- rowMeans(y4[, 1:3])
vrow <- function(series, sample, T, v, hv) {
  tibble(series = series, sample = sample, T = T, sd = sqrt(v), var = v,
         ratio_median = v / hv, ratio_q05 = v / hv, ratio_q95 = v / hv)
}
draw_ratio <- function(fit_file, vname, hv, label, sample) {
  if (!file.exists(fit_file)) return(NULL)
  d <- readRDS(fit_file)$draws(vname, format = "draws_matrix")
  vs <- apply(unclass(d), 1, stats::var)
  r <- vs / hv
  tibble(series = label, sample = sample, T = ncol(d), sd = median(sqrt(vs)),
         var = median(vs), ratio_median = median(r),
         ratio_q05 = unname(quantile(r, 0.05)), ratio_q95 = unname(quantile(r, 0.95)))
}
s1 <- sprintf("full sample %s-%s", format(gdp$date[1], "%Y"), format(max(gdp$date), "%Y"))
s2 <- sprintf("hours subsample %s-%s", format(gdp4$date[1], "%Y"), format(max(gdp4$date), "%Y"))
hv1 <- var(headline_full); hv2 <- var(headline_sub)
vt <- bind_rows(
  vrow("Headline (mean of E, I, P)", s1, nrow(y3), hv1, hv1),
  vrow("GDP(E)", s1, nrow(y3), var(y3[, 1]), hv1),
  vrow("GDP(I)", s1, nrow(y3), var(y3[, 2]), hv1),
  vrow("GDP(P)", s1, nrow(y3), var(y3[, 3]), hv1),
  bind_rows(lapply(group_models$`1`, function(m)
    draw_ratio(sprintf("output/fit_%s.rds", m), "mu", hv1, sprintf("Latent mu (%s)", m), s1))),
  vrow("Headline (mean of E, I, P)", s2, nrow(y4), hv2, hv2),
  bind_rows(lapply(group_models$`2`, function(m)
    draw_ratio(sprintf("output/fit_%s.rds", m), "mu", hv2, sprintf("Latent mu (%s)", m), s2)))
)
# m11: quarterly aggregate of the monthly latent, against headline on m11's quarters
if (file.exists("output/fit_m11.rds")) {
  f11 <- readRDS("output/fit_m11.rds")
  T_q <- f11$metadata()$stan_variable_sizes$mu_q
  emp <- readRDS("data/processed/employment_growth_monthly.rds")
  emp <- emp[emp$date <= max(gdp$date), ]
  q_ok <- which(!is.na(match(gdp$date, emp$date)) & match(gdp$date, emp$date) >= 3)
  if (length(q_ok) == T_q) {
    hv3 <- var(headline_full[q_ok])
    s3 <- sprintf("mixed-frequency quarters %s-%s", format(gdp$date[q_ok[1]], "%Y"),
                  format(gdp$date[max(q_ok)], "%Y"))
    vt <- bind_rows(vt, vrow("Headline (mean of E, I, P)", s3, T_q, hv3, hv3),
                    draw_ratio("output/fit_m11.rds", "mu_q", hv3, "Latent mu_q (m11)", s3))
  } else {
    log_msg("variance decomposition: fit_m11.rds has %d quarters but current data give %d; m11 row skipped (re-run R/fit_m11.R)",
            T_q, length(q_ok))
  }
}
cat("\n=== Variance of growth series relative to headline (per-draw ratios) ===\n")
print(as.data.frame(vt), digits = 3)
write.csv(vt, "output/tables/variance_decomposition.csv", row.names = FALSE)

# ---- remove the invalid LOO objects once the new comparison exists ------------------
if (all(names(specs) %in% done)) {
  stale <- c("output/loos.rds", "output/loos_g1.rds", "output/loos_g2.rds", "output/loo_compare.rds")
  stale <- stale[file.exists(stale)]
  if (length(stale)) { unlink(stale); log_msg("removed stale LOO files: %s", paste(stale, collapse = ", ")) }
}
