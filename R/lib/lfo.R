# PSIS-LFO-CV for 1-step-ahead prediction (Burkner, Gabry & Vehtari 2020,
# "Approximate leave-future-out cross-validation for Bayesian time series models",
# J. Stat. Comput. Simul. 90(14)), on the marginal likelihood: theta = the
# hyperparameters, latent states (mu, h) integrated out by R/lib/filters.R.
#
# For i = L, ..., T-1 the target is elpd_{i+1} = log p(y_{i+1} | y_{1:i}).
#   * A fit on y_{1:i*} (i* <= i) gives draws theta_s ~ p(theta | y_{1:i*}).
#   * Importance ratios  log r_s = sum_{j=i*+1}^{i} log p(y_j | y_{1:j-1}, theta_s)
#     turn them into draws from p(theta | y_{1:i}); PSIS smooths the ratios.
#   * If the Pareto k of the ratios exceeds k_thres (0.7), refit on y_{1:i}
#     (unless the refit cap is reached; then the high-k step is used and logged).
#   * elpd_{i+1} = log sum_s w_s p(y_{i+1} | y_{1:i}, theta_s).
# The per-period densities p(y_j | y_{1:j-1}, theta_s) come from one filter pass of
# the full series per fit (Kalman: exact; RBPF: unbiased estimate).
#
# Particle counts (SV models). With `pf` given, each fit's filter pass uses a per-period
# schedule N_t fixed before the pass (pf_schedule): a pilot pass on pf$pilot_draws
# evenly spaced draws at pf$pilot_n particles measures the per-period ESS fraction
# f_t = ESS_t / N (its pf$q quantile over the pilot draws), and
#   N_t = min(pf$n_max, max(pf$n_base, ceiling(pf$ess_target / f_t)))   (rounded up to 500),
# so N is raised only where the weights collapse (e.g. the COVID quarters) and is capped
# at n_max. The schedule depends on the pilot's random numbers only, not on the pass that
# uses it, so each pass's estimate stays unbiased. The schedule is stored with the fit
# and reused by every re-run pass (lfo_rerun_filter).
#
# Optional "pre" component (used by m11): information that arrives before the target
# within period i+1 (monthly employment inside the quarter). Then the secondary
# target `sub` is log p(GDP_{i+1} | y_{1:i}, pre_{i+1}) and its weights also include
# log p(pre_{i+1} | y_{1:i}, theta_s) in the ratios.

suppressPackageStartupMessages({
  library(loo)
  library(parallel)
})

# Run the filter for all draws, in parallel over `cores` forked workers.
filter_par <- function(model, dm, y, ctx, N, seed, cores) {
  S <- nrow(dm)
  if (cores <= 1 || model_defs[[model]]$type == "kf") {
    return(filter_loglik(model, dm, y, ctx, N = N, seed = seed))
  }
  chunks <- split(seq_len(S), cut(seq_len(S), cores, labels = FALSE))
  res <- mclapply(seq_along(chunks), function(j) {
    filter_loglik(model, dm[chunks[[j]], , drop = FALSE], y, ctx, N = N,
                  seed = seed * 100L + j)
  }, mc.cores = cores, mc.preschedule = FALSE)
  bad <- vapply(res, function(r) inherits(r, "try-error") || is.null(r$full), TRUE)
  if (any(bad)) stop("filter worker failed: ", paste(sapply(res[bad], as.character), collapse = "; "))
  list(full = do.call(rbind, lapply(res, `[[`, "full")),
       sub  = if (is.null(res[[1]]$sub)) NULL else do.call(rbind, lapply(res, `[[`, "sub")),
       ess_min = unlist(lapply(res, `[[`, "ess_min")),
       ess_frac = if (is.null(res[[1]]$ess_frac)) NULL else do.call(rbind, lapply(res, `[[`, "ess_frac")))
}

# Per-period particle-count schedule from a pilot pass (see the header). Returns an
# integer vector of length nrow(y) with attribute "f" (the pilot ESS-fraction quantile).
pf_schedule <- function(model, dm, y, ctx, pf, seed, cores) {
  S <- nrow(dm)
  idx <- unique(round(seq(1, S, length.out = min(S, pf$pilot_draws))))
  r <- filter_par(model, dm[idx, , drop = FALSE], y, ctx, pf$pilot_n, seed, cores)
  if (is.null(r$ess_frac)) stop("filter for ", model, " does not report ess_frac")
  f <- apply(r$ess_frac, 2, quantile, probs = pf$q, na.rm = TRUE, names = FALSE)
  f[!is.finite(f) | f <= 0] <- 1e-12
  Nt <- pmin(pf$n_max, pmax(pf$n_base, ceiling(pf$ess_target / f / 500) * 500))
  Nt <- as.integer(Nt)
  attr(Nt, "f") <- f
  Nt
}
sched_summary <- function(Nt) {
  sprintf("N_t: base %d in %d periods, raised in %d (max %d, in %d periods), cost x%.2f of N=%d",
          min(Nt), sum(Nt == min(Nt)), sum(Nt > min(Nt)), max(Nt), sum(Nt == max(Nt) & max(Nt) > min(Nt)),
          mean(Nt) / min(Nt), min(Nt))
}

psis_w <- function(lr, chain_id) {
  r_eff <- relative_eff(exp(lr - max(lr)), chain_id = chain_id)
  ps <- suppressWarnings(psis(lr, r_eff = r_eff))
  list(k = unname(pareto_k_values(ps)), lw = as.vector(weights(ps, log = TRUE, normalize = TRUE)))
}

# Weights and pointwise elpd for target i+1 from a fit at cur_i with filter output ll.
lfo_step <- function(ll, cur_i, i, chain_id) {
  S <- nrow(ll$full)
  lr <- if (i > cur_i) rowSums(ll$full[, (cur_i + 1):i, drop = FALSE]) else NULL
  w_joint <- if (is.null(lr)) list(k = NA_real_, lw = rep(-log(S), S)) else psis_w(lr, chain_id)
  w_tgt <- w_joint
  if (!is.null(ll$pre)) w_tgt <- psis_w((if (is.null(lr)) 0 else lr) + ll$pre[, i + 1], chain_id)
  list(k = w_joint$k, k_tgt = w_tgt$k,
       ess_w = 1 / sum(exp(2 * w_joint$lw)),
       elpd = lse(w_joint$lw + ll$full[, i + 1]),
       elpd_sub = if (is.null(ll$sub)) NA_real_ else lse(w_tgt$lw + ll$sub[, i + 1]))
}

# psis_lfo: run (or resume) PSIS-LFO for one model.
#   refit_fun(i) -> list(dm = draws matrix of theta from a fit on y_{1:i},
#                        chain_id = integer vector, diag = one-row data.frame,
#                        attempts = optional data.frame, one row per sampling attempt
#                        (failed ones included); stored in st$attempts)
#   ctx_fun(i)   -> filter context for a fit on y_{1:i} (e.g. m7 covid indices)
#   y            -> observation matrix passed to the filter (all periods)
#   dates        -> length-T vector of period labels (T = number of LFO periods)
#   ll_post      -> maps the filter output to LFO periods (identity except for m11,
#                   where months are aggregated to quarters)
psis_lfo <- function(model, y, dates, L, refit_fun, ctx_fun, k_thres = 0.7,
                     max_refits = 15L, N = 1000L, cores = 4L, ckpt = NULL,
                     log_fun = message, ll_post = identity, pf = NULL) {
  T <- length(dates)
  st <- if (!is.null(ckpt) && file.exists(ckpt)) readRDS(ckpt) else NULL
  if (!is.null(st) && isTRUE(st$done)) {
    log_fun(sprintf("[%s] complete checkpoint found, skipping", model)); return(st)
  }
  do_fit <- function(i) {
    t0 <- Sys.time()
    f <- refit_fun(i)
    seed <- 7000L + i
    Nt <- N
    if (!is.null(pf) && model_defs[[model]]$type == "sv") {
      t1 <- Sys.time()
      Nt <- pf_schedule(model, f$dm, y, ctx_fun(i), pf, seed + 900000L, cores)
      log_fun(sprintf("[%s] fit 1:%d particle schedule (pilot %d draws, N=%d, %.0f s): %s", model, i,
                      min(nrow(f$dm), pf$pilot_draws), pf$pilot_n,
                      as.numeric(difftime(Sys.time(), t1, units = "secs")), sched_summary(Nt)))
    }
    ll <- ll_post(filter_par(model, f$dm, y, ctx_fun(i), Nt, seed, cores))
    if (anyNA(ll$full) || any(ll$full > 50, na.rm = TRUE))
      stop(sprintf("[%s] filter output on fit 1:%d has NaN or implausible values (max %.3g)",
                   model, i, max(ll$full, na.rm = TRUE)))
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    log_fun(sprintf("[%s] fit on periods 1:%d (to %s): %d draws, %.0f s; %s", model, i,
                    format(dates[i]), nrow(f$dm), secs,
                    paste(sprintf("%s=%s", names(f$diag),
                                  sapply(f$diag, function(v)
                                    if (is.numeric(v) && !isTRUE(all.equal(v, round(v)))) format(signif(v, 4))
                                    else if (is.numeric(v)) sprintf("%.0f", v) else format(v))),
                          collapse = " ")))
    list(i = i, dm = f$dm, chain_id = f$chain_id, ll = ll, seed = seed, Nt = Nt,
         attempts = if (is.null(f$attempts)) NULL else
           cbind(data.frame(model = model, i_fit = i, date_fit = format(dates[i])), f$attempts),
         diag = cbind(data.frame(model = model, i_fit = i, date_fit = format(dates[i]),
                                 secs = round(secs)), f$diag))
  }
  if (is.null(st)) {
    cur <- do_fit(L)
    st <- list(model = model, L = L, T = T, k_thres = k_thres, max_refits = max_refits,
               N = N, pf = pf, cur = cur, fits = list(cur[c("i", "dm", "chain_id", "seed", "Nt")]),
               fit_diag = cur$diag, attempts = cur$attempts, rows = list(), next_i = L,
               n_refits = 0L, done = FALSE)
    if (!is.null(ckpt)) saveRDS(st, ckpt)
  }
  for (i in st$next_i:(T - 1)) {
    s1 <- lfo_step(st$cur$ll, st$cur$i, i, st$cur$chain_id)
    k_dec <- suppressWarnings(max(s1$k, s1$k_tgt, na.rm = TRUE))  # -Inf if both NA
    refit <- FALSE; capped <- FALSE
    if (i > st$cur$i && k_dec > k_thres) {
      if (st$n_refits < max_refits) {
        cur <- do_fit(i)
        st$cur <- cur
        st$fits[[length(st$fits) + 1]] <- cur[c("i", "dm", "chain_id", "seed", "Nt")]
        st$fit_diag <- dplyr::bind_rows(st$fit_diag, cur$diag)
        st$attempts <- dplyr::bind_rows(st$attempts, cur$attempts)   # every refit attempt
        st$n_refits <- st$n_refits + 1L
        refit <- TRUE
        s2 <- lfo_step(cur$ll, cur$i, i, cur$chain_id)
        s1$elpd <- s2$elpd; s1$elpd_sub <- s2$elpd_sub; s1$ess_w <- s2$ess_w
        s1$k_tgt_after <- s2$k_tgt
      } else {
        capped <- TRUE
      }
    }
    st$rows[[length(st$rows) + 1]] <- data.frame(
      model = model, i_train = i, t = i + 1, date = format(dates[i + 1]),
      i_fit = st$cur$i, pareto_k = if (is.finite(k_dec)) k_dec else NA_real_,
      refit = refit, k_high_unrefit = capped, ess_w = s1$ess_w,
      elpd = s1$elpd, elpd_sub = s1$elpd_sub)
    st$next_i <- i + 1
    if (refit && !is.null(ckpt)) saveRDS(st, ckpt)
  }
  st$done <- TRUE
  st$pointwise <- do.call(rbind, st$rows)
  st$cur$ll <- NULL           # large; everything needed is in st$fits / st$pointwise
  if (!is.null(ckpt)) saveRDS(st, ckpt)
  log_fun(sprintf("[%s] done: elpd_lfo = %.2f over %d predictions, %d refits, %d high-k steps used",
                  model, sum(st$pointwise$elpd), nrow(st$pointwise), st$n_refits,
                  sum(st$pointwise$k_high_unrefit)))
  st
}

# Recompute the LFO pointwise elpd with a fresh filter seed, holding the fitted draws,
# the refit schedule and each fit's particle schedule (unless N is given) fixed. For the RBPF models the spread between passes is the
# Monte Carlo error of the particle filter.
lfo_rerun_filter <- function(st, y, ctx_fun, seed_offset, N = st$N, cores = 4L,
                             ll_post = identity) {
  pw <- st$pointwise
  out <- pw
  for (f in st$fits) {
    Nf <- if (!is.null(f$Nt) && missing(N)) f$Nt else N
    ll <- ll_post(filter_par(st$model, f$dm, y, ctx_fun(f$i), Nf, f$seed + seed_offset, cores))
    if (anyNA(ll$full) || any(ll$full > 50, na.rm = TRUE)) stop("re-run filter output has NaN or implausible values")
    for (r in which(pw$i_fit == f$i)) {
      s <- lfo_step(ll, f$i, pw$i_train[r], f$chain_id)
      out$elpd[r] <- s$elpd
      out$elpd_sub[r] <- s$elpd_sub
    }
  }
  out
}
