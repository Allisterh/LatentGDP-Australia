# Fit m12 (stan/m12_struct_sv_serial.stan: m9 + AR(1) measurement errors) to the full
# sample and evaluate it with the production quality gate of R/fit_production.R.
#
#   Rscript R/fit_m12.R                      # seed 1, saves output/fit_m12.rds
#   Rscript R/fit_m12.R --seed=2 --out=PATH  # another seed, saved to PATH
#   Rscript R/fit_m12.R --parallel=2         # parallel chains (default 4)
#   Rscript R/fit_m12.R stall --chains=100 --init=2 --timeout=180 --parallel=2 --dir=DIR
#
# Fit mode: the production settings (4 chains x (2000 warmup + 2000 draws), adapt_delta
# 0.99, max_treedepth 12). The gate is fit_production.R's, unchanged (same thresholds,
# same NA handling); its hyperparameter list is extended with psi. The fit is saved before
# the gate; a failing gate ends the script with an error (non-zero exit status).
#
# Stall mode: the warmup-stall harness used for m9 (see the header of
# stan/m9_struct_sv.stan). Many independent single-chain runs of 150 warmup iterations
# from Stan's random inits (init = 2 by default); a chain counts as stalled if it has not
# finished within --timeout seconds or its last step size is below 1e-6 (step-size
# collapse). Writes DIR/stall_summary.rds and prints the stalled chains.

#
# Sample window and output root: R/sample_window.R (LATENTGDP_START, LATENTGDP_OUT). The
# default --out and gate-record paths are under the output root, e.g.
#   LATENTGDP_START=1980-03-01 LATENTGDP_OUT=output/sample_1980 Rscript R/fit_m12.R --seed=2
# writes output/sample_1980/fit_m12_seed2.rds. Both unset = the full sample and the
# paths above, exactly as before.

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
})
source("R/sample_window.R")
cat(sw_describe(), "\n")

args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[length(hit)]) else default
}
mode <- if (length(args) && !startsWith(args[1], "--")) args[1] else "fit"
stopifnot(mode %in% c("fit", "stall"))

gdp_growth <- sw_trim_q(readRDS("data/processed/gdp_growth.rds"))
stan_data <- list(
  T = nrow(gdp_growth),
  K = 3,
  y = as.matrix(gdp_growth[, c("g_E", "g_I", "g_P")])
)
stan_file <- "stan/m12_struct_sv_serial.stan"
mod <- cmdstan_model(stan_file)

# ---------------------------------------------------------------------------
# Stall harness
# ---------------------------------------------------------------------------
if (mode == "stall") {
  n_chains <- as.integer(opt("chains", "100"))
  init     <- opt("init", "2")
  timeout  <- as.numeric(opt("timeout", "180"))
  par      <- as.integer(opt("parallel", "2"))
  seed0    <- as.integer(opt("seed0", "1"))
  od       <- opt("dir", file.path(tempdir(), "m12_stall"))
  dir.create(od, showWarnings = FALSE, recursive = TRUE)
  jf <- file.path(od, "data.json")
  write_stan_json(stan_data, jf)
  exe <- mod$exe_file()
  run_one <- function(s) {
    csv <- file.path(od, sprintf("s%d.csv", s))
    processx::process$new(exe, c(
      "sample", "num_warmup=150", "num_samples=1", "save_warmup=1",
      "adapt", "delta=0.99", "algorithm=hmc", "engine=nuts", "max_depth=12",
      "data", paste0("file=", jf), paste0("init=", init), "random", paste0("seed=", s),
      "output", paste0("file=", csv), "refresh=0"), stdout = NULL, stderr = NULL)
  }
  read_last <- function(s) {
    x <- tryCatch(suppressWarnings(read.csv(file.path(od, sprintf("s%d.csv", s)),
                                            comment.char = "#")), error = function(e) NULL)
    if (is.null(x) || nrow(x) == 0) return(list(step = NA_real_, rho_EI = NA_real_,
                                               rho = NA_real_, n = 0L, td_max = NA_real_))
    list(step = tail(x$stepsize__, 1), rho_EI = tail(x$rho_EI, 1), rho = tail(x$rho, 1),
         n = nrow(x), td_max = max(x$treedepth__))
  }
  seeds <- seed0 + seq_len(n_chains) - 1L
  res <- vector("list", n_chains)
  running <- list()
  nxt <- 1L
  while (nxt <= n_chains || length(running)) {
    while (length(running) < par && nxt <= n_chains) {
      running[[as.character(seeds[nxt])]] <- list(p = run_one(seeds[nxt]), t0 = Sys.time(),
                                                  i = nxt)
      nxt <- nxt + 1L
    }
    Sys.sleep(1)
    for (key in names(running)) {
      r <- running[[key]]
      el <- as.numeric(difftime(Sys.time(), r$t0, units = "secs"))
      if (!r$p$is_alive() || el > timeout) {
        timed_out <- r$p$is_alive()
        if (timed_out) r$p$kill()
        s <- as.integer(key)
        l <- read_last(s)
        res[[r$i]] <- data.frame(seed = s, secs = round(el, 1), timed_out = timed_out,
                                 exit = if (timed_out) NA_integer_ else r$p$get_exit_status(),
                                 rows = l$n, last_step = l$step, max_treedepth = l$td_max,
                                 last_rho_EI = l$rho_EI, last_rho = l$rho,
                                 stalled = timed_out || (!is.na(l$step) && l$step < 1e-6))
        cat(sprintf("seed %d: %.0f s, last step size %s, stalled %s\n", s, el,
                    format(signif(l$step, 3)), res[[r$i]]$stalled))
        running[[key]] <- NULL
      }
    }
  }
  res <- do.call(rbind, res)
  saveRDS(res, file.path(od, "stall_summary.rds"))
  cat(sprintf("\nTOTAL init=%s: %d stalled of %d chains (timeout %.0f s); run time median %.0f s, max %.0f s\n",
              init, sum(res$stalled), nrow(res), timeout, median(res$secs), max(res$secs)))
  if (any(res$stalled)) print(res[res$stalled, ])
  quit(save = "no", status = 0)
}

# ---------------------------------------------------------------------------
# Full fit
# ---------------------------------------------------------------------------
seed <- as.integer(opt("seed", "1"))
out  <- opt("out", if (seed == 1L) sw_out("fit_m12.rds") else
  sw_out(sprintf("fit_m12_seed%d.rds", seed)))
par  <- as.integer(opt("parallel", "4"))
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)

fit <- mod$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = par,
  iter_warmup = 2000,
  iter_sampling = 2000,
  seed = seed,
  refresh = 500,
  adapt_delta = 0.99,
  max_treedepth = 12
)
fit$save_object(file = out)
cat(sprintf("\nSaved %s\n", out))

cat(sprintf("\n=== m12 fit (seed %d) ===\n", seed))
print(fit$summary(variables = c("mu_bar", "rho", "tau", "rho_EI", "psi", "phi", "omega",
                                "mu_h", "eps_sd_bar")), n = 40)
tm <- fit$time()
cat(sprintf("Wall time: %.0f s (parallel_chains = %d); per-chain total s: %s\n", tm$total, par,
            paste(round(tm$chains$total), collapse = ", ")))

# ---------------------------------------------------------------------------
# Convergence gate: R/fit_production.R's, verbatim except that hyper_vars has psi.
# ---------------------------------------------------------------------------
GATE <- list(
  max_div_frac   = 0.001,
  min_ebfmi      = 0.3,
  ebfmi_floor    = 0.2,
  max_rhat_hyper = 1.01,
  max_rhat_state = 1.05,
  min_ess        = 400
)
hyper_vars <- c("c", "rho", "tau", "mu_h", "phi", "omega", "rho_EI", "mu_bar", "psi")

n_chains <- fit$num_chains()
sdiag    <- fit$sampler_diagnostics(format = "draws_df")
n_draws  <- nrow(sdiag)
n_div    <- sum(sdiag$divergent__)
n_tdmax  <- sum(sdiag$treedepth__ >= fit$metadata()$max_treedepth)
ebfmi    <- fit$diagnostic_summary(diagnostics = "ebfmi", quiet = TRUE)$ebfmi

summ_conv <- function(vars) {
  posterior::summarise_draws(fit$draws(variables = vars), "rhat", "ess_bulk", "ess_tail")
}
s_hyper <- summ_conv(hyper_vars)
s_mu    <- summ_conv("mu")
s_h     <- summ_conv("h")
stopifnot(nrow(s_mu) == stan_data$T, nrow(s_h) == stan_data$T * stan_data$K)

worst_max <- function(x) if (anyNA(x)) NA_real_ else max(x)
worst_min <- function(x) if (anyNA(x)) NA_real_ else min(x)
which_at  <- function(s, col, fn) {
  v <- s[[col]]
  if (anyNA(v)) return(paste(s$variable[is.na(v)][1], "(NA)"))
  s$variable[fn(v)]
}
div_frac <- n_div / n_draws
min_eb   <- worst_min(ebfmi)

checks <- data.frame(
  metric = c(
    sprintf("divergences (%d of %d draws, %d chains)", n_div, n_draws, n_chains),
    "min E-BFMI over chains",
    "max R-hat, hyperparameters",
    "max R-hat, mu[t]",
    "max R-hat, h[t,k]",
    "min bulk ESS, hyperparameters",
    "min tail ESS, hyperparameters",
    "min bulk ESS, mu[t]",
    "min tail ESS, mu[t]"
  ),
  value = c(div_frac, min_eb,
            worst_max(s_hyper$rhat), worst_max(s_mu$rhat), worst_max(s_h$rhat),
            worst_min(s_hyper$ess_bulk), worst_min(s_hyper$ess_tail),
            worst_min(s_mu$ess_bulk), worst_min(s_mu$ess_tail)),
  op = c("<", ">=", "<=", "<=", "<=", ">=", ">=", ">=", ">="),
  threshold = c(GATE$max_div_frac, GATE$min_ebfmi,
                GATE$max_rhat_hyper, GATE$max_rhat_state, GATE$max_rhat_state,
                GATE$min_ess, GATE$min_ess, GATE$min_ess, GATE$min_ess),
  worst = c(
    "",
    sprintf("chain %d", if (anyNA(ebfmi)) which(is.na(ebfmi))[1] else which.min(ebfmi)),
    which_at(s_hyper, "rhat", which.max),
    which_at(s_mu, "rhat", which.max),
    which_at(s_h, "rhat", which.max),
    which_at(s_hyper, "ess_bulk", which.min),
    which_at(s_hyper, "ess_tail", which.min),
    which_at(s_mu, "ess_bulk", which.min),
    which_at(s_mu, "ess_tail", which.min)
  ),
  stringsAsFactors = FALSE
)
checks$pass <- mapply(function(v, op, th) {
  if (is.na(v)) return(FALSE)
  switch(op, "<" = v < th, "<=" = v <= th, ">=" = v >= th)
}, checks$value, checks$op, checks$threshold)

fmt_val <- function(v, i) {
  if (is.na(v)) return("NA")
  if (i == 1) return(sprintf("%.3f%%", 100 * v))
  if (checks$op[i] == "<=") return(sprintf("%.4f", v))
  if (grepl("ESS", checks$metric[i])) return(sprintf("%.0f", v))
  sprintf("%.3f", v)
}
fmt_th <- function(i) {
  th <- checks$threshold[i]
  if (i == 1) return(sprintf("%s %.1f%%", checks$op[i], 100 * th))
  sprintf("%s %s", checks$op[i], format(th))
}

cat("\n=== Convergence gate (R/fit_production.R thresholds; hyperparameters + psi) ===\n")
cat(sprintf("Per-chain E-BFMI: %s\n", paste(sprintf("%.3f", ebfmi), collapse = ", ")))
cat(sprintf("Max-treedepth hits: %d of %d draws (informational)\n", n_tdmax, n_draws))
cat(sprintf("%-46s %10s  %-10s %-6s %s\n", "metric", "value", "threshold", "status", "worst"))
for (i in seq_len(nrow(checks))) {
  cat(sprintf("%-46s %10s  %-10s %-6s %s\n",
              checks$metric[i], fmt_val(checks$value[i], i), fmt_th(i),
              if (checks$pass[i]) "PASS" else "FAIL", checks$worst[i]))
}
saveRDS(list(seed = seed, checks = checks, ebfmi = ebfmi, n_tdmax = n_tdmax,
             time = tm), sub("\\.rds$", "_gate.rds", out))

failures <- character(0)
for (i in which(!checks$pass)) {
  failures <- c(failures, sprintf("%s = %s (required %s; worst: %s)",
                                  checks$metric[i], fmt_val(checks$value[i], i), fmt_th(i),
                                  checks$worst[i]))
}
if (!is.na(min_eb) && min_eb < GATE$ebfmi_floor) {
  failures <- c(failures, sprintf("min E-BFMI %.3f is below the hard floor %.1f",
                                  min_eb, GATE$ebfmi_floor))
}
if (length(failures) > 0) {
  stop(paste0("Fit quality gate FAILED (", length(failures), " problem",
              if (length(failures) > 1) "s" else "", "):\n  - ",
              paste(failures, collapse = "\n  - ")), call. = FALSE)
}
cat("\nFit quality gate passed: all checks within thresholds.\n")
