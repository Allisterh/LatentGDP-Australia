# Fit m13 (stan/m13_serial_latent_sv.stan: m12 + stochastic volatility on the latent
# innovation) to the full sample and evaluate it with the production quality gate of
# R/fit_production.R.
#
#   Rscript R/fit_m13.R                      # seed 1, saves output/fit_m13.rds
#   Rscript R/fit_m13.R --seed=2             # saves output/fit_m13_seed2.rds
#   Rscript R/fit_m13.R --seed=2 --out=PATH  # another seed, saved to PATH
#   Rscript R/fit_m13.R --parallel=2         # parallel chains (default 4)
#   Rscript R/fit_m13.R --omega_g_sd=1 --out=output/fit_m13_omegag1.rds
#                                            # prior sensitivity: omega_g ~ half-N(0, 1)
#   Rscript R/fit_m13.R stall --chains=100 --init=2 --timeout=300 --parallel=4 --dir=DIR
#
# Fit mode: the production settings (4 chains x (2000 warmup + 2000 draws), adapt_delta
# 0.99, max_treedepth 12). The gate is fit_production.R's (same thresholds, same NA
# handling); its hyperparameter list is extended with psi, mu_g, phi_g and omega_g (tau is
# the generated quantity exp(mu_g / 2)), and R-hat <= 1.05 is also required on every g[t].
# The fit is saved before the gate; the gate record is saved to <out>_gate.rds and, for
# the standard prior, its rows replace this fit's rows in output/tables/m13_gate_seeds.csv
# (--gate_csv=PATH to change, --gate_csv=none to skip); a failing gate then ends the script
# with an error (non-zero exit status).
#
# --omega_g_sd=S (S != 0.5) fits a copy of the Stan program, written to --stan_dir (default
# tempdir()), whose only change is the prior line omega_g ~ normal(0, S); its gate record is
# not written to m13_gate_seeds.csv unless --gate_csv is given.
#
# Stall mode: the warmup-stall harness of R/fit_m12.R (from m9). Many independent
# single-chain runs of 150 warmup iterations from Stan's random inits (init = 2 by
# default); a chain counts as stalled if it has not finished within --timeout seconds or
# its last step size is below 1e-6. Writes DIR/stall_summary.rds and prints the stalled
# chains.

#
# Sample window and output root: R/sample_window.R (LATENTGDP_START, LATENTGDP_OUT). The
# default --out and gate-record paths are under the output root, e.g.
#   LATENTGDP_START=1980-03-01 LATENTGDP_OUT=output/sample_1980 Rscript R/fit_m13.R --seed=2
# writes output/sample_1980/fit_m13_seed2.rds and output/sample_1980/tables/m13_gate_seeds.csv. Both unset = the full sample and the
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
stan_file <- "stan/m13_serial_latent_sv.stan"
omega_g_sd <- as.numeric(opt("omega_g_sd", "0.5"))
std_prior <- isTRUE(all.equal(omega_g_sd, 0.5))
if (!std_prior) {
  # Prior-sensitivity copy: one line changed, checked to occur exactly once.
  src <- readLines(stan_file)
  pat <- "^(\\s*)omega_g ~ normal\\(0, 0\\.5\\);"
  hit <- grep(pat, src)
  if (length(hit) != 1) stop("omega_g prior line not found exactly once in ", stan_file)
  src[hit] <- sub(pat, sprintf("\\1omega_g ~ normal(0, %s);", format(omega_g_sd)), src[hit])
  sdir <- opt("stan_dir", file.path(tempdir(), "m13_prior"))
  dir.create(sdir, showWarnings = FALSE, recursive = TRUE)
  stan_file <- file.path(sdir, sprintf("m13_omegag_sd%s.stan", gsub("\\.", "p", format(omega_g_sd))))
  writeLines(src, stan_file)
  cat(sprintf("Prior sensitivity: omega_g ~ half-N(0, %s); Stan file %s\n  changed line: %s\n",
              format(omega_g_sd), stan_file, trimws(src[hit])))
}
mod <- cmdstan_model(stan_file)

# ---------------------------------------------------------------------------
# Stall harness (as R/fit_m12.R)
# ---------------------------------------------------------------------------
if (mode == "stall") {
  n_chains <- as.integer(opt("chains", "100"))
  init     <- opt("init", "2")
  timeout  <- as.numeric(opt("timeout", "300"))
  par      <- as.integer(opt("parallel", "2"))
  seed0    <- as.integer(opt("seed0", "1"))
  od       <- opt("dir", file.path(tempdir(), "m13_stall"))
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
                                               rho = NA_real_, omega_g = NA_real_, n = 0L,
                                               td_max = NA_real_))
    list(step = tail(x$stepsize__, 1), rho_EI = tail(x$rho_EI, 1), rho = tail(x$rho, 1),
         omega_g = tail(x$omega_g, 1), n = nrow(x), td_max = max(x$treedepth__))
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
                                 last_omega_g = l$omega_g,
                                 stalled = timed_out || (!is.na(l$step) && l$step < 1e-6))
        cat(sprintf("seed %d: %.0f s, last step size %s, omega_g %s, stalled %s\n", s, el,
                    format(signif(l$step, 3)), format(signif(l$omega_g, 3)), res[[r$i]]$stalled))
        running[[key]] <- NULL
        saveRDS(do.call(rbind, res[!vapply(res, is.null, TRUE)]),
                file.path(od, "stall_partial.rds"))
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
out  <- opt("out", if (seed == 1L) sw_out("fit_m13.rds") else
  sw_out(sprintf("fit_m13_seed%d.rds", seed)))
par  <- as.integer(opt("parallel", "4"))
gate_csv <- opt("gate_csv", if (std_prior) sw_out("tables", "m13_gate_seeds.csv") else "none")
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)

t_wall0 <- Sys.time()
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
wall_secs <- as.numeric(difftime(Sys.time(), t_wall0, units = "secs"))
fit$save_object(file = out)
cat(sprintf("\nSaved %s\n", out))

cat(sprintf("\n=== m13 fit (seed %d, omega_g ~ half-N(0, %s)) ===\n", seed, format(omega_g_sd)))
print(fit$summary(variables = c("mu_bar", "rho", "tau", "mu_g", "phi_g", "omega_g", "rho_EI",
                                "psi", "phi", "omega", "mu_h", "eps_sd_bar")), n = 40)
tm <- fit$time()
cat(sprintf("Sampling time (fit$time()$total): %.0f s; wall time of $sample(): %.0f s (parallel_chains = %d); per-chain total s: %s\n",
            tm$total, wall_secs, par, paste(round(tm$chains$total), collapse = ", ")))

# ---------------------------------------------------------------------------
# Convergence gate: R/fit_production.R's thresholds and logic; hyper_vars extended with
# psi, mu_g, phi_g, omega_g, and a state R-hat row for g[t].
# ---------------------------------------------------------------------------
GATE <- list(
  max_div_frac   = 0.001,
  min_ebfmi      = 0.3,
  ebfmi_floor    = 0.2,
  max_rhat_hyper = 1.01,
  max_rhat_state = 1.05,
  min_ess        = 400
)
hyper_vars <- c("c", "rho", "tau", "mu_h", "phi", "omega", "rho_EI", "psi", "mu_bar",
                "mu_g", "phi_g", "omega_g")

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
s_g     <- summ_conv("g")
stopifnot(nrow(s_mu) == stan_data$T, nrow(s_h) == stan_data$T * stan_data$K,
          nrow(s_g) == stan_data$T)

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
    "max R-hat, g[t]",
    "min bulk ESS, hyperparameters",
    "min tail ESS, hyperparameters",
    "min bulk ESS, mu[t]",
    "min tail ESS, mu[t]"
  ),
  value = c(div_frac, min_eb,
            worst_max(s_hyper$rhat), worst_max(s_mu$rhat), worst_max(s_h$rhat),
            worst_max(s_g$rhat),
            worst_min(s_hyper$ess_bulk), worst_min(s_hyper$ess_tail),
            worst_min(s_mu$ess_bulk), worst_min(s_mu$ess_tail)),
  op = c("<", ">=", "<=", "<=", "<=", "<=", ">=", ">=", ">=", ">="),
  threshold = c(GATE$max_div_frac, GATE$min_ebfmi,
                GATE$max_rhat_hyper, GATE$max_rhat_state, GATE$max_rhat_state,
                GATE$max_rhat_state,
                GATE$min_ess, GATE$min_ess, GATE$min_ess, GATE$min_ess),
  worst = c(
    "",
    sprintf("chain %d", if (anyNA(ebfmi)) which(is.na(ebfmi))[1] else which.min(ebfmi)),
    which_at(s_hyper, "rhat", which.max),
    which_at(s_mu, "rhat", which.max),
    which_at(s_h, "rhat", which.max),
    which_at(s_g, "rhat", which.max),
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

cat("\n=== Convergence gate (R/fit_production.R thresholds; hyperparameters + psi, mu_g, phi_g, omega_g; R-hat on g[t]) ===\n")
cat(sprintf("Per-chain E-BFMI: %s\n", paste(sprintf("%.3f", ebfmi), collapse = ", ")))
cat(sprintf("Max-treedepth hits: %d of %d draws (informational)\n", n_tdmax, n_draws))
cat(sprintf("%-46s %10s  %-10s %-6s %s\n", "metric", "value", "threshold", "status", "worst"))
for (i in seq_len(nrow(checks))) {
  cat(sprintf("%-46s %10s  %-10s %-6s %s\n",
              checks$metric[i], fmt_val(checks$value[i], i), fmt_th(i),
              if (checks$pass[i]) "PASS" else "FAIL", checks$worst[i]))
}

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
gate_pass <- length(failures) == 0

saveRDS(list(seed = seed, omega_g_sd = omega_g_sd, checks = checks, ebfmi = ebfmi,
             n_tdmax = n_tdmax, time = tm, wall_secs = wall_secs, parallel_chains = par,
             gate_pass = gate_pass, failures = failures),
        sub("\\.rds$", "_gate.rds", out))

if (gate_csv != "none") {
  label <- sprintf("%s (R/fit_m13.R%s%s)", out, if (seed != 1L) sprintf(" --seed=%d", seed) else "",
                   if (par != 4L) sprintf(" --parallel=%d", par) else "")
  rows <- data.frame(
    fit = label, seed = seed, metric = checks$metric, value = checks$value, op = checks$op,
    threshold = checks$threshold, worst = checks$worst, pass = checks$pass,
    gate_pass = gate_pass,
    ebfmi_by_chain = paste(sprintf("%.3f", ebfmi), collapse = " "),
    treedepth_hits = n_tdmax, parallel_chains = par,
    sampling_secs_total = round(tm$total), wall_secs = round(wall_secs),
    chain_secs = paste(round(tm$chains$total), collapse = " "),
    stringsAsFactors = FALSE)
  if (file.exists(gate_csv)) {
    old <- read.csv(gate_csv, stringsAsFactors = FALSE)
    old <- old[old$seed != seed, , drop = FALSE]        # one record per seed
    rows <- rbind(old[, names(rows)], rows)
    rows <- rows[order(rows$seed), , drop = FALSE]
  }
  dir.create(dirname(gate_csv), showWarnings = FALSE, recursive = TRUE)
  write.csv(rows, gate_csv, row.names = FALSE)
  cat(sprintf("Gate record written to %s\n", gate_csv))
}

if (!gate_pass) {
  stop(paste0("Fit quality gate FAILED (", length(failures), " problem",
              if (length(failures) > 1) "s" else "", "):\n  - ",
              paste(failures, collapse = "\n  - ")), call. = FALSE)
}
cat("\nFit quality gate passed: all checks within thresholds.\n")
