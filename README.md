# LatentGDP

This project produces a Bayesian state-space estimate of Australian quarterly
GDP growth. It filters measurement error out of the three ABS measures of GDP
(expenditure-, income- and production-based). It replicates and extends RBA
RDP 2014-12 (Rees, Lancaster & Finlay).

The latest estimate is published at **https://davidastephan.github.io/LatentGDP/**.
It is rebuilt automatically after each ABS 5206.0 release.

## What this measures

The headline ABS GDP figure is the *average* of three separately constructed
quarterly measures: GDP(E), GDP(I) and GDP(P). They disagree quarter by
quarter. This project treats *true* GDP growth as a latent state, observed
through three noisy signals, and infers it with Bayesian methods.

The latent series is less volatile than headline GDP:
- **Production model (m12).** On the current vintage (1959 Q4 - 2026 Q2,
  T = 267), the variance of the latent path is **78% of the variance of
  headline growth**. That is the posterior median of var(mu draw) / var(headline),
  with a 90% interval of 71% to 85%.
- **Other models.** The other full-sample specifications give 78-83% (m9, the
  previous production model, 79%). The exception is m7, the COVID-scale
  model, at 61%. See NOTES.md for the full table.

So, according to m12, roughly a fifth of the quarter-to-quarter variance of
headline growth is measurement error rather than movement in true activity.
The model comparison does not test this split (see "What it does not test"
below). Earlier versions of
this project claimed a 30-35% reduction. That figure is too large: it was the
variance of the posterior-*mean* path, which is smoother than any single
plausible path.

## Production model

The production model is `stan/m12_struct_sv_serial.stan` (since 2026-10-17;
before that it was m9, `stan/m9_struct_sv.stan`). It has five parts:
- **Latent state.** Latent growth follows a stationary AR(1),
  mu_t = c + rho·mu_{t-1} + tau·eta_t, with mu_1 drawn from the stationary
  distribution.
- **Observations.** Each measure equals mu_t plus a measurement error.
- **Serially correlated errors.** Each error follows its own AR(1),
  eps_t^k = psi_k·eps_{t-1}^k + u_t^k, as in m5. The period-1 errors start
  from the stationary distribution at the period-1 variances.
- **Structural Sigma.** The E and I innovations have a constant correlation
  rho_EI. P's is independent of both. The ABS builds Production from separate
  industry surveys, while E and I share national-accounts balancing items.
- **Stochastic volatility (SV).** Each innovation log-variance h_{t,k}
  follows its own stationary AR(1), as in m9.

With psi = 0 the model is m9. The estimated psi are negative,
(-0.47, -0.34, -0.16) for (E, I, P), with 90% intervals [-0.57, -0.38],
[-0.46, -0.22] and [-0.36, 0.04]: a measurement error is partly reversed the
next quarter.

**How it is sampled.**
- Given h, the model is linear-Gaussian in mu. Quasi-differencing
  (y_t - psi·y_{t-1}) removes the error autocorrelation, so mu is integrated
  out exactly by a Kalman filter on the state (mu_t, mu_{t-1}) in the model
  block, and HMC only explores the hyperparameters and h.
- mu is drawn exactly from p(mu | y, h, theta) by forward-filtering
  backward-sampling in generated quantities.
- rho_EI is sampled as atanh(rho_EI), with the Jacobian added explicitly.
- The filter and FFBS are checked against the dense multivariate-normal
  likelihood and smoother in `tests/test_m12.R`.

The model and priors are in [priors.md](priors.md): m9's priors plus
psi_k ~ N(0, 0.5) on (-1, 1), as in m5.

**Sampler settings and diagnostics.**
- Settings: 4 chains × (2000 warmup + 2000 draws), adapt_delta 0.99,
  max_treedepth 12, seed 1.
- Diagnostics on the current data (production fit, seed 1):
  - 3 divergences out of 8000 draws;
  - per-chain E-BFMI 0.86-0.99;
  - max R-hat 1.002 on hyperparameters;
  - minimum bulk / tail ESS 2118 / 3520 on hyperparameters.
- A seed-2 fit also passes the gate (4 divergences, E-BFMI 0.84, R-hat 1.002,
  bulk ESS 1911). Both seeds' gate records are in
  `output/tables/m12_gate_seeds.csv`.
- A fit takes about 2.5 times as long as m9's. The production fit took
  1260 s (21 minutes) on 4 local cores, with chains of 952-1260 s. On a
  heavily loaded machine (load average 17-24 on 6 cores) the same fit took
  1544 s (26 minutes). With only 2 parallel chains it takes 32-40 minutes.

**Why m12.** m12 replaced m9 under a rule fixed before m12's results were
known. It had to meet all four conditions:
1. pass the production gate (below) on two seeds: it did (seeds 1 and 2);
2. have a reliable LFO score: the particle filter's Monte Carlo error on the
   m12 - m9 difference is 0.37 nats and its measured bias about 0.1;
3. beat m9 on PSIS-LFO-CV by at least 2 Newey-West SEs: +16.77 nats,
   SE 6.12, z = 2.74;
4. fit locally in at most 30 minutes: it does with 4 parallel chains (21
   minutes; 26 on a heavily loaded machine), so the margin is thin.

The gain over m9 is broad rather than driven by a few quarters: m12 is ahead
in 99 of 146 quarters and by +19.9 (Newey-West SE 3.8) before 2020. It is
behind in 2020, by 4.2 nats over the four quarters (mostly 2020 Q2 and Q4).

The structural restriction is still not supported by prediction: m9 and m4
(SV without it) are tied, +1.5 nats (SE 1.5). It is kept because it reflects
how the ABS builds the series. Earlier versions of this README said m9 "won
the LOO comparison". That claim was withdrawn because the LOO computation was
invalid (see NOTES.md).

## Model comparison

`R/03_compare.R` scores every model on PSIS-LFO-CV (Bürkner, Gabry & Vehtari
2020) of the one-step-ahead predictive density log p(y_t | y_{1:t-1}).
- The latent states are integrated out: an exact Kalman filter for the
  linear-Gaussian specs, and a Rao-Blackwellised particle filter for the SV
  specs m4, m9, m10 and m12. The particle count starts at 2000 per quarter
  and rises to up to 100000 where the effective sample size collapses (the
  COVID quarters). Each SV score is the mean of 4 independent filter runs.
- The window is 1990 Q1 - 2026 Q2: 146 predictions, the same for every model.
- Every stored fit meets the refit rule: 75 of 75 (12 initial fits and 63
  refits). The rule is divergences <= 1%, R-hat <= 1.05, E-BFMI >= 0.2 and
  bulk ESS >= 100. Three fits met it only after retries with new seeds (m9
  on data to 2006 Q3 and to 2020 Q4, m12 to 2020 Q4; 6 failed attempts,
  listed in `output/tables/lfo_refit_attempts.csv`). Keeping the seed that
  passes is seed selection, so the targets those fits serve were rescored
  without it: m12 - m9 is +16.09 (z 2.47) with each fit's first attempt and
  +16.62 (z 2.64) with all attempts pooled, against +16.77 (z 2.74) as
  reported. Ten fits that meet the rule still have bulk ESS below 400
  (five on tau, in the SV fits ending in 2020 Q2 or Q4; the other five on
  phi, omega or m3's Omega); see NOTES.md.

Results on the 3-measure models (E, I, P, full sample):

| model | elpd_lfo | diff vs m12 | SE (iid) | SE (Newey-West) | MCSE |
|---|---|---|---|---|---|
| m12 structural Sigma + SV + AR(1) errors | -1105.5 | 0 | - | - | 0.20 |
| m9 structural Sigma + SV | -1122.2 | -16.8 | 6.7 | 6.1 | 0.31 |
| m4 SV | -1123.7 | -18.3 | 6.9 | 7.5 | 0.02 |
| m5 serially correlated errors | -1168.8 | -63.4 | 14.3 | 20.7 | exact |
| m7 COVID variance scale | -1177.5 | -72.1 | 31.0 | 34.4 | exact |
| m6, m1, m3, m2 | -1195.6 to -1196.2 | about -90 | about 15 | about 21 | exact |

The iid SE, sqrt(n var(d)), treats the 146 quarters as independent. The
pointwise differences are autocorrelated, so the table also gives a
Newey-West SE with an 8-quarter lag. Two features each help, and their gains
add up:
- **Stochastic volatility.** m9 - m5 = +46.6 (SE 12.5 iid, 19.7 Newey-West).
- **Serially correlated errors.** m5 - m2 = +27.4 (SE 2.9 iid, 2.0
  Newey-West), and m12 - m9 = +16.8 on top of SV (SE 6.7 iid, 6.1
  Newey-West, z 2.7).
- **Both.** m12 - m5 = +63.4 (Newey-West SE 20.7, z 3.1).

Further results:
- **Hours.** m10 was re-specified after its old version's results were seen,
  so its gains below are post hoc (a second specification scored on the
  same window) and its z-values overstate the evidence. On the hours
  subsample, the re-specified m10 (4 signals, SV on
  the measurement errors *and* on latent growth) is ahead of m8 (4 signals,
  constant variances) by 124.2 (SE 53.5 iid, 64.7 Newey-West, z 1.9). On the
  GDP-only target, neither hours nor monthly employment does detectably
  better than m12:
  - m10 - m12: +47.0 (SE 28.5 iid, 30.1 Newey-West). Against m9 the gap is
    larger, +63.8 (Newey-West SE 25.6, z 2.5), so SV on latent growth may be
    worth adding to the 3-measure models;
  - m11 - m12: -7.9 (SE 18.7 iid, 27.0 Newey-West). m11 also sees employment
    inside the predicted quarter.
- **Particle-filter accuracy.** Against high-particle reference runs, the
  filter's bias summed over the COVID quarters is +0.16 (m12), +0.02 (m9),
  -0.27 (m4) and -0.15 (m10; -0.01 on its GDP-only target), and under 0.02
  in absolute value over 2015-2019 (checked for m12 and m9). None of these
  changes a comparison by more than 0.5 nats. The Monte Carlo error comes
  from only 4 filter runs per model; see NOTES.md.
- **What it does not test.** LFO scores predictions of the observed measures,
  not the split between latent growth and measurement error. Models that
  predict GDP about equally well can imply very different splits: over the
  same 1978-2026 quarters m10's latent path has 75% of headline variance and
  m12's 91% (56% and 88% excluding 2020 Q1 - 2021 Q1). The "about a fifth is
  measurement error" figure is a property of m12, not something the
  comparison confirms.

## Pipeline

Production (run by CI):

```
R/01_fetch_data.R      ABS 5206.0 Table 24 (GDP E/I/P) -> data/processed/gdp_growth.rds
R/fit_production.R     fits m12, saves output/fit_production.rds, applies the quality gate
R/output_estimate.R    writes docs/estimates.csv, estimates_latest.json, diagnostics.json
R/build_site.R         renders docs/index.html and charts
R/input_hash.R         prints the SHA-256 fingerprint of the inputs (used by CI)
```

**Quality gate.** `R/fit_production.R` stops with an error, so CI publishes
nothing, unless all of these hold:
- divergent transitions are under 0.1% of post-warmup draws;
- the minimum E-BFMI over chains is at least 0.3;
- R-hat is at most 1.01 on every hyperparameter and at most 1.05 on every
  mu[t] and h[t,k];
- bulk and tail ESS are at least 400 on every hyperparameter and every mu[t].

The hyperparameters are c, rho, tau, mu_h, phi, omega, rho_EI, psi and mu_bar.
Chain and draw counts are read from the fit. Any NA diagnostic counts as a
failure.

Research suite (not run by CI):

```
R/01b_fetch_research_data.R  5206.0 Table 1 hours worked + 6202.0 monthly employment
                             (run after 01_fetch_data.R; needed by m8, m10, m11)
R/02_fit_models.R            fits m1-m10 and m12 (Rscript R/02_fit_models.R m9 m12 fits a subset)
                             -> output/fit_m*.rds, output/tables/fit_diagnostics.csv
R/fit_m11.R                  fits m11 (mixed frequency)
R/fit_m12.R                  fits m12 with the production gate (any seed); warmup-stall harness
R/03_compare.R               PSIS-LFO-CV comparison + variance decomposition
                             -> output/tables/model_comparison*.csv, lfo_*.csv,
                                variance_decomposition.csv; checkpoints in output/lfo/
R/04_figures.R               output/figures/*.png
tests/test_filters.R         checks the Kalman / particle filters used by 03_compare.R
tests/test_m12.R             checks m12's Stan filter, FFBS and particle filter (R/lib/filter_m12.*)
tests/test_m3_m10.R          checks the m3 and m10 Stan filters and m10's particle filter (R/lib/filter_m10.*)
R/diag_m10_modes.R           m10 mode check -> output/tables/m10_mode_check.csv
R/diag_m3_ident.R            m3 identification checks -> output/tables/m3_identification.csv
```

## Automation

The GitHub Actions workflow is `.github/workflows/recompute.yml`.

- **Schedule:** Wednesdays at 02:00 UTC. This covers the usual
  first-Wednesday release window for 5206.0.
- **Input check:** each run fetches the data, then `R/input_hash.R`
  fingerprints two things: the GDP growth data, and the production code (the
  m12 Stan file, the four production scripts, `R/utils_abs.R`,
  `R/input_hash.R` and the workflow file). If the hash matches the
  `input_hash` in the committed `docs/diagnostics.json`, every later step is
  skipped and the job still succeeds. Later steps are the cmdstanr and CmdStan
  install, the fit, the site build, the commit and the Pages deploy.
- **Manual trigger:** from the Actions tab ("Recompute latent GDP estimate" ->
  Run workflow). Tick **force** to refit and redeploy even when the inputs are
  unchanged.
- **Side effects:** when the inputs changed, the workflow commits the updated
  `docs/` to `master` and redeploys GitHub Pages.
- **Quality gate:** the job fails, and nothing is published, if the m12 fit
  fails any check listed under Pipeline.
- **Caching:** the CmdStan install (CmdStan 2.39.0) and the compiled m12
  executable are cached. The cache keys are the CmdStan version, the runner
  image and the hash of `stan/m12_struct_sv_serial.stan`. The job timeout is
  120 minutes (raised from 90 for m12). The fit's run time on GitHub's
  runners has not been measured; see "Production model" for the local time.

## Local reproduction

```bash
# One-time
Rscript -e 'install.packages(c("readabs","posterior","loo","jsonlite","dplyr","tidyr","ggplot2","openssl","Rcpp","testthat","mvtnorm")); install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev","https://cloud.r-project.org")); cmdstanr::install_cmdstan()'

# Production estimate
Rscript R/01_fetch_data.R
Rscript R/fit_production.R
Rscript R/output_estimate.R
Rscript R/build_site.R

# Research suite (about 1 h for 02; several hours for 03 on 4-6 cores)
Rscript R/01b_fetch_research_data.R
Rscript R/02_fit_models.R
Rscript R/fit_m11.R
Rscript tests/test_filters.R
Rscript tests/test_m12.R
Rscript tests/test_m3_m10.R
Rscript R/03_compare.R
Rscript R/04_figures.R
```

Two environment variables redirect the production scripts' input and output,
for smoke tests that must not touch `docs/`:
- `LATENTGDP_FIT` (default `output/fit_production.rds`), read by
  `output_estimate.R`;
- `LATENTGDP_DOCS` (default `docs`), used by `output_estimate.R` and
  `build_site.R`.

## References

- Rees, Lancaster & Finlay (2014), "Estimating the Australian Output Gap
  and Other State-Space Refinements to GDP Measurement", RBA RDP 2014-12.
  https://www.rba.gov.au/publications/rdp/2014/2014-12.html
- Aruoba, Diebold, Nalewaik, Schorfheide & Song (2016), "Improving GDP
  Measurement: A Measurement-Error Perspective", J. Econometrics 191(2).
- Bürkner, Gabry & Vehtari (2020), "Approximate leave-future-out
  cross-validation for Bayesian time series models", J. Statistical
  Computation and Simulation 90(14).
- ABS Cat. 5206.0, Australian National Accounts: National Income,
  Expenditure and Product.
