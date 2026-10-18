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
- **Production model (m13).** On the current vintage (1959 Q4 - 2026 Q2,
  T = 267), the variance of the latent path is **79% of the variance of
  headline growth**. That is the posterior median of var(mu draw) / var(headline),
  with a 90% interval of 71% to 89%.
- **Other models.** The other full-sample specifications give 78-83% (m12,
  production from 2026-10-17 to 2026-10-18, 78%; m9, production before that,
  79%). The exception is m7, the COVID-scale model, at 61%. See NOTES.md for
  the full table.

So, according to m13, roughly a fifth of the quarter-to-quarter variance of
headline growth is measurement error rather than movement in true activity.
That share depends on the period: excluding the COVID quarters (2020 Q1 -
2021 Q1) it is about a quarter (73% of headline variance on the full sample,
75% on 1978 Q4 - 2026 Q2, against 88% for m12 on those quarters). The model
comparison does not test this split (see "What it does not test" below).
Earlier versions of
this project claimed a 30-35% reduction. That figure is too large: it was the
variance of the posterior-*mean* path, which is smoother than any single
plausible path.

## Production model

The production model is `stan/m13_serial_latent_sv.stan` (since 2026-10-18;
from 2026-10-17 it was m12, `stan/m12_struct_sv_serial.stan`, and before that
m9, `stan/m9_struct_sv.stan`). It is m12 with stochastic volatility added to
latent growth. It has six parts:
- **Latent state.** Latent growth follows a stationary AR(1) with a
  time-varying innovation sd, mu_t = c + rho·mu_{t-1} + tau_t·eta_t, with
  tau_t = exp(g_t / 2). mu_1 ~ N(c/(1-rho), tau_1²/(1-rho²)).
- **Latent stochastic volatility.** The latent log-variance g_t follows its
  own stationary AR(1), g_t = mu_g + phi_g (g_{t-1} - mu_g) + omega_g zeta_t,
  as in the re-specified m10.
- **Observations.** Each measure equals mu_t plus a measurement error.
- **Serially correlated errors.** Each error follows its own AR(1),
  eps_t^k = psi_k·eps_{t-1}^k + u_t^k, as in m5 and m12. The period-1 errors
  start from the stationary distribution at the period-1 variances.
- **Structural Sigma.** The E and I innovations have a constant correlation
  rho_EI. P's is independent of both. The ABS builds Production from separate
  industry surveys, while E and I share national-accounts balancing items.
- **Measurement stochastic volatility.** Each innovation log-variance h_{t,k}
  follows its own stationary AR(1), as in m9 and m12.

With g_t held constant (omega_g = 0) the model is m12; with psi = 0 as well
it is m9. Posterior (seed 1; medians and 90% intervals; annualised units):
- **Latent volatility.** tau_t falls from about 4.9 in the 1960s to 1.6 in
  the 2000s and 1.4 in the 2010s (posterior median of each decade's average),
  and is 12.0 [8.3, 19.1] in 2020 Q2. m12's constant tau is 3.99. The
  long-run value exp(mu_g / 2) is 2.47 [1.71, 3.61]; phi_g 0.91, omega_g
  0.61 [0.42, 0.89].
- **Measurement errors.** psi = (-0.46, -0.36, -0.19) for (E, I, P), with
  90% intervals [-0.56, -0.37], [-0.47, -0.23] and [-0.35, -0.02]: a
  measurement error is partly reversed the next quarter. rho_EI is 0.12
  [-0.03, 0.27].
- **Latent path.** The median path differs from m12's by 0.38 pp
  (annualised) on average; 2020 Q2 is -27.9 (m12 -27.4, headline -28.1).

**How it is sampled.**
- Given g and h, the model is linear-Gaussian in mu. Quasi-differencing
  (y_t - psi·y_{t-1}) removes the error autocorrelation, so mu is integrated
  out exactly by a Kalman filter on the state (mu_t, mu_{t-1}) in the model
  block (m12's filter with tau replaced by tau_t), and HMC only explores the
  hyperparameters, g and h. g and h are non-centred.
- mu is drawn exactly from p(mu | y, g, h, theta) by forward-filtering
  backward-sampling in generated quantities.
- rho_EI is sampled as atanh(rho_EI), with the Jacobian added explicitly.
- The filter and FFBS are checked against the dense multivariate-normal
  likelihood and smoother in `tests/test_m13.R` (133 checks). With tau_t
  constant the filter is bit-identical to m12's.

The model and priors are in [priors.md](priors.md): m12's priors, with
tau ~ half-N(0, 2) replaced by m10's latent SV priors mu_g ~ N(2, 1),
phi_g ~ U(-1, 1), omega_g ~ half-N(0, 0.5). Widening the omega_g prior to
half-N(0, 1) changes omega_g from 0.61 to 0.66, the variance ratio from
0.790 to 0.786 and the latent path by 0.03 pp on average (NOTES.md).

**Sampler settings and diagnostics.**
- Settings: 4 chains × (2000 warmup + 2000 draws), adapt_delta 0.99,
  max_treedepth 12, seed 1.
- Diagnostics on the current data (production fit, seed 1):
  - 4 divergences out of 8000 draws;
  - per-chain E-BFMI 0.85-0.97;
  - max R-hat 1.002 on hyperparameters, 1.001 on mu, 1.002 on h and g;
  - minimum bulk / tail ESS 1553 / 3213 on hyperparameters and 4344 / 4866
    on mu.
- A seed-2 fit also passes the gate (0 divergences, E-BFMI 0.84, max R-hat
  1.003 on hyperparameters, bulk ESS 1784). Both seeds' gate records are in
  `output/tables/m13_gate_seeds.csv`.
- A fit takes about as long as m12's: 1384 s (23 minutes) wall on 4 local
  cores for the seed-1 research fit, chains 989-1375 s; seed 2 took 1664 s
  (28 minutes) at load average 10-11. The production run on 2026-10-18 (load average about 4)
  took 1173 s of sampling (chains 856-1173 s), about 20 minutes wall; it
  reproduced the seed-1 research fit draw for draw.

**Why m13.** m13 replaced m12 under a rule fixed before m13's LFO results
were known. The rule added one condition to the one used for m12: the z
condition also had to hold with each LFO refit's first sampling attempt
only, so that no retry with a new seed could create the margin. That
condition guards only against seed selection; it does not correct for m13
having been proposed after the re-specified m10's LFO results were seen
(see the end of this section). m13 had to meet all five conditions:
1. pass the production gate (below) on two seeds: it did (seeds 1 and 2);
2. have a reliable LFO score: the particle filter's Monte Carlo error on the
   m13 - m12 difference is 0.29 nats and its measured bias about -0.1; all
   6 of m13's LFO fits met the refit rule on the first attempt;
3. beat m12 on PSIS-LFO-CV by at least 2 Newey-West SEs: +59.65 nats,
   SE 28.89, z = 2.06;
4. the same with first attempts only: +59.79, SE 29.00, z = 2.06;
5. fit locally in at most 30 minutes: 23 minutes (28 under load).

The z and fit-time margins are thin. One quarter dominates the SE: 2020 Q2
alone is worth +23.7 nats. Without the COVID quarters the gain is +28.5
(Newey-West SE 6.0, z 4.7); before 2020 it is +20.1 (SE 4.9, z 4.1), and m13
is ahead in 122 of 146 quarters. From a 1980 Q1 start the pre-2020 gain is
only +6.4 (Newey-West SE 2.75, z 2.3), and m13 is ahead in 107 quarters. So
much of this pre-2020 support also comes from the pre-1980 training data
(see "1980 Q1 sample" below). The comparison is post hoc in the sense
above: m13 is a second latent-SV specification scored on the same window
after the first (m10) was seen, so its z overstates the evidence somewhat.
Nothing in the rule adjusts for that: a Bonferroni bar for two tries
(z about 2.24) would not be met.

**Why m12 (superseded 2026-10-18).** m12 replaced m9 on 2026-10-17 under a
rule fixed before m12's results were known: (1) the production gate on two
seeds, (2) a reliable LFO score (MCSE 0.37 on m12 - m9, bias about 0.1),
(3) beat m9 by z_NW >= 2 (+16.77 nats, SE 6.12, z = 2.74), (4) a local fit of
at most 30 minutes (21 minutes; 26 under load). Its gain over m9 is broad
(ahead in 99 of 146 quarters, +19.9 before 2020) and it is behind m9 in 2020
by 4.2 nats.

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
  specs m4, m9, m10, m12 and m13. The particle count starts at 2000 per quarter
  and rises to up to 100000 where the effective sample size collapses (the
  COVID quarters). Each SV score is the mean of 4 independent filter runs.
- The window is 1990 Q1 - 2026 Q2: 146 predictions, the same for every model.
- Every stored fit meets the refit rule: 81 of 81 (13 initial fits and 68
  refits). The rule is divergences <= 1%, R-hat <= 1.05, E-BFMI >= 0.2 and
  bulk ESS >= 100. Three fits met it only after retries with new seeds (m9
  on data to 2006 Q3 and to 2020 Q4, m12 to 2020 Q4; 6 failed attempts,
  listed in `output/tables/lfo_refit_attempts.csv`). Keeping the seed that
  passes is seed selection, so the targets those fits serve were rescored
  without it: m12 - m9 is +16.09 (z 2.47) with each fit's first attempt and
  +16.62 (z 2.64) with all attempts pooled, against +16.77 (z 2.74) as
  reported. m13 needed no retry; m13 - m12 is +59.79 (z 2.06) with first
  attempts only and +59.60 (z 2.06) pooled, against +59.65 (z 2.06) as
  reported. Twelve fits that meet the rule still have bulk ESS below 400
  (six of them on data ending in 2020 Q2 or Q4; two are m13's, on data to
  2002 Q2 and 2010 Q3). The diagnostics table records each fit's lowest bulk
  ESS but not which variable has it; see NOTES.md.

Results on the 3-measure models (E, I, P, full sample):

| model | elpd_lfo | diff vs m13 | SE (iid) | SE (Newey-West) | MCSE |
|---|---|---|---|---|---|
| m13 m12 + latent SV | -1045.8 | 0 | - | - | 0.22 |
| m12 structural Sigma + SV + AR(1) errors | -1105.5 | -59.6 | 24.7 | 28.9 | 0.20 |
| m9 structural Sigma + SV | -1122.2 | -76.4 | 21.8 | 24.8 | 0.31 |
| m4 SV | -1123.7 | -77.9 | 21.9 | 23.7 | 0.02 |
| m5 serially correlated errors | -1168.8 | -123.0 | 15.3 | 25.7 | exact |
| m7 COVID variance scale | -1177.5 | -131.7 | 12.2 | 20.7 | exact |
| m6, m1, m3, m2 | -1195.6 to -1196.2 | about -150 | about 15.5 | about 25.5 | exact |

The iid SE, sqrt(n var(d)), treats the 146 quarters as independent. The
pointwise differences are autocorrelated, so the table also gives a
Newey-West SE with an 8-quarter lag. Three features each help, each
measured as one step added to a model that has the ones before (only serial
correlation was measured both with and without SV, +27.4 and +16.8; whether
the three gains add up was not tested):
- **Stochastic volatility in the measurement errors.** m9 - m5 = +46.6
  (SE 12.5 iid, 19.7 Newey-West). m9 - m5 is not a nested step: m9 has SV
  but not m5's AR(1) errors. The comparisons that add only measurement SV
  are m4 - m2 = +72.5 (Newey-West SE 19.1, z 3.8) and m9 - m6 = +73.3
  (Newey-West SE 19.9, z 3.7). From a 1980 Q1 start they roughly halve but
  stay above z 2: +33.4 (z 2.1) and +38.9 (z 2.6). See "1980 Q1 sample"
  below.
- **Serially correlated errors.** m5 - m2 = +27.4 (SE 2.9 iid, 2.0
  Newey-West), and m12 - m9 = +16.8 on top of SV (SE 6.7 iid, 6.1
  Newey-West, z 2.7).
- **Stochastic volatility in latent growth.** m13 - m12 = +59.6 (SE 24.7
  iid, 28.9 Newey-West, z 2.06). 2020 Q2 alone contributes +23.7; without
  the COVID quarters the gain is +28.5 (SE 6.0 Newey-West). From a 1980 Q1
  start it is +57.4 but z 1.57.
- **All three.** m13 - m5 = +123.0 (Newey-West SE 25.7, z 4.8).

Further results:
- **Hours.** m10 was re-specified after its old version's results were seen,
  so its gains below are post hoc (a second specification scored on the
  same window) and its z-values overstate the evidence. On the hours
  subsample, the re-specified m10 (4 signals, SV on
  the measurement errors *and* on latent growth) is ahead of m8 (4 signals,
  constant variances) by 124.2 (SE 53.5 iid, 64.7 Newey-West, z 1.9). On the
  GDP-only target, neither the hours model nor the employment model beats
  m13, but these models differ from m13 in more than the labour data (m10
  has no AR(1) measurement errors or structural E/I-only correlation; m11
  has no SV), and no m13-plus-hours model has been fitted, so the gaps do
  not measure the value of hours or employment:
  - m13 - m10: +12.6 (SE 6.9 iid, 5.5 Newey-West; z 1.8 iid, 2.3
    Newey-West). m10 also sees past hours and starts in 1978;
  - m13 - m11: +67.5 (SE 12.7 iid, 18.4 Newey-West). m11 also sees
    employment inside the predicted quarter;
  - *superseded by m13:* m10 - m12 was +47.0 (SE 28.5 iid, 30.1 Newey-West)
    and m10 - m9 +63.8 (Newey-West SE 25.6, z 2.5), which suggested adding
    latent SV to the 3-measure models; m13 is that model.
- **Particle-filter accuracy.** Against high-particle reference runs, the
  filter's bias summed over the COVID quarters is +0.07 (m13), +0.16 (m12),
  +0.02 (m9), -0.27 (m4) and -0.15 (m10; -0.01 on its GDP-only target), and
  under 0.03 in absolute value over 2015-2019 (checked for m13, m12 and m9). None of these
  changes a comparison by more than 0.5 nats. The Monte Carlo error comes
  from only 4 filter runs per model; see NOTES.md.
- **What it does not test.** LFO scores predictions of the observed measures,
  not the split between latent growth and measurement error. Models that
  predict GDP about equally well can imply very different splits: over the
  same 1978-2026 quarters m10's latent path has 75% of headline variance,
  m13's 86% and m12's 91% (56%, 75% and 88% excluding 2020 Q1 - 2021 Q1).
  On the full sample m13 and m12 nearly agree (79% and 78%). The "about a
  fifth is measurement error" figure is a property of m13, not something the
  comparison confirms.

## 1980 Q1 sample (research only; production still uses the full sample)

RDP 2014-12 starts its data in 1980 Q1 (section 3.2). Its reason is that,
although quarterly national accounts go back to 1959 Q3, "the pattern of
measurement errors in the early years of each GDP series may be
unrepresentative of their current performance". The whole research suite was
re-run from that start on 2026-10-18/19:
- **Data.** The first growth observation is 1980 Q1 and the last 2026 Q2,
  T = 186 for every model.
- **LFO.** The LFO window is unchanged (1990 Q1 - 2026 Q2, 146 targets),
  with 40 quarters of initial training.
- **Outputs.** Everything is in `output/sample_1980/`.
- **Details.** Full results, reproduction commands and diagnostics are in
  NOTES.md, section "1980Q1 sample (RDP 2014-12 alignment)".

Only the start date follows the paper. **The RDP's fourth signal is the
unemployment rate; our m8 and m10 use hours worked.**

**Production is unchanged.** It still fits m13 on 1959 Q4 - 2026 Q2. On the
1980 sample m13 still passes the production gate on two seeds and is still
the best model on every LFO target. Its latent path barely moves: the
correlation with the full-sample path is 0.9985, and the paths differ by
0.13 pp on average. The comparison author's recommendation is to **keep
production on m13 with the full sample**, not to move it to 1980 Q1. This
is a recommendation; nothing was changed. For reference, m13 trained from
1980 predicts the same 146 targets 6.3 nats better than m13 trained from
1959 (NW z 1.69). That is below the z >= 2 bar used for every production
switch. Most of that gain is in 1990-99 (+4.5). From 2000 on it is +1.8
(z 0.7), so for current production the start date hardly matters for
prediction.

What the later start does to the full-sample conclusions (Newey-West z in
brackets):
- **Survives: serially correlated measurement errors.** m5 - m2 is +20.4
  (9.4) and m12 - m9 +20.6 (3.7).
- **Survives: m13 is first** in group 1 and on the GDP-only target.
- **Survives: the structural restriction is still unsupported.**
- **Survives: the hours model's lead over m8.** m10 - m8 is +119.6 (1.8).
- **Weakens: latent SV.** m13 - m12 is +57.4, but z falls to 1.57, and it is
  only +6.4 (2.3) before 2020. Excluding all of 2020-22 it is +12.7 (2.6),
  against +27.1 (4.6) on the full sample. The z >= 2 condition used to
  switch production to m13 would not be met on this sample.
- **Weakens: SV in the measurement errors.** The gain roughly halves but is
  still supported. The nested comparisons, which add only measurement SV,
  are:
  - m4 - m2: +33.4 (2.1), against +72.5 (3.8) on the full sample;
  - m9 - m6: +38.9 (2.6), against +73.3 (3.7).

  Before 2020 they are +19.4 (3.0) and +21.4 (3.2). About half of the
  full-sample gain came from the constant-variance models' error variances
  being inflated by the noisy pre-1980 training data. An earlier version
  said this gain "largely disappears". That used m9 - m5 (+16.7, 1.13) and
  m4 - m1 (+28.2, 1.46), and neither is nested: m9 lacks m5's AR(1) errors,
  and m4 adds a latent AR(1) to m1.
- **Fails: m4, m9 and m12 on COVID.** These three (measurement SV with a
  constant latent sd) treat the 2020 Q2 collapse as common measurement
  error on this sample: latent growth about 0 instead of -28, and a
  variance ratio of 0.32-0.37. That is a specification failure, not a
  finding. m13 and m10, which have latent SV, keep COVID in latent growth.
  So, mostly, do the constant-variance models (mu about -24 in 2020 Q2).
- **Changes: m13's variance ratio** is 0.87 [0.81, 0.94] on 1980-2026
  (0.78 excluding 2020 Q1 - 2021 Q1), against 0.79 on the full sample.
  This is entirely because the quarters covered differ. Over exactly
  1980 Q1 - 2026 Q2 the full-sample fit gives 0.88 (0.79).

## Pipeline

Production (run by CI):

```
R/01_fetch_data.R      ABS 5206.0 Table 24 (GDP E/I/P) -> data/processed/gdp_growth.rds
R/fit_production.R     fits m13, saves output/fit_production.rds, applies the quality gate
R/output_estimate.R    writes docs/estimates.csv, estimates_latest.json, diagnostics.json, model_summary.json
R/build_site.R         renders docs/index.html and charts
R/input_hash.R         prints the SHA-256 fingerprint of the inputs (used by CI)
```

**Quality gate.** `R/fit_production.R` stops with an error, so CI publishes
nothing, unless all of these hold:
- divergent transitions are under 0.1% of post-warmup draws;
- the minimum E-BFMI over chains is at least 0.3;
- R-hat is at most 1.01 on every hyperparameter and at most 1.05 on every
  mu[t], h[t,k] and g[t];
- bulk and tail ESS are at least 400 on every hyperparameter and every mu[t].

The hyperparameters are c, rho, tau, mu_h, phi, omega, rho_EI, psi, mu_bar,
mu_g, phi_g and omega_g (tau = exp(mu_g / 2) in m13).
Chain and draw counts are read from the fit. Any NA diagnostic counts as a
failure.

Research suite (not run by CI):

```
R/01b_fetch_research_data.R  5206.0 Table 1 hours worked + 6202.0 monthly employment
                             (run after 01_fetch_data.R; needed by m8, m10, m11)
R/02_fit_models.R            fits m1-m10, m12 and m13 (Rscript R/02_fit_models.R m9 m12 fits a subset)
                             -> output/fit_m*.rds, output/tables/fit_diagnostics.csv
R/fit_m11.R                  fits m11 (mixed frequency)
R/fit_m12.R                  fits m12 with the production gate (any seed); warmup-stall harness
R/fit_m13.R                  the same for m13 (also --omega_g_sd=S for the prior-sensitivity fit)
R/03_compare.R               PSIS-LFO-CV comparison + variance decomposition
                             -> output/tables/model_comparison*.csv, lfo_*.csv,
                                variance_decomposition.csv; checkpoints in output/lfo/
R/04_figures.R               output/figures/*.png
R/sample_window.R            optional sample start / output root for 02, fit_m11-13, 03, 04
                             (LATENTGDP_START, LATENTGDP_OUT; unset = full sample, unchanged)
tests/test_filters.R         checks the Kalman / particle filters used by 03_compare.R
tests/test_m12.R             checks m12's Stan filter, FFBS and particle filter (R/lib/filter_m12.*)
tests/test_m13.R             checks m13's Stan filter, FFBS and particle filter (R/lib/filter_m13.*)
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
  m13 Stan file, the four production scripts, `R/utils_abs.R`,
  `R/input_hash.R` and the workflow file). If the hash matches the
  `input_hash` in the committed `docs/diagnostics.json`, every later step is
  skipped and the job still succeeds. Later steps are the cmdstanr and CmdStan
  install, the fit, the site build, the commit and the Pages deploy.
- **Manual trigger:** from the Actions tab ("Recompute latent GDP estimate" ->
  Run workflow). Tick **force** to refit and redeploy even when the inputs are
  unchanged.
- **Side effects:** when the inputs changed, the workflow commits the updated
  `docs/` to `master` and redeploys GitHub Pages.
- **Quality gate:** the job fails, and nothing is published, if the m13 fit
  fails any check listed under Pipeline.
- **Caching:** the CmdStan install (CmdStan 2.39.0) and the compiled m13
  executable are cached. The cache keys are the CmdStan version, the runner
  image and the hash of `stan/m13_serial_latent_sv.stan`. The job timeout is
  120 minutes (raised from 90 for m12, unchanged for m13). The fit's run time
  on GitHub's runners has not been measured for m12 or m13; see "Production
  model" for the local time.

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
Rscript tests/test_m13.R
Rscript tests/test_m3_m10.R
Rscript R/03_compare.R
Rscript R/04_figures.R
```

To run the research suite on the 1980 Q1 sample instead (outputs in
`output/sample_1980/`, nothing in `output/` touched), use the drivers. Each
is sequential and resumable:

```bash
nohup bash output/sample_1980/run_1980.sh >> output/sample_1980/logs/driver.log 2>&1 &
nohup bash output/sample_1980/lfo/run_lfo_1980.sh >> output/sample_1980/lfo/driver.log 2>&1 &
```

Alternatively, set `LATENTGDP_START=1980-03-01 LATENTGDP_OUT=output/sample_1980`
before the research commands above; NOTES.md gives the step-by-step list.

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
