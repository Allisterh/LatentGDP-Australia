# NOTES

This is the research record. Sections dated 2026-10 replace earlier results
that turned out to be wrong; the superseded claims are noted where they were
made so the history stays readable. The newest round is "Round of
2026-10-17" below: it adds m12, which is now the production model, and
supersedes the 2026-10-16 comparison tables and the m3 and m10 sections.

## Paper

RBA RDP 2014-12, "Estimating the Australian Output Gap and Other State-Space
Refinements to GDP Measurement", by Daniel Rees, David Lancaster & Richard Finlay.
https://www.rba.gov.au/publications/rdp/2014/2014-12.html

## Data (refreshed 2026-10-16, ABS release for 2026 Q2)

ABS Cat. 5206.0, Table 24 (Selected Analytical Series): GDP E-, I- and P-based,
chain volume measures, seasonally adjusted (`R/01_fetch_data.R`). The levels
run from 1959 Q3 to 2026 Q2, giving T = 267 quarters of 400 x dlog growth
rates (1959 Q4 to 2026 Q2).

Research-only series (`R/01b_fetch_research_data.R`):
- **Hours worked** (5206.0 Table 1). Joined to the GDP growth rates, this gives
  1978 Q4 to 2026 Q2, T = 191 (m8, m10).
- **Employed persons** (6202.0 Table 1, SA, monthly). Growth runs from 1978-03
  to 2026-09 (583 months). m11 uses months up to the last GDP quarter (T_m = 580)
  and the T_q = 193 quarters whose three months are all present.

Sample statistics (annualised %, T = 267):
- means: g_E 3.25, g_I 3.19, g_P 3.26
- sd: g_E 6.0, g_I 5.1, g_P 4.7
- pairwise correlations 0.58-0.62

The 2020 Q2 and Q3 COVID quarters (about -28% and +13% annualised) are extreme
outliers. In 2020 Q2 the three measures agree to within 0.5 pp (-28.0, -28.4,
-27.9).

The previous run (2026-05-21) used 1959 Q4 - 2025 Q4, T = 265, with the hours
sample described as "post-1979, T = 189". All numbers below are from the
current vintage.

## Fits (current data, `R/02_fit_models.R` and `R/fit_m11.R`)

Settings for every fit: 4 chains x (2000 warmup + 2000 draws), seed 1. The
table shows posterior means. The diagnostics for m1-m10 and m12 come from
`output/tables/fit_diagnostics.csv`; m11's come from the `R/fit_m11.R` run
log, because that script does not write to the table:
- **div** is divergent transitions out of 8000 draws.
- **E-BFMI** is the minimum over chains.
- **R-hat** is the maximum over the model's hyperparameters.

| Spec | T | mu_bar | rho | tau | sig_E | sig_I | sig_P | div | E-BFMI | R-hat | Notes |
|------|---|--------|-----|-----|-------|-------|-------|-----|--------|-------|-------|
| m1  | 267 | 3.23 | -    | 4.01 | 4.08 | 3.25 | 2.65 | 0 | 0.73 | 1.001 | clean |
| m2  | 267 | 3.12 | 0.04 | 4.02 | 4.10 | 3.26 | 2.64 | 0 | 0.54 | 1.004 | rho indistinguishable from 0 |
| m3  | 267 | 3.12 | 0.06 | 4.04 | 4.35 | 3.13 | 2.43 | 0 | 0.89 | 1.007 | mu integrated out (2026-10-17); Omega weakly identified; min bulk / tail ESS 882 / 438 |
| m3 *(superseded, mu sampled)* | 267 | 3.12 | 0.06 | 4.06 | 4.35 | 3.11 | 2.43 | 0 | **0.03** | 1.029 | min bulk ESS 132 |
| m4  | 267 | 3.08 | 0.05 | 4.05 | -    | -    | -    | 5 | 0.89 | 1.002 | SV; clean |
| m5  | 267 | 3.13 | 0.08 | 4.00 | 3.53 | 2.98 | 2.52 | 0 | 0.50 | 1.003 | sig = innovation sd; clean |
| m6  | 267 | 3.13 | 0.05 | 3.99 | 4.20 | 3.36 | 2.56 | 0 | 0.48 | 1.002 | rho_EI = 0.09 |
| m7  | 267 | 3.26 | 0.10 | 3.54 | 4.14 | 3.26 | 2.64 | 0 | 0.54 | 1.002 | tau drops: COVID quarters absorbed by kappa |
| m8  | 191 | 2.86 | 0.12 | 3.62 | 2.86 | 3.08 | 1.83 | 0 | 0.49 | 1.003 | +hours; 1978 Q4 onward |
| m9  | 267 | 3.08 | 0.06 | 4.01 | -    | -    | -    | 1 | 0.90 | 1.002 | production until 2026-10-17; rho_EI = 0.14 |
| m10 | 191 | 2.76 | 0.44 | 1.59 | -    | -    | -    | 4 | 0.86 | 1.002 | re-specified 2026-10-17: latent SV; tau = exp(mu_g/2) |
| m10 *(superseded, constant tau)* | 191 | 2.98 | 0.39 | 2.19 | -    | -    | -    | 1 | 0.84 | 1.011 | multimodal posterior; min tail ESS 88 |
| m12 | 267 | 3.10 | 0.09 | 3.99 | -    | -    | -    | 3 | 0.86 | 1.002 | production from 2026-10-17; m9 + AR(1) errors; psi = (-0.47, -0.34, -0.16) |
| m11 | 193 q / 580 m | - | 0.49 | 4.48 | 2.87 | 3.09 | 1.86 | 0 | 0.45 | 1.002 | monthly latent; sig = sigma_q |

**m11.** Its parameters are at monthly frequency and do not compare directly
with the other models.
- Long-run monthly latent: mu_bar_m = 2.65.
- rho = 0.49 per month, which is about 0.12 at quarterly aggregation
  (0.49^3), consistent with m8.
- Monthly employment loads on the latent with beta_emp = 0.61.

**m4 and m9.** Their measurement-error sd varies over time. The long-run sd
exp(mu_h / 2), using the posterior mean of mu_h, is:

| | E | I | P |
|---|---|---|---|
| m4 | 2.97 | 2.31 | 1.94 |
| m9 | 3.05 | 2.41 | 1.95 |

The log-variances are highly persistent: m4 phi = (0.98, 0.94, 0.96) and
m9 phi = (0.98, 0.94, 0.97).

**m5.** The AR(1) measurement-error coefficients are negative:
psi_E = -0.51, psi_I = -0.43, psi_P = -0.18, with 90% intervals all below 0.
Negative serial correlation means the measurement error mean-reverts. That is
consistent with statistical-agency overshoot and undershoot, and with
revisions.

**Production is the cleanest signal.** Production has the smallest
measurement-error sd in every spec, so production-side GDP is the cleanest
single signal of true growth. Expenditure has the largest sd on the full
sample. On the post-1978 samples (m8, m11), Income's is slightly larger than
Expenditure's.

## Round of 2026-10-17: m12, m3 rewrite, m10 re-specification, stricter LFO, production switch

This section supersedes the 2026-10-16 comparison tables below (they are kept,
marked superseded, so the history stays readable). The data vintage is
unchanged (2026 Q2, T = 267).

### Summary

- **m12** (m9 + AR(1) measurement errors) is new. It wins group 1 on
  PSIS-LFO-CV: -1105.47 against -1122.24 for m9 (+16.77, Newey-West SE 6.12,
  z 2.74).
- **Production switched from m9 to m12**, under a rule fixed before m12's
  results were known (below).
- **m3** now samples cleanly (mu integrated out). Its Omega is still weakly
  identified; that is a property of the model.
- **m10** had a genuine two-mode posterior. It was re-specified with
  stochastic volatility on the latent innovation, which removes the second
  mode. Its old LFO score (-1507.1) is invalid. The re-specification was
  chosen after seeing the data, so the new m10's LFO gains are post hoc.
- **The LFO refit rule** now also checks E-BFMI and bulk ESS. Every one of
  the 75 stored fits (12 initial fits and 63 refits) meets it, but three
  of them met it only after 6 failed attempts with other seeds (see "LFO
  method changes").
- **The particle filter** now uses an adaptive particle count, four runs per
  model, and was checked against high-particle reference runs.

### m12: specification

`stan/m12_struct_sv_serial.stan`:
- y_t^k = mu_t + eps_t^k, k = E, I, P;
- eps_t^k = psi_k eps_{t-1}^k + u_t^k;
- u_t ~ N(0, R_t), R_t = D_t C D_t, D_t = diag(exp(h_t / 2)),
  C = [1 r 0; r 1 0; 0 0 1], r = rho_EI = tanh(z_EI): m9's covariance;
- h_{t,k}: stationary AR(1), non-centred, exactly as in m9;
- mu_t = c + rho mu_{t-1} + tau eta_t, mu_1 stationary.

sigma_t = exp(h / 2) is the sd of the innovation u_t, as m5's sigma is. With
psi = 0 the model is m9.

- **Initial errors.** eps_1 ~ N(0, V_1), V_1[j,k] = R_1[j,k] / (1 - psi_j psi_k):
  the stationary covariance if the variance were held at its period-1 level.
  This matches m5 and m9, reduces to m9 at psi = 0, and keeps log p(y_1) in
  the likelihood. The LFO filter uses the same V_1.
- **Priors.** m9's, unchanged, plus psi_k ~ N(0, 0.5) on (-1, 1) as in m5
  (`priors.md`).
- **Sampling.** Given h, the model block integrates mu out exactly on the
  quasi-differenced data y_t - psi∘y_{t-1}, which load on (mu_t, mu_{t-1}).
  The log-determinant and quadratic form are written as sums of squares
  (Cauchy-Binet), extending m9's stable formula to a 2-dimensional state.
  mu is drawn by FFBS; the backward step uses the joint filtered law of
  (mu_{t+1}, mu_t) given y_{1:t+1}, because y_{t+1} also depends on mu_t.
- **Tests** (`tests/test_m12.R`, 57 checks, all pass):
  - the Stan filter against the dense multivariate-normal log-density, to
    1e-9, for 5 psi settings, near rho_EI = -1 and at h = -15;
  - at psi = 0, equality with m9's Stan filter and m9's particle filter step;
  - FFBS against the dense smoother (20000 draws). Reusing m9's backward
    step here would give a worst error of 51.6 Monte Carlo SEs; the correct
    step gives 2.2;
  - the C++ filter (`R/lib/filter_m12.cpp`) against Stan, to 1e-11;
  - the particle filter against the exact Kalman filter at omega -> 0, against
    m9's filter at psi = 0, and against a brute-force average over 100000
    simulated h paths;
  - the C++ filter against the fit's `log_lik_pred` on 5 posterior draws.
- **Warmup stalls.** 0 of 100 short warmup chains stalled
  (`R/fit_m12.R stall`, Stan's default random inits). An init = 0.5 batch,
  as run for m9, has not been run.

### m12: fits and posterior

Both seeds use the production settings (4 x (2000 + 2000), adapt_delta 0.99,
max_treedepth 12) and pass the unchanged production gate, with psi added to
the hyperparameter list. The gate records of both seeds are in
`output/tables/m12_gate_seeds.csv`; the fits (`output/fit_m12.rds`,
`output/fit_m12_seed2.rds`, with `output/fit_m12_gate.rds` and
`output/fit_m12_seed2_gate.rds`) are gitignored and are regenerated by
`Rscript R/fit_m12.R` and `Rscript R/fit_m12.R --seed=2 --out=output/fit_m12_seed2.rds`:

| Check | seed 1 | seed 2 | Gate |
|---|---|---|---|
| divergences | 3 of 8000 (0.037%) | 4 of 8000 (0.050%) | < 0.1% |
| E-BFMI by chain | 0.862, 0.899, 0.923, 0.989 | 0.925, 0.899, 0.843, 0.923 | min >= 0.3 |
| max R-hat, hyperparameters | 1.0020 (omega[3]) | 1.0017 (omega[3]) | <= 1.01 |
| max R-hat, mu / h | 1.0015 / 1.0019 | 1.0012 / 1.0021 | <= 1.05 |
| min bulk / tail ESS, hyperparameters | 2118 (omega[2]) / 3520 | 1911 (omega[2]) / 3078 (phi[1]) | >= 400 |
| min bulk / tail ESS, mu | 2852 / 6476 | 3055 / 6061 | >= 400 |
| max-treedepth hits | 0 | 0 | - |

The two seeds agree on every summary to about 0.01. Posterior means, with
90% intervals where useful (seed 1, `output/fit_m12.rds`; m5 and m9 for
comparison):

| | m12 | m9 | m5 |
|---|---|---|---|
| psi (E, I, P) | -0.47, -0.34, -0.16 | - | -0.51, -0.43, -0.18 |
| psi 90% intervals | [-0.57, -0.38], [-0.46, -0.22], [-0.36, 0.04] | - | all below 0 |
| rho_EI | 0.175 [0.003, 0.323] | 0.138 [-0.02, 0.29] | - |
| phi | 0.976, 0.947, 0.975 | 0.981, 0.940, 0.970 | - |
| omega | 0.243, 0.369, 0.436 | 0.226, 0.449, 0.459 | - |
| rho, tau, mu_bar | 0.089, 3.99, 3.10 | 0.06, 4.01, 3.08 | 0.08, 4.00, 3.13 |
| variance ratio (median, 90%) | 0.779 [0.707, 0.850] | 0.785 [0.709, 0.866] | 0.783 [0.710, 0.857] |

- **The m5 finding survives SV.** The measurement errors still mean-revert
  once their variance is allowed to change, slightly less strongly than in m5.
  Only psi_P's interval includes 0.
- **Nothing else moves much.** The SV parameters, the latent AR(1) and the
  variance ratio are close to m9's.
- **Speed.** A fit is about 2.5 times slower per chain than m9's. Measured
  run times (CmdStan's total sampling time, `fit$time()$total`):
  - production fit, 4 parallel chains, lightly loaded machine: 1260 s
    (chains 952-1260 s);
  - the same fit rerun by a verifier on a heavily loaded machine (MATLAB and
    other R jobs, load average 17-24 on 6 cores): 1544 s (chains
    1229-1544 s), 26.3 minutes wall including the compile check and the
    gate; it reproduced the seed-1 fit exactly;
  - the research fits, run with 2 parallel chains (two batches of two):
    seed 1 2396 s, seed 2 1916 s.

  The 30-minute condition holds only with 4 parallel chains, and with little
  margin on a busy machine. `output/tables/fit_diagnostics.csv` gives the
  research fits' times (wall_secs; for saved fits it is `fit$time()$total`;
  m12's 2396 s is the 2-parallel-chain seed-1 fit).

### m3: rewrite and identification

`stan/m3_corr.stan` now integrates mu out with a constant-R Kalman filter
(R = diag(sigma) Omega diag(sigma), stable quadratic form) and draws mu by
FFBS. The model, priors and outputs are unchanged; `log_lik_pred` is added.
Tests: `tests/test_m3_m10.R` (dense MVN to 1e-8, FFBS against the dense
smoother, equality with the LFO filter on 5 posterior draws).

| | old (mu sampled) | new (mu integrated out) |
|---|---|---|
| divergences | 0 | 0 |
| min E-BFMI | 0.03 | 0.89 |
| max R-hat | 1.029 | 1.007 (Omega[3,1]) |
| min bulk / tail ESS | 132 / - | 882 / 438 |
| mu: max R-hat / min bulk ESS | - | 1.005 / 944 |
| run time (CmdStan sampling time) | - | 66 s |

The posterior is the same: sigma 4.35 / 3.13 / 2.43, tau 4.04, Omega
0.07 / 0.04 / -0.23.

Omega is still weakly identified, now as a clear property of the model rather
than a sampler failure (`Rscript R/diag_m3_ident.R`, which writes the numbers
below to `output/tables/m3_identification.csv`). With rho near 0 the data pin down only
v·11' + Sigma, where v is the latent variance, plus a lag-1 autocovariance
rho·v of 0.97 with sd 1.20, which is essentially uninformative. Evidence:
- corr(v, Sigma_jk) is -0.56 to -0.78, and corr(tau, mean off-diagonal
  Omega) is -0.70;
- in the lowest-tau decile (tau 3.16) Omega is (0.29, 0.31, 0.21); in the
  highest (tau 4.55) it is (-0.02, -0.09, -0.43);
- the posterior sd of Omega is 46-83% of the LKJ(4) prior sd (0.316);
- the variance ratio's 90% interval is [0.59, 0.92], wider than m2's.

m3's tail ESS for tau, sigma and Omega is 438-522, just above 400, because of
the heavy low-tau tail along this ridge.

### m10: two modes, diagnosis and re-specification

**Diagnosis.** 16 chains of the old specification (constant tau) from
dispersed starting values, 1000 + 1000 each, found two modes:

| | mode A (9 chains) | mode B (3 chains) |
|---|---|---|
| rho, tau | 0.39, 2.1 | 0.12, 3.6 |
| phi_E, omega_E | 0.66, 0.67 | 0.97, 0.22 |
| lp__ | about -2287 | about -2330 |
| 2020 Q2 mu (headline -28.1) | -1.3 | -26.9 |
| 2020 Q2 measurement sds | about 9-11 for all four | about 1.5 |
| variance ratio | 0.36 | 0.86 |

The other 4 chains were stranded between the two (tau 2.4-3.2). Mode A says
four separately built series (E, I, P and hours) were each wrong by about -27
in the same quarter. Mode B inflates tau for every quarter. Comparing the
modes' mass with h integrated out was not possible: the particle filter's
estimate for mode A kept rising with the particle count (about -2020 at
N = 2000, -2006 at 20000, -2000 at 100000, against a stable -1982 for B),
because mode A needs four volatilities to jump at once.

Both modes react to one mis-specification. With a constant latent sd on the
post-1978 sample (pre-2020 common covariance implies tau of about 2.3), COVID
is a roughly 13-sigma latent shock. The model must either inflate tau
everywhere (B) or call a common shock measurement error (A). Changing the
sampling coordinates cannot fix that.

**Decision: re-specified.** The latent innovation sd is now exp(g_t / 2),
with g_t a stationary AR(1) that has the same priors as each measurement SV
block (mu_g ~ N(2, 1), phi_g ~ U(-1, 1), omega_g ~ half-N(0, 0.5)). This
replaces tau ~ half-N(0, 2). SV on all four measurement errors is kept. The
justification is economic and does not depend on the fit: the volatility of
real activity has changed over time (lower since the 1980s), and COVID was a
real shock on which E, I and P agree within 0.5 pp and hours within about
6 pp. mu is integrated out by a Kalman filter with time-varying state
variance; g and h are non-centred.

**Caveat: chosen after seeing the data.** The re-specification was chosen
after the old model's two modes, its fit and its LFO score had been seen.
The new m10's LFO gains (+124 over m8 on the joint target; +47.0 over m12
and +63.8 over m9 on GDP-only) are therefore post hoc: they come from a
second specification scored on the same evaluation window, so their
z-values overstate the evidence. m10 stays a research model and plays no
part in the production decision.

**Result.**
- 16 of 16 dispersed chains land in one mode (R-hat 1.00 across all 16;
  rho 0.43-0.45, tau 1.51-1.57; 8 divergences in 16000). Per-chain table:
  `output/tables/m10_mode_check.csv` (summarised from the saved
  diagnosis fit). `Rscript R/diag_m10_modes.R dispersed` repeats the design
  (same init distributions, seeds and settings); its inits are drawn by a
  tidied copy of the original scratch script, so they are not bit-identical. These chains start g_raw and h_raw near N(0, 0.5),
  so none starts in old mode A; the "modeA" run below tests that case.
- 4 chains started at an old-mode-A configuration (all four measurement sds
  about 10 in 2020 Q2 and 8 in 2020 Q3, flat latent volatility, tau 2.1,
  rho 0.39, phi_E 0.66; `Rscript R/diag_m10_modes.R modeA`, 700 + 500
  iterations, seed 2024) all leave it and end in the same single mode:
  lp -2333 to -2339, rho 0.42-0.45, tau 1.53-1.56, 2020 Q2 mu -28.1 to
  -28.2 with sd_E about 1.4, 0 divergences, R-hat 1.01 across the 4 chains
  (run "modeA" in `output/tables/m10_mode_check.csv`). This repeats a
  reviewer's independent test with the same result.
- The production-settings fit (seed 1) has 4 divergences in 8000, E-BFMI
  0.86-0.93, max hyperparameter R-hat 1.0024 (phi[1]), min bulk / tail ESS
  1426 / 1758, and meets the production gate's thresholds.
- Posterior means: rho 0.44, phi_g 0.81, omega_g 0.94, beta 0.95;
  tau = exp(mu_g / 2) 1.59 (median 1.54). tau_t is about 0.8-1.9 by decade
  and about 14 in 2020 Q2-Q3. In 2020 Q2 mu is -28.2 and the measurement sds
  are about 1.4.
- Variance ratio (1978 Q4 - 2026 Q2) 0.747 [0.678, 0.824]. Excluding
  2020 Q1 - 2021 Q1 it is 0.56, against 0.84 for m8, 0.86 for m9 and 0.88 for
  m12 on the same quarters.
- omega_g's posterior (90% interval [0.63, 1.33]) sits in the upper tail of
  its prior, which puts 6.6% of its mass above 0.92. No prior-sensitivity
  check has been run.

**LFO filter.** `R/lib/filter_m10.cpp` / `filter_m10.R`: particles over
(g, h), mu by per-particle Kalman filter, a locally adapted proposal for g
with exact weight correction. Tested in `tests/test_m3_m10.R` against the
dense MVN, the m8 Kalman filter at omega -> 0, brute-force Monte Carlo over
200000 g paths, and the Stan `log_lik_pred`.

### LFO method changes

- **Refit rule (tightened 2026-10-16).** A refit is accepted when divergences
  <= 1% of draws, max hyperparameter R-hat <= 1.05, min E-BFMI >= 0.2 and min
  bulk ESS >= 100. "Hyperparameters" are the parameters the filter uses. A
  failing refit is retried with a new seed up to 3 times. The earlier rule
  checked only divergences and R-hat. The rule is deliberately looser than
  the production gate, because refits have 4000 draws, not 8000.
- **Result.** All 75 stored fits (12 initial fits and 63 refits) meet it
  (`n_fits_below_rule` = 0 for every model).
- **Failed attempts.** Three fits were kept only after earlier attempts
  failed the rule: 6 failed attempts in all
  (`output/tables/lfo_refit_attempts.csv`; `n_attempts_failed` and
  `n_fits_retried` in `model_comparison.csv`):

  | fit | failed attempts | kept |
  |---|---|---|
  | m9, data to 2006 Q3 (1:188) | 1: 106 divergences, bulk ESS 63 | attempt 2: 0 divergences, ESS 1071 |
  | m9, data to 2020 Q4 (1:245) | 1-3: tau R-hat 1.069 / 1.053 / 1.141, ESS 53 / 58 / 20 | attempt 4 (seed 34002): R-hat 1.040, ESS 106 |
  | m12, data to 2020 Q4 (1:245) | 1: ESS 99; 2: R-hat 1.068, ESS 59 | attempt 3: R-hat 1.013, ESS 125 |

  Retrying until a seed passes is seed selection. At the 2020 Q4 data end
  point both posteriors are not sampled reliably by 4 x (1000 + 1000) chains
  (tau mixes slowly after the COVID quarters), so the kept attempt is the one
  that happened to pass, not evidence that the posterior is well sampled.
  The effect is measured below ("Seed-selection check").
- **Seed-selection check** (`Rscript R/03_compare.R retrycheck m9 m12`;
  `output/tables/lfo_retry_sensitivity.csv` and
  `lfo_retry_sensitivity_summary.csv`). The failed attempts were re-created
  from their seeds; every one reproduced its logged diagnostics exactly. The
  53 targets served by the three retried fits (m9: 2006 Q4 - 2014 Q3 and
  2021 Q1 - 2026 Q1; m12: 2021 Q1 - 2025 Q1) were then rescored with
  attempt 1's draws alone (no retry) and with all attempts' chains pooled
  (equal draws per attempt, no choice among seeds), holding the refit
  schedule and particle schedules fixed:

  | m12 - m9 | diff | SE (NW) | z (NW) |
  |---|---|---|---|
  | as reported (kept attempts) | +16.77 | 6.12 | 2.74 |
  | attempt 1 only | +16.09 | 6.52 | 2.47 |
  | all attempts pooled | +16.62 | 6.29 | 2.64 |
  | the 53 targets dropped (93 left) | +11.41 | 6.31 | 1.81 |

  Selecting the passing seed moved m12 - m9 by +0.7 nats at most, and the
  z >= 2 condition holds without any selection. Pooling does not make the
  2020 Q4 posteriors well sampled: over all 16 m9 chains tau's R-hat is
  1.080 (min bulk ESS 147); over m12's 12 chains it is 1.037 (ESS 270). The
  "dropped" row removes evidence rather than measuring selection (most of
  the 53 targets are the 2006-2014 quarters served by m9's 1:188 fit, whose
  attempt choice changes m9's score by 0.1 nats), and with z 1.81 it falls
  short of the decision threshold; it is shown for completeness. Other key
  differences change by at most 0.5 nats (m9 - m4: +1.48 reported, +2.02
  attempt 1 only, +1.67 pooled).
- **Weak mixing.** Ten fits that meet the rule have min bulk ESS below 400
  (`weak_mixing` in `lfo_fit_diagnostics.csv`; none has E-BFMI below 0.3).
  They would not pass the production gate:
  - m12: data to 1994 Q2 (ESS 336, phi[3]), 2020 Q2 (210, tau), 2020 Q4
    (125, tau);
  - m9: data to 2000 Q3 (231, phi[3]), 2020 Q2 (210, tau), 2020 Q4 (106,
    tau);
  - m4: data to 2003 Q2 (332, omega[2]), 2008 Q2 (399.9, phi[3]), 2020 Q2
    (117, tau; R-hat 1.033);
  - m3: data to 2020 Q2 (391, Omega[2,1]).

  The weak spots are tau on the fits just after COVID (all three SV models)
  and an SV persistence or volatility parameter on some pre-COVID fits.
- **Reuse.** m4 and m9 refits with the same Stan file, data end point and seed
  as an earlier checkpoint were taken from it (5 m4 and 4 m9; the m4 1:175 cache hit is
  logged twice because of a restart) and checked
  against the new rule. The linear models' stored refits also meet the new
  rule, and recomputing their scores on the current data reproduced them
  exactly, so they were not rerun.
- **Particle schedule.** Per refit, a pilot pass sets the particle count per
  quarter: 2000 by default, raised where the ESS fraction collapses, capped
  at 100000. The filter runs on every 2nd posterior draw (every 4th for m10).
- **Monte Carlo error.** Each SV model's score is the mean of 4 independent
  filter runs; MCSE = sd(run totals) / 2. m12 runs -1105.88 / -1105.72 /
  -1105.24 / -1105.04 (MCSE 0.20); m9 -1121.50 / -1121.97 / -1122.66 /
  -1122.83 (MCSE 0.31); m4 MCSE 0.02; m10 0.08. The spread comes almost
  entirely from 2020 Q2 (per-run sd 0.41 for m12, 0.65 for m9).
  - With 4 runs each sd has 3 degrees of freedom, so the MCSE is itself
    very uncertain: the 95% upper bound of an sd estimated from 4 runs is
    about 3.7 times the estimate (sqrt(3 / qchisq(0.025, 3))). For the
    m12 - m9 difference (MCSE 0.37) that bound is about 1.4 nats, still far
    below the 16.8-nat difference.
  - The run totals drift in opposite directions for the two models (m12
    rises, m9 falls across the runs). The main run alone, which made the
    refit decisions, gives m12 - m9 = +15.62 rather than the 4-run mean of
    +16.77, with Newey-West SE 6.98 (iid 7.31), so z_NW = 2.24: still above
    2, but closer to it. Averaging the runs removes particle noise from the
    pointwise differences (mostly at 2020 Q2), which is why the 4-run SE is
    smaller.

### Particle-filter accuracy (bias against high-particle references)

`Rscript R/03_compare.R pfcheck` / `pfref`; tables
`output/tables/lfo_pf_covid_check.csv` and `lfo_pf_precovid_check.csv`. The
reference uses at least 8000 particles, up to 400000. 2020 Q2 has 3 reference
runs (m12 -47.03, se 0.10; m9 -43.84, se 0.27; m10 -26.67, se 0.02; m4
-43.82, se 0.04); every other target has 1. m10 and m4 were added on
2026-10-17 afternoon (m10 also on its GDP-only target); the 2015-19 check
was run for m12 and m9 only.

| model | COVID quarters, schedule | COVID quarters, fixed 2000 | 2015-19, schedule | 2015-19, fixed 2000 |
|---|---|---|---|---|
| m12 | +0.16 | -0.37 | -0.017 | -0.022 |
| m9 | +0.02 | -0.47 | +0.005 | -0.019 |
| m4 | -0.27 | -0.28 | not run | not run |
| m10, joint target | -0.15 | +0.11 | not run | not run |
| m10, GDP-only target | -0.01 | -0.07 (against the first reference run) | not run | not run |

m4's and m10's COVID bias is negative (their scores are slightly too low),
almost all of it at 2020 Q2 (m4 -0.27, m10 -0.14). For m10 this is far
smaller than the filter's finite-N bias in `tests/test_filters.R` (vii),
where the mean log p-hat on simulated data is 1.9 nats below the brute-force
value, because the schedule raises N where the weights collapse. Effects on the reported differences: m12 - m4 is overstated by
about 0.4 nats (out of 18.3); m10 - m8 is understated by about 0.15, and
m10 - m12 on GDP-only by about 0.17, so the m10 "ahead" figures are, if
anything, slightly conservative.

The COVID quarters are 2020 Q1 - 2021 Q1 (bias summed over them). The
schedule's net bias on m12 - m9 is about +0.1 nats, far below the 16.8-nat
difference. One caveat: in 2020 Q3 - 2021 Q1 the per-quarter estimates sit up
to 0.27 nats from the reference (in both directions), although repeated runs
agree to 0.01, so the run-to-run spread understates the error after a
volatility jump. These deviations are almost the same for m9 and m12 and
cancel in their difference. A comparison between an SV model and a
constant-variance model carries about 0.1-0.5 nats of uncorrected bias over
the COVID quarters. For m4 the check gives -0.27 nats over the COVID quarters
(table above), consistent with its 4 runs agreeing to 0.03 nats; no 2015-2019
check was run for m4 or m10.

The 2026-10-16 statement that the finite-N bias is "about 0.3-0.4 nats
against the SV models" referred to the fixed N = 2000 filter; it is
consistent with the -0.37 / -0.47 above and is superseded by the schedule.

### Comparison results (current)

Group 1: E, I, P, full sample; 146 predictions, 1990 Q1 - 2026 Q2
(`output/tables/model_comparison.csv`, `model_comparison_pairwise.csv`):

| model | elpd_lfo | SE | diff vs m12 | SE diff | SE diff (NW) | refits | below rule | MCSE |
|---|---|---|---|---|---|---|---|---|
| m12 m9 + AR(1) errors     | -1105.47 | 42.7 | 0 | - | - | 7 | 0 | 0.20 |
| m9 structural Sigma + SV  | -1122.24 | 40.0 | -16.77 | 6.67 | 6.12 | 7 | 0 | 0.31 |
| m4 SV                     | -1123.72 | 40.0 | -18.25 | 6.92 | 7.53 | 6 | 0 | 0.02 |
| m5 serial errors          | -1168.84 | 29.3 | -63.37 | 14.33 | 20.71 | 5 | 0 | exact |
| m7 COVID scale            | -1177.54 | 14.1 | -72.07 | 31.05 | 34.41 | 4 | 0 | exact |
| m6 structural Sigma       | -1195.59 | 29.1 | -90.12 | 15.08 | 21.39 | 3 | 0 | exact |
| m1 baseline               | -1195.66 | 29.6 | -90.19 | 14.59 | 21.72 | 4 | 0 | exact |
| m3 free correlation       | -1195.93 | 29.2 | -90.46 | 15.09 | 21.02 | 3 | 0 | exact |
| m2 AR(1)                  | -1196.25 | 29.3 | -90.78 | 14.95 | 21.25 | 4 | 0 | exact |

Group 2: E, I, P, hours, post-1978 (joint 4-vector target):

| model | elpd_lfo | SE | diff | SE diff | SE diff (NW) | refits | MCSE |
|---|---|---|---|---|---|---|---|
| m10 4-signal + SV + latent SV (re-specified) | -1417.72 | 26.4 | 0 | - | - | 7 | 0.08 |
| m8 4-signal | -1541.96 | 72.2 | -124.25 | 53.51 | 64.72 | 6 | exact |

GDP-only target, same window: m10 -1058.44, m12 -1105.47, m11 -1113.36,
m9 -1122.24, m4 -1123.72, m8 -1160.69, m5 -1168.84. m8 and m10 also see past
hours; m11 sees monthly employment up to the end of the target quarter. The
m11 joint score is -2469.88, not comparable with either group.

Key pairwise differences (a - b):

| a | b | target | diff | SE | SE (NW) | z (NW) |
|---|---|---|---|---|---|---|
| m12 | m9 | group 1 | +16.77 | 6.67 | 6.12 | 2.74 |
| m12 | m4 | group 1 | +18.25 | 6.92 | 7.53 | 2.42 |
| m12 | m5 | group 1 | +63.37 | 14.33 | 20.71 | 3.06 |
| m12 | m7 | group 1 | +72.07 | 31.05 | 34.41 | 2.09 |
| m9 | m5 | group 1 | +46.59 | 12.46 | 19.66 | 2.37 |
| m9 | m4 | group 1 | +1.48 | 1.53 | 2.18 | 0.68 |
| m5 | m2 | group 1 | +27.41 | 2.95 | 1.95 | 14.05 |
| m6 | m2 | group 1 | +0.66 | 0.34 | 0.38 | 1.74 |
| m10 | m8 | group 2 joint | +124.25 | 53.51 | 64.72 | 1.92 |
| m10 | m8 | GDP-only | +102.24 | 35.40 | 45.65 | 2.24 |
| m10 | m12 | GDP-only | +47.03 | 28.49 | 30.11 | 1.56 |
| m10 | m9 | GDP-only | +63.80 | 24.66 | 25.61 | 2.49 |
| m12 | m11 | GDP-only | +7.89 | 18.73 | 27.05 | 0.29 |

m12 - m9 by sub-period (`model_comparison_subperiod.csv`):

| period | n | diff | SE iid | SE NW |
|---|---|---|---|---|
| 1990-99 | 40 | +6.69 | 3.06 | 2.05 |
| 2000-07 | 32 | +9.23 | 2.05 | 1.06 |
| 2008-19 | 48 | +3.94 | 2.22 | 1.52 |
| 2020 | 4 | -4.18 | 4.94 | 1.74 |
| 2021-22 | 8 | +0.41 | 0.62 | 0.52 |
| 2023-26 | 14 | +0.68 | 1.00 | 0.87 |
| pre-2020 | 120 | +19.86 | 4.36 | 3.82 |
| 2020 on | 26 | -3.09 | 4.95 | 3.58 |
| all excl. COVID | 142 | +20.53 | 4.53 | 4.25 |
| **all** | **146** | **+16.77** | **6.67** | **6.12** |

m12 is ahead in 99 of 146 quarters (median gain +0.15 per quarter). Its
largest losses are 2020 Q4 (-3.1) and 2020 Q2 (-3.0).

**Conclusions.**
1. **Serial correlation and SV are separate gains that add up.** m5 - m2 =
   +27.4 (constant variance), m12 - m9 = +16.8 (with SV), m9 - m5 = +46.6
   and m12 - m5 = +63.4 (z 3.1). Conclusion 4 of 2026-10-16 ("m9 + AR(1)
   measurement errors is the most promising next extension") is borne out.
2. **m9 against m4 is still a tie** (+1.5, SE 1.5 iid, 2.2 NW). The
   structural restriction is kept on construction grounds only.
3. **Hours.** The re-specified m10 now clearly beats m8 on the joint target
   (+124, z 1.9 NW, 2.3 iid) and on GDP-only (+102, z 2.2 NW). On GDP-only it
   is ahead of m9 by +63.8 (z 2.5 NW) and of m12 by +47.0 (z 1.6 NW): the
   latent-volatility feature, not just hours, may be what helps. It uses a
   shorter sample and past hours, so it is a research model.
4. **The comparison still does not test the variance decomposition.** Over
   the same 1978 Q4 - 2026 Q2 quarters, m10's latent variance ratio is 0.75
   [0.68, 0.82] and m12's 0.91 [0.84, 0.97] (0.56 and 0.88 excluding
   2020 Q1 - 2021 Q1); m10 predicts GDP somewhat better. "About a fifth of
   headline variance is measurement error" is m12's answer, not a result the
   comparison validates. A model with latent SV (m10) puts the
   measurement-error share much higher outside COVID.

**Caveats.**
- z of 2.5-2.7 for m12 - m9 comes from a single evaluation window; the
  pointwise differences are autocorrelated (lag-1 autocorrelation -0.39, so
  here the iid SE is the larger), and the Newey-West lag of 8 is a choice.
  The z >= 2 result holds at every lag tried: Newey-West SE 5.61 / 6.12 /
  6.57 / 6.99 at lags 4 / 8 / 12 / 16 (z 2.99 / 2.74 / 2.55 / 2.40), and a
  moving-block bootstrap of the sum gives SE 5.45 / 6.13 / 6.65 for blocks
  of 4 / 8 / 12 quarters.
- The 2020 Q2 reference has 3 runs and every other reference 1, so the
  reference measures the bias with error.
- The tightened rule is applied to cached refits from the earlier run (same
  Stan file and seed, so the same draws) and to the earlier linear-model
  refits; it is not a fresh run.
- m11's joint score includes employment inside the target quarter; only its
  GDP-only column is comparable, and with a larger information set.
- Three times on 2026-10-17 a background LFO job and its parent shell
  disappeared without an exit status (11:05-11:19, 12:07-12:25 and, in the
  m4 particle-count check, 14:44-15:09), at about the same time as other
  heavy jobs on the machine ended. The cause is unknown. Saved work was kept
  and only the missing parts were rerun; all three events are in
  `output/lfo/lfo_log.txt`.
- Provenance of m12's initial LFO fit (data to 1989 Q4): it was made by a
  process started before a 2026-10-16 change to `R/lib/lfo.R` (a further
  bookkeeping-only change on 2026-10-17 records every refit attempt), and
  the restarted run resumed from it. Re-creating that fit and its main
  filter pass with the current code reproduced the stored draws, particle
  schedule and pointwise scores exactly (logged in `lfo_log.txt`).
- Pseudo-real-time only, as before (current vintage, one-step-ahead).

### Production decision: m9 -> m12

The rule was fixed before m12's LFO results were known: switch production to
m12 only if (1) it passes the production gate on two seeds, (2) its LFO score
is reliable, (3) it beats m9 by z_NW >= 2, and (4) a local fit takes at most
30 minutes.

| Condition | Observed | Met |
|---|---|---|
| gate on two seeds | seeds 1 and 2 both pass (table above) | yes |
| reliable LFO | MCSE of m12 - m9 0.37 (from 4 runs, 95% upper bound about 1.4); measured bias about +0.1; all 75 fits meet the refit rule, 3 of them only after retries (rescored without seed selection: +16.09 / +16.62, z 2.47 / 2.64) | yes |
| z_NW >= 2 | +16.77 / 6.12 = 2.74 (main filter run alone: +15.62 / 6.98 = 2.24) | yes |
| local fit <= 30 min | 21 min (1260 s, 4 parallel chains; chains 952-1260 s); 25.7 min (1544 s) on a heavily loaded machine; 32-40 min with 2 parallel chains | yes, with 4 parallel chains; thin margin under load |

Production therefore uses `stan/m12_struct_sv_serial.stan` from 2026-10-17:
`R/fit_production.R` (psi added to the gate's hyperparameters),
`R/output_estimate.R`, `R/input_hash.R`, `R/build_site.R` and the CI
workflow (executable cache, job timeout raised from 90 to 120 minutes). The
production fit on the current data (`Rscript R/fit_production.R`, seed 1,
2026-10-17) reproduces the seed-1 research fit exactly: 3 divergences in 8000,
E-BFMI 0.862-0.989, max R-hat 1.0020 (omega[3]) on hyperparameters and
1.0015 / 1.0019 on mu / h, min bulk / tail ESS 2118 / 3520 on
hyperparameters and 2852 / 6476 on mu, 0 treedepth hits; gate passed.

Open items:
- the warmup-stall harness for m12 has not been run at init = 0.5;
- the fit's run time on GitHub's runners has not been measured;
- m12 does worse than m9 in 2020 (-4.2 over the four quarters) in one-step
  prediction. Its published 2020 path is close to m9's: median QoQ growth
  -0.33 / -6.84 / +3.38 / +3.52 in 2020 Q1-Q4, against m9's last published
  -0.27 / -6.84 / +3.40 / +3.41;
- m9 + latent SV (as in the new m10) is the natural next candidate.

## Model comparison: PSIS-LFO-CV (2026-10-16; results superseded 2026-10-17)

This section replaces the PSIS-LOO comparison in earlier versions of this file,
which was invalid (see the next section). **Its result tables are superseded**
by "Round of 2026-10-17" above, which adds m12, re-scores m3, m4, m9 and m10
with the tightened refit rule and the adaptive particle filter, and replaces
the old m10. The method below still applies except where that section says
otherwise. The old numbers are kept for the record.

### Method

`R/03_compare.R` with `R/lib/filters.R` and `R/lib/lfo.R`.
- **Score.** The one-step-ahead predictive density
  log p(y_t | y_{1:t-1}), computed with PSIS-LFO-CV (Bürkner, Gabry & Vehtari
  2020). The hyperparameters are integrated over the posterior fitted to
  y_{1:t-1}. The latent states are integrated out by a filter, so each score
  conditions only on the past:
  - linear-Gaussian specs (m1-m3, m5-m8, m11): an exact Kalman filter on all
    4000 draws;
  - SV specs (m4, m9, m10): a Rao-Blackwellised particle filter (particles over
    h, Kalman over mu) with N = 2000 particles, on 2000 draws.
- **Window.** 1990 Q1 - 2026 Q2, 146 one-step predictions. It is the same
  calendar window for every model.
- **Initial fits.**
  - Group 1 (m1-m7, m9): 121 quarters, to 1989 Q4.
  - Group 2 (m8, m10): 45 quarters.
  - m11: 47 quarters plus the matching months.
- **Refits.**
  - A refit is triggered when the Pareto k exceeds 0.7. The cap is 12 refits
    per model; none reached it (3-8 refits per model). No step with k > 0.7
    was used without a refit.
  - Each refit runs the unmodified Stan file with 4 x (1000 + 1000) iterations.
  - Quality rule: a refit with more than 1% divergences or an R-hat above 1.05
    is retried with a new seed, up to 3 attempts. If all attempts fail, the
    best is kept and logged. The rule was fixed before the SV runs and does
    not test E-BFMI or ESS. The summary flags refits with min E-BFMI < 0.3 or
    min bulk ESS < 400 (the production gate's thresholds) as `weak_mixing` in
    `lfo_fit_diagnostics.csv`: m3 5 of its 5 fits, m4 3 of 8, m9 3 of 8, m10 2 of 9, all other
    models 0. They were used as they are.
  - `output/lfo/lfo_log.txt` from this run prints numbers to 4 significant
    figures, seeds included, so its seeds are wrong (10121 is logged as
    10120). The seeds in `lfo_fit_diagnostics.csv` are correct; the logging is
    fixed for future runs.
- **Uncertainty.**
  - The SE of a difference is sqrt(n var(pointwise difference)). This ignores
    serial correlation in the pointwise differences, which is positive
    (lag-1 autocorrelation 0.1-0.6 for the SV-versus-constant-variance
    pairs). `R/03_compare.R` therefore also reports a Newey-West SE
    (Bartlett kernel, 8-quarter lag; columns `se_diff_nw`, `z_nw`,
    `se_diff_vs_best_nw`). Both are shown below.
  - Monte Carlo error for the particle filter comes from re-running every
    filter pass with new seeds: m9 gave -1122.57 against -1123.36, m4 -1123.54
    against -1123.54, and m10 -1504.90 against -1507.11. This has two limits.
    With one rerun per model the MCSE has 1 degree of freedom. And a rerun at
    the same N = 2000 measures spread, not the downward bias of log p-hat at
    finite N. That bias is concentrated in the COVID quarters: for m9, the
    2020 Q3 prediction from the fit to 2020 Q2 is -13.11 / -13.12 at N = 2000
    (two seeds) against -12.79 at N = 20000. Most of the m9 and m10 rerun
    spread also comes from one fit each (m9: the fit predicting 2020 Q2). The
    bias is about 0.3-0.4 nats against the SV models in total and changes no
    conclusion.
  - The particle filter draws the initial h_1 from a proposal with sd capped
    at 3, with an exact importance correction. When the stationary sd of h
    exceeds 3 sqrt(2) = 4.24 (up to 8.8% of m9's LFO draws, 5.1% of m4's), the
    weights have infinite variance and the period-1 estimate is biased low in
    practice. Period 1 is 1959 Q4, outside the window: recomputing m9's LFO
    with the plain bootstrap initialisation gives -1123.50 against -1123.36,
    within Monte Carlo error.
- **Tests.** `tests/test_filters.R` checks the filters:
  - the Kalman filter against the dense multivariate-normal likelihood for
    every linear-Gaussian spec;
  - the particle filter at omega -> 0 against the Kalman filter;
  - the particle filter against the Stan programs' `log_lik_pred`;
  - Monte Carlo variability across seeds;
  - numerical robustness at extreme h.
- **Outputs.** `output/tables/model_comparison.csv`,
  `model_comparison_pairwise.csv`, `lfo_pointwise.csv`,
  `lfo_fit_diagnostics.csv`, `output/lfo/lfo_log.txt`, and figure
  `output/figures/07_lfo_cumulative_elpd.png`.

### Group 1: E, I, P, full sample (target: the 3-vector of GDP measures) *(superseded 2026-10-17)*

| model | elpd_lfo | SE | diff vs m9 | SE diff | SE diff (NW) | refits | MCSE |
|---|---|---|---|---|---|---|---|
| m9 structural Sigma + SV | -1123.36 | 41.2 | 0 | - | - | 7 | 0.56 |
| m4 SV                    | -1123.54 | 40.6 | -0.18 | 1.59 | 2.41 | 7 | 0.004 |
| m5 serial errors         | -1168.84 | 29.3 | -45.48 | 13.48 | 20.24 | 5 | exact |
| m7 COVID scale           | -1177.54 | 14.1 | -54.18 | 29.21 | 32.45 | 4 | exact |
| m6 structural Sigma      | -1195.59 | 29.1 | -72.23 | 13.28 | 20.51 | 3 | exact |
| m3 free correlation      | -1195.63 | 29.4 | -72.27 | 13.07 | 20.27 | 4 | exact |
| m1 baseline              | -1195.66 | 29.6 | -72.30 | 13.07 | 20.78 | 4 | exact |
| m2 AR(1)                 | -1196.25 | 29.3 | -72.89 | 13.12 | 20.40 | 4 | exact |

### Group 2: E, I, P, hours, post-1978 (target: the joint 4-vector) *(superseded: old m10)*

| model | elpd_lfo | SE | diff | SE diff | SE diff (NW) | refits | MCSE |
|---|---|---|---|---|---|---|---|
| m10 4-signal + SV | -1507.11 | 72.0 | 0 | - | - | 8 | 1.57 |
| m8 4-signal       | -1541.96 | 72.2 | -34.85 | 11.39 | 18.83 | 6 | exact |

The m11 joint score (employment plus GDP) is -2469.88. It is not comparable
with either group.

### GDP-only target, same window, all models *(superseded 2026-10-17)*

| model | m11 | m9 | m4 | m10 | m8 | m5 | m7 | m1, m2, m3, m6 |
|---|---|---|---|---|---|---|---|---|
| GDP-only elpd | -1113.36 | -1123.36 | -1123.54 | -1124.54 | -1160.69 | -1168.84 | -1177.54 | -1195.6 to -1196.2 |

The information sets differ. m8 and m10 also see past hours. m11 also sees
monthly employment up to the end of the target quarter, because employment is
released before GDP.

### Key pairwise differences (a - b, SE) *(superseded 2026-10-17)*

| a | b | target | diff | SE | SE (NW) | z (NW) |
|---|---|---|---|---|---|---|
| m9 | m4 | group 1 | +0.18 | 1.59 | 2.41 | 0.08 |
| m9 | m5 | group 1 | +45.48 | 13.48 | 20.24 | 2.25 |
| m4 | m5 | group 1 | +45.29 | 13.05 | 19.40 | 2.33 |
| m9 | m6 | group 1 | +72.23 | 13.28 | 20.51 | 3.52 |
| m9 | m2 | group 1 | +72.89 | 13.12 | 20.40 | 3.57 |
| m9 | m7 | group 1 | +54.18 | 29.21 | 32.45 | 1.67 |
| m5 | m2 | group 1 | +27.41 | 2.95 | 1.95 | 14.05 |
| m7 | m2 | group 1 | +18.71 | 17.87 | 17.32 | 1.08 |
| m6 | m2 | group 1 | +0.66 | 0.34 | 0.38 | 1.74 |
| m10 | m8 | group 2 joint | +34.85 | 11.39 | 18.83 | 1.85 |
| m10 | m8 | GDP-only | +36.14 | 9.25 | 16.34 | 2.21 |
| m9 | m10 | GDP-only | +1.18 | 11.88 | 22.60 | 0.05 |
| m11 | m9 | GDP-only | +10.00 | 16.58 | 24.66 | 0.41 |
| m9 | m8 | GDP-only | +37.33 | 12.62 | 21.63 | 1.73 |

m9 against m4 depends on the period: +2.43 (SE 0.75) before 2020 and -2.25
from 2020 on.

### Where the SV advantage comes from (elpd summed by period, 2026-10-16 scores)

| period | m2 | m5 | m7 | m4 | m9 |
|---|---|---|---|---|---|
| 1990-99 | -340.1 | -332.5 | -340.1 | -334.2 | -334.6 |
| 2000-07 | -254.6 | -247.5 | -254.6 | -242.6 | -241.2 |
| 2008-19 | -367.5 | -359.1 | -367.0 | -335.2 | -333.7 |
| 2020    |  -64.9 |  -64.6 |  -46.1 |  -73.0 |  -74.9 |
| 2021-22 |  -66.2 |  -64.9 |  -67.8 |  -60.9 |  -61.5 |
| 2023-26 | -102.9 | -100.2 | -101.9 |  -77.7 |  -77.5 |

The SV gain builds from 2000 onward and is largest in 2008-19 and 2023-26, as
measurement-error variance falls. The SV models lose in 2020: m9 loses 9.7
nats to m5 in 2020 Q2 alone. m7's whole advantage over m2 comes from 2020,
because the model knows in advance which four quarters are the COVID quarters.

### Conclusions (2026-10-16; see "Round of 2026-10-17" for the current ones)

1. **SV versus constant variance: consistent evidence, of moderate
   strength once serial correlation is allowed for.** m9 beats the best
   constant-variance 3-measure model (m5) by 45.5 nats over 146 quarters,
   about 0.3 nats per quarter: SE 13.5 iid, 20.2 Newey-West (z about 2.2).
   It beats m1, m2, m3 and m6 by about 72 (z about 3.5 Newey-West). Against
   m7 the margin is only +54.2 (SE 32.5 Newey-West, z 1.7), but m7 knows the
   COVID quarters in advance. The direction is robust: SV is ahead in every
   sub-period except 2020 (table above), apart from the 1990s, when m5 is
   about 2 nats ahead. m9 - m2 is +52.7 (SE 5.8) before
   2020 and +84.1 excluding the four COVID quarters.
2. **m9 versus m4: no evidence either way.** The difference is +0.18
   (SE 1.59): m9 is ahead before 2020 and behind after. Keeping m9 rather
   than m4 in production rests on the a-priori construction argument:
   Production's errors are independent by how the ABS builds the series.
   Diagnostics do not separate them either: both full fits are clean (m4 5
   divergences, E-BFMI 0.89, bulk ESS 1512; m9 1, 0.90, 1692).

   Predictive accuracy is not part of the argument. The earlier "m9 best by LOO, +6
   over m4" claim is withdrawn.
3. **The structural restriction alone is neutral.** m6 against m2 is +0.66
   (SE 0.34).
4. **Serially correlated errors help.** m5 against m2 is +27.4 (SE 3.0). No SV
   model has this feature, so m9 + AR(1) measurement errors is the most
   promising next extension. *(Done 2026-10-17: that model is m12, which beats
   m9 by +16.8, NW z 2.7, and is now the production model.)*
5. **Hours and monthly employment add nothing detectable for GDP.**
   - m10 is ahead of m8 by 34.9 (SE 11.4 iid, 18.8 Newey-West, z 1.85),
     suggestive rather than decisive. On the GDP-only target it is level
     with m9: m9 - m10 = +1.2 (SE 11.9). *(Superseded: these are the old,
     invalid m10 scores. The re-specified m10 is ahead of m8 by +124 and of
     m9 by +63.8 on GDP-only.)*
   - m11 - m9 is +10.0 (SE 16.6), even though m11 has constant variance and
     sees same-quarter employment.

   An SV version of m11 is worth testing as a monthly nowcast.

### Caveats (2026-10-16)

- **Pseudo-real-time only.** The evaluation uses the current data vintage for
  every past quarter, so revisions are ignored. It is one-step-ahead only, and
  quarters before 1990 are not evaluated.
- **Serial correlation.** It makes the iid SEs of differences too small. With
  Newey-West SEs the SV-versus-constant-variance margins fall from z of 3.4-5.6
  to 2.2-3.6 (1.7 against m7), and m10 - m8 from 3.1 to 1.85.
- **The comparison does not test the variance decomposition.** LFO scores the
  observables. Models tied on GDP-only prediction can imply very different
  latent paths: m9 - m10 is +1.2 (SE 11.9), yet over the same 1978-2026
  quarters m10's latent variance ratio is 0.35 [0.26, 0.49] and m9's 0.90
  [0.83, 0.97]. "About a fifth of headline variance is measurement error" is
  m9's answer, not a result the comparison validates.
- **Refits below the quality rule** *(superseded 2026-10-17: under the
  tightened rule all 75 stored fits, 12 initial fits and 63 refits, pass; the
  m9 refit to 2020 Q4 was retried and its attempt 4 kept)*. Two refits were kept below it:
  - m9 on data to 2020 Q4: tau R-hat 1.053, bulk ESS 58. It covers 20
    predictions, 2021 Q1 - 2025 Q4.
  - m10 on data to 2025 Q1: tau R-hat 1.073. It covers 5 predictions.
- **m3.** Every m3 refit mixes poorly (E-BFMI 0.034-0.062), so its score is
  less reliable. It is level with m1, m2 and m6 in any case. *(Superseded:
  the rewritten m3's refits have E-BFMI 0.93-1.01; its score is -1195.93.)*
- **m10.** m10 is multimodal (see below), and its particle-filter Monte Carlo
  error is the largest (about 1.6 nats). *(Superseded: the old score was
  also biased low by roughly 15-20 nats per draw over COVID, because the
  fit sat in mode A; the re-specified m10 has MCSE 0.08.)*
- **Shorter refits.** SV refits use 1000 + 1000 iterations per chain, not the
  production 2000 + 2000.

## Why the earlier PSIS-LOO comparison was invalid (note added 2026-10-16)

Earlier versions of this file ranked the models by PSIS-LOO on the saved
pointwise `log_lik[t]`. That comparison claimed:
- m9 was best, +6.3 over m4;
- SV gained about 240 elpd;
- m10 beat m8 by 98.

None of these numbers is evidence, for three reasons:

1. **The latent states were not integrated out.** Every spec has a
   per-quarter latent mu_t, and the SV specs also have a per-quarter h_t.
   `log_lik[t]` = log p(y_t | mu_t, h_t, theta) conditions on states that
   were themselves fitted to the left-out observation y_t. Leaving y_t out
   then requires a large change in the posterior of mu_t (and h_t), which
   importance sampling cannot capture.
2. **p_loo was close to T.** The old objects had p_loo of 208-336 for
   T = 265. Recomputed on the current fits, p_loo is 205-390 for T = 267
   (m1-m7, m9), 159 and 316 for T = 191 (m8, m10), and 144 for m11's
   T_q = 193. Each quarter effectively has its own parameter.
3. **Most Pareto k were above 0.7.** On the current fits, 187-242 of 267 k
   values exceed 0.7 for m1-m7 and m9, 155 and 149 of 191 for m8 and m10, and
   131 of 193 for m11. The PSIS estimates are unreliable for most
   observations.

The LFO comparison above fixes all three problems: it scores
p(y_t | y_{1:t-1}) with the states integrated out by a filter. The old LOO
objects (`output/loos*.rds`, `output/loo_compare.rds`) were deleted once the
LFO results existed.

## Headline replication result (per-draw variance ratios, 2026-10-16; m3, m10, m12 updated 2026-10-17)

For each posterior draw s, the ratio var_t(mu_t^(s)) / var_t(headline) is
computed. The headline is the mean of E, I and P. The table shows the median
and 5-95% quantiles across draws (`output/tables/variance_decomposition.csv`).

Full sample (T = 267, 1959 Q4 - 2026 Q2; headline sd 4.54):

| Series | var / var(headline): median | 5% | 95% |
|---|---|---|---|
| GDP(E)              | 1.75  | - | - |
| GDP(I)              | 1.26  | - | - |
| GDP(P)              | 1.09  | - | - |
| Latent mu (m1)      | 0.786 | 0.708 | 0.864 |
| Latent mu (m2)      | 0.785 | 0.712 | 0.861 |
| Latent mu (m3)      | 0.831 | 0.593 | 0.924 |
| Latent mu (m3) *(superseded, mu-sampled fit)* | 0.826 | 0.625 | 0.921 |
| Latent mu (m4)      | 0.798 | 0.727 | 0.879 |
| Latent mu (m5)      | 0.783 | 0.710 | 0.857 |
| Latent mu (m6)      | 0.779 | 0.696 | 0.861 |
| Latent mu (m7, COVID scale) | **0.614** | 0.544 | 0.685 |
| Latent mu (m9, production until 2026-10-17)  | 0.785 | 0.709 | 0.866 |
| Latent mu (m12, production from 2026-10-17) | 0.779 | 0.707 | 0.850 |

Hours subsample (T = 191, 1978 Q4 - 2026 Q2; headline sd 3.95):

| Series | median | 5% | 95% |
|---|---|---|---|
| Latent mu (m8)  | 0.848 | 0.767 | 0.930 |
| Latent mu (m10, re-specified with latent SV) | 0.747 | 0.678 | 0.824 |
| Latent mu (m10) *(superseded: old spec, mode A)* | 0.347 | 0.263 | 0.487 |
| Latent mu (m9), same quarters  | 0.899 | 0.829 | 0.968 |
| Latent mu (m12), same quarters | 0.909 | 0.840 | 0.973 |

The m9 and m12 rows are their full-sample fits restricted to the
1978 Q4 - 2026 Q2 quarters (computed 2026-10-17, not in the CSV). Excluding
2020 Q1 - 2021 Q1 the medians are m10 0.558, m9 0.861, m12 0.875.

Mixed-frequency quarters (T_q = 193, 1978 - 2026; headline sd 3.93):

| Series | median | 5% | 95% |
|---|---|---|---|
| Latent mu_q (m11) | 0.846 | 0.764 | 0.929 |

**What changed.** Earlier versions of this table reported the variance of the
posterior-*mean* path: about 0.65 for most models, 0.48 for m7 and 0.68 for
m9. They concluded that latent growth is "30-35% less variable than
headline". The posterior-mean path is smoother than any single plausible path,
so that statistic overstates the reduction. Per draw, the latent series has
about **20% less variance than headline** (m9: 0.785 [0.709, 0.866]; m12:
0.779 [0.707, 0.850]) and about 55% less than GDP(E).

This still supports the paper's qualitative claim that part of the
quarter-to-quarter volatility in Australian GDP growth is measurement error.
The share, however, is about a fifth of headline variance, not a third.

## Identification issue in m3

In `m3_corr.stan` the model has weak identification between two things:
(i) movements in mu_t, and (ii) positively correlated measurement errors
across the approaches. A common shift to all three observed series at time t
can be explained either way.

- **LKJ(2), the earlier run.** The sampler produced about 2% divergences and
  an E-BFMI of 0.02.
- **LKJ(4), the current code.** With adapt_delta = 0.99, the current-data fit
  has 0 divergences, but E-BFMI is still 0.03 and R-hat reaches 1.029
  (Omega[3,2]), with a minimum bulk ESS of 132. Every LFO refit also has
  E-BFMI below 0.07.
- **Posterior.** The correlations all straddle zero: Omega[2,1] = 0.07,
  Omega[3,1] = 0.05 and Omega[3,2] = -0.22, with all 90% intervals including 0.

m3 should not be relied on. Two fixes are possible (neither implemented):
- a structural restriction, which is what m6 does;
- marginalising mu, as m9 does.

**Update 2026-10-17.** mu is now marginalised. The bullets above describe the
superseded mu-sampled fit. The new fit samples cleanly (E-BFMI 0.89-1.01,
R-hat 1.007), the posterior is unchanged, and the weak identification of
Omega is now shown to be a property of the model (see "Round of 2026-10-17",
"m3: rewrite and identification").

## Identification issue in m4

The initial random-walk SV parameterisation was unidentified, with R-hat about
1.6 across the board. It had h_t = h_{t-1} + omega*xi_t, with h_1 anchored at
2 log(sigma).

The model was switched to a stationary AR(1) on log-variance
(Kim-Shephard-Chib). With mu sampled jointly with h, that left phi[3]
(Production) looking bimodal, with R-hat 1.30.

**Update 2026-10.** m4 is now sampled with mu integrated out, the same
Rao-Blackwellisation as m9 (see the next section). In the current fit phi[3]
mixes cleanly:
- R-hat at most 1.002 on every hyperparameter;
- the four chain medians of phi[3] are 0.968-0.970;
- 5 divergences in 8000 draws.

The earlier "bimodality" was a sampler artefact of the mu-h funnel, not a
feature of the posterior.

## m9 reparameterisation and production diagnostics (2026-10; m9 was production until 2026-10-17)

**The problem.** With mu sampled jointly with h, m9 got worse as T grew. The
last CI fit with that parameterisation (published 2026-10-14, T = 267,
`docs/diagnostics.json` at that commit) had 44 divergences in 8000 draws and a
minimum E-BFMI of 0.087.

**The fix.** It changes coordinates only. The observation equation, the
structural zero, the AR(1) latent, the stationary SV and every prior are
unchanged.
- **mu integrated out.** Given h and the hyperparameters, m9 is
  linear-Gaussian in mu. The model block integrates mu out with a scalar
  Kalman filter, so HMC explores only (theta, h). Generated quantities draw mu
  exactly from p(mu | y, h, theta) by forward-filtering backward-sampling
  (FFBS).
- **Checks.** The filter and FFBS were checked against the dense
  multivariate-normal likelihood and smoother.
- **rho_EI.** It is sampled as z_EI = atanh(rho_EI), with the Jacobian added
  explicitly. 1 +- rho_EI and the filter's quadratic form are computed
  stably from z_EI. With Stan's built-in (-1, 1) bound, a few chains in a
  stall test collapsed their step size in early warmup at rho_EI = -1. The
  counts are recorded in the header of `stan/m9_struct_sv.stan`.
- **h** stays non-centred.
- **Data guard.** The data block declares `int<lower=3, upper=3> K`. m6 has
  the same guard.

m4 and m10 use the same Kalman / FFBS scheme, with a diagonal measurement
covariance.

**Production fit on the current data** (`R/fit_production.R`; seed 1;
4 x (2000 + 2000); adapt_delta 0.99; max_treedepth 12):

| Check | Value | Gate |
|---|---|---|
| divergences | 1 of 8000 (0.013%) | < 0.1% |
| E-BFMI by chain | 0.920, 0.942, 0.950, 0.900 | min >= 0.3 |
| max R-hat, hyperparameters | 1.0024 (phi[3]) | <= 1.01 |
| max R-hat, mu[t] / h[t,k] | 1.0022 / 1.0025 | <= 1.05 |
| min bulk / tail ESS, hyperparameters | 1692 (omega[2]) / 2807 | >= 400 |
| min bulk / tail ESS, mu[t] | 2872 / 5934 | >= 400 |
| max-treedepth hits | 0 | - |

The fit takes about 8-9 minutes on 4 local cores. A refit with seed 2 also
passed the gate; that fit is not kept in `output/`.

Posterior means:
- rho = 0.06, tau = 4.01, mu_bar = 3.08;
- rho_EI = 0.14, with 90% interval [-0.02, 0.29];
- phi = (0.98, 0.94, 0.97), omega = (0.23, 0.45, 0.46).

## m10 multimodality (superseded 2026-10-17)

*This describes the old constant-tau m10. A 16-chain diagnosis confirmed the
two modes, traced them to the constant latent sd, and m10 was re-specified
with latent SV; see "Round of 2026-10-17", "m10: two modes, diagnosis and
re-specification".*

The m10 posterior has at least two modes:
- **Mode A:** rho about 0.40, tau about 2.1, phi_E about 0.65-0.70.
- **Mode B:** rho about 0.12, tau about 3.7, phi_E about 0.97.

Short runs put chains in different modes. The two-chain smoke run that showed
this is recorded in the header of `stan/m10_labour_sv.stan`.

In the current full fit, all four chains found mode A:
- chain medians: rho 0.39-0.40, tau 2.13-2.16;
- R-hat 1.011 on phi[1];
- minimum tail ESS 88.

The earlier (mu-sampled) version had R-hat 1.4-1.5 on rho, tau, phi[1] and
omega. The post-COVID LFO refits needed retries, and one (data to 2025 Q1)
could not get tau R-hat under 1.05.

m10's very low variance ratio (0.35) comes from mode A's small tau.
Re-specifying the model or making an identification decision would be needed
to fix this (on beta, tau or the SV priors). Reparameterising cannot remove
genuine multimodality.

## Substantive surprises

1. **ABS data quality has improved enormously over 65 years.** The 1960s and
   1970s had very large measurement-error sd. In m9 the E-side sd peaks at
   about 8.2 pp annualised in 1973 Q3. By 2026 Q2 the sd is about 1.4 (E),
   0.6 (I) and 0.4 (P). The measurement-error variance falls over the sample,
   and that fall is exactly where the SV models gain in the LFO comparison.
2. **The SV models attribute COVID to genuine movement in mu_t.** They do not
   treat it as inflated measurement error. In m9 the median measurement-error
   sd barely moves through 2020:
   - E is 1.38-1.43 from 2020 Q1 to 2021 Q1, against 1.39 in 2019 Q4;
   - the median mu in 2020 Q2 is -27.3, against a headline of -28.1.

   That is the right call, because the three measures agreed to within 0.5 pp.
   It costs the SV models in one-step prediction, though: m9 loses about 10
   nats to m5 at 2020 Q2, because an AR(1) latent with tau about 4 does not
   foresee a 7-sigma quarter. In m4 the 1970s and the mid-1990s (I-side sd
   peaks at 5.7 in 1994 Q4) are much bigger spikes in measurement variance
   than COVID.
3. **m7 disagrees with the SV models about COVID.** m7 is free to inflate the
   variance for just the four pandemic quarters.
   - Only 2020 Q2 needs a large scale-up: kappa = 8.4, 90% interval
     [4.7, 14.3].
   - 2020 Q3 and Q4 have posterior-mean kappa below 1 (0.80, 0.86).
   - 2021 Q1 is about 1 (1.10, interval [0.47, 2.30]).

   m7 therefore treats part of the Q2 collapse as noise, which smooths its
   latent path (variance ratio 0.61 against about 0.78 for the others). In the
   LFO comparison this gives m7 a large gain in 2020, but it relies on knowing
   the COVID quarters in advance.
4. **m6's restriction is well identified but not supported by prediction.**
   m6 restricts Cov(eps_P, .) = 0, which is well identified: rho_EI = 0.09,
   90% interval [-0.06, 0.23]. The data are consistent with the
   construction-based assumption that Production has independent errors. The
   restriction is predictively neutral: m6 against m2 is +0.66 (SE 0.34), and
   m9 against m4 is +0.18 (SE 1.59). It is principled, not data-driven.
   - *Superseded:* "m6 loses ~3 elpd to m1/m2 by LOO".
5. **Hours worked sharpen the Production signal.** Adding hours (m8, on the
   1978 Q4+ sample) gives sigma_P = 1.83, against 2.64 for m2 on the full
   sample. The samples differ, so part of that gap is the lower post-1978
   noise. m8 also has alpha = -0.76 and beta = 0.83. That implies
   productivity growth of about 1.3 pp/yr at the sample mean, in line with
   historical estimates.
6. **SV is the gain; the structural restriction adds no predictive accuracy.**
   SV is a consistent gain (+45 over the best constant-variance model, z
   about 2.2 with a Newey-West SE; ahead of every constant-variance model in every sub-period but 2020, except m5 in the 1990s). The
   structural restriction is predictively neutral, both alone (m6 against m2)
   and combined with SV (m9 against m4).
   - *Superseded:* "m9 (m6 + m4) wins another 6 elpd over m4" and "fixes m4's
     bimodal phi[3]". Both are withdrawn: the first came from invalid LOO, the
     second was a sampler artefact (see "Identification issue in m4").
7. *(Superseded 2026-10-17: the multimodality came from the constant latent
   sd, and m10 now has latent SV; see "Round of 2026-10-17".)*
   **m10 (4-signal + SV) is multimodal on T = 191.** Stacking extensions does
   not always help, because SV needs sample size to identify its T x K
   log-variance trajectory. Even so, m10 is ahead of m8 on prediction (+34.9
   joint, +36.1 GDP-only; z 1.85 and 2.2 with Newey-West SEs), so SV probably
   helps on the hours sample too.
8. **m11 concentrates the COVID collapse in April 2020.** Its monthly latent
   path puts the Q2 2020 collapse mostly into April 2020 (median about -43%
   annualised), with May at about -31%. The model uses monthly employment to
   identify when within the quarter the collapse happened. This is the
   value-add of mixed-frequency models for real-time and nowcasting use. On
   GDP-only prediction, m11 is level with m9 (+10.0, SE 16.6).

## Sampler settings

- **Full fits** (`R/02_fit_models.R`, `R/fit_m11.R`, `R/fit_production.R`):
  4 chains x (2000 warmup + 2000 sampling) = 8000 post-warmup draws per
  model, seed 1.
- **adapt_delta and max_treedepth** (m1-m10 as recorded in
  `output/tables/fit_diagnostics.csv`; m11 as set in `R/fit_m11.R`, which
  does not write to that table):
  - m1, m2: default adapt_delta (0.8), max_treedepth 10;
  - m3, m5, m6, m7, m8: adapt_delta 0.99, max_treedepth 10;
  - m4, m9, m10, m11, m12: adapt_delta 0.99, max_treedepth 12 (m12's are
    the production settings, as m9's were).
- **LFO refits:** 4 x (1000 + 1000), with each model's own adapt_delta and
  max_treedepth, and the refit rule described in "Round of 2026-10-17".
- **Diagnostics.** The research suite reports diagnostics and never stops on
  them. Only `R/fit_production.R` enforces a gate:
  - divergences < 0.1%;
  - min E-BFMI >= 0.3;
  - R-hat <= 1.01 on hyperparameters and <= 1.05 on mu[t] and h[t,k];
  - bulk and tail ESS >= 400 on hyperparameters and mu[t].

  The hyperparameters are c, rho, tau, mu_h, phi, omega, rho_EI, psi and
  mu_bar (psi added with the switch to m12).

## Files

**Data**
- `R/01_fetch_data.R`: pulls ABS 5206.0 Table 24 (production data).
- `R/01b_fetch_research_data.R`: pulls 5206.0 Table 1 hours and 6202.0
  employment (research only; run after 01).
- `R/utils_abs.R`: retry wrapper and ABS connection-check workaround, shared
  by both fetch scripts.

**Research suite**
- `R/02_fit_models.R`: fits m1-m10 and m12 (`Rscript R/02_fit_models.R m9 m12`
  for a subset). Writes `output/tables/fit_diagnostics.csv`.
- `R/fit_m11.R`: fits m11.
- `R/fit_m12.R`: fits m12 with the production gate (any seed), and the
  warmup-stall harness.
- `R/03_compare.R`: PSIS-LFO-CV comparison and per-draw variance
  decomposition.
- `R/lib/filters.R`, `R/lib/filters.cpp`, `R/lib/lfo.R`: Kalman filter,
  Rao-Blackwellised particle filter and the LFO algorithm.
- `R/lib/filter_m12.R`, `R/lib/filter_m12.cpp`: m12's Kalman filter and
  particle filter (registered by `register_m12()` in `R/03_compare.R`).
- `R/lib/filter_m10.R`, `R/lib/filter_m10.cpp`: the re-specified m10's
  particle filter over (g, h) (registered by `register_m10_filter()`).
- `tests/test_filters.R`, `tests/test_m12.R`, `tests/test_m3_m10.R`: filter
  tests, each against a dense multivariate-normal likelihood.
- `R/diag_m10_modes.R`: m10 mode check (16 dispersed chains; 4 chains started
  in old mode A) -> `output/tables/m10_mode_check.csv`.
- `R/diag_m3_ident.R`: m3 identification checks ->
  `output/tables/m3_identification.csv`.
- `R/04_figures.R`: figures 01-07 in `output/figures/`.

**Production**
- `R/fit_production.R`, `R/output_estimate.R`, `R/build_site.R`: production
  fit (m12 since 2026-10-17; m9 before), estimate writer and site.
- `R/input_hash.R`: CI input fingerprint.

**Models**
- `stan/m{1..12}_*.stan`: the 12 specifications.
  - m1-m5: original replication of Rees-Lancaster-Finlay.
  - m6: structural Sigma restriction (a cheap novel extension).
  - m7: Lenza-Primiceri COVID variance scale.
  - m8: 4-signal (Rees-Lancaster-Finlay Model 4, adding hours worked).
  - m9: m6 + SV (production until 2026-10-17; tied with m4 on LFO).
  - m10: m8 + SV on the measurement errors and on latent growth
    (re-specified 2026-10-17; the old constant-tau version was multimodal).
  - m11: mixed-frequency, Banbura-Modugno style, with monthly employment.
  - m12: m9 + AR(1) measurement errors (production from 2026-10-17).

**Outputs**
- `output/fit_m{1..12}.rds`, `output/fit_production.rds`: saved CmdStanFit
  objects (gitignored). `output/fit_m12_gate.rds` and
  `output/fit_m12_seed2_gate.rds` are the seed-1 and seed-2 gate records
  (also gitignored; their content is in `output/tables/m12_gate_seeds.csv`).
- `output/lfo/`: LFO log (`lfo_log.txt`, committed) and checkpoints (`*.rds`,
  gitignored, regenerable by `R/03_compare.R`). `output/lfo/superseded/`
  holds the checkpoints of the 2026-10-16 run; they are local only and also
  gitignored. `R/03_compare.R` reads two of them (m4, m9) as a refit cache
  (`REFIT_CACHE`); without them the same refits are recomputed from their
  seeds, with identical draws.
- `output/tables/`: comparison, diagnostics and variance-decomposition tables.
- `output/figures/`: generated PNGs.

**Documentation**
- `priors.md`: every prior, as coded.
