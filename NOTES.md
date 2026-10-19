# NOTES

This is the research record. Sections dated 2026-10 replace earlier results
that turned out to be wrong; the superseded claims are noted where they were
made so the history stays readable. The newest section is "Production
sample moved to 1980 Q1, 2026-10-19": production m13 is now fitted on
1980 Q1 onward. Before it, "1980Q1 sample (RDP 2014-12 alignment)" re-ran
the whole research suite on data starting in 1980 Q1 (at the time leaving
production on the full sample). It does not supersede the full-sample
research results: those stay as recorded. It does say which
full-sample conclusions survive the later start. The round before it, "Round
of 2026-10-18", added m13 (m12 + latent SV), which is now the production
model. "Round of 2026-10-17" added m12 (production from 2026-10-17 to
2026-10-18) and superseded the 2026-10-16 comparison tables and the m3 and
m10 sections.

## Paper

RBA RDP 2014-12, "A State-space Approach to Australian GDP Measurement",
by Daniel Rees, David Lancaster & Richard Finlay.
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
table shows posterior means. The diagnostics for m1-m10, m12 and m13 come from
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
| m12 | 267 | 3.10 | 0.09 | 3.99 | -    | -    | -    | 3 | 0.86 | 1.002 | production 2026-10-17 to 2026-10-18; m9 + AR(1) errors; psi = (-0.47, -0.34, -0.16) |
| m13 | 267 | 2.77 | 0.27 | 2.54 | -    | -    | -    | 4 | 0.85 | 1.002 | production from 2026-10-18; m12 + latent SV; tau = exp(mu_g/2); psi = (-0.46, -0.36, -0.19) |
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

## Production sample moved to 1980 Q1, 2026-10-19

At the user's request the production estimate is now fitted on 1980 Q1 - 2026 Q2
(T = 186), the sample start of RBA RDP 2014-12, instead of 1959 Q4 - 2026 Q2
(T = 267). The model is unchanged (m13).

- **Why.** The paper's reason (section 3.2): the pattern of measurement errors in
  the early years of each GDP series may be unrepresentative of current
  performance. The prediction evidence is roughly neutral: on the same 146 LFO
  targets, m13 trained from 1980 scores 6.3 nats better than m13 trained from
  1959 (NW z 1.69), mostly in 1990-99, and +1.8 (z 0.7) from 2000 on. The
  comparison author had recommended keeping the full sample (section
  "Production: unchanged, and the recommendation" below, now superseded).
- **Code.** `R/production_config.R` defines `PRODUCTION_START = 1980-03-01` and
  `load_production_growth()`, used by `R/fit_production.R`,
  `R/output_estimate.R` and `R/input_hash.R`. The input hash now covers only
  the fitted quarters plus the config file. The research suite still defaults
  to the full sample (`LATENTGDP_START` unset).
- **First production fit on 1980 Q1 (local, seed 1).** Gate passed: 0 of 8000
  divergences, E-BFMI 0.882-0.948, max R-hat 1.0013 (hyperparameters), 1.0013
  (mu), 1.0015 (h), 1.0017 (g), min bulk/tail ESS 1315/882 (hyperparameters)
  and 5531/6268 (mu); 18 min wall. The draws reproduce the research fit
  `output/sample_1980/fit_m13.rds` exactly.
- **Local vs published run.** The figures above are from the local seed-1
  fit (`output/fit_production.rds`), which reproduces the 1980 research fit
  draw for draw. The published `docs/` come from the forced CI run on
  2026-10-19 (commit d623f3d), which has different draws: 2 divergences,
  min E-BFMI 0.846, max R-hat 1.003 on hyperparameters, min bulk / tail ESS
  1285 / 915. Its median latent path differs from the local fit's by at most
  0.017 pp quarter-on-quarter.
- **Published numbers.** 2026 Q2 latent growth +0.48% QoQ [0.31, 0.64] (full
  sample: +0.48% [0.30, 0.63]). Variance ratio 0.87 [0.81, 0.94] over
  1980 Q1 - 2026 Q2 (full-sample fit: 0.79 [0.70, 0.89] over 1959 Q4 - 2026 Q2).
  Latent volatility tau_t by decade (annualised, posterior median of decade
  means): 1980s 3.1, 1990s 2.4, 2000s 1.5, 2010s 1.3, 2020s 4.0; 2020 Q2 13.7.
- **Site.** The historical chart now starts in 1980 Q1. The "Why this
  specification" text notes that the LFO comparison starts every model in
  1959 Q4 and that on the 1980 sample m13's lead over m12 is weaker
  (+57.4, z 1.6).

## 1980Q1 sample (RDP 2014-12 alignment), 2026-10-18/19

The whole research suite was re-run on data starting in 1980 Q1:
- every model's fit (m1-m13);
- the m13 production gate on two seeds;
- the per-draw variance decomposition;
- PSIS-LFO-CV.

Fits ran on 2026-10-18 and the LFO from 21:16 on 2026-10-18 to 04:21 on
2026-10-19. The data vintage is the 2026 Q2 release used everywhere else.
**Production was unchanged at the time** (it fitted m13 on the full sample,
1959 Q4 - 2026 Q2; *superseded 2026-10-19*, see "Production sample moved to
1980 Q1"). Every 1980-sample output is under `output/sample_1980/`, and no
full-sample output (`output/fit_*.rds`, `output/lfo/`, `output/tables/`,
`output/figures/`, `output/fit_production.rds`, `docs/`) was changed. All
numbers below come from `output/sample_1980/tables/` unless another source
is named.

### Why 1980 Q1

RBA RDP 2014-12 (Rees, Lancaster & Finlay 2014, "A State-space Approach to
Australian GDP Measurement"), section
3.2:

> "Our data span 1980:Q1-2013:Q2. The starting date reflects the fact that,
> while Australian national accounts data are available on a quarterly basis
> from 1959:Q3, the quality of the underlying data sources has changed over
> time, so that the pattern of measurement errors in the early years of each
> GDP series may be unrepresentative of their current performance."

The paper uses that one sample for all its models. Our full-sample runs
start in 1959 Q4, so their measurement-error parameters are partly estimated
from the period the paper deliberately leaves out. Only the start date is
aligned with the paper. Two differences remain:
- **The fourth signal.** The RDP's fourth signal is the **unemployment
  rate**. Our m8 and m10 use **hours worked** (ABS 5206.0 Table 1). This was
  not changed. m11 (monthly employment, 6202.0) is our own extension.
- **The end date.** The paper's data end in 2013 Q2 (its vintage). Ours end
  in 2026 Q2 (the current vintage).

### Definitions

- **Sample.** The first growth observation is 1980 Q1 (dated 1980-03-01),
  400 × (log level 1980 Q1 - log level 1979 Q4). The last is 2026 Q2. T = 186
  for every model.
  - The 3-measure models lose 1959 Q4 - 1979 Q4 (81 quarters; T 267 -> 186).
  - m8 and m10 lose 1978 Q4 - 1979 Q4 (5 quarters; T 191 -> 186).
  - m11 uses months from 1980-01, so 1980 Q1 is complete: T_m = 558 months,
    T_q = 186 quarters.
  - In this sample 2020 Q2 is t = 162.
- **Data on this sample** (annualised %):
  - means: E 2.87, I 2.91, P 2.91;
  - sd: E 4.46, I 4.66, P 4.09; headline sd 3.88 (full sample 4.54);
  - pairwise correlations 0.64-0.71.

  GDP(I), not GDP(E), is now the most variable measure in raw terms: its
  variance is 1.44 times the headline's, against 1.32 for E and 1.11 for P.
  On the full sample the ratios are 1.75, 1.26 and 1.09.
- **LFO.**
  - The evaluation window is unchanged: 1990 Q1 - 2026 Q2, 146 one-step
    predictions, the same targets as the full-sample run.
  - Initial training is the 40 quarters 1980 Q1 - 1989 Q4 (full sample:
    121 quarters, 1959 Q4 - 1989 Q4).
  - The method is unchanged:
    - the refit rule: divergences <= 1%, R-hat <= 1.05, E-BFMI >= 0.2, bulk
      ESS >= 100;
    - a refit when Pareto k exceeds 0.7;
    - up to 4 attempts;
    - the particle schedule of 2000 to 100000 particles;
    - 4 filter passes for each SV model;
    - a Newey-West lag of 8;
    - the retry (seed-selection) sensitivity.
  - The refit cache (`REFIT_CACHE`, which holds full-sample fits) is
    disabled whenever `LATENTGDP_START` is set.
  - In this run "1:i" counts rows of the 1980 sample: 1:146 is data to
    2016 Q2, not the full-sample 1:146 (data to 1996 Q2).

### How to reproduce

`R/sample_window.R` reads two optional environment variables, which the
research scripts (`R/02_fit_models.R`, `R/fit_m11.R`, `R/fit_m12.R`,
`R/fit_m13.R`, `R/03_compare.R`, `R/04_figures.R`) use:
- `LATENTGDP_START` is the date of the first growth observation
  (`1980-03-01` = 1980 Q1);
- `LATENTGDP_OUT` is the output root.

With both unset the scripts behave exactly as before. That was checked: the
Stan data lists are `identical()` and serialize-identical to the old code
path. A default `R/03_compare.R summary` run in a sandbox regenerated the 8
tables that summary mode writes (`lfo_fit_diagnostics`, `lfo_pointwise`,
`lfo_refit_attempts`, `lfo_retry_sensitivity_summary`, `model_comparison`,
`model_comparison_pairwise`, `model_comparison_subperiod`,
`variance_decomposition`) byte-identical. The other 8 files of
`output/tables/` are written by other modes or scripts (`fit_diagnostics`,
the two particle-filter checks, `lfo_retry_sensitivity`, `m10_mode_check`,
`m12_gate_seeds`, `m13_gate_seeds`, `m3_identification`). That run did not
reproduce them; it only left them unchanged.

**Guard.** `R/sample_window.R` stops as soon as it is sourced if
`LATENTGDP_START` is set and the output root is `output/`. That covers
`LATENTGDP_OUT` unset or set to `output`, `output/`, `./output` or the
absolute path. Every research script that uses it (`R/02_fit_models.R`,
`R/fit_m11.R`, `R/fit_m12.R`, `R/fit_m13.R`, `R/03_compare.R`,
`R/04_figures.R`) therefore refuses to write a truncated-sample fit or
table into `output/`, and it stops before creating any directory.
- **Before 2026-10-19 only `R/03_compare.R` had this check.** In a sandbox,
  `LATENTGDP_START=1980-03-01 Rscript R/02_fit_models.R m1` without
  `LATENTGDP_OUT` wrote `output/fit_m1.rds` with T = 186 and replaced m1's
  row in `output/tables/fit_diagnostics.csv`.
- **That never happened in the real tree.** The 87 full-sample files under
  `output/` (fits, `lfo/`, `tables/`, `figures/`) and the production files
  still match their checksums from before the 1980 run.
- **Tested in a sandbox.** All six scripts exit with the error for all five
  spellings of the root, and no file is written. With both variables unset,
  `sw_out()`, `sw_is_default()` and the trimmed data are unchanged. A
  separate root without `LATENTGDP_START` is still allowed.

The production scripts never source `R/sample_window.R`.

```bash
cd /Users/davidstephan/Documents/latentgdp/LatentGDP
# Fits, gate, variance decomposition, figures 01-06 (sequential, resumable;
# logs in output/sample_1980/logs/):
nohup bash output/sample_1980/run_1980.sh >> output/sample_1980/logs/driver.log 2>&1 &
# LFO (sequential, resumable from checkpoints in output/sample_1980/lfo/):
nohup bash output/sample_1980/lfo/run_lfo_1980.sh >> output/sample_1980/lfo/driver.log 2>&1 &
bash output/sample_1980/lfo/run_lfo_1980_extra.sh     # COVID particle checks for m9, m4
# or step by step:
export LATENTGDP_START=1980-03-01 LATENTGDP_OUT=output/sample_1980
Rscript R/02_fit_models.R m1 m2 m3 m4 m5 m6 m7 m8 m9 m10 m12
Rscript R/fit_m13.R && Rscript R/fit_m13.R --seed=2 && Rscript R/02_fit_models.R --diag-only m13
Rscript R/fit_m11.R && Rscript R/02_fit_models.R --diag-only m11
Rscript R/02_fit_models.R --vardecomp
Rscript R/03_compare.R m1 ...; Rscript R/03_compare.R summary   # see run_lfo_1980.sh for the order
Rscript R/04_figures.R
unset LATENTGDP_START LATENTGDP_OUT        # back to the full sample in output/
```

Compute: the 13 fits' wall times sum to about 42 minutes, or about 55 with
m13's seed-2 fit (809 s), excluding compilation. m13 (1042 s) and m12
(629 s) are the longest. The LFO took about 7 h 05 min, including all
re-runs and particle checks.

### Fits and diagnostics

Every fit uses 4 chains × (2000 warmup + 2000 draws) and seed 1
(`tables/fit_diagnostics.csv`). The columns are:
- **div**: divergent transitions out of 8000;
- **E-BFMI**: the minimum over chains;
- **R-hat**: the maximum over hyperparameters, with the variable that has it;
- **ESS**: the minimum bulk / tail ESS over hyperparameters;
- **s**: wall seconds (CmdStan's sampling total for m13).

| model | adapt_delta / treedepth | div | E-BFMI | R-hat | ESS | s |
|---|---|---|---|---|---|---|
| m1 | 0.8 / 10 | 0 | 0.761 | 1.0004 (sigma[1]) | 4141 / 4635 | 3 |
| m2 | 0.8 / 10 | 0 | 0.543 | 1.0008 (tau) | 2145 / 3255 | 5 |
| m3 | 0.99 / 10 | 0 | 0.949 | 1.0021 (Omega[3,1]) | 2353 / 2841 | 29 |
| m4 | 0.99 / 12 | **9** | 0.830 | 1.0015 (phi[2]) | 1742 / 2497 | 140 |
| m5 | 0.99 / 10 | 0 | 0.557 | 1.0020 (tau) | 1675 / 2599 | 28 |
| m6 | 0.99 / 10 | 0 | 0.435 | 1.0048 (rho) | 1214 / 1918 | 53 |
| m7 | 0.99 / 10 | 0 | 0.476 | 1.0027 (tau) | 1558 / 3207 | 16 |
| m8 | 0.99 / 10 | 0 | 0.565 | 1.0030 (rho) | 1370 / 2401 | 11 |
| m9 | 0.99 / 12 | 1 | 0.848 | 1.0063 (phi[1]) | 1273 / 2609 | 196 |
| m10 | 0.99 / 12 | 2 | 0.878 | 1.0033 (omega[3]) | 1658 / 749 | 298 |
| m11 | 0.99 / 12 | 0 | **0.380** | 1.0007 (c) | 1312 / 2055 | 50 |
| m12 | 0.99 / 12 | 1 | 0.845 | 1.0038 (phi[1]) | 513 / 1687 | 629 |
| m13 | 0.99 / 12 | 0 | 0.882 | 1.0013 (rho_EI) | 1315 / 882 | 1042 |

- **mu and treedepth.** Max R-hat on mu is at most 1.0033 in every model,
  and min bulk ESS on mu is at least 2484 (m3). There are no max-treedepth
  hits.
- **m4** has 9 divergences out of 8000 (0.11%; 5 on the full sample). That
  is above the production gate's 0.1%, but m4 is not gated.
- **m11**'s E-BFMI of 0.38 (0.45 on the full sample) is above 0.3 but the
  lowest of the suite.
- **m12**'s min bulk ESS of 513 (phi[1]) is the lowest; it is above the
  gate's 400.

### m13 gate on the 1980 sample

`R/fit_m13.R` applies exactly `R/fit_production.R`'s gate. The code was
compared line by line: the same thresholds, NA handling and E-BFMI floor
(`tables/m13_gate_seeds.csv`). **Both seeds pass.**

| Check | seed 1 | seed 2 | Gate |
|---|---|---|---|
| divergences | 0 of 8000 | 1 of 8000 (0.0125%) | < 0.1% |
| E-BFMI by chain | 0.948, 0.926, 0.900, 0.882 | 0.903, 0.957, 0.862, 0.897 | min >= 0.3 |
| max R-hat, hyperparameters | 1.0013 (rho_EI) | 1.0056 (phi[1]) | <= 1.01 |
| max R-hat, mu / h / g | 1.0013 / 1.0015 / 1.0017 | 1.0009 / 1.0031 / 1.0019 | <= 1.05 |
| min bulk / tail ESS, hyperparameters | 1315 / 882 (phi[1]) | 1094 / 782 (phi[1]) | >= 400 |
| min bulk / tail ESS, mu | 5531 / 6268 | 5059 / 6040 | >= 400 |
| max-treedepth hits | 0 | 0 | - |
| sampling time (4 parallel chains) | 1042 s (wall 1048 s; chains 739-1042 s) | 802 s (wall 809 s; chains 621-802 s) | <= 30 min |

### Posterior: what changes with the later start

Posterior means, full sample -> 1980 sample
(`tables/posterior_comparison_full_vs_1980.csv`; m13's tau_t path is from
`logs/compare_full_vs_1980.log`).

**1. Constant-variance models (m1, m2, m3, m5, m6, m7).** The GDP(E)
measurement-error sd falls sharply. That is what the paper's rationale
predicts: the early decades' measurement errors inflated it.

| | m1 | m2 | m3 | m5 (innovation sd) | m6 | m7 |
|---|---|---|---|---|---|---|
| sig_E | 4.08 -> 2.53 | 4.10 -> 2.55 | 4.35 -> 2.62 | 3.53 -> 2.27 | 4.20 -> 2.63 | 4.14 -> 2.62 |
| sig_I | 3.25 -> 3.00 | 3.26 -> 3.02 | 3.13 -> 3.00 | 2.98 -> 2.76 | 3.36 -> 3.10 | 3.26 -> 3.05 |
| sig_P | 2.65 -> 2.07 | 2.64 -> 2.04 | 2.43 -> 1.85 | 2.52 -> 2.02 | 2.56 -> 1.96 | 2.64 -> 1.98 |
| tau | 4.01 -> 3.57 | 4.02 -> 3.56 | 4.04 -> 3.61 | 4.00 -> 3.55 | 3.99 -> 3.56 | 3.54 -> 2.62 |
| rho | - | 0.04 -> 0.14 | 0.06 -> 0.14 | 0.08 -> 0.16 | 0.05 -> 0.15 | 0.10 -> 0.34 |
| mu_bar | 3.23 -> 2.90 | 3.12 -> 2.82 | 3.13 -> 2.82 | 3.13 -> 2.81 | 3.13 -> 2.83 | 3.26 -> 3.02 |

- **The noisiest measure.** Income, not Expenditure, now has the largest
  error sd, as on the post-1978 samples (m8, m11). Production is still the
  cleanest signal in every specification.
- **m5.** psi goes from (-0.51, -0.43, -0.18) to (-0.44, -0.40, -0.17).
- **m7.** kappa_1 (2020 Q2) rises from 8.4 to 11.3.
- **m6.** rho_EI goes from 0.09 to 0.08.

**2. Measurement SV with a constant latent sd (m4, m9, m12): COVID becomes
common measurement error.** On the 1980 sample these three fits put the
2020 Q2 collapse into the measurement errors:
- **2020 Q2.** The posterior mean of mu in 2020 Q2 is +0.20 (m4), +0.25
  (m9) and -0.99 (m12), with 90% intervals about [-4, +4] (m12
  [-5.0, 2.8]). The full-sample fits put it near the headline -28: the
  full-sample and 1980 mean paths differ by up to 27.4 pp (m4, m9) and
  26.2 pp (m12) on the common quarters, and the paths' correlation is only
  0.77-0.80. The three measurement sds in that quarter are 9-13 (posterior
  means).
- **Parameters.** tau falls from about 4.0 to 2.05 / 2.03 / 2.16, and rho
  rises from 0.05 / 0.06 / 0.09 to 0.40 / 0.41 / 0.43.
- **m12.** phi_E falls from 0.98 to 0.01 (90% [-0.65, 0.80]), omega_E rises
  from 0.24 to 0.80, rho_EI goes from 0.175 to 0.27 and psi_E from -0.47 to
  -0.37. m4's and m9's phi_E fall to 0.58 and 0.55.
- **phi_E is prior-dominated on this sample** (added 2026-10-19; from the
  `fits$dm` draws in `lfo/lfo_{m4,m9,m12,m13}.rds` and the final fits).
  The prior on phi is Uniform(-1, 1), from its bounds: mean 0, 90%
  [-0.9, 0.9].
  - **Up to 2020 Q2.** In every 1980 LFO refit of m4, m9, m12 and m13
    trained up to 2020 Q2 (i = 40 to 162), the posterior mean of phi[1] is
    between -0.28 and +0.12, with 90% intervals of about [-0.8, 0.8]. That
    is essentially the prior. Before 2020 omega_E is also small (0.20-0.30),
    so GDP(E)'s log-variance barely moves and its persistence is not
    identified.
  - **When 2020 Q2 enters the training data** (i = 162), omega_E jumps to
    0.64-0.78 in m4, m9 and m12, while phi_E stays near 0. E's "SV" then
    acts as i.i.d. heavy-tailed noise, not persistent volatility. That is
    the mechanism that lets these fits absorb 2020 Q2 as common measurement
    error. m13's omega_E stays at 0.29-0.33.
  - **Later refits.** phi[1] rises to 0.41-0.53 only in the later refits
    (m9 1:168 and 1:185, m4 1:181, m13 1:177). m12's last refit (1:168) is
    still -0.15.
  - **Full sample.** Every LFO refit has phi[1] 0.70-0.98.
  - **Final fits, phi_E mean [90%] / omega_E, full -> 1980:**
    - m4: 0.98 / 0.24 -> 0.58 [0.12, 0.86] / 0.67;
    - m9: 0.98 / 0.23 -> 0.55 [0.04, 0.86] / 0.63;
    - m12: 0.98 / 0.24 -> 0.015 [-0.65, 0.80] / 0.80;
    - m13: 0.98 / 0.25 -> 0.84 [0.11, 0.99] / 0.27.

    m13's phi_E has a long left tail, and phi[1] has the lowest bulk / tail
    ESS of m13's hyperparameters on both gate seeds (1315 / 882 and
    1094 / 782).
  - **What this means for the comparisons.** The 1980-sample measurement-SV
    comparisons (m4 - m2, m9 - m6, m12 - m5) mostly measure the value of a
    heavier-tailed, roughly i.i.d. GDP(E) error plus persistent SV in I and
    P (phi_I and phi_P about 0.9 in the 1980 fits of m4 and m12). They do not measure persistent GDP(E)
    volatility.
- **Not a stranded mode.** This is the same failure as the old constant-tau
  m10 ("mode A", below). For each of m9, m4 and m12, four extra chains were
  run on the 1980 data. They started from full-sample posterior draws, where
  COVID is in mu (mu in 2020 Q2 about -25 to -28, tau about 4), with h_raw
  sliced to the 1980 quarters. The sampler settings were the production
  ones (adapt_delta 0.99, max_treedepth 12, 2000 + 2000, seed 11). In all
  three models all 4 chains moved to the measurement-error solution:

  | model | chains' mu (2020 Q2) | tau | lp__ by chain | reference 1980 fit lp__, tau, mu | div | min E-BFMI |
  |---|---|---|---|---|---|---|
  | m9 | 0.2 to 0.3 | 2.02-2.04 | -1664 to -1666 | -1665.4, 2.03, 0.25 | 5 | 0.83 |
  | m4 | 0.2 to 0.3 | 2.04-2.06 | -1665 to -1667 | -1665.9, 2.05, 0.20 | 0 | 0.80 |
  | m12 | -0.9 to -1.0 | 2.15-2.18 | -1643 to -1644 | -1643.3, 2.16, -0.99 | 1 | 0.81 |

  - **Warmup draws** were saved for m4 and m12. Every chain left the
    COVID-in-mu region within the first 8-52 warmup iterations (first
    iteration with mu in 2020 Q2 above -10: m4 11, 14, 8, 52; m12 21, 10,
    10, 34). The chains never came back. So that configuration is not a
    competing local mode at this posterior. It is a region the sampler
    leaves at once.
  - **Mixing.** m12's phi[1] again mixes worst (bulk ESS 409, R-hat 1.01;
    reference fit 513).
  - **Limits.** This is still one alternative starting region per model,
    and it is evidence about the sampler, not a proof that no other mode
    exists. The m9 check (`wf4_scratch/fits/m9_mode_check.R`) and the m4 /
    m12 checks (`wf4_scratch/repair/mode_check.R`, logs `mode_m4.log`,
    `mode_m12.log`) are in the session scratch directory, not the repo.
- **Checked here.** The 2020 Q2 values were recomputed from
  `output/sample_1980/fit_{m4,m9,m12,m13,m10}.rds` for this note.

**3. m13 (production model): COVID stays latent.**
- **2020 Q2.** mu in 2020 Q2 is -27.9 [-29.0, -26.8] (full sample -27.9).
  The measurement sds in that quarter are 0.9-1.6.
- **The path.** The mean path has correlation 0.9985 with the full-sample
  path, a mean absolute difference of 0.13 pp (annualised) and at most
  1.03 pp.

| m13 (posterior means) | full | 1980 |
|---|---|---|
| tau = exp(mu_g / 2) | 2.54 | 1.87 |
| mu_g, phi_g, omega_g | 1.81, 0.90, 0.63 | 1.20, 0.82, 0.89 |
| rho, mu_bar | 0.27, 2.77 | 0.39, 2.63 |
| psi (E, I, P) | -0.46, -0.35, -0.19 | -0.43, -0.35, -0.18 |
| rho_EI | 0.12 | 0.09 |
| long-run measurement sd eps_sd_bar (E, I, P) | 3.36, 2.59, 2.20 | 2.50, 2.61, 2.01 |
| mu_h[E]; phi[E] (90%) | 2.08; 0.975 | 1.58; 0.84 [0.11, 0.99] |

**tau_t** (posterior median [90%], full -> 1980):

| quarter | full | 1980 |
|---|---|---|
| 1982 Q4 | 5.08 [3.03, 8.89] | 5.12 [2.88, 10.10] |
| 2008 Q4 | 1.37 | 1.32 |
| 2015 Q1 | 0.77 | 0.67 |
| 2019 Q4 | 4.53 | 3.11 |
| 2020 Q2 | 11.95 [8.32, 19.07] | 13.67 [9.09, 24.06] |
| 2026 Q2 | 1.25 | 1.03 |

The average of the posterior-median tau_t over 1993 Q1 - 2019 Q4 is 1.64 on
the full sample and 1.46 on the 1980 sample. In 1980 Q1 tau_t is 3.01
against 1.62, but 1980 Q1 is period 1 of the 1980 sample, set by the initial
condition. The latent-volatility profile is the same: high in the early
1980s, low in the 2000s-2010s, a large real shock at COVID. The 1980 fit
has a lower long-run level and less persistent, more variable latent
volatility.

**4. m8, m10, m11: essentially unchanged.** The hours sample loses only
1978 Q4 - 1979 Q4.
- **m10.** tau 1.59 -> 1.56, rho 0.43 -> 0.45, omega_g 0.94 -> 0.96, mu in
  2020 Q2 -28.2. Its 4 chains agree (lp__ -2260 to -2262), with no sign of
  two modes. The mean path has correlation 0.9998 with the full-sample path.
- **m8.** sig_E 2.86 -> 2.66; path correlation 0.9996.
- **m11.** beta_emp 0.61 -> 0.65; path correlation 0.9994.

**5. m3: still weakly identified, but less so**
(`logs/compare_full_vs_1980.log`, computed as in `R/diag_m3_ident.R`):
- corr(tau, mean off-diagonal Omega) goes from -0.70 to -0.54;
- Omega means go from (0.07, 0.04, -0.23) to (0.03, -0.07, -0.14);
- the posterior sd of Omega is 49 / 71 / 67% of the prior sd (full sample
  46 / 57 / 83%);
- min bulk / tail ESS goes from 882 / 438 to 2353 / 2841, and E-BFMI is
  0.95.

### Variance decomposition (1980 sample)

`tables/variance_decomposition.csv` gives the per-draw var(mu) /
var(headline): the median and 5-95%. Note that "excl. COVID" here drops
2020 Q1 - 2021 Q1 (5 quarters; the LFO tables' "excl. COVID" drops
2020 Q2 - 2021 Q1).

| Series | 1980 Q1 - 2026 Q2 (T 186) | excl. 2020 Q1 - 2021 Q1 (T 181) | full-sample fit, for reference |
|---|---|---|---|
| m1 | 0.850 [0.768, 0.936] | 0.925 | 0.786 (1959-2026) |
| m2 | 0.850 [0.765, 0.938] | 0.929 | 0.785 |
| m3 | 0.893 [0.716, 0.982] | 0.950 | 0.831 |
| m4 | **0.328** [0.254, 0.416] | 0.552 | 0.798 |
| m5 | 0.852 [0.773, 0.930] | 0.926 | 0.783 |
| m6 | 0.851 [0.762, 0.940] | 0.934 | 0.779 |
| m7 | 0.509 [0.422, 0.584] | 0.781 | 0.614 |
| m9 | **0.323** [0.240, 0.414] | 0.546 | 0.785 |
| m12 | **0.373** [0.300, 0.462] | 0.611 | 0.779 (excl. COVID 0.727) |
| **m13** | **0.872 [0.811, 0.938]** | **0.776 [0.675, 0.883]** | 0.790 [0.705, 0.889] (excl. 0.731); on exactly 1980 Q1 - 2026 Q2 0.883 [0.819, 0.950] (excl. 0.793) |
| m8 | 0.853 [0.769, 0.936] | 0.850 | 0.848 (1978 Q4 - 2026 Q2) |
| m10 | 0.784 [0.719, 0.854] | 0.607 | 0.747 (excl. 0.558) |
| m11 | 0.856 [0.772, 0.941] | 0.883 | 0.846 (1978-2026) |

- **m13: the quarters, not the fit.** m13's ratio rises from 0.79 to 0.87.
  This is entirely a change of quarters, not of estimates. The full-sample
  m13 fit restricted to exactly the same quarters (1980 Q1 - 2026 Q2) gives
  0.883 [0.819, 0.950], or 0.793 excluding COVID. That is slightly above the
  1980 fit's 0.872 and 0.776. Over 1978 Q4 - 2026 Q2 it gives 0.855 (0.755).
  (Recomputed 2026-10-19 from `output/fit_m13.rds` mu draws; not in the
  tables.)
  On the 1980 sample m13 attributes about 13% of headline variance to
  measurement error, or about 22% excluding COVID.
- **m4, m9 and m12** fall to 0.32-0.37. This is not a finding about
  measurement error. It reflects the failure in section 2 above: COVID is
  counted as measurement error.

### LFO on the 1980 sample

Group 1 (E, I, P; `tables/model_comparison.csv`; full-sample scores from
`tables/lfo_full_vs_sample_models.csv`):

| model | elpd_lfo 1980 | full sample | change | diff vs m13 | SE iid / NW | refits | weak mixing | MCSE |
|---|---|---|---|---|---|---|---|---|
| **m13** | **-1039.50** | -1045.82 | +6.3 | 0 | - | 7 | 1 | 0.02 |
| m12 | -1096.90 | -1105.47 | +8.6 | -57.40 | 27.2 / 36.6 | 6 | 0 | 0.75 |
| m9 | -1117.50 | -1122.24 | +4.7 | -78.01 | 30.6 / 40.4 | 6 | 2 | 0.60 |
| m7 | -1120.42 | -1177.54 | +57.1 | -80.93 | 12.2 / 19.0 | 6 | 0 | exact |
| m4 | -1121.23 | -1123.72 | +2.5 | -81.74 | 31.1 / 42.3 | 5 (1 retried) | 2 | 0.91 |
| m5 | -1134.18 | -1168.84 | +34.7 | -94.68 | 33.8 / 42.9 | 5 | 0 | exact |
| m1 | -1149.40 | -1195.66 | +46.3 | -109.90 | 31.8 / 35.7 | 4 | 0 | exact |
| m2 | -1154.59 | -1196.25 | +41.7 | -115.09 | 33.8 / 42.9 | 5 | 0 | exact |
| m6 | -1156.36 | -1195.59 | +39.2 | -116.87 | 34.3 / 43.7 | 4 | 0 | exact |
| m3 | -1158.29 | -1195.93 | +37.6 | -118.79 | 36.5 / 45.0 | 5 | 0 | exact |

The group 1 order is m13, m12, m9, m7, m4, m5, m1, m2, m6, m3. On the full
sample it was m13, m12, m9, m4, m5, m7, m6, m1, m3, m2.

Group 2 (with hours) and the other targets:
- **Group 2.** m10 -1415.42 (full -1417.72) and m8 -1535.03 (full
  -1541.96). m10 - m8 is +119.62 (SE 53.0 iid, 65.0 NW, z 1.84).
- **Group 3.** m11 joint -2469.50 (full -2469.88).
- **GDP-only target.**
  - m10 -1057.97 (full -1058.44), m11 -1120.04 (full -1113.36), m8 -1159.28
    (full -1160.69).
  - The top five are m13, m10, m12, m9, m11. On the full sample they were
    m13, m10, m12, m11, m9.

**Key differences, whole window** (`tables/lfo_full_vs_sample.csv`; SE iid /
Newey-West lag 8; "change" is 1980 minus full, with its NW SE). The nested
measurement-SV pairs m4 - m2 and m9 - m6 were added to the 1980-only tables
on 2026-10-19 (`R/03_compare.R summary` rerun from the existing checkpoints).
That run only added rows: no row of any 1980 table changed, and no
full-sample table was touched.

| a - b | target | full sample | 1980 sample | change (NW SE) |
|---|---|---|---|---|
| m4 - m2 (measurement SV, nested) | group 1 | +72.5 (12.1 / 19.1), z 3.79 | **+33.4** (9.4 / 15.8), **z 2.11** | -39.2 (14.2) |
| m9 - m6 (measurement SV, nested) | group 1 | +73.3 (12.2 / 19.9), z 3.69 | **+38.9** (9.7 / 15.1), **z 2.58** | -34.5 (15.3) |
| m12 - m5 (measurement SV + structural E/I, with AR(1) errors) | group 1 | +63.4 (14.3 / 20.7), z 3.06 | +37.3 (13.3 / 15.1), z 2.47 | -26.1 (20.5) |
| m9 - m5 (not nested: adds SV, drops AR(1) errors) | group 1 | +46.6 (12.5 / 19.7), z 2.37 | +16.7 (9.3 / 14.7), z 1.13 | -29.9 (14.3) |
| m4 - m1 (not nested: adds SV and latent AR(1)) | group 1 | +71.9 (12.0 / 19.4), z 3.71 | +28.2 (10.4 / 19.4), z 1.46 | -43.8 (13.8) |
| m5 - m2 (serial errors) | group 1 | +27.4 (2.9 / 2.0), z 14.0 | +20.4 (3.6 / 2.2), z 9.40 | -7.0 (1.0) |
| m12 - m9 (serial errors, with SV) | group 1 | +16.8 (6.7 / 6.1), z 2.74 | +20.6 (7.9 / 5.5), z 3.72 | +3.8 (8.8) |
| m13 - m12 (latent SV) | group 1 | +59.6 (24.7 / 28.9), z 2.06 | +57.4 (27.2 / 36.6), **z 1.57** | -2.2 (11.1) |
| m9 - m4 (structural, SV) | group 1 | +1.5 (1.5 / 2.2), z 0.68 | +3.7 (1.6 / 2.4), z 1.58 | +2.2 (4.2) |
| m6 - m2 (structural, constant) | group 1 | +0.66 (0.34 / 0.38), z 1.74 | -1.77 (0.89 / 1.01), z -1.76 | -2.4 (1.2) |
| m7 - m1 (COVID scale) | group 1 | +18.1 (18.4 / 16.7), z 1.08 | +29.0 (33.8 / 28.4), z 1.02 | +10.9 (12.6) |
| m13 - m5 (all three features) | group 1 | +123.0 (15.3 / 25.7), z 4.78 | +94.7 (33.8 / 42.9), z 2.20 | -28.3 (25.8) |
| m13 - m7 | group 1 | +131.7 (12.2 / 20.7), z 6.35 | +80.9 (12.2 / 19.0), z 4.25 | -50.8 (13.3) |
| m10 - m8 | group 2 | +124.2 (53.5 / 64.7), z 1.92 | +119.6 (53.0 / 65.0), z 1.84 | -4.6 (2.0) |
| m10 - m8 | GDP-only | +102.2 (35.4 / 45.6), z 2.24 | +101.3 (37.9 / 48.6), z 2.08 | -0.9 (4.0) |
| m13 - m10 | GDP-only | +12.6 (6.9 / 5.5), z 2.30 | +18.5 (4.6 / 3.3), z 5.66 | +5.9 (3.9) |
| m13 - m11 | GDP-only | +67.5 (12.7 / 18.4), z 3.67 | +80.6 (25.0 / 24.6), z 3.27 | +13.0 (13.4) |
| m10 - m12 | GDP-only | +47.0 (28.5 / 30.1), z 1.56 | +38.9 (29.0 / 36.9), z 1.05 | -8.1 (10.4) |

**By sub-period**, 1980 sample, with full sample in brackets (NW SE for the
long periods; with n = 4 in 2020 the NW SE is not meaningful):

| a - b | 1990-99 | 2000-07 | 2008-19 | 2020 | 2021-22 | 2023-26 | pre-2020 (n 120) | excl. 2020 Q2 - 2021 Q1 (n 142) |
|---|---|---|---|---|---|---|---|---|
| m13 - m12 | -0.4 [+1.3] | +0.4 [+4.0] | +6.4 [+14.8] | +32.4 [+30.1] | +12.3 [+2.5] | +6.4 [+7.0] | +6.4, SE 2.75, z 2.31 [+20.1, z 4.14] | +22.6, SE 10.5 [+28.5, SE 6.0] |
| m9 - m5 (not nested) | -3.7 [-2.1] | -1.2 [+6.2] | +8.1 [+25.4] | +5.1 [-8.7] | -8.3 [+3.1] | +16.5 [+22.7] | +3.3, SE 6.2 [+29.6, SE 11.2] | +12.5, SE 14.6 [+56.8, SE 17.7] |
| m4 - m2 (measurement SV) | +1.6 [+5.9] | +3.4 [+12.0] | +14.5 [+31.5] | +4.1 [-7.5] | -8.8 [+5.3] | +18.6 [+25.2] | +19.4, SE 6.4, z 3.02 [+49.5, z 4.68] | +30.2, SE 15.6, z 1.93 [+81.6] |
| m9 - m6 (measurement SV) | +2.0 [+5.6] | +3.7 [+13.3] | +15.7 [+33.6] | +6.4 [-8.8] | -7.4 [+4.3] | +18.4 [+25.4] | +21.4, SE 6.6, z 3.22 [+52.5, z 4.78] | +33.4, SE 14.9, z 2.24 [+83.6] |
| m4 - m1 (not nested) | +2.8 [+6.1] | +4.3 [+12.2] | +16.3 [+31.5] | -4.4 [-8.3] | -9.5 [+5.1] | +18.7 [+25.3] | +23.4, SE 7.1 [+49.7] | +33.1, SE 16.4, z 2.02 [+81.7] |
| m5 - m2 | +5.3 [+7.6] | +5.0 [+7.1] | +6.8 [+8.4] | +0.3 [+0.2] | +1.0 [+1.3] | +1.9 [+2.7] | +17.2, SE 2.1 [+23.1] | +20.0, SE 2.1 [+27.0] |
| m12 - m9 | +5.7 [+6.7] | +6.0 [+9.2] | +2.9 [+3.9] | +5.2 [-4.2] | +0.5 [+0.4] | +0.4 [+0.7] | +14.5, SE 2.8 [+19.9] | +16.4, SE 3.0 [+20.5] |
| m13 - m10 (GDP-only) | +3.5 [-0.5] | +6.4 [+7.0] | +4.7 [+6.3] | -1.1 [-2.3] | +2.6 [+1.8] | +2.5 [+0.4] | +14.6, SE 2.8 [+12.7] | +19.2, SE 3.2 [+14.9] |
| m10 - m12 (GDP-only) | -3.9 [+1.9] | -6.0 [-3.0] | +1.7 [+8.5] | +33.5 [+32.4] | +9.7 [+0.7] | +3.9 [+6.6] | -8.2, SE 3.9 [+7.4] | +3.4, SE 9.5 [+13.6] |

**m13 - m12 in detail.**
- **Where the margin comes from.** 2020 Q2 alone is worth +25.6. m13 is
  ahead in 107 of 146 quarters (full sample 122), with a median gain per
  quarter of +0.04.
- **Largest single quarters.** Gains: 2020 Q2 +25.6, 2020 Q4 +7.4,
  2021 Q4 +4.7. Losses: 2016 Q3 -1.7, 2020 Q1 -0.9, 2008 Q4 -0.2.
- **Lag sensitivity.** The NW z is 1.78 / 1.57 / 1.49 / 1.45 at lags
  4 / 8 / 12 / 16 (iid z 2.11).
- **MCSE.** The MCSE of the difference is 0.75 (m12's 4 passes disagree,
  almost all in 2020 Q2).

**Refits and fit quality** (`tables/lfo_fit_diagnostics.csv`,
`tables/lfo_refit_attempts.csv`).
- **All 86 stored fits meet the refit rule** (13 initial fits and 73
  refits). Every initial fit passed on its first attempt, so no training
  window was lengthened. No high-k step was used without a refit (largest
  k used 0.699).
- **One failed attempt.** m4's refit on 1:146 (data to 2016 Q2) failed
  attempt 1: 12 divergences, E-BFMI 0.011, R-hat 1.28, bulk ESS 11. It
  passed attempt 2 (1 divergence, E-BFMI 0.88, R-hat 1.006, ESS 790).
- **Weak mixing.** Five fits meet the rule but have bulk ESS below 400:
  m4 1:81 (226) and 1:129 (271), m9 1:65 (376) and 1:106 (160), and m13
  1:122 (351). None has E-BFMI below 0.3.
- **m13's fits.** m13 has 8 fits: 1:40, 64, 93, 122, 140, 160, 162 and 177
  (data to 1989 Q4, 1995 Q4, 2003 Q1, 2010 Q2, 2014 Q4, 2019 Q4, 2020 Q2
  and 2024 Q1).

**Retry (seed-selection) sensitivity**
(`tables/lfo_retry_sensitivity_summary.csv`). Only m4's 1:146 fit was
retried, so only the 16 targets it serves are rescored:
- m9 - m4 is +3.73 as reported, +0.46 with first attempts only and +1.28
  pooled;
- m13 - m4 is +81.7, +78.5 and +79.3 (NW z 1.93, 2.00, 1.98);
- m13 - m12 is unaffected (neither model was retried).

The retried attempt was re-created and reproduced its logged diagnostics
exactly.

**Particle filter.**
- **Pass totals.** Every SV model has 4 passes:
  - m13: -1039.45 / -1039.54 / -1039.52 / -1039.47 (MCSE 0.02);
  - m12: -1095.45 / -1097.74 / -1098.58 / -1095.84 (MCSE 0.75);
  - m9: -1118.38 / -1117.68 / -1118.20 / -1115.75 (MCSE 0.60);
  - m4: -1120.14 / -1119.31 / -1123.24 / -1122.24 (MCSE 0.91);
  - m10: MCSE 0.18.
- **Bias over COVID.** Bias against high-particle references
  (`tables/lfo_pf_covid_check.csv`), summed over 2020 Q1 - 2021 Q1, for the
  schedule against the mean of 3 reference passes at 2020 Q2:

  | | m13 | m12 | m9 | m4 |
  |---|---|---|---|---|
  | schedule | -0.02 | +0.08 | -1.55 | +0.55 |
  | fixed N = 2000 | -0.03 | -3.27 | -3.71 | -3.31 |
- **Bias over 2015-19.** The schedule's bias over the 20 targets
  2015-2019 (`tables/lfo_pf_precovid_check.csv`) is -0.030 (m13) and -0.020
  (m12).
- **The 2020 Q2 problem.** For the three models that treat COVID as
  measurement error, the 2020 Q2 predictive density is hard to estimate.
  The run-to-run sd at N = 100000 is 1.5 (m12), 1.2 (m9) and 1.8 (m4) nats.
  The 400000-particle reference passes have se 0.56, 1.11 and 0.39, against
  0.003 for m13. Comparisons that involve m4, m9 or m12 therefore carry
  about 1-2 nats of particle uncertainty from that one quarter. That is
  small against m13 - m12 (57 nats) but not against m9 - m4 (3.7 nats).
- **Not run.** The 2015-19 checks were run for m13 and m12 only. No COVID
  check was run for m10 on this sample.

**Process.**
- **The job vanished once.** The LFO driver and its R process disappeared
  without an exit status between 00:46 and 01:31 on 2026-10-19, during the
  m12 COVID particle check. The system log shows no jetsam or memorystatus
  kill, so the cause is unknown. This is logged in
  `output/sample_1980/lfo/lfo_log.txt`. The driver was restarted, resumed
  from its checkpoints, and the m12 check then completed in 8 minutes.
- **The budget was overrun.** Compute took about 7 h 05 min against a
  planned 5 h. Nothing in the protocol was cut, and COVID checks for m9 and
  m4 were added (`run_lfo_1980_extra.sh`).
- **The fits driver was interrupted once.** `run_1980.sh` ran m1-m9 and
  both m13 seeds from 07:42 on 2026-10-18. At 08:27:09 the m10 fit received
  SIGTERM (exit 143; `logs/driver.log`). The driver then logged "start
  fit_m12", and that fit and the driver vanished with no end line. The
  cause is not recorded. The driver was restarted at 20:37. It skipped the
  finished steps (the fit file and its diagnostics row both exist) and then
  ran m10, m12, m11, the m11 diagnostics row, the variance decomposition and
  figures 01-06 (DRIVER DONE 20:54:58). So m10, m12 and m11 were fitted
  about 12 hours after the others. `R/02_fit_models.R` was edited
  in between and again after the restart (mtime 20:55). Results are not
  affected:
  - the m1 refit with the current code is bit-identical;
  - the Stan data lists match;
  - a sandbox re-run of the 1980 variance decomposition with the current
    `R/02_fit_models.R` is byte-identical;
  - a sandbox re-run of m1's 1980 LFO with the current `R/03_compare.R`
    gives identical pointwise values (max |diff| 0). That file was edited
    at 21:16:25, after m1's checkpoint at 21:16:19.
- **Figure 07** (`figures/07_lfo_cumulative_elpd.png`) was first drawn in a
  sandbox output root so that figures 01-06 were not touched, and its
  subtitle did not say "1980 sample". On 2026-10-19 `R/04_figures.R` was
  changed to append the sample note to the subtitles of figures 04, 06 and
  07 too (01 and 02 already had it; 05 already names its sample start). The
  note is empty on the full sample, so the full-sample figures are
  unchanged. All seven 1980 figures were redrawn in
  `output/sample_1980/figures/` from the same fits and tables.

### Which full-sample conclusions survive the 1980 start

The targets are identical (146 quarters, 1990 Q1 - 2026 Q2). What differs is
the training data: 40 against 121 quarters of initial training, with no
1960s-70s. So each "change" mixes the effect of dropping the early data with
the effect of estimating parameters from less history. All of this is
pseudo-real-time on the current vintage, one step ahead.

**Survive**
- **m13 is the best model** on every target: group 1, and GDP-only, ahead
  of m10 by +18.5 (z 5.66; full sample +12.6, z 2.30). Its posterior path
  barely moves (correlation 0.9985). On the 1980 sample it is the only model
  with measurement SV (m4, m9, m12, m13) that keeps COVID in the latent
  series. It is not the only 3-measure model that does:
  - the constant-variance models m1, m2, m3, m5 and m6 also keep most of
    2020 Q2 in latent growth on both samples. Their posterior mean of mu in
    2020 Q2 is -23.5 to -24.4 on the 1980 sample (-22.5 to -23.4 on the full
    sample), against -27.9 for m13 and a headline of about -28;
  - m7 moves 2020 Q2 into its COVID variance scale by design (mu +3.6; full
    sample +1.7).

  (Corrected 2026-10-19. An earlier version said m13 is "the only 3-measure
  model that keeps COVID in the latent series on both samples".)
- **Serially correlated measurement errors**, the most robust feature:
  - m5 - m2 is +20.4 (z 9.4; full sample +27.4, z 14.0);
  - m12 - m9 is +20.6 (z 3.72; full sample +16.8, z 2.74), +14.5 before
    2020;
  - m5 - m2 is positive in every sub-period on both samples. m12 - m9 is
    positive in every sub-period except 2020 on the full sample (-4.2);
  - m13's psi is nearly unchanged at (-0.43, -0.35, -0.18).
- **The structural E/I restriction is still not supported.** m9 - m4 is
  +3.7 (z 1.58), but +0.5 with first attempts only and +1.3 pooled, within
  about 2 nats of particle uncertainty. m6 - m2 turns slightly negative
  (-1.77, z -1.76).
- **Hours model.** m10 - m8 is +119.6 (z 1.84; full sample +124.2,
  z 1.92), and m10's posterior is essentially unchanged.
- **Production is the cleanest single signal** (smallest error sd) in
  every specification.

**Weaken**
- **Latent SV (m13 - m12).** The point estimate holds: +57.4 against +59.6.
  The evidence is weaker:
  - the NW SE rises to 36.6, so z is 1.57 (1.45-1.78 over lags 4-16);
  - before 2020 the gain falls from +20.1 (z 4.1) to +6.4 (z 2.3);
  - outside 2020 Q2 - 2021 Q1 it falls from +28.5 to +22.6 (SE 10.5).
    Outside all of 2020-22, which also removes m12's post-COVID
    misspecification on this sample, it falls from +27.1 (NW SE 5.9, z 4.59)
    to +12.7 (NW SE 4.9, z 2.59);
  - m13 is ahead in 107 of 146 quarters (full sample 122);
  - the 1980 margin leans on 2020 (+32.4) and 2021-22 (+12.3, against +2.5
    on the full sample), where m12's 1980 fit mis-handles COVID, so part of
    the margin scores that specification failure rather than a general
    gain;
  - **the production switch rule's z_NW >= 2 condition over m12 would not
    be met on this sample**.
- **All three features together (m13 - m5).** +94.7 (z 2.20) against
  +123.0 (z 4.78).
- **SV in the measurement errors: the gain roughly halves but is still
  supported.** Two comparisons add only measurement SV to an otherwise
  identical model (nested): m4 - m2 (independent errors) and m9 - m6
  (structural E/I correlation). A third, m12 - m5, adds SV and the
  structural restriction to a model with AR(1) errors. The restriction
  alone is worth little (m9 - m4 +3.7, m6 - m2 -1.8). Group 1, NW lag 8:

  | a - b | full sample | 1980 sample | 1980, pre-2020 (full) | 1980, excl. 2020-22 |
  |---|---|---|---|---|
  | m4 - m2 | +72.5, z 3.79 | +33.4, z 2.11 | +19.4, z 3.02 (+49.5) | +38.0, z 2.88 |
  | m9 - m6 | +73.3, z 3.69 | +38.9, z 2.58 | +21.4, z 3.22 (+52.5) | +39.9, z 3.05 |
  | m12 - m5 | +63.4, z 3.06 | +37.3, z 2.47 | +17.8, z 3.54 (+49.4) | +34.8, z 3.00 |

  - **Seed selection.** With first attempts only, m4 - m2 is +36.6
    (z 2.28), and +35.8 (z 2.25) pooled
    (`tables/lfo_retry_sensitivity_summary.csv`). m9 and m6 had no retried
    fit.
  - **Correction.** An earlier version of this note said the gain "largely
    disappears". That rested on m9 - m5 (+46.6, z 2.37 -> +16.7, z 1.13)
    and m4 - m1 (+71.9, z 3.71 -> +28.2, z 1.46). Neither is nested. m9 has
    SV but not m5's AR(1) errors, so m9 - m5 nets the SV gain against losing
    serial correlation, which is worth +20.6 on this sample (m12 - m9). m4 -
    m1 also adds a latent AR(1), worth -5.2 on this sample (m2 - m1; -0.6
    on the full sample).
  - **Caveats.** On the 1980 sample m4, m9 and m12 treat COVID as
    measurement error ("Posterior", section 2), so their 2020 and 2021-22
    scores partly reflect that failure. The nested gains are +4.1 / +6.4 in
    2020 and -8.8 / -7.4 in 2021-22. GDP(E)'s measurement SV is also not
    persistent on this sample (phi_E is prior-dominated; same section), so
    what SV adds for E there is mostly a heavier-tailed error, not
    persistent volatility.
  - **Why it halves.** Before 2020 every model gains from the later start,
    but the constant-variance models gain about twice as much: m1 +51.6, m2
    +55.2 and m5 +49.3, against m4 +25.2, m9 +23.1 and m12 +17.7 (1980
    minus full, pre-2020 sums). Their error variances, GDP(E)'s above all,
    are no longer inflated by the early decades.
  - **The whole-window totals exaggerate this asymmetry.** They are +35 to
    +46 for the constant-variance models and +2.5 to +8.6 for the SV models,
    but they also contain two other periods:
    - 2020, where both kinds of model lose from the later start: m1 -16.5,
      m5 -24.1, m4 -12.6, m9 -10.3, m12 -0.9;
    - 2021-22, where only the SV models lose: m4 -12.3, m9 -9.9, m12 -9.9,
      against m1 +2.3 and m5 +1.4.
  - So about half of the full-sample measurement-SV gain came from pre-1980
    training data. That is the paper's point about unrepresentative early
    measurement errors. The other half does not depend on the start date.
    It is largest in 2008-19 and 2023-26 (m9 - m6 +15.7 and +18.4).

**Change or reverse**
- **m7 (constant variance plus a COVID scale) rises from 6th to 4th** in
  group 1, level with m9 and m4 (-1120.4 against -1117.5 and -1121.2).
- **The variance decompositions of m4, m9 and m12 break down.** They fall
  to 0.32-0.37 because these fits count COVID as measurement error. The
  full-sample statement that "the other specifications give 78-83%" does
  not hold on this sample.
- **m13's variance ratio.** It is 0.87 [0.81, 0.94] (0.78 excluding
  2020 Q1 - 2021 Q1) against 0.79 on the full sample. The difference is
  entirely which quarters are covered. Over exactly 1980 Q1 - 2026 Q2 the
  full-sample fit gives 0.883 [0.819, 0.950] (0.793 excluding COVID),
  slightly above the 1980 fit (see "Variance decomposition"). The
  "about a fifth is measurement error" figure is specific to the full
  sample; on 1980-2026 it is about an eighth, or about a fifth excluding
  COVID.
- **m10 against m12 (GDP-only).** m10 - m12 falls from +47.0 to +38.9
  (z 1.05), and before 2020 m10 is now behind m12 (-8.2). As before, none of
  the labour-model gaps measures the value of the labour data.

### Production: unchanged, and the recommendation *(superseded 2026-10-19: production moved to 1980 Q1 at the user's request; see "Production sample moved to 1980 Q1")*

Production still fits m13 on the full sample. This round changed no
production file.

**Recommendation (comparison author): keep production on m13, unchanged
(full sample; do not move to 1980 Q1).** Reasons:
- On the 1980 sample m13 is still the best model on every target. It is
  also the only model with measurement SV that keeps COVID in the latent
  series. The constant-variance models also keep most of it there, but
  they score 95-119 nats worse in group 1. (Corrected 2026-10-19: this
  reason first said "the only 3-measure model", which is false.)
- If anyone wants to firm up the latent-SV evidence, the useful next step
  is a test that does not reuse this evaluation window, such as real-time
  scoring of the next ABS releases. Re-scoring on yet another sample would
  not settle it.

That recommendation does not separately argue the sample-start question.
Two facts from the tables bear on it (added by the documentation editor):
- **Prediction.** m13 trained from 1980 Q1 scores +6.3 nats better on the
  same 146 targets than m13 trained from 1959 Q4 (-1039.50 against
  -1045.82; NW SE 3.74, so z 1.69; `tables/lfo_full_vs_sample_models.csv`).
  That is below the z_NW >= 2 bar that every production switch has used.
  The gain is mostly in 1990-99 (+4.5, NW SE 2.31, z 1.96). Those are the
  targets nearest the dropped 1960s-70s data, and the 1980 run's training
  sets are shortest there. From 2000 on the gain is about zero: +1.8 over
  2000-2026 (NW SE 2.51, z 0.72) and +1.3 over 2008-2026 (z 0.63). For
  current production the start date hardly matters for prediction.
- **The estimate.** The two fits' latent paths differ by 0.13 pp
  (annualised) on average, at most 1.03 pp.
- **The variance ratio.** The published ratio would move from 0.79 to 0.87
  only because it would be reported over 1980-2026 instead of 1959-2026.
  Over exactly 1980 Q1 - 2026 Q2 the full-sample fit gives 0.883
  [0.819, 0.950] (0.793 excluding 2020 Q1 - 2021 Q1). The 1980 fit gives
  0.872 (0.776). The estimation start changes neither the path nor the
  ratio materially. (An earlier version said moving the start "would
  change the published variance ratio (0.79 -> 0.87) more than it changes
  the path". That mixed up the estimation start with the reporting window.)

### Open items (1980 sample)

- **Version control (fixed 2026-10-19).** `.gitignore` now has
  `output/sample_*/**/*.rds`. That covers the 45 fit objects, LFO
  checkpoints and particle-check files under `output/sample_1980/` (about
  1.7 GB; several are over GitHub's 100 MB limit). `git check-ignore`
  matches `output/sample_1980/fit_m1.rds` and
  `output/sample_1980/lfo/lfo_m1.rds`. The 72 files that are meant to be
  kept (tables, logs, figures, driver scripts; 3.1 MB) are still untracked
  but not ignored.
- **Tests.** The unit tests (283 expectations in 4 files, all pass) test
  only the full-sample fits in `output/`. There are no 1980-sample checks.
- **Multimodality.** One alternative starting region (COVID in mu, from
  full-sample draws) has now been tried for each of m4, m9 and m12. In all
  three, every chain left it within about 50 warmup iterations and matched
  the reference fit ("Posterior", section 2). Other regions, such as tempered
  or optimisation-based searches, were not tried.
- **Particle checks not run:** the 2015-19 particle checks for m9, m4 and
  m10, and a COVID check for m10.
- **m11 retries.** `run_retrycheck()` in `R/03_compare.R` scores m11
  without `ll_post_m11`. This is pre-existing code, not exercised here (m11
  had no retried fit), but it would be wrong if an m11 refit is ever
  retried.
- **Post hoc.** m13 was proposed after the re-specified m10's results were
  seen on this same evaluation window. Re-scoring on the 1980 sample reuses
  that window, so it is not independent evidence.

## Round of 2026-10-18: m13 (m12 + latent SV) and the production switch m12 -> m13

The data vintage is unchanged (2026 Q2, T = 267). The scores of every
earlier model are unchanged; m13 is added to group 1. The 2026-10-17 round
below stays valid except where this section marks it superseded.

### Summary

- **m13** = m12 + stochastic volatility on the latent innovation (the
  feature that the re-specified m10 added). It is the best model on every
  comparison: group 1 -1045.82 against -1105.47 for m12 (+59.65, Newey-West
  SE 28.89, z 2.06), and on the GDP-only target it is ahead of m10 by +12.62
  (NW SE 5.48, z 2.30).
- **Production switched from m12 to m13**, under a rule fixed before m13's
  LFO results were known and stricter than the m9 -> m12 rule (below).
- **Post hoc.** m13 was proposed after the re-specified m10's LFO results
  had been seen (m10 - m9 +63.8 and m10 - m12 +47.0 on GDP-only, suggesting
  that latent SV is what helps). m13 is therefore a second latent-SV
  specification scored on the same evaluation window, and its z-values
  overstate the evidence somewhat, and nothing in the decision rule corrects
  for that. The rule's extra condition (z >= 2 with first attempts only)
  guards only against choosing a lucky refit seed, and it changed nothing
  here (m13 needed no retry; z 2.06 either way). Counting the re-specified
  m10 and m13 as two latent-SV tries on the same window, a Bonferroni bar
  for two tries (z about 2.24) would not be met by z 2.06.
- **The gain is not only COVID.** 2020 Q2 alone is worth +23.7 nats and
  dominates the SE; without the COVID quarters m13 - m12 is +28.5 (NW SE 6.0,
  z 4.7), before 2020 +20.1 (SE 4.9, z 4.1), and m13 is ahead in 122 of 146
  quarters.
- **The variance decomposition barely moves on the full sample** (0.790
  against m12's 0.779), but outside COVID m13 puts more of headline variance
  into measurement error than m12 (see "m13: fits and posterior").

### m13: specification

`stan/m13_serial_latent_sv.stan`:
- y_t^k = mu_t + eps_t^k, k = E, I, P; eps_t^k = psi_k eps_{t-1}^k + u_t^k;
  u_t ~ N(0, R_t) with m12's structural covariance (rho_EI = tanh(z_EI), P
  independent) and stationary AR(1) SV h_{t,k}, all exactly as m12;
- mu_t = c + rho mu_{t-1} + exp(g_t / 2) eta_t;
- g_t = mu_g + phi_g (g_{t-1} - mu_g) + omega_g zeta_t, stationary and
  non-centred, g_1 ~ N(mu_g, omega_g² / (1 - phi_g²));
- mu_1 ~ N(c / (1 - rho), exp(g_1) / (1 - rho²)): m10's convention (the
  stationary law at the period-1 variance; the marginal stationary law under
  latent SV has no closed form, and period 1 is outside the LFO window).

With g_t = log(tau²) for all t the model is m12 exactly.
- **Priors.** mu_g ~ N(2, 1), phi_g ~ U(-1, 1), omega_g ~ half-N(0, 0.5)
  (m10's) replace tau ~ half-N(0, 2); every other prior is m12's
  (`priors.md`). `tau = exp(mu_g / 2)` and `tau_t = exp(g_t / 2)` are
  generated quantities.
- **Sampling.** Given (g, h) the model block integrates mu out with m12's
  quasi-differenced Kalman filter on (mu_t, mu_{t-1}), with tau replaced by
  tau_t in the predicted-state factor and in the backward-sampling variance;
  mu is drawn by FFBS. The data block guards K = 3.
- **Tests** (`tests/test_m13.R`, 133 checks, all pass, none skipped):
  - the Stan filter's total and per-period log-likelihood against the dense
    multivariate-normal density, to 1e-9 / 1e-8, for 5 psi settings
    (including psi = 0, psi_E = psi_I and |psi| up to 0.95) with non-constant
    g and h, a COVID-like jump of +5 in g, g = -12 in some periods,
    z_EI = -5 and h = -15;
  - with constant tau_t, bit-identical (`expect_identical`) to m12's
    `kf_m12`, in Stan and in C++;
  - FFBS against the dense smoother with non-constant g, 20000 draws: worst
    mean error 1.8-2.0 Monte Carlo SEs, worst covariance error 2.3;
  - the C++ filter (`R/lib/filter_m13.cpp`) against Stan, to 1e-11;
  - the particle filter (particles over (g, h), Kalman over
    (mu_t, mu_{t-1})) against the exact Kalman filter at omega, omega_g ->
    0, and against m12's particle filter at omega_g = 0 (to 1e-10, same
    seed); against brute-force Monte Carlo over 200000 (g, h) paths with a
    latent-volatility jump: -37.488 (bootstrap proposal) and -37.497
    (adapted proposal) against -37.496 (se 0.005);
  - the C++ conditional filter against the fit's `log_lik_pred` on 5 draws.
  The LFO dispatch (`filter_loglik("m13")`, with a particle schedule and
  through `filter_par`'s workers) is checked against the dense MVN in
  `tests/test_filters.R` (x), 9 checks. The whole of `tests/test_filters.R`
  was later run and passes (63 checks, none failed or skipped).
- **Warmup stalls** (`R/fit_m13.R stall`, 150 warmup iterations,
  adapt_delta 0.99, treedepth 12, timeout 600 s): 0 of 100 chains stalled at
  init = 2 (median 116 s, max 265 s) and 0 of 40 at init = 0.5 (median 57 s).
  Process note: seeds 5 and 6 of the init = 2 batch were killed by mistake
  mid-run and rerun separately (neither stalled); the harness output
  (`stall_init2/`, `stall_init2_redo56/`, `stall_init05/`, each with a
  `stall_summary.rds`) was written to the session scratch directory, not the
  repo, and `stall_init2/stall_summary.rds` still holds the killed rows
  (exit -15). Rerun with `Rscript R/fit_m13.R stall` to reproduce.

### m13: fits and posterior

Both seeds use the production settings and pass the production gate, with
mu_g, phi_g and omega_g added to the hyperparameters and R-hat <= 1.05
required on every g[t] as well (`output/tables/m13_gate_seeds.csv`; fits
`output/fit_m13.rds`, `output/fit_m13_seed2.rds`, gitignored, regenerated by
`Rscript R/fit_m13.R` and `Rscript R/fit_m13.R --seed=2`):

| Check | seed 1 | seed 2 | Gate |
|---|---|---|---|
| divergences | 4 of 8000 (0.050%) | 0 of 8000 | < 0.1% |
| E-BFMI by chain | 0.935, 0.854, 0.969, 0.920 | 0.843, 0.844, 0.869, 0.855 | min >= 0.3 |
| max R-hat, hyperparameters | 1.0020 (omega[1]) | 1.0027 (omega[2]) | <= 1.01 |
| max R-hat, mu / h / g | 1.0011 / 1.0018 / 1.0020 | 1.0020 / 1.0019 / 1.0015 | <= 1.05 |
| min bulk / tail ESS, hyperparameters | 1553 (omega[2]) / 3213 | 1784 (omega[2]) / 2618 | >= 400 |
| min bulk / tail ESS, mu | 4344 / 4866 | 4550 / 5191 | >= 400 |
| max-treedepth hits | 0 | 0 | - |
| sampling time (4 parallel chains) | 1375 s (wall 1384 s; chains 989-1375 s) | 1664 s (wall 1674 s) at load average 10-11 | <= 30 min |

`output/tables/fit_diagnostics.csv` tabulates m13 from the saved fit
(`R/02_fit_models.R --diag-only m13`), so its `wall_secs` is CmdStan's
sampling total, 1375 s, not the 1384 s wall time.

The two seeds' posterior means agree to within 0.007. Seed 1, posterior
medians with 90% intervals (annualised units):

| | m13 | m12 |
|---|---|---|
| latent innovation sd | tau_t, see below; exp(mu_g / 2) 2.47 [1.71, 3.61] | tau 3.99 [3.66, 4.34] |
| phi_g, omega_g | 0.91 [0.80, 0.97], 0.61 [0.42, 0.89] | - |
| psi (E, I, P) | -0.46, -0.36, -0.19 | -0.47, -0.34, -0.16 |
| psi_P 90% interval | [-0.35, -0.02] (now excludes 0) | [-0.36, 0.04] |
| rho_EI | 0.12 [-0.03, 0.27] | 0.175 [0.00, 0.32] |
| rho, mu_bar (means) | 0.27, 2.77 | 0.09, 3.10 |
| variance ratio, full sample | 0.790 [0.705, 0.889] | 0.779 [0.707, 0.850] |
| variance ratio, full sample excl. 2020 Q1 - 2021 Q1 | 0.731 | 0.727 |
| variance ratio, 1978 Q4 - 2026 Q2 | 0.855 [0.790, 0.923] | 0.909 |
| the same excl. 2020 Q1 - 2021 Q1 | 0.755 | 0.875 (m10 0.558) |

- **Latent volatility.** tau_t (posterior median of each decade's average):
  1960s 4.9, 1970s 4.6, 1980s 3.4, 1990s 2.5, 2000s 1.6, 2010s 1.4,
  2020s 4.0; 2020 Q2 12.0 [8.3, 19.1], 2020 Q3 11.6. Real activity was more
  volatile before the mid-1980s and much less volatile in the 2000s-2010s
  than m12's constant tau of 4 allows, and COVID is a large real shock.
- **Latent path.** The median path differs from m12's by 0.38 pp
  (annualised) on average, at most 1.7 (1963 Q3). 2020 Q2 mu is -27.9
  (m12 -27.4, headline -28.1).
- **Decomposition.** With a time-varying latent sd, a quiet-period quarter's
  disagreement between the measures is attributed more to measurement error:
  on 1978 Q4 - 2026 Q2 excluding COVID, m13's latent path has 75.5% of
  headline variance against m12's 87.5% (m10, with hours, 55.8%). On the
  full sample the two nearly agree because the 1960s-70s volatility and
  COVID dominate.
- **Prior sensitivity.** Seed 1 refitted with omega_g ~ half-N(0, 1)
  (`Rscript R/fit_m13.R --omega_g_sd=1 --out=output/fit_m13_omegag1.rds`)
  passes the gate (4 divergences, E-BFMI 0.839, max R-hat 1.0019 on
  hyperparameters and 1.0018 on g, min bulk ESS 1577; 1441 s). omega_g goes
  from 0.61 [0.42, 0.89] to 0.66 [0.44, 1.02], phi_g from 0.91 to 0.90, tau
  from 2.47 to 2.43. The median tau_t path changes by 1.9% on average (at
  most 7.4%, 2019 Q3; relative change of each quarter's median); 2020 Q2 tau_t 12.0 -> 12.3; decade means change by at
  most 0.1. The latent median path moves by 0.033 pp on average (at most
  0.18, 1969 Q1); the variance ratio goes from 0.790 to 0.786 (full sample)
  and from 0.755 to 0.753 (1978 Q4 - 2026 Q2 excluding COVID). The results
  are insensitive to this widening. (omega_g's posterior median sits at the
  78th percentile of the half-N(0, 0.5) prior; m10's 0.94 sat higher.)

### m13: LFO

`Rscript R/03_compare.R m13`: same window (1990 Q1 - 2026 Q2, 146 targets,
initial training 121 quarters), refit rule, 4 attempts, particle schedule,
4 filter passes and Newey-West lag 8 as every other model.

- **Fits.** 6 (the initial fit and refits at 1:146, 1:171, 1:204, 1:223 and
  1:243, the last on data to 2020 Q2). All 6 meet the refit rule on the first
  attempt; 0 divergences, min E-BFMI 0.81, max R-hat 1.013, min bulk ESS 222;
  max Pareto k used 0.696, no high-k step without a refit. Two fits are weak
  mixing (bulk ESS below 400): 1:171 (309; highest R-hat on phi_g) and 1:204
  (222; highest R-hat on phi[3]). The fit diagnostics record the lowest bulk
  ESS but not the variable that has it.
  m13 needs one refit after COVID; m12 also needed a 2020 Q4 refit, which
  passed only on its third attempt.
- **Filter.** The RBPF in `R/lib/filter_m13.cpp` with the locally adapted g
  proposal (defensive share 0.1) and the pf_schedule particle schedule. On
  12 posterior draws with 4 runs each, at fixed N = 2000 the adapted proposal
  beats the bootstrap one (COVID run-sd 0.65 against 1.38 nats, bias -0.27
  against -0.96; min pilot ESS fraction at 2020 Q2 0.0028 against 0.0006);
  under the schedule they are close. The schedule raised N in 23 of 267
  quarters (100000 only at 2020 Q2; cost x1.31). The same 12-draw check
  against a 2-run N = 50000 reference gives a whole-series bias of about
  -0.6 to -0.7 nats per draw over 267 quarters (COVID quarters -0.19); that
  is on the full-sample likelihood, not the LFO targets (below).
- **Monte Carlo error.** Pass totals -1046.05 (main), -1045.17, -1046.10,
  -1045.96: mean -1045.82, MCSE 0.22. MCSE of m13 - m12 0.29. The main pass
  alone gives m13 - m12 = +59.83 (iid SE 25.08, NW SE 29.24).
- **Particle bias against high-particle references** (`pfcheck`, m12's
  reference settings): over the COVID quarters 2020 Q1 - 2021 Q1 the summed
  bias is +0.069 under the schedule and -0.173 at N = 2000 (m12: +0.160 and
  -0.374); the 2020 Q2 reference is -23.205 (se 0.109, 3 passes) against
  m12's -47.03. Over the 20 targets 2015-2019 it is -0.028 under the
  schedule (per target -0.015 to +0.017) and -0.035 at N = 2000 (m12
  -0.017). Net bias on m13 - m12 is about -0.09 (COVID) and -0.011
  (2015-19): negligible against a 59.6-nat difference. Only 2020 Q2 has 3
  reference passes; every other target has 1.
- **Process.** The LFO process (20:21-22:25 on 2026-10-17) exited with
  status 1: "Error while shutting down parallel: unable to terminate some
  child processes", printed at R shutdown after every output had been
  written. The summary was rerun and exited cleanly; logged in
  `output/lfo/lfo_log.txt`. Compute used about 3 h 55 min including the
  particle checks.
- **Existing models untouched.** Every elpd_lfo and GDP-only score in
  `model_comparison.csv` and all 206 earlier rows of
  `model_comparison_pairwise.csv` are identical to commit 277f609;
  `lfo_retry_sensitivity.csv` is byte-identical.

Group 1 (E, I, P):

| model | elpd_lfo | SE | refits | retried | below rule | weak mixing | MCSE |
|---|---|---|---|---|---|---|---|
| **m13** | **-1045.82** | 22.7 | 5 | 0 | 0 | 2 | 0.22 |
| m12 | -1105.47 | 42.7 | 7 | 1 | 0 | 3 | 0.20 |
| m9 | -1122.24 | 40.0 | 7 | 2 | 0 | 3 | 0.31 |
| m4 | -1123.72 | 40.0 | 6 | 0 | 0 | 3 | 0.02 |

(m5, m7, m6, m1, m3, m2 are unchanged; m13 leads each by 123.0-150.4,
NW SE 20.7-25.8.)

Key differences (`model_comparison_pairwise.csv`):

| a - b | target | diff | SE iid | z iid | SE NW(8) | z NW | MCSE |
|---|---|---|---|---|---|---|---|
| m13 - m12 | group 1 | **+59.65** | 24.74 | 2.41 | 28.89 | **2.06** | 0.29 |
| m13 - m9 | group 1 | +76.42 | 21.80 | 3.51 | 24.78 | 3.08 | 0.38 |
| m13 - m4 | group 1 | +77.90 | 21.89 | 3.56 | 23.74 | 3.28 | 0.22 |
| m13 - m5 | group 1 | +123.01 | 15.29 | 8.04 | 25.72 | 4.78 | 0.22 |
| m13 - m10 | GDP-only | +12.62 | 6.87 | 1.84 | 5.48 | 2.30 | 0.22 |
| m13 - m11 | GDP-only | +67.54 | 12.67 | 5.33 | 18.42 | 3.67 | |
| m13 - m8 | GDP-only | +114.86 | 32.98 | 3.48 | 44.91 | 2.56 | |

Newey-West lags 4 / 8 / 12 / 16 give SE 27.5 / 28.9 / 29.1 / 29.2 for
m13 - m12, so z 2.17 / 2.06 / 2.05 / 2.05: the z >= 2 condition holds at
every lag tried, but only just. m13 - m10 compares models with different
information (m10 also sees past hours and its fits start in 1978) and
different structure (m10 has independent, non-AR measurement errors; m11 has
no SV), so neither gap measures the value of the labour data: no m13 + hours
model exists.

Seed-selection (retry) sensitivity for m13 - m12
(`lfo_retry_sensitivity_summary.csv`). m13 had no retried refit, so only
the 17 targets served by m12's retried 1:245 fit change:

| variant | diff | SE iid | SE NW | z NW |
|---|---|---|---|---|
| as reported (kept attempts) | +59.65 | 24.74 | 28.89 | 2.06 |
| first attempt only | +59.79 | 24.75 | 29.00 | 2.06 |
| all attempts pooled | +59.60 | 24.74 | 28.90 | 2.06 |
| the 53 retried targets dropped (93 left) | +42.17 | 24.53 | 28.78 | 1.47 |

The "dropped" row removes evidence (including 2020 Q2's +23.7) rather than
measuring selection, as for m12 - m9. (In this file `n_targets_rescored`
counts all 53 retried targets for every pair, including pairs where neither
model was retried; that was already so before, and the values are right.)

By sub-period (`model_comparison_subperiod.csv`, which now has a `target`
column: earlier rows are "group 1 joint", the new m13 - m10, m12 - m10 and
m10 - m9 rows "GDP-only"):

| period | n | m13 - m12 (SE iid / NW) | m13 - m9 (SE iid / NW) | m13 - m10, GDP-only (SE iid / NW) |
|---|---|---|---|---|
| 1990-99 | 40 | +1.31 (0.81 / 0.54) | +8.00 (3.16 / 2.34) | -0.54 (3.86 / 3.56) |
| 2000-07 | 32 | +3.99 (1.11 / 0.86) | +13.23 (2.22 / 1.29) | +6.96 (2.18 / 1.21) |
| 2008-19 | 48 | +14.78 (3.35 / 2.76) | +18.72 (3.65 / 3.74) | +6.32 (2.69 / 2.18) |
| 2020 | 4 | +30.07 (22.39 / 8.28*) | +25.89 (19.63 / 7.72*) | -2.32 (4.09 / 1.84*) |
| 2021-22 | 8 | +2.51 (2.97 / 1.58) | +2.92 (2.92 / 1.83) | +1.84 (2.02 / 1.51) |
| 2023-26 | 14 | +6.99 (0.71 / 0.89) | +7.66 (1.00 / 1.09) | +0.37 (1.29 / 0.82) |
| pre-2020 | 120 | +20.09 (3.83 / 4.86) | +39.95 (5.37 / 5.48) | +12.74 (5.24 / 5.02) |
| all excl. COVID (2020 Q2 - 2021 Q1) | 142 | +28.53 (5.19 / 6.02) | +49.06 (6.44 / 6.26) | +14.92 (5.73 / 5.58) |
| **all** | 146 | **+59.65 (24.74 / 28.89)** | +76.42 (21.80 / 24.78) | +12.62 (6.87 / 5.48) |

\* With n = 4 the Newey-West SE is not meaningful.

m13 is ahead of m12 in 122 of 146 quarters (median gain +0.17 per quarter).
Largest gains: 2020 Q2 +23.7, 2020 Q4 +5.5, 2020 Q3 +2.3. Largest losses:
2016 Q3 -2.0, 2020 Q1 -1.5, 2008 Q4 -1.0. Outside 2020 the gain per
quarter grows over the window as the latent volatility falls further below
m12's constant tau: about +0.03 per quarter in 1990-99, +0.12 in 2000-07,
+0.31 in 2008-19 and +0.50 in 2023-26.

### Production decision: m12 -> m13

The rule was fixed before m13's LFO results were known. Because m13 is a
post-hoc candidate (above), it adds a fifth condition to the four used for
m9 -> m12: switch production to m13 only if (1) it passes the production
gate on two seeds, (2) its LFO score is reliable, (3) it beats m12 by
z_NW >= 2, (4) it still does so using each LFO refit's first attempt only,
and (5) a local fit takes at most 30 minutes.

| Condition | Observed | Met |
|---|---|---|
| gate on two seeds | seeds 1 and 2 both pass (table above) | yes |
| reliable LFO | MCSE of m13 - m12 0.29 (4 passes each); measured bias about -0.1; all 6 m13 fits meet the refit rule on attempt 1 | yes |
| z_NW >= 2 | +59.65 / 28.89 = 2.06 (lags 4-16: 2.05-2.17; main pass alone 59.83 / 29.24 = 2.05) | yes, narrowly |
| z_NW >= 2, first attempts only | +59.79 / 29.00 = 2.06 | yes, narrowly |
| local fit <= 30 min | 23.1 min (1384 s wall, seed 1, 4 parallel chains); 27.9 min (1674 s wall, 1664 s sampling) for seed 2 at load average 10-11 | yes; thin margin under load |

The comparison author recommended making m13 the leading candidate rather
than switching automatically, given the thin margins and the post-hoc
origin. The pre-set rule decides, and all its conditions are met, so
production uses `stan/m13_serial_latent_sv.stan` from 2026-10-18:
`R/fit_production.R` (gate hyperparameters add mu_g, phi_g, omega_g; R-hat
<= 1.05 also on every g[t]), `R/output_estimate.R` (model string,
hyperparameter list), `R/input_hash.R` (the m13 Stan file replaces m12's),
`R/build_site.R` (model description, equations, priors, comparison text)
and the CI workflow (cached executable `stan/m13_serial_latent_sv`, key on
its .stan hash; timeout unchanged at 120 minutes).

The production fit on the current data (`Rscript R/fit_production.R`, seed 1,
2026-10-18, load average about 4) reproduces the seed-1 research fit draw
for draw (identical mu, g and omega_g draws): 4 divergences in 8000, E-BFMI
0.854-0.969, max R-hat 1.0020 (omega[1]) on hyperparameters and 1.0011 /
1.0018 / 1.0020 on mu / h / g, min bulk / tail ESS 1553 / 3213 on
hyperparameters and 4344 / 4866 on mu, 0 treedepth hits; gate passed.
Sampling took 1173 s (chains 856-1173 s), about 20 minutes wall. A second
production run (another process, about 00:39-00:59 on 2026-10-18) then
rewrote `output/fit_production.rds` and `docs/`: the same seed gave the same
draws (mu and omega_g identical to the seed-1 research fit), sampling took
1177 s (chains 859-1177 s), and the input hash is unchanged.
`R/output_estimate.R` and `R/build_site.R` were then rerun. A later review
fix made the site's posterior numbers in the model description (tau_t by
decade and at 2020 Q2, the psi medians) come from the fit:
`R/output_estimate.R` now also writes `docs/model_summary.json`, which
`R/build_site.R` fills into the text. The other models' variance splits in
"Why this specification" remain research-suite numbers, labelled as such.
The estimates are unchanged (`estimates.csv` and `estimates_latest.json`
byte-identical); because both scripts changed, the input hash changed, and
the hash in `docs/diagnostics.json` equals `Rscript R/input_hash.R`
(37ae6eaa...). The next weekly CI run will therefore refit once even
without new data.

Published numbers, m13 against the last m12 publication (quarter-on-quarter
%, posterior median and 90% interval):
- 2026 Q2: 0.48 [0.30, 0.63] against 0.47 [0.29, 0.62]; headline 0.42;
- 2020 Q1-Q4: -0.26 / -6.98 / +3.43 / +3.50 against -0.33 / -6.84 / +3.38 /
  +3.52;
- the median path differs by 0.10 pp on average (at most 0.43);
- variance ratio 0.790 [0.705, 0.889] against 0.779 [0.707, 0.850].

What does and does not change for users: the full-sample variance ratio
moves from 0.779 to 0.790, and the published path moves by 0.38 pp
(annualised; about 0.1 pp quarter-on-quarter) on average.

Open items:
- the fit's run time on GitHub's runners has not been measured (m12's was
  not either); at 23-28 minutes locally the 120-minute job limit leaves
  room, but the 30-minute local condition is close under load;
- the 2015-2019 particle check has 1 reference pass per target; only
  2020 Q2 has 3;
- two m13 LFO fits (1:171, 1:204) have bulk ESS of 309 and 222, below the
  production level;
- the post-hoc caveat cannot be removed on this evaluation window; a real
  test is the out-of-sample record from the next ABS releases.

## Round of 2026-10-17: m12, m3 rewrite, m10 re-specification, stricter LFO, production switch

This section supersedes the 2026-10-16 comparison tables below (they are kept,
marked superseded, so the history stays readable). The data vintage is
unchanged (2026 Q2, T = 267).

### Summary

- **m12** (m9 + AR(1) measurement errors) is new. It wins group 1 on
  PSIS-LFO-CV: -1105.47 against -1122.24 for m9 (+16.77, Newey-West SE 6.12,
  z 2.74).
- **Production switched from m9 to m12**, under a rule fixed before m12's
  results were known (below). *(Superseded 2026-10-18: production is now
  m13; see "Round of 2026-10-18".)*
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
  They would not pass the production gate. In brackets: the lowest bulk ESS
  and the variable with the highest R-hat (`worst_rhat`); the table does not
  record which variable has the lowest ESS:
  - m12: data to 1994 Q2 (ESS 336, phi[3]), 2020 Q2 (210, tau), 2020 Q4
    (125, tau);
  - m9: data to 2000 Q3 (231, phi[3]), 2020 Q2 (210, tau), 2020 Q4 (106,
    tau);
  - m4: data to 2003 Q2 (332, omega[2]), 2008 Q2 (399.9, phi[3]), 2020 Q2
    (117, tau; R-hat 1.033);
  - m3: data to 2020 Q2 (391, Omega[2,1]).

  The highest R-hat is on tau on the fits just after COVID (all three SV
  models) and on an SV persistence or volatility parameter on some pre-COVID
  fits.
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

### Comparison results (2026-10-17; m13 added 2026-10-18)

*These scores are unchanged, but m13 (-1045.82) now leads group 1 and the
GDP-only target; "diff vs m12" below is no longer the difference from the
best model. See "Round of 2026-10-18" for m13's rows.*

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
   *(2026-10-18: m13 tested that suggestion; it is ahead of m10 by +12.6 on
   GDP-only, z_NW 2.3.)*
4. **The comparison still does not test the variance decomposition.** Over
   the same 1978 Q4 - 2026 Q2 quarters, m10's latent variance ratio is 0.75
   [0.68, 0.82] and m12's 0.91 [0.84, 0.97] (0.56 and 0.88 excluding
   2020 Q1 - 2021 Q1); m10 predicts GDP somewhat better. "About a fifth of
   headline variance is measurement error" is m12's answer, not a result the
   comparison validates. A model with latent SV (m10) puts the
   measurement-error share much higher outside COVID. *(2026-10-18: m13,
   with latent SV on the full sample, gives 0.755 on those quarters
   excluding COVID, between m10 and m12.)*

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

### Production decision: m9 -> m12 *(superseded 2026-10-18 by m12 -> m13)*

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
- m9 + latent SV (as in the new m10) is the natural next candidate. *(Done
  2026-10-18 as m12 + latent SV: m13, now production.)*

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
   m9 by +16.8, NW z 2.7, and was the production model from 2026-10-17 to
   2026-10-18, when m13 replaced it.)*
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

## Headline replication result (per-draw variance ratios, 2026-10-16; m3, m10, m12 updated 2026-10-17; m13 added 2026-10-18)

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
| Latent mu (m12, production 2026-10-17 to 2026-10-18) | 0.779 | 0.707 | 0.850 |
| Latent mu (m13, production from 2026-10-18) | 0.790 | 0.705 | 0.889 |

Hours subsample (T = 191, 1978 Q4 - 2026 Q2; headline sd 3.95):

| Series | median | 5% | 95% |
|---|---|---|---|
| Latent mu (m8)  | 0.848 | 0.767 | 0.930 |
| Latent mu (m10, re-specified with latent SV) | 0.747 | 0.678 | 0.824 |
| Latent mu (m10) *(superseded: old spec, mode A)* | 0.347 | 0.263 | 0.487 |
| Latent mu (m9), same quarters  | 0.899 | 0.829 | 0.968 |
| Latent mu (m12), same quarters | 0.909 | 0.840 | 0.973 |
| Latent mu (m13), same quarters | 0.855 | 0.790 | 0.923 |

The m9, m12 and m13 rows are their full-sample fits restricted to the
1978 Q4 - 2026 Q2 quarters (computed 2026-10-17, m13 2026-10-18; not in the
CSV). Excluding 2020 Q1 - 2021 Q1 the medians are m10 0.558, m9 0.861, m12
0.875, m13 0.755.

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
0.779 [0.707, 0.850]; m13: 0.790 [0.705, 0.889]) and about 55% less than
GDP(E). Excluding the COVID quarters the share differs more between models
(full sample: m13 0.731, m12 0.727; 1978 Q4 onward: m13 0.755, m12 0.875).

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
  - m4, m9, m10, m11, m12, m13: adapt_delta 0.99, max_treedepth 12 (m13's
    are the production settings, as m12's and m9's were).
- **LFO refits:** 4 x (1000 + 1000), with each model's own adapt_delta and
  max_treedepth, and the refit rule described in "Round of 2026-10-17".
- **Diagnostics.** The research suite reports diagnostics and never stops on
  them. Only `R/fit_production.R` enforces a gate:
  - divergences < 0.1%;
  - min E-BFMI >= 0.3;
  - R-hat <= 1.01 on hyperparameters and <= 1.05 on mu[t], h[t,k] and g[t];
  - bulk and tail ESS >= 400 on hyperparameters and mu[t].

  The hyperparameters are c, rho, tau, mu_h, phi, omega, rho_EI, psi,
  mu_bar, mu_g, phi_g and omega_g (psi added with the switch to m12; mu_g,
  phi_g, omega_g and the g[t] R-hat check with the switch to m13, where tau
  is the generated quantity exp(mu_g / 2)).

## Files

**Data**
- `R/01_fetch_data.R`: pulls ABS 5206.0 Table 24 (production data).
- `R/01b_fetch_research_data.R`: pulls 5206.0 Table 1 hours and 6202.0
  employment (research only; run after 01).
- `R/utils_abs.R`: retry wrapper and ABS connection-check workaround, shared
  by both fetch scripts.

**Research suite**
- `R/02_fit_models.R`: fits m1-m10, m12 and m13 (`Rscript R/02_fit_models.R m9 m12`
  for a subset; `--diag-only m13` tabulates a saved fit). Writes
  `output/tables/fit_diagnostics.csv`.
- `R/fit_m11.R`: fits m11.
- `R/fit_m12.R`: fits m12 with the production gate (any seed), and the
  warmup-stall harness.
- `R/fit_m13.R`: the same for m13, plus `--omega_g_sd=S` for the
  prior-sensitivity fit; writes `output/tables/m13_gate_seeds.csv`.
- `R/03_compare.R`: PSIS-LFO-CV comparison and per-draw variance
  decomposition.
- `R/lib/filters.R`, `R/lib/filters.cpp`, `R/lib/lfo.R`: Kalman filter,
  Rao-Blackwellised particle filter and the LFO algorithm.
- `R/lib/filter_m12.R`, `R/lib/filter_m12.cpp`: m12's Kalman filter and
  particle filter (registered by `register_m12()` in `R/03_compare.R`).
- `R/lib/filter_m10.R`, `R/lib/filter_m10.cpp`: the re-specified m10's
  particle filter over (g, h) (registered by `register_m10_filter()`).
- `R/lib/filter_m13.R`, `R/lib/filter_m13.cpp`: m13's Kalman filter and
  particle filter over (g, h) (registered by `register_m13()`).
- `tests/test_filters.R`, `tests/test_m12.R`, `tests/test_m13.R`,
  `tests/test_m3_m10.R`: filter tests, each against a dense
  multivariate-normal likelihood.
- `R/diag_m10_modes.R`: m10 mode check (16 dispersed chains; 4 chains started
  in old mode A) -> `output/tables/m10_mode_check.csv`.
- `R/diag_m3_ident.R`: m3 identification checks ->
  `output/tables/m3_identification.csv`.
- `R/04_figures.R`: figures 01-07 in `output/figures/`.
- `R/sample_window.R`: the optional sample start (`LATENTGDP_START`) and
  output root (`LATENTGDP_OUT`) for the research scripts. With both unset,
  behaviour is unchanged. Production never sources it. See "1980Q1 sample
  (RDP 2014-12 alignment)".

**Production**
- `R/fit_production.R`, `R/output_estimate.R`, `R/build_site.R`: production
  fit (m13 since 2026-10-18; m12 from 2026-10-17; m9 before), estimate
  writer and site.
- `R/input_hash.R`: CI input fingerprint.

**Models**
- `stan/m{1..13}_*.stan`: the 13 specifications.
  - m1-m5: baseline family motivated by Rees-Lancaster-Finlay. RDP Model 1
    (AR(1) latent, independent errors) corresponds to m2; m1 has an iid
    latent state. RDP Model 2 (latent innovation correlated with the errors,
    zeta = 0.5 restriction) has no counterpart: m3 frees the error
    correlations instead. RDP Model 3 adds the unemployment rate.
  - m6: structural Sigma restriction (a cheap novel extension).
  - m7: Lenza-Primiceri COVID variance scale.
  - m8: 4-signal, in the spirit of RDP Model 3 (the paper has three
    models), but with hours worked as the fourth signal where the paper
    uses the unemployment rate.
  - m9: m6 + SV (production until 2026-10-17; tied with m4 on LFO).
  - m10: m8 + SV on the measurement errors and on latent growth
    (re-specified 2026-10-17; the old constant-tau version was multimodal).
  - m11: mixed-frequency, Banbura-Modugno style, with monthly employment.
  - m12: m9 + AR(1) measurement errors (production 2026-10-17 to
    2026-10-18).
  - m13: m12 + SV on latent growth (production from 2026-10-18).

**Outputs**
- `output/fit_m{1..13}.rds`, `output/fit_production.rds`: saved CmdStanFit
  objects (gitignored). `output/fit_m12_gate.rds` and
  `output/fit_m12_seed2_gate.rds` are the seed-1 and seed-2 gate records
  (also gitignored; their content is in `output/tables/m12_gate_seeds.csv`);
  likewise `output/fit_m13_gate.rds`, `output/fit_m13_seed2_gate.rds` and
  `output/tables/m13_gate_seeds.csv`. `output/fit_m13_omegag1.rds` is the
  prior-sensitivity fit.
- `output/lfo/`: LFO log (`lfo_log.txt`, committed) and checkpoints (`*.rds`,
  gitignored, regenerable by `R/03_compare.R`). `output/lfo/superseded/`
  holds the checkpoints of the 2026-10-16 run; they are local only and also
  gitignored. `R/03_compare.R` reads two of them (m4, m9) as a refit cache
  (`REFIT_CACHE`); without them the same refits are recomputed from their
  seeds, with identical draws.
- `output/tables/`: comparison, diagnostics and variance-decomposition tables.
- `output/figures/`: generated PNGs.
- `output/sample_1980/`: the 1980 Q1 research run. It holds:
  - the fits (`fit_m*.rds`, m13's two seeds and gate records);
  - `tables/` (the same tables as `output/tables/`, plus
    `posterior_comparison_full_vs_1980.csv`, `lfo_full_vs_sample.csv` and
    `lfo_full_vs_sample_models.csv`);
  - `figures/` (01-07);
  - `lfo/` (checkpoints, `lfo_log.txt`, `run_lfo_1980.sh`,
    `run_lfo_1980_extra.sh`);
  - `logs/`;
  - the fits driver `run_1980.sh`.

  Its `*.rds` files (about 1.7 GB) are gitignored (`output/sample_*/**/*.rds`).

**Documentation**
- `priors.md`: every prior, as coded.
