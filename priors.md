# Priors

Every prior below is transcribed from the `model` block of the Stan file named
in each section (checked against the files on 2026-10-17). If a Stan file and
this document disagree, the Stan file is right and this document needs fixing.

**Units.** Every model works in annualised growth, 400 × Δlog(GDP), in percent
(m11's monthly series are 1200 × Δlog, also annualised). Since 1959, Australian
real GDP growth has averaged about 3% annualised, with quarter-to-quarter
swings of a few percentage points. The priors are weakly informative on that
scale.

**Notation.**
- `N(m, s)` is a normal with mean m and **standard deviation** s.
- `half-N(0, s)` is a `normal(0, s)` statement on a parameter declared
  `<lower=0>`.
- A normal on a parameter declared `<lower=-1, upper=1>` is that normal
  truncated to (-1, 1). Stan applies it without a normalising constant, which
  is fine because the bounds are fixed.
- `std_normal` on `mu_raw` / `h_raw` is the non-centred form of the state
  equation. The priors it implies on the latent states are written out
  wherever they apply.

**Prior predictive checks.** None have been run. There is no prior predictive
simulation anywhere in this repository. The scales were chosen by hand from
the sample moments above, not checked by simulation.

## Shared AR(1) latent block (m2-m9, m12; m10 in part)

m2, m3, m4, m5, m6, m7, m8, m9 and m12 all use the same latent-state equation
and the same priors on it. m10 uses the same `c` and `rho` but, since
2026-10-17, replaces the constant `tau` with stochastic volatility on the
latent innovation (see its section):

| Parameter | Declared | Prior | Note |
|---|---|---|---|
| `c`   | real | N(1.5, 1) | intercept. `mu_bar = c / (1 - rho)` is a generated quantity; at rho = 0.5 the prior centre corresponds to mu_bar = 3 |
| `rho` | (-1, 1) | N(0.5, 0.3) | AR(1) persistence of latent growth |
| `tau` | > 0 | half-N(0, 2) | sd of latent-growth innovations |

The implied prior on the latent states is mu_1 ~ N(c/(1-rho), tau²/(1-rho²))
(the stationary distribution) and mu_t | mu_{t-1} ~ N(c + rho·mu_{t-1}, tau²).
- m2, m5, m6, m7 and m8 implement it with `mu_raw ~ std_normal()`.
- m3, m4, m9, m10 and m12 integrate mu out with a Kalman filter that starts
  from the same stationary initial distribution. This is the same prior (for
  m10, with tau² replaced by exp(g_t)).

## m1 baseline (`stan/m1_baseline.stan`)

| Parameter | Declared | Prior |
|---|---|---|
| `mu_bar`   | real | N(3, 2) |
| `tau`      | > 0  | half-N(0, 2) |
| `sigma[k]`, k = E, I, P | > 0 | half-N(0, 2) |
| `mu[t]`    | real | N(mu_bar, tau), iid over t |

## m2 AR(1) latent (`stan/m2_ar1.stan`)

| Parameter | Declared | Prior |
|---|---|---|
| `c`, `rho`, `tau` | | shared AR(1) block |
| `sigma[k]` | > 0 | half-N(0, 2) |

## m3 correlated measurement errors (`stan/m3_corr.stan`)

| Parameter | Declared | Prior |
|---|---|---|
| `c`, `rho`, `tau` | | shared AR(1) block |
| `sigma[k]` | > 0 | half-N(0, 2) |
| `L_Omega`  | cholesky_factor_corr[3] | `lkj_corr_cholesky(4)`, i.e. Omega ~ LKJ(4) |

LKJ(2) was tried first. It sampled badly, so the prior was tightened to LKJ(4)
(see NOTES.md, "Identification issue in m3"). Since 2026-10-17 mu is
integrated out by a Kalman filter instead of being sampled through `mu_raw`;
the priors are unchanged.

## m4 stochastic volatility (`stan/m4_sv.stan`)

| Parameter | Declared | Prior |
|---|---|---|
| `c`, `rho`, `tau` | | shared AR(1) block |
| `mu_h[k]`  | real | N(2, 1). Long-run log-variance; exp(2/2) ≈ 2.7 pp sd |
| `phi[k]`   | (-1, 1) | Uniform(-1, 1) |
| `omega[k]` | > 0 | half-N(0, 0.5) |
| `h_raw[t,k]` | real | std_normal |

The implied prior on the log-variances is a stationary AR(1):
- h_{1,k} ~ N(mu_h_k, omega_k²/(1-phi_k²));
- h_{t,k} | h_{t-1,k} ~ N(mu_h_k + phi_k(h_{t-1,k} - mu_h_k), omega_k²).

The measurement-error sd is exp(h_{t,k}/2).

## m5 serially correlated measurement errors (`stan/m5_serial.stan`)

| Parameter | Declared | Prior |
|---|---|---|
| `c`, `rho`, `tau` | | shared AR(1) block |
| `sigma[k]` | > 0 | half-N(0, 2). sd of the measurement-error innovation u |
| `psi[k]`   | (-1, 1) | N(0, 0.5). AR(1) coefficient of the measurement error |

## m6 structural Sigma (`stan/m6_structural_corr.stan`)

| Parameter | Declared | Prior |
|---|---|---|
| `c`, `rho`, `tau` | | shared AR(1) block |
| `sigma[k]` | > 0 | half-N(0, 2) |
| `rho_EI`   | (-1, 1) | N(0, 0.5) |

Cov(eps_P, eps_E) = Cov(eps_P, eps_I) = 0 is a restriction in the model
structure, not a prior.

## m7 COVID variance scale (`stan/m7_covid_scale.stan`)

| Parameter | Declared | Prior |
|---|---|---|
| `c`, `rho`, `tau` | | shared AR(1) block |
| `sigma[k]` | > 0 | half-N(0, 2) |
| `kappa[q]`, q = 2020Q2, 2020Q3, 2020Q4, 2021Q1 | > 0 | LogNormal(0, 1). Median 1 |

In every other quarter, kappa_t = 1.

## m8 four signals: E, I, P and hours (`stan/m8_labour.stan`)

| Parameter | Declared | Prior |
|---|---|---|
| `c`, `rho`, `tau` | | shared AR(1) block |
| `sigma_gdp[k]` | > 0 | half-N(0, 2) |
| `sigma_U` | > 0 | half-N(0, 2) |
| `alpha`   | real | N(-1, 1). Hours growth = alpha + beta·mu, so alpha absorbs trend productivity growth |
| `beta`    | > 0 | N(1, 0.5) truncated to beta > 0 (sign restriction) |

## m9 structural Sigma + SV (`stan/m9_struct_sv.stan`; production until 2026-10-17)

| Parameter | Declared | Prior |
|---|---|---|
| `c`, `rho`, `tau` | | shared AR(1) block: N(1.5, 1), N(0.5, 0.3) on (-1, 1), half-N(0, 2) |
| `mu_h[k]`  | real | N(2, 1) |
| `phi[k]`   | (-1, 1) | Uniform(-1, 1) |
| `omega[k]` | > 0 | half-N(0, 0.5) |
| `rho_EI`   | (-1, 1) | N(0, 0.5), the same prior as m6 |
| `h_raw[t,k]` | real | std_normal. Implies the same stationary AR(1) prior on h as m4 |

- **rho_EI.** It is sampled as `z_EI = atanh(rho_EI)` and written as
  `target += normal_lpdf(tanh(z_EI) | 0, 0.5) + log(1 - rho_EI^2)`, with the
  Jacobian computed stably from `z_EI`. This is a change of coordinates only:
  the prior on rho_EI is the same N(0, 0.5) on (-1, 1).
- **mu.** mu is integrated out by the Kalman filter and drawn by FFBS. Its
  prior is the shared stationary AR(1) prior.

## m10 four signals + SV, with latent SV (`stan/m10_labour_sv.stan`)

Re-specified on 2026-10-17 (NOTES.md, "m10: two modes, diagnosis and
re-specification"). The latent innovation sd is now time-varying,
exp(g_t / 2), with its own stationary AR(1) log-variance g_t. The new latent
SV block takes the same priors as each measurement SV block, and replaces the
earlier `tau ~ half-N(0, 2)`.

| Parameter | Declared | Prior |
|---|---|---|
| `c`, `rho` | | shared AR(1) block: N(1.5, 1), N(0.5, 0.3) on (-1, 1) |
| `mu_g`    | real | N(2, 1). Long-run log-variance of the latent innovation |
| `phi_g`   | (-1, 1) | Uniform(-1, 1) |
| `omega_g` | > 0 | half-N(0, 0.5) |
| `g_raw[t]` | real | std_normal. Implies g_1 ~ N(mu_g, omega_g²/(1-phi_g²)), g_t \| g_{t-1} ~ N(mu_g + phi_g(g_{t-1} - mu_g), omega_g²) |
| `alpha` | real | N(-1, 1) |
| `beta`  | > 0 | N(1, 0.5) truncated to beta > 0 |
| `mu_h[k]`, k = E, I, P, hours | real | N(2, 1) |
| `phi[k]`   | (-1, 1) | Uniform(-1, 1) |
| `omega[k]` | > 0 | half-N(0, 0.5) |
| `h_raw[t,k]` | real | std_normal. Stationary AR(1) prior on h, as in m4 |

- **Latent states.** mu_1 ~ N(c/(1-rho), exp(g_1)/(1-rho²)) and
  mu_t | mu_{t-1} ~ N(c + rho·mu_{t-1}, exp(g_t)). mu is integrated out by a
  Kalman filter and drawn by FFBS.
- **tau.** `tau = exp(mu_g / 2)` and `tau_t = exp(g_t / 2)` are generated
  quantities, not parameters, and have no prior of their own.
- **omega_g.** Its posterior (mean 0.94, 90% interval [0.63, 1.33]) sits in
  the upper tail of the half-N(0, 0.5) prior, which puts 6.6% of its mass
  above 0.92. No prior-sensitivity check has been run.
- **Superseded.** Before 2026-10-17, m10 had a constant `tau ~ half-N(0, 2)`
  from the shared block.

## m11 mixed frequency (`stan/m11_mixed_freq.stan`)

The latent state is monthly. These priors are at monthly frequency and differ
from the shared block:

| Parameter | Declared | Prior |
|---|---|---|
| `c`   | real | N(0.3, 0.5). Monthly intercept; `mu_bar_m = c / (1 - rho)` |
| `rho` | (-1, 1) | N(0.9, 0.1). Monthly persistence, high a priori |
| `tau` | > 0 | half-N(0, 5). Monthly innovation sd (annualised %) |
| `sigma_q[k]` | > 0 | half-N(0, 2). Quarterly GDP measurement errors |
| `alpha_emp` | real | N(-1.5, 1.5) |
| `beta_emp`  | > 0 | N(1, 0.5) truncated to beta_emp > 0 |
| `sigma_emp` | > 0 | half-N(0, 5) |
| `mu_raw[m]` | real | std_normal, giving a stationary AR(1) monthly latent |

## m12 structural Sigma + SV + AR(1) measurement errors: production (`stan/m12_struct_sv_serial.stan`)

Production model since 2026-10-17. Every prior is m9's, unchanged, plus m5's
prior on the measurement-error AR(1) coefficients.

| Parameter | Declared | Prior |
|---|---|---|
| `c`, `rho`, `tau` | | shared AR(1) block: N(1.5, 1), N(0.5, 0.3) on (-1, 1), half-N(0, 2) |
| `mu_h[k]`  | real | N(2, 1) |
| `phi[k]`   | (-1, 1) | Uniform(-1, 1) |
| `omega[k]` | > 0 | half-N(0, 0.5) |
| `rho_EI`   | (-1, 1) | N(0, 0.5), written on `z_EI = atanh(rho_EI)` with the Jacobian, as in m9 |
| `psi[k]`   | (-1, 1) | N(0, 0.5), the same prior as m5 |
| `h_raw[t,k]` | real | std_normal. Stationary AR(1) prior on h, as in m4 and m9 |

- **What h is.** sigma_{t,k} = exp(h_{t,k}/2) is the sd of the measurement-error
  *innovation* u_t^k (eps_t^k = psi_k eps_{t-1}^k + u_t^k), as `sigma` is in m5.
  The error sd at the long-run log-variance is the generated quantity
  `eps_sd_bar = exp(mu_h/2) / sqrt(1 - psi²)`.
- **Initial errors.** eps_1 ~ N(0, V_1) with
  V_1[j,k] = R_1[j,k] / (1 - psi_j psi_k), the stationary covariance at the
  period-1 innovation covariance R_1. This is part of the likelihood, not a
  separate prior; at psi = 0 it is m9's period-1 density.
- **mu.** mu is integrated out by a Kalman filter on the quasi-differenced
  data and drawn by FFBS. Its prior is the shared stationary AR(1) prior.
