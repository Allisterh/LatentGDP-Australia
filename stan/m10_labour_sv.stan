// m10: 4-signal (E, I, P, hours) model with stochastic volatility on every
// measurement error AND on the latent innovation (re-specified 2026-10-16).
//   mu_t  = c + rho * mu_{t-1} + exp(g_t / 2) * eta_t
//   g_t   = mu_g + phi_g (g_{t-1} - mu_g) + omega_g zeta_t      (stationary AR(1))
//   mu_1  ~ N(c / (1 - rho), exp(g_1) / (1 - rho^2))
//   y_t^E = mu_t + eps_t^E,             eps ~ N(0, exp(h_{t,E}/2)^2)
//   y_t^I = mu_t + eps_t^I,             eps ~ N(0, exp(h_{t,I}/2)^2)
//   y_t^P = mu_t + eps_t^P,             eps ~ N(0, exp(h_{t,P}/2)^2)
//   y_t^U = alpha + beta * mu_t + eps_t^U,  eps ~ N(0, exp(h_{t,U}/2)^2)
//   h_{t,k} = mu_h_k + phi_k (h_{t-1,k} - mu_h_k) + omega_k xi_{t,k}  (stationary AR(1))
// Priors: c, rho, alpha, beta, mu_h, phi, omega as before; the latent SV block
// uses the same priors as each measurement SV block: mu_g ~ N(2, 1),
// phi_g ~ U(-1, 1), omega_g ~ half-N(0, 0.5). It replaces tau ~ half-N(0, 2).
// Output `tau` = exp(mu_g / 2) is the latent innovation sd at the long-run
// log-variance (a generated quantity, kept so summaries that ask for tau work);
// `tau_t` = exp(g_t / 2) is its path.
// Sample: post-1978 hours subsample (T = 191 at the 2026Q2 vintage).
//
// Why the latent SV (previous m10 = m8 + measurement SV only, constant tau).
// That posterior had two modes (16 chains from dispersed inits, 1000 + 1000):
//   A (9 of 16 chains): rho ~0.39, tau ~2.1, phi_E ~0.66, omega_E ~0.67,
//     lp__ ~ -2287. In 2020Q2 the latent mu is -1.3 (headline -28.1) and all
//     four measurement-error sds jump to ~10: COVID is read as four coincident
//     errors of about -27 in independently compiled series (expenditure, income,
//     industry production and Labour Force hours). Latent variance ratio 0.36.
//   B (3 chains): rho ~0.12, tau ~3.6, phi_E ~0.97,
//     omega_E ~0.22, lp__ ~ -2330. mu tracks COVID (-26.9) with measurement
//     sds ~1.5, but tau is inflated for every other quarter. Ratio 0.86.
//   The other 4 chains sat between them (tau 2.4-3.2, lp__ -2301 to -2310).
// Both are responses to one misspecification: with a constant tau the
// post-1978 sample (the Great Moderation) makes COVID a ~13-sigma latent shock,
// so the model must either inflate tau everywhere (B) or call a common shock
// measurement error (A). Time-varying volatility of true activity is standard
// (Great Moderation; COVID lockdowns were a real shock - the four measures agree
// within 0.5 pp in 2020Q2 and within ~6 pp including hours) and lets the model
// put a one-quarter jump in the latent volatility instead. With it, 16 of 16
// dispersed-init chains land in one mode (mu ~ -28.2 in 2020Q2, tau_t ~14,
// measurement sds ~1.4); see NOTES.md.
//
// Sampler parameterisation: conditional on g, h, alpha and beta the model is
// linear-Gaussian in mu, so mu is integrated out exactly with a Kalman filter in
// the model block (time-varying state variance exp(g_t)) and drawn exactly from
// p(mu | y, g, h, theta) by forward-filtering backward-sampling in generated
// quantities. g and h are non-centred. The filter and FFBS are checked against
// the dense multivariate-normal likelihood and smoother in tests/test_m3_m10.R.
// LFO filter: R/lib/filter_m10.R / filter_m10.cpp (particles over (g, h)).

functions {
  // Scalar-state Kalman filter with time-varying innovation variance qv[t]
  // (qv[1] enters only through P_1 = qv[1] / (1 - rho^2)) and independent
  // measurement errors with log-variances h[t, k]. Stable quadratic form
  // (Lagrange identity, as in m9). Returns [log p(y_t | y_{1:t-1}), m_{t|t}, P_{t|t}].
  matrix kf_tv(matrix y, matrix h, vector Z, vector d, real c, real rho, vector qv) {
    int T = rows(y);
    int K = cols(y);
    matrix[T, 3] out;
    real a = c / (1 - rho);
    real P = qv[1] / (1 - square(rho));
    real rho2 = square(rho);
    vector[K] z;
    vector[K] w;
    for (t in 1:T) {
      for (k in 1:K) {
        real u = exp(-0.5 * h[t, k]);
        z[k] = Z[k] * u;
        w[k] = (y[t, k] - d[k] - a * Z[k]) * u;
      }
      real s = dot_self(z);
      real q = dot_product(z, w);
      real Q = dot_self(w);
      real cross = 0;
      for (i in 1:(K - 1))
        for (j in (i + 1):K)
          cross += square(z[i] * w[j] - z[j] * w[i]);
      real g = 1 + P * s;
      out[t, 1] = -0.5 * (K * log(2 * pi()) + sum(h[t]) + log(g) + (Q + P * cross) / g);
      out[t, 2] = a + P * q / g;
      out[t, 3] = P / g;
      if (t < T) {
        a = c + rho * out[t, 2];
        P = rho2 * out[t, 3] + qv[t + 1];
      }
    }
    return out;
  }

  vector ffbs_tv_rng(matrix f, real c, real rho, vector qv) {
    int T = rows(f);
    vector[T] mu;
    mu[T] = normal_rng(f[T, 2], sqrt(f[T, 3]));
    for (tt in 1:(T - 1)) {
      int t = T - tt;
      real Pp = square(rho) * f[t, 3] + qv[t + 1];
      real J = rho * f[t, 3] / Pp;
      mu[t] = normal_rng(f[t, 2] + J * (mu[t + 1] - c - rho * f[t, 2]),
                         sqrt(f[t, 3] * qv[t + 1] / Pp));
    }
    return mu;
  }
}
data {
  int<lower=1> T;
  matrix[T, 3] y_gdp;
  vector[T] y_U;
}
transformed data {
  matrix[T, 4] y_all = append_col(y_gdp, y_U);
}
parameters {
  real c;
  real<lower=-1, upper=1> rho;
  real alpha;
  real<lower=0> beta;
  real mu_g;
  real<lower=-1, upper=1> phi_g;
  real<lower=0> omega_g;
  vector[4] mu_h;
  vector<lower=-1, upper=1>[4] phi;
  vector<lower=0>[4] omega;
  vector[T] g_raw;
  matrix[T, 4] h_raw;
}
transformed parameters {
  vector[T] g;
  matrix[T, 4] h;
  g[1] = mu_g + omega_g / sqrt(1 - square(phi_g)) * g_raw[1];
  for (t in 2:T)
    g[t] = mu_g + phi_g * (g[t - 1] - mu_g) + omega_g * g_raw[t];
  for (k in 1:4) {
    real h_sd_uncond = omega[k] / sqrt(1 - phi[k] * phi[k]);
    h[1, k] = mu_h[k] + h_sd_uncond * h_raw[1, k];
    for (t in 2:T)
      h[t, k] = mu_h[k] + phi[k] * (h[t-1, k] - mu_h[k]) + omega[k] * h_raw[t, k];
  }
}
model {
  c       ~ normal(1.5, 1);
  rho     ~ normal(0.5, 0.3);
  alpha   ~ normal(-1, 1);
  beta    ~ normal(1, 0.5);
  mu_g    ~ normal(2, 1);
  phi_g   ~ uniform(-1, 1);
  omega_g ~ normal(0, 0.5);
  mu_h    ~ normal(2, 1);
  phi     ~ uniform(-1, 1);
  omega   ~ normal(0, 0.5);
  g_raw   ~ std_normal();
  to_vector(h_raw) ~ std_normal();
  target += sum(kf_tv(y_all, h, [1, 1, 1, beta]', [0, 0, 0, alpha]', c, rho, exp(g))[, 1]);
}
generated quantities {
  real mu_bar = c / (1 - rho);
  real tau = exp(0.5 * mu_g);
  vector[T] tau_t = exp(0.5 * g);
  vector[T] mu;
  matrix[T, 4] sigma_t = exp(0.5 * h);
  matrix[T, 4] y_rep;
  vector[T] log_lik;
  vector[T] log_lik_pred;
  {
    vector[T] qv = exp(g);
    matrix[T, 3] f = kf_tv(y_all, h, [1, 1, 1, beta]', [0, 0, 0, alpha]', c, rho, qv);
    log_lik_pred = f[, 1];
    mu = ffbs_tv_rng(f, c, rho, qv);
  }
  for (t in 1:T) {
    log_lik[t] = 0;
    for (k in 1:3) {
      y_rep[t, k] = normal_rng(mu[t], sigma_t[t, k]);
      log_lik[t] += normal_lpdf(y_gdp[t, k] | mu[t], sigma_t[t, k]);
    }
    y_rep[t, 4] = normal_rng(alpha + beta * mu[t], sigma_t[t, 4]);
    log_lik[t] += normal_lpdf(y_U[t] | alpha + beta * mu[t], sigma_t[t, 4]);
  }
}
