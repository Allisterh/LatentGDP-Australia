// m4: AR(1) latent growth with stationary-SV measurement errors.
//   mu_t  = c + rho * mu_{t-1} + tau * eta_t,  mu_1 ~ stationary AR(1)
//   y_t^k = mu_t + eps_t^k,           eps_t^k ~ N(0, exp(h_{t,k}/2)^2)
//   h_{t,k} = mu_h_k + phi_k * (h_{t-1,k} - mu_h_k) + omega_k * xi_{t,k}
//   h_{1,k} ~ N(mu_h_k, omega_k^2 / (1 - phi_k^2))    (stationary init)
// Stationary AR(1) on log-variance (Kim-Shephard-Chib). Drops the unidentified
// random-walk parameterisation tried earlier - see NOTES.md.
//
// Sampler parameterisation (same model, same priors, different coordinates;
// identical to m9's so that model comparison is like-for-like):
// conditional on the log-variance paths h, the model is linear-Gaussian in mu,
// so mu is integrated out exactly with a Kalman filter in the model block
// (Rao-Blackwellisation) and drawn exactly from p(mu | y, h, theta) by
// forward-filtering backward-sampling in generated quantities. HMC only
// explores (theta, h). h stays non-centred. Outputs (mu, sigma_t, y_rep,
// log_lik, mu_bar) keep their previous meaning; log_lik_pred is new.

functions {
  // Scalar-state Kalman filter for
  //   mu_t = c + rho * mu_{t-1} + tau * eta_t,  mu_1 ~ N(c/(1-rho), tau^2/(1-rho^2))
  //   y_{t,k} = d_k + Z_k * mu_t + e_{t,k},  e_{t,k} ~ N(0, exp(h_{t,k}))  independent over k.
  // Returns T x 3 matrix: [log p(y_t | y_{1:t-1}, h, theta), m_{t|t}, P_{t|t}].
  matrix kf_diag(matrix y, matrix h, vector Z, vector d, real c, real rho, real tau) {
    int T = rows(y);
    int K = cols(y);
    matrix[T, 3] out;
    matrix[T, K] ih = exp(-h);                          // 1 / sigma_{t,k}^2
    matrix[T, K] yd = y - rep_matrix(d', T);
    vector[T] s = ih * square(Z);                       // Z' R_t^{-1} Z
    vector[T] base = K * log(2 * pi()) + h * rep_vector(1, K);
    real a = c / (1 - rho);
    real P = square(tau) / (1 - square(rho));
    real tau2 = square(tau);
    real rho2 = square(rho);
    for (t in 1:T) {
      row_vector[K] e = yd[t] - a * Z';                 // innovation
      row_vector[K] w = e .* ih[t];
      real q = w * Z;                                   // Z' R^{-1} e
      real Q = w * e';                                  // e' R^{-1} e
      real g = 1 + P * s[t];
      real Pf = P / g;
      real m = a + Pf * q;
      out[t, 1] = -0.5 * (base[t] + log(g) + Q - P * square(q) / g);
      out[t, 2] = m;
      out[t, 3] = Pf;
      a = c + rho * m;
      P = rho2 * Pf + tau2;
    }
    return out;
  }

  // Backward sampling of mu_{1:T} | y, h, theta from the filter output f = kf_diag(...):
  //   mu_T ~ N(m_T, P_{T|T});  mu_t | mu_{t+1} ~ N(m_t + J (mu_{t+1} - c - rho m_t),
  //   P_{t|t} tau^2 / P_{t+1|t}),  J = rho P_{t|t} / P_{t+1|t}.
  vector ffbs_rng(matrix f, real c, real rho, real tau) {
    int T = rows(f);
    vector[T] mu;
    mu[T] = normal_rng(f[T, 2], sqrt(f[T, 3]));
    for (tt in 1:(T - 1)) {
      int t = T - tt;
      real Pp = square(rho) * f[t, 3] + square(tau);
      real J = rho * f[t, 3] / Pp;
      mu[t] = normal_rng(f[t, 2] + J * (mu[t + 1] - c - rho * f[t, 2]),
                         sqrt(f[t, 3] * square(tau) / Pp));
    }
    return mu;
  }
}
data {
  int<lower=1> T;
  int<lower=1> K;
  matrix[T, K] y;
}
transformed data {
  vector[K] Z = rep_vector(1, K);   // loadings of each measure on mu
  vector[K] d = rep_vector(0, K);   // intercepts
}
parameters {
  real c;
  real<lower=-1, upper=1> rho;
  real<lower=0> tau;
  vector[K] mu_h;                    // long-run log-variance per series
  vector<lower=-1, upper=1>[K] phi;  // persistence of log-variance
  vector<lower=0>[K] omega;          // sd of log-variance innovations
  matrix[T, K] h_raw;
}
transformed parameters {
  matrix[T, K] h;
  for (k in 1:K) {
    real h_sd_uncond = omega[k] / sqrt(1 - phi[k] * phi[k]);
    h[1, k] = mu_h[k] + h_sd_uncond * h_raw[1, k];
    for (t in 2:T)
      h[t, k] = mu_h[k] + phi[k] * (h[t-1, k] - mu_h[k]) + omega[k] * h_raw[t, k];
  }
}
model {
  c     ~ normal(1.5, 1);
  rho   ~ normal(0.5, 0.3);
  tau   ~ normal(0, 2);
  mu_h  ~ normal(2, 1);              // log-var ~ 2 -> sd ~ exp(1) ~ 2.7
  phi   ~ uniform(-1, 1);            // weakly informative; prior on a half-life via posterior
  omega ~ normal(0, 0.5);
  to_vector(h_raw) ~ std_normal();

  // log p(y | h, theta) with mu (AR(1) prior as before) integrated out.
  target += sum(kf_diag(y, h, Z, d, c, rho, tau)[, 1]);
}
generated quantities {
  real mu_bar = c / (1 - rho);
  vector[T] mu;                      // exact draw from p(mu | y, h, theta) (FFBS)
  matrix[T, K] sigma_t = exp(0.5 * h);
  matrix[T, K] y_rep;
  vector[T] log_lik;                 // log p(y_t | mu_t, h_t, theta) (pointwise conditional)
  vector[T] log_lik_pred;            // log p(y_t | y_{1:t-1}, h, theta), mu integrated out
  {
    matrix[T, 3] f = kf_diag(y, h, Z, d, c, rho, tau);
    log_lik_pred = f[, 1];
    mu = ffbs_rng(f, c, rho, tau);
  }
  for (t in 1:T) {
    log_lik[t] = 0;
    for (k in 1:K) {
      y_rep[t, k] = normal_rng(mu[t], sigma_t[t, k]);
      log_lik[t] += normal_lpdf(y[t, k] | mu[t], sigma_t[t, k]);
    }
  }
}
