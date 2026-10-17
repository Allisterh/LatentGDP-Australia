// m9: m6 (structural Sigma) + m4 (stochastic volatility) hybrid.
//   mu_t        = c + rho * mu_{t-1} + tau * eta_t,     mu_1 ~ stationary AR(1)
//   y_t^E, y_t^I | mu_t ~ MVN([mu_t, mu_t], Sigma_EI_t)
//   y_t^P       | mu_t ~ N(mu_t, sigma_{t,P}^2)
//   Sigma_EI_t = [s_E^2, rho_EI s_E s_I; rho_EI s_E s_I, s_I^2],  s_k = exp(h_{t,k}/2)
//   h_{t,k} = mu_h_k + phi_k (h_{t-1,k} - mu_h_k) + omega_k xi_{t,k}  (stationary AR(1))
// Correlation rho_EI between E and I errors is constant; P is independent
// (justified by ABS construction: P comes from independent industry surveys),
// i.e. Cov(eps_P, eps_E) = Cov(eps_P, eps_I) = 0.
//
// Sampler parameterisation (same model, same priors, different coordinates):
// conditional on the log-variance paths h and the hyperparameters, the model is
// linear-Gaussian in mu (scalar AR(1) state, 3-dim observation with covariance
// R_t = blockdiag(Sigma_EI_t, s_P^2)). mu is therefore integrated out exactly
// with a Kalman filter in the model block (Rao-Blackwellisation), and HMC only
// explores (theta, h). In generated quantities mu is drawn exactly from
// p(mu | y, h, theta) by forward-filtering backward-sampling, so every output
// (mu, sigma_t, y_rep, log_lik, mu_bar) has the same meaning as before.
// The previous version sampled mu_raw jointly with h_raw; the mu-h coupling
// produced a funnel whose E-BFMI and divergences worsened as T grew
// (T=267: 26/8000 divergences, chain E-BFMI 0.11, lp__ R-hat 1.07).
// h stays non-centred. The filter and FFBS were checked against the dense
// multivariate-normal likelihood and smoother (exact to machine precision /
// within Monte Carlo error).
//
// rho_EI is sampled as z_EI = atanh(rho_EI) (Stan's own (-1,1) transform up to a
// factor 2) with the Jacobian added explicitly, so 1 +- rho_EI can be computed
// from z_EI without rounding. With Stan's built-in bound, 3 of 100 chains started
// from Stan's default random inits (2 of 100 with init = 0.5) stalled in the
// first warmup iterations: the step size collapsed to ~1e-17 and the chain ran
// at max treedepth without moving. 4 of those 5 sat at rho_EI = -1, where
// 1 - rho_EI^2 rounds to 0 and the gradients become garbage. With this form:
// 0 of 200 chains stalled. Same model, same prior on rho_EI.

functions {
  // Scalar-state Kalman filter for
  //   mu_t = c + rho * mu_{t-1} + tau * eta_t,  mu_1 ~ N(c/(1-rho), tau^2/(1-rho^2))
  //   y_t  = 1_3 * mu_t + e_t,  e_t ~ N(0, R_t),
  //   R_t  = [ sE^2, r sE sI, 0 ; r sE sI, sI^2, 0 ; 0, 0, sP^2 ],  s_k = exp(h_{t,k}/2),
  //   r = rho_EI = tanh(z_EI).
  // Numerically stable form. Whiten each measure by s_k and rotate E/I onto the
  // eigenvectors (1,1)/sqrt2 and (1,-1)/sqrt2 of the correlation matrix, whose
  // eigenvalues e_+ = 1 + r = 2 inv_logit(2 z_EI) and e_- = 1 - r = 2 inv_logit(-2 z_EI)
  // are computed from z_EI directly, so they never round to 0 even when r rounds
  // to +-1. With whitened innovation w_j, loading u_j and variance e_j
  // (j = +, -, P; e_P = 1):
  //   s = sum u_j^2/e_j,  q = sum u_j w_j/e_j,  Q = sum w_j^2/e_j,  lambda = 1/P_{t|t-1},
  //   Q - q^2/(lambda + s) = [lambda Q + sum_{j<k} (u_j w_k - u_k w_j)^2/(e_j e_k)] / (lambda + s)
  // (Lagrange identity: every term is non-negative, so there is no cancellation).
  // Returns T x 3 matrix: [log p(y_t | y_{1:t-1}, h, theta), m_{t|t}, P_{t|t}].
  matrix kf_m9(matrix y, matrix h, real c, real rho, real tau, real z_EI) {
    int T = rows(y);
    matrix[T, 3] out;
    real ep = 2 * inv_logit(2 * z_EI);                 // 1 + r
    real em = 2 * inv_logit(-2 * z_EI);                // 1 - r
    real log_det_C = 2 * log2() + log_inv_logit(2 * z_EI) + log_inv_logit(-2 * z_EI);
    real a = c / (1 - rho);                            // predicted mean m_{t|t-1}
    real P = square(tau) / (1 - square(rho));          // predicted var P_{t|t-1}
    real tau2 = square(tau);
    real rho2 = square(rho);
    real isq2 = inv_sqrt(2);
    for (t in 1:T) {
      real uE = exp(-0.5 * h[t, 1]);
      real uI = exp(-0.5 * h[t, 2]);
      real uP = exp(-0.5 * h[t, 3]);
      real wE = (y[t, 1] - a) * uE;
      real wI = (y[t, 2] - a) * uI;
      real wP = (y[t, 3] - a) * uP;
      real up = (uE + uI) * isq2;
      real um = (uE - uI) * isq2;
      real wp = (wE + wI) * isq2;
      real wm = (wE - wI) * isq2;
      real lam = 1 / P;
      real s = square(up) / ep + square(um) / em + square(uP);
      real q = up * wp / ep + um * wm / em + uP * wP;
      real Q = square(wp) / ep + square(wm) / em + square(wP);
      real cross = square(up * wm - um * wp) / (ep * em)
                   + square(up * wP - uP * wp) / ep
                   + square(um * wP - uP * wm) / em;
      real ls = lam + s;
      real quad = (lam * Q + cross) / ls;              // = Q - P q^2 / (1 + P s)
      real Pf = 1 / ls;
      out[t, 1] = -0.5 * (3 * log(2 * pi()) + h[t, 1] + h[t, 2] + h[t, 3] + log_det_C
                          + log(P) + log(ls) + quad);
      out[t, 2] = a + Pf * q;
      out[t, 3] = Pf;
      a = c + rho * out[t, 2];
      P = rho2 * Pf + tau2;
    }
    return out;
  }

  // Backward sampling of mu_{1:T} | y, h, theta from the filter output f = kf_m9(...):
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
  int<lower=3, upper=3> K;           // E, I, P: the structural Sigma is 3 x 3
  matrix[T, K] y;
}
transformed data {
  array[T] vector[2] y_EI;
  vector[T] y_P;
  for (t in 1:T) {
    y_EI[t] = y[t, 1:2]';
    y_P[t]  = y[t, 3];
  }
}
parameters {
  real c;
  real<lower=-1, upper=1> rho;
  real<lower=0> tau;
  vector[K] mu_h;
  vector<lower=-1, upper=1>[K] phi;
  vector<lower=0>[K] omega;
  real z_EI;                         // atanh(rho_EI); see header
  matrix[T, K] h_raw;
}
transformed parameters {
  real<lower=-1, upper=1> rho_EI = tanh(z_EI);
  matrix[T, K] h;
  for (k in 1:K) {
    real h_sd_uncond = omega[k] / sqrt(1 - phi[k] * phi[k]);
    h[1, k] = mu_h[k] + h_sd_uncond * h_raw[1, k];
    for (t in 2:T)
      h[t, k] = mu_h[k] + phi[k] * (h[t-1, k] - mu_h[k]) + omega[k] * h_raw[t, k];
  }
}
model {
  c       ~ normal(1.5, 1);
  rho     ~ normal(0.5, 0.3);
  tau     ~ normal(0, 2);
  mu_h    ~ normal(2, 1);
  phi     ~ uniform(-1, 1);
  omega   ~ normal(0, 0.5);
  // rho_EI ~ normal(0, 0.5) on (-1, 1), written on z_EI = atanh(rho_EI) with the
  // Jacobian log(1 - rho_EI^2) = log(1 + r) + log(1 - r), computed stably.
  target += normal_lpdf(rho_EI | 0, 0.5)
            + 2 * log2() + log_inv_logit(2 * z_EI) + log_inv_logit(-2 * z_EI);
  to_vector(h_raw) ~ std_normal();

  // log p(y | h, theta) with mu (AR(1) prior as before) integrated out.
  target += sum(kf_m9(y, h, c, rho, tau, z_EI)[, 1]);
}
generated quantities {
  real mu_bar = c / (1 - rho);
  vector[T] mu;                      // exact draw from p(mu | y, h, theta) (FFBS)
  matrix[T, K] sigma_t = exp(0.5 * h);
  matrix[T, K] y_rep;
  vector[T] log_lik;                 // log p(y_t | mu_t, h_t, theta) (pointwise conditional)
  vector[T] log_lik_pred;            // log p(y_t | y_{1:t-1}, h, theta), mu integrated out
  {
    matrix[T, 3] f = kf_m9(y, h, c, rho, tau, z_EI);
    log_lik_pred = f[, 1];
    mu = ffbs_rng(f, c, rho, tau);
  }
  real sqrt_1m_r2 = sqrt(4 * inv_logit(2 * z_EI) * inv_logit(-2 * z_EI));  // sqrt(1 - rho_EI^2)
  for (t in 1:T) {
    real s_E = sigma_t[t, 1];
    real s_I = sigma_t[t, 2];
    real s_P = sigma_t[t, 3];
    matrix[2, 2] L_EI;
    L_EI[1, 1] = s_E;
    L_EI[2, 1] = s_I * rho_EI;
    L_EI[1, 2] = 0;
    L_EI[2, 2] = s_I * sqrt_1m_r2;
    vector[2] mu_vec = rep_vector(mu[t], 2);
    vector[2] yrep_EI = multi_normal_cholesky_rng(mu_vec, L_EI);
    y_rep[t, 1] = yrep_EI[1];
    y_rep[t, 2] = yrep_EI[2];
    y_rep[t, 3] = normal_rng(mu[t], s_P);
    log_lik[t]  = multi_normal_cholesky_lpdf(y_EI[t] | mu_vec, L_EI)
                + normal_lpdf(y_P[t] | mu[t], s_P);
  }
}
