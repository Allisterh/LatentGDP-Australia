// m3: AR(1) latent growth with correlated measurement errors.
//   y_t = mu_t * 1_K + eps_t,        eps_t ~ MVN(0, Sigma)
//   Sigma = diag(sigma) * Omega * diag(sigma),  Omega ~ LKJ(4)
//   (LKJ(4) on the Cholesky factor; LKJ(2) was tried first and sampled badly,
//    see NOTES.md "Identification issue in m3".)
//   mu_t = c + rho * mu_{t-1} + eta_t,  eta_t ~ N(0, tau^2),  mu_1 ~ stationary AR(1)
// Common shocks (e.g. balancing-item adjustments across approaches) get
// absorbed by Omega instead of leaking into mu_t.
//
// Sampler parameterisation (same model, same priors, different coordinates):
// the model is linear-Gaussian in mu with a constant measurement covariance
// Sigma, so mu is integrated out exactly by a scalar Kalman filter in the model
// block and HMC explores only (c, rho, tau, sigma, L_Omega). mu is drawn exactly
// from p(mu | y, theta) by forward-filtering backward-sampling in generated
// quantities, so mu, Omega, log_lik, y_rep and mu_bar keep their meaning.
// The previous version sampled mu_raw jointly with Omega: E-BFMI 0.03, R-hat 1.029
// (Omega[3,2]), min bulk ESS 132 on the 2026Q2 vintage, and every LFO refit had
// E-BFMI < 0.07. The filter is checked against the dense multivariate-normal
// likelihood in tests/test_m3_m10.R.
//
// Identification (not changed by the reparameterisation): with rho near 0 the
// latent state is close to iid, so the data identify the cross-sectional
// covariance v 1 1' + Sigma (v = tau^2 / (1 - rho^2)) and the lag-1
// autocovariance rho v 1 1' only. Adding d to every element of Sigma and
// subtracting d from v leaves the first unchanged, so (tau, Sigma) are
// identified along that ridge only by rho, the PD constraint and the priors.

functions {
  // Scalar-state Kalman filter for
  //   mu_t = c + rho * mu_{t-1} + tau * eta_t,  mu_1 ~ N(c/(1-rho), tau^2/(1-rho^2))
  //   y_t  = 1_K * mu_t + e_t,  e_t ~ N(0, L L'),  L = diag(sigma) * L_Omega.
  // Whitening by L: u = L^{-1} 1, w_t = L^{-1}(y_t - a_t 1). With lambda = 1/P_{t|t-1},
  //   s = u'u, q = u'w, Q = w'w,
  //   log p(y_t | y_{1:t-1}) = -0.5 [K log 2pi + log|R| + log(P) + log(lambda + s) + quad],
  //   quad = Q - q^2 / (lambda + s) = (lambda Q + (s Q - q^2)) / (lambda + s),
  // and s Q - q^2 = |u|^2 |w|^2 - (u'w)^2 >= 0 is computed as a sum of squares
  // (Lagrange identity), so the quadratic form has no cancellation.
  // Returns T x 3 matrix: [log p(y_t | y_{1:t-1}, theta), m_{t|t}, P_{t|t}].
  matrix kf_m3(matrix y, matrix L, real c, real rho, real tau) {
    int T = rows(y);
    int K = cols(y);
    matrix[T, 3] out;
    vector[K] u = mdivide_left_tri_low(L, rep_vector(1, K));
    matrix[K, T] Wy = mdivide_left_tri_low(L, y');     // L^{-1} y_t, column t
    real s = dot_self(u);
    real log_det_R = 2 * sum(log(diagonal(L)));
    real a = c / (1 - rho);
    real P = square(tau) / (1 - square(rho));
    real tau2 = square(tau);
    real rho2 = square(rho);
    for (t in 1:T) {
      vector[K] w = Wy[, t] - a * u;
      real q = dot_product(u, w);
      real Q = dot_self(w);
      real cross = 0;
      for (j in 1:(K - 1))
        for (k in (j + 1):K)
          cross += square(u[j] * w[k] - u[k] * w[j]);
      real lam = 1 / P;
      real ls = lam + s;
      real quad = (lam * Q + cross) / ls;
      real Pf = 1 / ls;
      out[t, 1] = -0.5 * (K * log(2 * pi()) + log_det_R + log(P) + log(ls) + quad);
      out[t, 2] = a + Pf * q;
      out[t, 3] = Pf;
      a = c + rho * out[t, 2];
      P = rho2 * Pf + tau2;
    }
    return out;
  }

  // Backward sampling of mu_{1:T} | y, theta from the filter output f = kf_m3(...):
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
  int<lower=2> K;
  matrix[T, K] y;
}
parameters {
  real c;
  real<lower=-1, upper=1> rho;
  real<lower=0> tau;
  vector<lower=0>[K] sigma;
  cholesky_factor_corr[K] L_Omega;
}
transformed parameters {
  matrix[K, K] L_Sigma = diag_pre_multiply(sigma, L_Omega);
}
model {
  c       ~ normal(1.5, 1);
  rho     ~ normal(0.5, 0.3);
  tau     ~ normal(0, 2);
  sigma   ~ normal(0, 2);
  L_Omega ~ lkj_corr_cholesky(4);   // tighter than LKJ(2): identification through
                                     // regularisation. mu_t and corr(eps) are
                                     // weakly identified - see NOTES.md.

  // log p(y | theta) with mu (stationary AR(1) prior as before) integrated out.
  target += sum(kf_m3(y, L_Sigma, c, rho, tau)[, 1]);
}
generated quantities {
  real mu_bar = c / (1 - rho);
  matrix[K, K] Omega = multiply_lower_tri_self_transpose(L_Omega);
  vector[T] mu;                      // exact draw from p(mu | y, theta) (FFBS)
  vector[T] log_lik;                 // log p(y_t | mu_t, theta) (pointwise conditional)
  vector[T] log_lik_pred;            // log p(y_t | y_{1:t-1}, theta), mu integrated out
  matrix[T, K] y_rep;
  {
    matrix[T, 3] f = kf_m3(y, L_Sigma, c, rho, tau);
    log_lik_pred = f[, 1];
    mu = ffbs_rng(f, c, rho, tau);
  }
  for (t in 1:T) {
    vector[K] mu_vec = rep_vector(mu[t], K);
    log_lik[t] = multi_normal_cholesky_lpdf(y[t]' | mu_vec, L_Sigma);
    y_rep[t] = multi_normal_cholesky_rng(mu_vec, L_Sigma)';
  }
}
