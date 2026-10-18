// m12: m9 (structural Sigma + stochastic volatility) + serially correlated
// measurement errors (as in m5).
//   mu_t    = c + rho * mu_{t-1} + tau * eta_t,     mu_1 ~ stationary AR(1)
//   y_t^k   = mu_t + eps_t^k,                       k = E, I, P
//   eps_t^k = psi_k * eps_{t-1}^k + u_t^k
//   (u_t^E, u_t^I) ~ MVN(0, Sigma_EI_t),  Sigma_EI_t = [s_E^2, r s_E s_I; r s_E s_I, s_I^2]
//   u_t^P ~ N(0, s_P^2), independent of (u^E, u^I);   s_k = s_{t,k} = exp(h_{t,k}/2)
//   h_{t,k} = mu_h_k + phi_k (h_{t-1,k} - mu_h_k) + omega_k xi_{t,k}   (stationary AR(1),
//             exactly as m9)
// so R_t = Cov(u_t) = blockdiag(Sigma_EI_t, s_P^2) is m9's measurement covariance, and
// sigma_t (generated quantity) is the sd of the innovation u_t, as sigma is in m5.
// psi = 0 gives m9 exactly.
//
// Initial condition for eps (t = 1). eps_1 ~ N(0, V_1) with
//   V_1[j,k] = R_1[j,k] / (1 - psi_j psi_k),
// i.e. the stationary covariance of the VAR(1) eps_t = diag(psi) eps_{t-1} + u_t if the
// innovation covariance were held at R_1 (the solution of V = Psi V Psi + R_1). Diagonal:
// exp(h_{1,k}) / (1 - psi_k^2), m5's stationary sd sigma_k / sqrt(1 - psi_k^2) with sigma_k
// replaced by s_{1,k}; off-diagonal: r s_E s_I / (1 - psi_E psi_I) (m5 has no correlation,
// so this is the only extension needed). Its correlation is
//   r_1 = r sqrt((1 - psi_E^2)(1 - psi_I^2)) / (1 - psi_E psi_I),   |r_1| <= |r|.
// Why this rather than conditioning on y_1: it is m5's convention and m9's (both model
// y_1 with a stationary start), it reduces to m9 at psi = 0, and it keeps log p(y_1) in the
// likelihood, so the total log-likelihood and every per-period predictive density are on
// the same footing as m4/m5/m9. Conditioning on y_1 would drop y_1 from the likelihood and
// leave mu_1's prior to be re-derived. The LFO particle filter (R/lib/filter_m12.cpp) uses
// exactly this V_1. Period 1 (1959 Q4) is outside the LFO window either way.
//
// Sampler parameterisation (as m9): conditional on h and the hyperparameters the model is
// linear-Gaussian in mu, and quasi-differencing removes eps exactly:
//   t >= 2:  ytil_t = y_t - psi .* y_{t-1} = mu_t * 1 - psi * mu_{t-1} + u_t,  u_t ~ N(0, R_t)
//   t  = 1:  y_1 = mu_1 * 1 + eps_1,  eps_1 ~ N(0, V_1)
// (unit Jacobian given y_{t-1}). The observation at t loads on (mu_t, mu_{t-1}), so the
// filter works on the 2-dim state (mu_t, mu_{t-1}); because mu_{t+1} depends on mu_t only,
// it carries just the scalar filtered moments (m_t, P_t) of mu_t between periods. mu is
// integrated out in the model block; HMC explores (theta, h). Generated quantities draw mu
// exactly from p(mu | y, h, theta) by forward-filtering backward-sampling: given mu_{t+1},
// mu_t depends on the data only through y_{1:t+1}, so the backward step uses the joint
// filtered law of (mu_{t+1}, mu_t) | y_{1:t+1}. The filter and FFBS are checked against the
// dense multivariate-normal likelihood and smoother in tests/test_m12.R.
//
// rho_EI is sampled as z_EI = atanh(rho_EI) with the Jacobian added explicitly, as in m9,
// and 1 +- rho_EI (and 1 +- r_1) are computed from z_EI without cancellation.

functions {
  // Whiten a 3-vector x (in measure units) for covariance D Corr(r) D with
  // 1 / D_kk = sc[k] and Corr(r) = [1, r, 0; r, 1, 0; 0, 0, 1]: scale by sc, rotate E/I onto
  // the eigenvectors (1, 1)/sqrt2 and (1, -1)/sqrt2 of Corr(r), whose eigenvalues are
  // ep = 1 + r and em = 1 - r, and divide by their square roots.
  vector wh3(vector x, vector sc, real ep, real em) {
    vector[3] o;
    real xE = x[1] * sc[1];
    real xI = x[2] * sc[2];
    o[1] = (xE + xI) / sqrt(2 * ep);
    o[2] = (xE - xI) / sqrt(2 * em);
    o[3] = x[3] * sc[3];
    return o;
  }
  vector cross3(vector a, vector b) {
    vector[3] o;
    o[1] = a[2] * b[3] - a[3] * b[2];
    o[2] = a[3] * b[1] - a[1] * b[3];
    o[3] = a[1] * b[2] - a[2] * b[1];
    return o;
  }

  // Kalman filter for m12 given h (see the header). Period t:
  //   predicted state  (mu_t, mu_{t-1}) = a + G z,  z ~ N(0, I_2),
  //     t >= 2: a = (c + rho m, m), G = [rho sqrt(P), tau; sqrt(P), 0],  (m, P) = filtered
  //             moments of mu_{t-1}
  //     t  = 1: a = (c/(1-rho), -), G row 1 = (0, tau/sqrt(1-rho^2))
  //   observation residual v and loadings on z, whitened by the measurement covariance H_t
  //   (R_t for t >= 2, V_1 for t = 1):  w = W v,  b_j = W Zt G[, j].
  // With B = [b1 b2] and M = I + B'B, the prediction covariance is F = W^-1 (I + B B') W^-T,
  //   log|F| = log|H_t| + log det M,
  //   det M  = 1 + |b1|^2 + |b2|^2 + |b1 x b2|^2,
  //   v'F^-1 v = w'(I + BB')^-1 w
  //            = (|w|^2 + |b1 x w|^2 + |b2 x w|^2 + det[b1 b2 w]^2) / det M
  // (from 1 + w'(I+BB')^-1 w = det(I + [B w]'[B w]) / det(I + B'B) and Cauchy-Binet: every
  // term is a sum of squares, so there is no cancellation when a measurement variance is
  // tiny; this is m9's Lagrange-identity form generalised from a scalar to a 2-dim state).
  // Posterior z | y_{1:t} ~ N(M^-1 B'w, M^-1); the components of M^-1 B'w are written with
  // the Binet-Cauchy identity (a x b).(c x d) = (a.c)(b.d) - (a.d)(b.c).
  // Returns T x 6: [log p(y_t | y_{1:t-1}, h, theta), m_t = E(mu_t | y_{1:t}),
  //   P_t = Var(mu_t | y_{1:t}), and, for t >= 2, the backward-sampling coefficients of
  //   mu_{t-1} | mu_t, y_{1:t} ~ N(col4 + col5 * (mu_t - m_t), col6)] (cols 4-6 are 0 at t = 1).
  matrix kf_m12(matrix y, matrix h, real c, real rho, real tau, real z_EI, vector psi) {
    int T = rows(y);
    matrix[T, 6] out;
    real ep = 2 * inv_logit(2 * z_EI);                 // 1 + r
    real em = 2 * inv_logit(-2 * z_EI);                // 1 - r
    real log_det_C = 2 * log2() + log_inv_logit(2 * z_EI) + log_inv_logit(-2 * z_EI);
    // t = 1: V_1 = D_1 Corr(r_1) D_1, D_1 = diag(s_{1,k} / sqrt(1 - psi_k^2)).
    // 1 +- r_1 = (dk + (1 +- r) k) / d with d = 1 - psi_E psi_I,
    // k = sqrt((1 - psi_E^2)(1 - psi_I^2)) and dk = d - k = (psi_E - psi_I)^2 / (d + k) >= 0.
    real d1 = 1 - psi[1] * psi[2];
    real k1 = sqrt((1 - square(psi[1])) * (1 - square(psi[2])));
    real dk1 = square(psi[1] - psi[2]) / (d1 + k1);
    real ep1 = (dk1 + ep * k1) / d1;
    real em1 = (dk1 + em * k1) / d1;
    vector[3] l1m = log1m(square(psi));                // log(1 - psi_k^2)
    vector[3] s1 = exp(0.5 * l1m);                     // sqrt(1 - psi_k^2)
    vector[3] ones = rep_vector(1, 3);
    for (t in 1:T) {
      vector[3] sc = exp(-0.5 * h[t]');
      vector[3] w;
      vector[3] b1;
      vector[3] b2;
      real a1;
      real g11;
      real g12;
      real sP = 0;
      real lH;
      if (t == 1) {
        vector[3] sc1 = sc .* s1;
        a1 = c / (1 - rho);
        g11 = 0;
        g12 = tau / sqrt(1 - square(rho));
        w = wh3(y[1]' - a1, sc1, ep1, em1);
        b1 = rep_vector(0, 3);
        b2 = wh3(rep_vector(g12, 3), sc1, ep1, em1);
        lH = sum(h[1]) - sum(l1m) + log(ep1) + log(em1);
      } else {
        real m = out[t - 1, 2];
        sP = sqrt(out[t - 1, 3]);
        a1 = c + rho * m;
        g11 = rho * sP;
        g12 = tau;
        w = wh3((y[t]' - a1) - psi .* (y[t - 1]' - m), sc, ep, em);
        b1 = wh3(sP * (rho - psi), sc, ep, em);
        b2 = wh3(rep_vector(tau, 3), sc, ep, em);
        lH = sum(h[t]) + log_det_C;
      }
      {
        vector[3] x12 = cross3(b1, b2);
        vector[3] x1w = cross3(b1, w);
        vector[3] x2w = cross3(b2, w);
        real d3 = dot_product(x12, w);
        real den = 1 + dot_self(b1) + dot_self(b2) + dot_self(x12);
        real num = dot_self(w) + dot_self(x1w) + dot_self(x2w) + square(d3);
        real zh1 = (dot_product(b1, w) - dot_product(x12, x2w)) / den;
        real zh2 = (dot_product(b2, w) + dot_product(x12, x1w)) / den;
        real g1aG1 = square(g11) + square(g12) + dot_self(g11 * b2 - g12 * b1);  // G1 adj(M) G1'
        out[t, 1] = -0.5 * (3 * log(2 * pi()) + lH + log(den) + num / den);
        out[t, 2] = a1 + g11 * zh1 + g12 * zh2;
        out[t, 3] = g1aG1 / den;
        if (t == 1) {
          out[t, 4] = 0;
          out[t, 5] = 0;
          out[t, 6] = 0;
        } else {
          // G2 = (sP, 0): E(mu_{t-1} | y_{1:t}) = m + sP zh1,
          // Cov(mu_t, mu_{t-1} | y_{1:t}) = sP (g11 (1 + |b2|^2) - g12 b1.b2) / den,
          // Var(mu_{t-1} | mu_t, y_{1:t}) = det(G)^2 / (den Var(mu_t | y_{1:t})) = tau^2 P / g1aG1.
          real cov12 = sP * (g11 + dot_product(b2, g11 * b2 - g12 * b1)) / den;
          out[t, 4] = out[t - 1, 2] + sP * zh1;
          out[t, 5] = cov12 / out[t, 3];
          out[t, 6] = square(tau * sP) / g1aG1;
        }
      }
    }
    return out;
  }

  // Backward sampling of mu_{1:T} | y, h, theta from f = kf_m12(...):
  //   mu_T ~ N(m_T, P_T);  mu_t | mu_{t+1} ~ N(f[t+1,4] + f[t+1,5] (mu_{t+1} - m_{t+1}), f[t+1,6]).
  vector ffbs_m12_rng(matrix f) {
    int T = rows(f);
    vector[T] mu;
    mu[T] = normal_rng(f[T, 2], sqrt(f[T, 3]));
    for (tt in 1:(T - 1)) {
      int t = T - tt;
      mu[t] = normal_rng(f[t + 1, 4] + f[t + 1, 5] * (mu[t + 1] - f[t + 1, 2]),
                         sqrt(f[t + 1, 6]));
    }
    return mu;
  }

  // log N(x | 0, D Corr(r) D) with 1/D_kk = sc[k], log|D Corr D| = lH.
  real lpdf_wh3(vector x, vector sc, real ep, real em, real lH) {
    return -0.5 * (3 * log(2 * pi()) + lH + dot_self(wh3(x, sc, ep, em)));
  }
}
data {
  int<lower=1> T;
  int<lower=3, upper=3> K;           // E, I, P: the structural Sigma is 3 x 3
  matrix[T, K] y;
}
parameters {
  real c;
  real<lower=-1, upper=1> rho;
  real<lower=0> tau;
  vector[K] mu_h;
  vector<lower=-1, upper=1>[K] phi;
  vector<lower=0>[K] omega;
  real z_EI;                         // atanh(rho_EI); see m9
  vector<lower=-1, upper=1>[K] psi;  // measurement-error persistence (as m5)
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
  psi     ~ normal(0, 0.5);
  // rho_EI ~ normal(0, 0.5) on (-1, 1), written on z_EI with the Jacobian (as m9).
  target += normal_lpdf(rho_EI | 0, 0.5)
            + 2 * log2() + log_inv_logit(2 * z_EI) + log_inv_logit(-2 * z_EI);
  to_vector(h_raw) ~ std_normal();

  // log p(y | h, theta) with mu (stationary AR(1) prior) and eps integrated out.
  target += sum(kf_m12(y, h, c, rho, tau, z_EI, psi)[, 1]);
}
generated quantities {
  real mu_bar = c / (1 - rho);
  vector[T] mu;                      // exact draw from p(mu | y, h, theta) (FFBS)
  matrix[T, K] sigma_t = exp(0.5 * h);   // sd of the innovation u_t (m5's sigma)
  vector[K] eps_sd_bar = exp(0.5 * mu_h) ./ sqrt(1 - square(psi));  // sd of eps at h = mu_h
  matrix[T, K] y_rep;                // y_t | mu_t, mu_{t-1}, y_{t-1}, h_t (m5's convention)
  vector[T] log_lik;                 // log p(y_t | mu_t, mu_{t-1}, y_{t-1}, h_t, theta)
  vector[T] log_lik_pred;            // log p(y_t | y_{1:t-1}, h, theta), mu integrated out
  {
    matrix[T, 6] f = kf_m12(y, h, c, rho, tau, z_EI, psi);
    log_lik_pred = f[, 1];
    mu = ffbs_m12_rng(f);
  }
  {
    real ep = 2 * inv_logit(2 * z_EI);
    real em = 2 * inv_logit(-2 * z_EI);
    real log_det_C = 2 * log2() + log_inv_logit(2 * z_EI) + log_inv_logit(-2 * z_EI);
    real sqrt_1m_r2 = sqrt(ep * em);                   // sqrt(1 - rho_EI^2)
    real d1 = 1 - psi[1] * psi[2];
    real k1 = sqrt((1 - square(psi[1])) * (1 - square(psi[2])));
    real dk1 = square(psi[1] - psi[2]) / (d1 + k1);
    real ep1 = (dk1 + ep * k1) / d1;                   // 1 + r_1
    real em1 = (dk1 + em * k1) / d1;                   // 1 - r_1
    real r1 = rho_EI * k1 / d1;
    vector[3] l1m = log1m(square(psi));
    for (t in 1:T) {
      vector[3] sc = exp(-0.5 * h[t]');
      vector[3] s = sigma_t[t]';
      vector[3] m_t;                                   // E(y_t | mu_t, mu_{t-1}, y_{t-1})
      vector[3] z3;
      real rr;
      real sq;
      if (t == 1) {
        s = s ./ exp(0.5 * l1m);                       // sd of eps_1
        m_t = rep_vector(mu[1], 3);
        rr = r1;
        sq = sqrt(ep1 * em1);
        log_lik[t] = lpdf_wh3(y[1]' - m_t, sc .* exp(0.5 * l1m), ep1, em1,
                              sum(h[1]) - sum(l1m) + log(ep1) + log(em1));
      } else {
        m_t = mu[t] + psi .* (y[t - 1]' - mu[t - 1]);
        rr = rho_EI;
        sq = sqrt_1m_r2;
        log_lik[t] = lpdf_wh3(y[t]' - m_t, sc, ep, em, sum(h[t]) + log_det_C);
      }
      for (k in 1:3) z3[k] = std_normal_rng();
      y_rep[t, 1] = m_t[1] + s[1] * z3[1];
      y_rep[t, 2] = m_t[2] + s[2] * (rr * z3[1] + sq * z3[2]);
      y_rep[t, 3] = m_t[3] + s[3] * z3[3];
    }
  }
}
