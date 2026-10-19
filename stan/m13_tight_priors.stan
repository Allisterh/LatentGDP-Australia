// m13_tight_priors: m13 (stan/m13_serial_latent_sv.stan) with tighter, substantively
// justified priors chosen by prior predictive simulation (NOTES.md, "Tighter priors for m13
// (prior predictive check and refit), 2026-10-19"). RESEARCH VARIANT ONLY: production
// (R/fit_production.R) is unchanged and still fits stan/m13_serial_latent_sv.stan.
//
// The likelihood is m13's, unchanged: the functions block (kf_m13, ffbs_m13_rng, wh3, ...)
// and the data block are byte-identical copies of m13's, the latent and measurement
// log-variance recursions produce the same g and h for the same (mu_g, phi_g, omega_g,
// g_raw) and (mu_h, phi, omega, h_raw), and the filter receives c = mu_bar (1 - rho).
// tests/test_m13_tight.R checks that kf_m13 here returns m13's output for the same inputs.
//
// What changes is the prior, through three things (m13's priors in brackets):
//   1. Latent long-run mean. mu_bar ~ N(3, 2) is placed on the long-run mean of latent growth
//      (the prior m1 uses), and c = mu_bar (1 - rho) is derived [c ~ N(1.5, 1)]. mu_bar is
//      now a sampled parameter and the prior is a density on the sampled (mu_bar, rho), so no
//      Jacobian is needed. (Keeping c sampled and writing mu_bar = c/(1-rho) ~ N(3, 2) would
//      need target += -log1m(rho), the log |d mu_bar / d c|; the two are the same model.)
//      Under m13's prior, rho near 1 made c/(1-rho) explode.
//   2. Latent persistence. rho ~ Beta(2.625, 2.625) on (0, 1): mean 0.5, sd 0.2, the prior
//      RBA RDP 2014-12 (Rees, Lancaster & Finlay) use for their latent persistence
//      [rho ~ N(0.5, 0.3) on (-1, 1)]. It excludes negative persistence and puts 1.3% of its
//      mass above 0.9 (m13's prior: 4.7%).
//   3. Log-variance scale. For each SV block (latent g and measurement h_k) the prior is put
//      on the STATIONARY sd of the log-variance, s = omega / sqrt(1 - phi^2), instead of on
//      the innovation sd omega: s_g, s_h[k] ~ half-N(0, 1), with phi_g, phi[k] ~ U(-1, 1)
//      kept, and omega = s sqrt(1 - phi^2) derived [omega ~ half-N(0, 0.5)]. Under m13's
//      prior, |phi| near 1 with a moderate omega made the stationary sd of h (and g)
//      explode; this was the main source of implausible prior predictive tails. half-N(0, 1)
//      on s equals m13's half-N(0, 0.5) on omega at |phi| = 0.87. s and phi are sampled
//      directly, so again no Jacobian. Initial log-variances g_1 ~ N(mu_g, s_g^2) and
//      h_{1,k} ~ N(mu_h_k, s_h_k^2) are the same stationary laws as in m13.
// Unchanged: mu_g ~ N(2, 1), mu_h ~ N(2, 1), psi ~ N(0, 0.5) on (-1, 1),
// rho_EI ~ N(0, 0.5) on (-1, 1) (sampled as z_EI = atanh(rho_EI) with the Jacobian).
// Kim-Shephard-Chib's (phi + 1)/2 ~ Beta(20, 1.5) was considered and rejected: with m13's
// omega prior it makes the prior predictive tails worse (NOTES.md).
//
// Outputs keep m13's names: c, omega_g and omega are transformed parameters, mu_bar a
// parameter, so the production gate's variable list applies unchanged. The new sampled
// scale parameters are s_g and s_h.
//
// m13's header follows (model equations, initial conditions, sampler parameterisation;
// its "Priors" paragraph is superseded by the list above).
//  | m13: m12 (structural Sigma + SV on the measurement innovations + AR(1) measurement
//  | errors) + stochastic volatility on the LATENT innovation, as in the re-specified m10.
//  |   mu_t    = c + rho * mu_{t-1} + exp(g_t / 2) * eta_t
//  |   g_t     = mu_g + phi_g (g_{t-1} - mu_g) + omega_g zeta_t      (stationary AR(1))
//  |   g_1     ~ N(mu_g, omega_g^2 / (1 - phi_g^2))
//  |   mu_1    ~ N(c / (1 - rho), exp(g_1) / (1 - rho^2))             (m10's convention)
//  |   y_t^k   = mu_t + eps_t^k,                                       k = E, I, P
//  |   eps_t^k = psi_k * eps_{t-1}^k + u_t^k
//  |   (u_t^E, u_t^I) ~ MVN(0, Sigma_EI_t), u_t^P ~ N(0, s_P^2) independent; s_{t,k} = exp(h_{t,k}/2)
//  |   h_{t,k} = mu_h_k + phi_k (h_{t-1,k} - mu_h_k) + omega_k xi_{t,k}  (stationary AR(1), as m12)
//  |   eps_1 ~ N(0, V_1), V_1[j,k] = R_1[j,k] / (1 - psi_j psi_k)      (as m12)
//  | Everything except the latent innovation sd is m12 unchanged. With g_t = log(tau^2) for all
//  | t (omega_g = 0, mu_g = log tau^2) the model is m12 exactly.
//  | 
//  | Initial latent state. mu_1 is drawn from the stationary AR(1) law the latent process would
//  | have if its innovation variance were held at its period-1 value exp(g_1): mean c/(1-rho),
//  | variance exp(g_1)/(1-rho^2). This is m10's convention (stan/m10_labour_sv.stan) and the
//  | analogue of m12's V_1 for the measurement errors; with g constant it is m12's stationary
//  | start. It is not the marginal stationary law of mu under latent SV (which has no closed
//  | form); period 1 (1959 Q4) is outside the LFO window either way.
//  | 
//  | Priors: m12's for c, rho, mu_h, phi, omega, psi and rho_EI; the latent SV block uses m10's
//  | (the same as each measurement SV block): mu_g ~ N(2, 1), phi_g ~ U(-1, 1),
//  | omega_g ~ half-N(0, 0.5). They replace m12's tau ~ half-N(0, 2).
//  | Output `tau` = exp(mu_g / 2) is the latent innovation sd at the long-run log-variance (a
//  | generated quantity, so summaries and gates that ask for tau work, as in m10); `tau_t` =
//  | exp(g_t / 2) is its path.
//  | 
//  | Sampler parameterisation (as m12): conditional on (g, h) and the hyperparameters the model
//  | is linear-Gaussian in mu; quasi-differencing (y_t - psi .* y_{t-1}) removes eps exactly and
//  | mu is integrated out by m12's Kalman filter on the state (mu_t, mu_{t-1}), with the state
//  | innovation sd tau replaced by tau_t = exp(g_t / 2) in period t. kf_m13 is kf_m12 with that
//  | single change (tau -> tau_t[t] in the predicted-state factor G and in the backward-sampling
//  | variance), so with tau_t constant it returns m12's kf_m12 bit for bit. mu is drawn exactly
//  | from p(mu | y, g, h, theta) by FFBS. g and h are non-centred. rho_EI is sampled as
//  | z_EI = atanh(rho_EI) with the Jacobian, and 1 +- rho_EI, 1 +- r_1 are computed from z_EI
//  | without cancellation (m12's handling, unchanged). The filter and FFBS are checked against
//  | the dense multivariate-normal likelihood and smoother in tests/test_m13.R. LFO filter:
//  | R/lib/filter_m13.R / filter_m13.cpp (particles over (g, h), Kalman over (mu_t, mu_{t-1})).

functions {
  // m12's whitening of a 3-vector for covariance D Corr(r) D (see stan/m12_struct_sv_serial.stan).
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

  // Kalman filter for m13 given (g, h): m12's kf_m12 with the state innovation sd tau
  // replaced by tau_t[t] = exp(g_t / 2). Period t:
  //   predicted state  (mu_t, mu_{t-1}) = a + G z,  z ~ N(0, I_2),
  //     t >= 2: a = (c + rho m, m), G = [rho sqrt(P), tau_t[t]; sqrt(P), 0]
  //     t  = 1: a = (c/(1-rho), -), G row 1 = (0, tau_t[1]/sqrt(1-rho^2))
  // and the Cauchy-Binet forms of log det F and v'F^-1 v exactly as in kf_m12 (see the m12
  // header). Returns T x 6: [log p(y_t | y_{1:t-1}, g, h, theta), m_t, P_t, and, for t >= 2,
  // the backward-sampling coefficients of mu_{t-1} | mu_t, y_{1:t} ~
  // N(col4 + col5 * (mu_t - m_t), col6)] (cols 4-6 are 0 at t = 1).
  matrix kf_m13(matrix y, matrix h, vector tau_t, real c, real rho, real z_EI, vector psi) {
    int T = rows(y);
    matrix[T, 6] out;
    real ep = 2 * inv_logit(2 * z_EI);                 // 1 + r
    real em = 2 * inv_logit(-2 * z_EI);                // 1 - r
    real log_det_C = 2 * log2() + log_inv_logit(2 * z_EI) + log_inv_logit(-2 * z_EI);
    real d1 = 1 - psi[1] * psi[2];
    real k1 = sqrt((1 - square(psi[1])) * (1 - square(psi[2])));
    real dk1 = square(psi[1] - psi[2]) / (d1 + k1);
    real ep1 = (dk1 + ep * k1) / d1;                   // 1 + r_1
    real em1 = (dk1 + em * k1) / d1;                   // 1 - r_1
    vector[3] l1m = log1m(square(psi));                // log(1 - psi_k^2)
    vector[3] s1 = exp(0.5 * l1m);                     // sqrt(1 - psi_k^2)
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
      real tt = tau_t[t];
      if (t == 1) {
        vector[3] sc1 = sc .* s1;
        a1 = c / (1 - rho);
        g11 = 0;
        g12 = tt / sqrt(1 - square(rho));
        w = wh3(y[1]' - a1, sc1, ep1, em1);
        b1 = rep_vector(0, 3);
        b2 = wh3(rep_vector(g12, 3), sc1, ep1, em1);
        lH = sum(h[1]) - sum(l1m) + log(ep1) + log(em1);
      } else {
        real m = out[t - 1, 2];
        sP = sqrt(out[t - 1, 3]);
        a1 = c + rho * m;
        g11 = rho * sP;
        g12 = tt;
        w = wh3((y[t]' - a1) - psi .* (y[t - 1]' - m), sc, ep, em);
        b1 = wh3(sP * (rho - psi), sc, ep, em);
        b2 = wh3(rep_vector(tt, 3), sc, ep, em);
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
        real g1aG1 = square(g11) + square(g12) + dot_self(g11 * b2 - g12 * b1);
        out[t, 1] = -0.5 * (3 * log(2 * pi()) + lH + log(den) + num / den);
        out[t, 2] = a1 + g11 * zh1 + g12 * zh2;
        out[t, 3] = g1aG1 / den;
        if (t == 1) {
          out[t, 4] = 0;
          out[t, 5] = 0;
          out[t, 6] = 0;
        } else {
          // Var(mu_{t-1} | mu_t, y_{1:t}) = det(G)^2 / (den Var(mu_t | y_{1:t})) = tau_t^2 P / g1aG1.
          real cov12 = sP * (g11 + dot_product(b2, g11 * b2 - g12 * b1)) / den;
          out[t, 4] = out[t - 1, 2] + sP * zh1;
          out[t, 5] = cov12 / out[t, 3];
          out[t, 6] = square(tt * sP) / g1aG1;
        }
      }
    }
    return out;
  }

  // Backward sampling of mu_{1:T} | y, g, h, theta from f = kf_m13(...) (m12's ffbs_m12_rng):
  //   mu_T ~ N(m_T, P_T);  mu_t | mu_{t+1} ~ N(f[t+1,4] + f[t+1,5] (mu_{t+1} - m_{t+1}), f[t+1,6]).
  vector ffbs_m13_rng(matrix f) {
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
  real mu_bar;                       // long-run mean of latent growth (c = mu_bar (1 - rho))
  real<lower=0, upper=1> rho;
  real mu_g;                         // latent log-variance: long-run mean
  real<lower=-1, upper=1> phi_g;     //   persistence
  real<lower=0> s_g;                 //   stationary sd (omega_g = s_g sqrt(1 - phi_g^2))
  vector[K] mu_h;
  vector<lower=-1, upper=1>[K] phi;
  vector<lower=0>[K] s_h;            // stationary sd of h_k (omega_k = s_h_k sqrt(1 - phi_k^2))
  real z_EI;                         // atanh(rho_EI); see m9
  vector<lower=-1, upper=1>[K] psi;  // measurement-error persistence (as m5, m12)
  vector[T] g_raw;
  matrix[T, K] h_raw;
}
transformed parameters {
  real c = mu_bar * (1 - rho);
  real<lower=0> omega_g = s_g * sqrt((1 - phi_g) * (1 + phi_g));
  vector<lower=0>[K] omega = s_h .* sqrt((1 - phi) .* (1 + phi));
  real<lower=-1, upper=1> rho_EI = tanh(z_EI);
  vector[T] g;
  matrix[T, K] h;
  g[1] = mu_g + s_g * g_raw[1];
  for (t in 2:T)
    g[t] = mu_g + phi_g * (g[t - 1] - mu_g) + omega_g * g_raw[t];
  for (k in 1:K) {
    h[1, k] = mu_h[k] + s_h[k] * h_raw[1, k];
    for (t in 2:T)
      h[t, k] = mu_h[k] + phi[k] * (h[t-1, k] - mu_h[k]) + omega[k] * h_raw[t, k];
  }
}
model {
  mu_bar  ~ normal(3, 2);
  rho     ~ beta(2.625, 2.625);
  mu_g    ~ normal(2, 1);
  phi_g   ~ uniform(-1, 1);
  s_g     ~ normal(0, 1);
  mu_h    ~ normal(2, 1);
  phi     ~ uniform(-1, 1);
  s_h     ~ normal(0, 1);
  psi     ~ normal(0, 0.5);
  // rho_EI ~ normal(0, 0.5) on (-1, 1), written on z_EI with the Jacobian (as m9, m12, m13).
  target += normal_lpdf(rho_EI | 0, 0.5)
            + 2 * log2() + log_inv_logit(2 * z_EI) + log_inv_logit(-2 * z_EI);
  g_raw ~ std_normal();
  to_vector(h_raw) ~ std_normal();

  // log p(y | g, h, theta) with mu and eps integrated out (m13's filter, unchanged).
  target += sum(kf_m13(y, h, exp(0.5 * g), c, rho, z_EI, psi)[, 1]);
}
generated quantities {
  real tau = exp(0.5 * mu_g);        // latent innovation sd at g = mu_g
  vector[T] tau_t = exp(0.5 * g);    // latent innovation sd path
  vector[T] mu;                      // exact draw from p(mu | y, g, h, theta) (FFBS)
  matrix[T, K] sigma_t = exp(0.5 * h);   // sd of the measurement innovation u_t (m12's)
  vector[K] eps_sd_bar = exp(0.5 * mu_h) ./ sqrt(1 - square(psi));  // sd of eps at h = mu_h
  matrix[T, K] y_rep;                // y_t | mu_t, mu_{t-1}, y_{t-1}, h_t (m12's convention)
  vector[T] log_lik;                 // log p(y_t | mu_t, mu_{t-1}, y_{t-1}, h_t, theta)
  vector[T] log_lik_pred;            // log p(y_t | y_{1:t-1}, g, h, theta), mu integrated out
  {
    matrix[T, 6] f = kf_m13(y, h, tau_t, c, rho, z_EI, psi);
    log_lik_pred = f[, 1];
    mu = ffbs_m13_rng(f);
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
