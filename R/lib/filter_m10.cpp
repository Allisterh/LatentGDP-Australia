// Rao-Blackwellised particle filter for m10 (stan/m10_labour_sv.stan): four signals
// (E, I, P, hours) with stochastic volatility on every measurement error AND on the
// latent innovation.
//
//   mu_1 ~ N(c/(1-rho), exp(g_1)/(1-rho^2)),  mu_t = c + rho mu_{t-1} + exp(g_t/2) eta_t
//   g_1 ~ N(mu_g, omega_g^2/(1-phi_g^2)),     g_t = mu_g + phi_g (g_{t-1} - mu_g) + omega_g zeta_t
//   h_{1,k} ~ N(mu_h_k, omega_k^2/(1-phi_k^2)), h_{t,k} = mu_h_k + phi_k (h_{t-1,k} - mu_h_k) + omega_k xi_{t,k}
//   y_{t,k} = d_k + Z_k mu_t + exp(h_{t,k}/2) e_{t,k},  e independent N(0, 1)
//
// Particles carry (g_t, h_t); mu is integrated out per particle by a scalar Kalman
// filter. A constant variance is the special case omega = 0 (then g or h_k stays at
// its mean), so the same code also gives the Kalman filter of m8 (all omegas 0) and
// the old m10 (omega_g = 0, mu_g = log tau^2); the tests use both.
//
// Proposal. h: the AR(1) transition (bootstrap), as in R/lib/filters.cpp. g: either
// the transition (prop = 0, bootstrap) or (prop = 1) a locally adapted proposal: for
// each particle, given its h_t and Kalman moments, the conditional density of y_t is
// a closed-form function of g_t, so the mode g^ of N(g; m_t, s_t^2) p(y_t | g, ...) is
// found by safeguarded Newton steps (analytic derivatives), and g_t is drawn from
//   q = (1 - eps) t_4(g^, 1.2 / sqrt(-f''(g^))) + eps N(m_t, s_t^2)
// (m_t, s_t = transition mean and sd; at t = 1 the stationary ones). The weights
// carry the exact correction N(g; m_t, s_t^2) / q(g), so the estimate is unbiased for
// either proposal; the adapted one matters at latent-volatility jumps such as 2020 Q2,
// where the bootstrap proposal leaves a handful of useful particles (min ESS 1-11 of
// 2000) and the estimate is several nats low. The defensive component bounds the
// weights by 1 / eps. Initial h_1 (and g_1 under prop = 0): proposal sd capped at
// h1_sd_max with the exact correction (same as filters.cpp, same caveat when the
// stationary sd exceeds sqrt(2) h1_sd_max).
//
// Output columns: (1) log p-hat(y_t | y_{1:t-1}, theta); (2) the same for the first
// nsub measures only (the GDP-only marginal, hours in the past only); (3) ESS of the
// weights at t before resampling; (4) the particle count N_t. Adaptive systematic
// resampling at ESS < thr * N. N is a scalar or a length-T schedule fixed in advance,
// handled as in R/lib/filters.cpp: when N_{t+1} != N_t the particles are resampled
// (systematic) to N_{t+1}, which keeps the estimate unbiased; a constant N leaves the
// random-number sequence unchanged.
//
// m10_kf_given_cpp: the conditional (given g and h) Kalman filter, the same
// computation as kf_tv in the Stan program; it is checked against log_lik_pred.

#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

static const double LOG2PI = std::log(2.0 * M_PI);

// Conditional step with independent measurement errors. On entry (a, P) are the
// predicted moments; on exit the filtered moments. Stable quadratic form:
//   e' F^{-1} e = Q - P q^2 / g = (Q + P sum_{i<j} (z_i w_j - z_j w_i)^2) / g.
static inline void step(int K, int nsub, const double* Z, const double* d,
                        const double* yt, const double* h, double& a, double& P,
                        double& lp_full, double& lp_sub, double* z, double* w) {
  double sumh = 0, sumh_sub = 0;
  for (int k = 0; k < K; ++k) {
    double u = std::exp(-0.5 * h[k]);
    z[k] = Z[k] * u;
    w[k] = (yt[k] - d[k] - a * Z[k]) * u;
    sumh += h[k];
    if (k < nsub) sumh_sub += h[k];
  }
  double s = 0, q = 0, Q = 0, cross = 0, ss = 0, qs = 0, Qs = 0, crs = 0;
  for (int i = 0; i < K; ++i) {
    s += z[i] * z[i]; q += z[i] * w[i]; Q += w[i] * w[i];
    if (i < nsub) { ss += z[i] * z[i]; qs += z[i] * w[i]; Qs += w[i] * w[i]; }
    for (int j = i + 1; j < K; ++j) {
      double c = z[i] * w[j] - z[j] * w[i];
      cross += c * c;
      if (j < nsub) crs += c * c;
    }
  }
  double g = 1.0 + P * s, gs = 1.0 + P * ss;
  lp_full = -0.5 * (K * LOG2PI + sumh + std::log(g) + (Q + P * cross) / g);
  lp_sub = -0.5 * (nsub * LOG2PI + sumh_sub + std::log(gs) + (Qs + P * crs) / gs);
  if (!std::isfinite(lp_full) || !std::isfinite(lp_sub)) {
    lp_full = -INFINITY; lp_sub = -INFINITY;
    return;
  }
  a = a + P * q / g;
  P = P / g;
}

// [[Rcpp::export]]
NumericMatrix m10_kf_given_cpp(NumericMatrix y, NumericMatrix h, NumericVector gpath,
                               NumericVector Z, NumericVector d, int nsub,
                               double c, double rho) {
  const int T = y.nrow(), K = y.ncol();
  NumericMatrix out(T, 2);
  double a = c / (1 - rho), P = std::exp(gpath[0]) / (1 - rho * rho);
  std::vector<double> z(K), w(K), yt(K), ht(K);
  for (int t = 0; t < T; ++t) {
    for (int k = 0; k < K; ++k) { yt[k] = y(t, k); ht[k] = h(t, k); }
    double lf, ls;
    step(K, nsub, Z.begin(), d.begin(), yt.data(), ht.data(), a, P, lf, ls, z.data(), w.data());
    out(t, 0) = lf; out(t, 1) = ls;
    if (t < T - 1) {
      a = c + rho * a;
      P = rho * rho * P + std::exp(gpath[t + 1]);
    }
  }
  return out;
}

static inline double lse_vec(const std::vector<double>& x, int n) {
  double mx = -INFINITY;
  for (int i = 0; i < n; ++i) if (x[i] > mx) mx = x[i];
  if (!std::isfinite(mx)) return -INFINITY;
  double s = 0;
  for (int i = 0; i < n; ++i) s += std::exp(x[i] - mx);
  return mx + std::log(s);
}
static inline int sched_max10(const IntegerVector& Nv, int T) {
  if (Nv.size() != 1 && Nv.size() != T) stop("N must have length 1 or nrow(y)");
  int mx = 0;
  for (int i = 0; i < Nv.size(); ++i) { if (Nv[i] < 1) stop("N must be >= 1"); if (Nv[i] > mx) mx = Nv[i]; }
  return mx;
}
static inline int sched_at10(const IntegerVector& Nv, int t) { return Nv.size() == 1 ? Nv[0] : Nv[t]; }

// Draw x ~ N(m, sd0^2) through a proposal with sd sp (sp <= sd0 allowed); returns the
// log importance correction log N(x; m, sd0^2) - log N(x; m, sp^2).
static inline double draw_capped(double m, double sd0, double sp, double& x) {
  double e = R::norm_rand();
  x = m + sp * e;
  if (sp < sd0) return std::log(sp / sd0) + 0.5 * e * e * (1.0 - (sp * sp) / (sd0 * sd0));
  return 0.0;
}

// Locally adapted proposal for g_t (see header). Given the particle's h_t, predicted
// mean a and P(g) = A + B exp(g), the conditional log density of y_t is
//   l(P) = const - 0.5 [log(1 + P s) + (Q + P C) / (1 + P s)]
// with s, q, Q, C the whitened sums of step(); dl/dP = -0.5 (s/G - q^2/G^2),
// d2l/dP2 = 0.5 s^2/G^2 - q^2 s/G^3 (G = 1 + P s), dP/dg = P - A.
// Draws g and returns log N(g; m, sg^2) - log q(g).
static inline double ldt4(double x) {              // log density of Student t_4
  return std::log(0.375) - 2.5 * std::log1p(x * x / 4.0);
}
static double draw_adapted(int K, const double* Z, const double* d, const double* yt,
                           const double* h, double a, double A, double B,
                           double m, double sg, double eps, double& g) {
  double s = 0, q = 0;
  for (int k = 0; k < K; ++k) {
    double u = std::exp(-0.5 * h[k]);
    double zk = Z[k] * u, wk = (yt[k] - d[k] - a * Z[k]) * u;
    s += zk * zk; q += zk * wk;
  }
  const double ip = 1.0 / (sg * sg);
  auto derivs = [&](double x, double& f1, double& f2) {
    double D = B * std::exp(x), P = A + D, G = 1 + P * s;
    double l1 = -0.5 * (s / G - q * q / (G * G));
    double l2 = 0.5 * s * s / (G * G) - q * q * s / (G * G * G);
    f1 = -(x - m) * ip + l1 * D;
    f2 = -ip + l2 * D * D + l1 * D;
  };
  double x = m, f1, f2;
  for (int it = 0; it < 30; ++it) {
    derivs(x, f1, f2);
    double H = (f2 < -0.01 * ip) ? f2 : -ip;          // safeguard: keep the step uphill
    double step = -f1 / H;
    if (step > 3) step = 3;
    if (step < -3) step = -3;
    x += step;
    if (std::fabs(step) < 1e-8) break;
  }
  derivs(x, f1, f2);
  double H = (f2 < -0.01 * ip) ? f2 : -ip;
  double scl = 1.2 / std::sqrt(-H);
  if (R::unif_rand() < eps) g = m + sg * R::norm_rand();
  else g = x + scl * R::rt(4.0);
  // log densities at g: prior N(m, sg^2), t_4 component, mixture proposal q
  double lp = -0.5 * (g - m) * (g - m) * ip - std::log(sg) - 0.5 * std::log(2 * M_PI);
  double lt = ldt4((g - x) / scl) - std::log(scl);
  double mx = std::max(lt, lp);
  double lq = mx + std::log((1 - eps) * std::exp(lt - mx) + eps * std::exp(lp - mx));
  return lp - lq;
}

// [[Rcpp::export]]
NumericMatrix m10_rbpf_cpp(NumericMatrix y, NumericVector Z, NumericVector d, int nsub,
                           double c, double rho,
                           double mu_g, double phi_g, double omega_g,
                           NumericVector mu_h, NumericVector phi, NumericVector omega,
                           IntegerVector Nv, double ess_thr, double h1_sd_max,
                           int prop, double eps_def) {
  const int T = y.nrow(), K = y.ncol();
  const int Nmax = sched_max10(Nv, T);
  int N = sched_at10(Nv, 0);
  NumericMatrix out(T, 4);
  std::vector<double> hp(Nmax * K), hp_new(Nmax * K), gp(Nmax), gp_new(Nmax);
  std::vector<double> ap(Nmax), Pp(Nmax), ap_new(Nmax), Pp_new(Nmax);
  std::vector<double> logW(Nmax, -std::log((double)N)), lw(Nmax), lws(Nmax), Wn(Nmax);
  std::vector<double> z(K), w(K), yt(K);
  const double a0 = c / (1 - rho);
  const double sd0g = omega_g / std::sqrt(1 - phi_g * phi_g);
  for (int i = 0; i < N; ++i) {
    if (prop == 0) logW[i] += draw_capped(mu_g, sd0g, std::min(sd0g, h1_sd_max), gp[i]);
    ap[i] = a0;
    Pp[i] = 0;                                   // filtered variance of "t = 0": unused
    for (int k = 0; k < K; ++k) {
      double s0 = omega[k] / std::sqrt(1 - phi[k] * phi[k]);
      logW[i] += draw_capped(mu_h[k], s0, std::min(s0, h1_sd_max), hp[i * K + k]);
    }
  }
  for (int t = 0; t < T; ++t) {
    for (int k = 0; k < K; ++k) yt[k] = y(t, k);
    for (int i = 0; i < N; ++i) {
      if (t > 0)
        for (int k = 0; k < K; ++k) {
          double& hk = hp[i * K + k];
          hk = mu_h[k] + phi[k] * (hk - mu_h[k]) + omega[k] * R::norm_rand();
        }
      // prior of g_t, and the predicted variance P(g) = A + B exp(g)
      double m, sg, A, B;
      if (t == 0) { m = mu_g; sg = sd0g; A = 0; B = 1 / (1 - rho * rho); }
      else {
        m = mu_g + phi_g * (gp[i] - mu_g); sg = omega_g;
        A = rho * rho * Pp[i]; B = 1;
        ap[i] = c + rho * ap[i];
      }
      if (t > 0 && prop == 0) gp[i] = m + sg * R::norm_rand();
      if (prop == 1) {
        if (sg > 0) logW[i] += draw_adapted(K, Z.begin(), d.begin(), yt.data(), &hp[i * K],
                                            ap[i], A, B, m, sg, eps_def, gp[i]);
        else gp[i] = m;
      }
      Pp[i] = A + B * std::exp(gp[i]);           // predicted variance of mu_t
    }
    for (int i = 0; i < N; ++i) {
      double lf, ls;
      step(K, nsub, Z.begin(), d.begin(), yt.data(), &hp[i * K], ap[i], Pp[i], lf, ls,
           z.data(), w.data());
      lw[i] = logW[i] + lf;
      lws[i] = logW[i] + ls;
    }
    // logW = normalised weights of t-1 (sum 1) plus the proposal correction, so
    // sum_i exp(logW_i + lf_i) is the unbiased SMC estimate of p(y_t | y_{1:t-1}).
    out(t, 0) = lse_vec(lw, N);
    out(t, 1) = lse_vec(lws, N);
    double ess;
    double tot = out(t, 0);
    if (!std::isfinite(tot)) {
      std::fill(logW.begin(), logW.end(), -std::log((double)N));
      std::fill(Wn.begin(), Wn.end(), 1.0 / N);
      ess = 0;
    } else {
      double s2 = 0;
      for (int i = 0; i < N; ++i) { logW[i] = lw[i] - tot; Wn[i] = std::exp(logW[i]); s2 += Wn[i] * Wn[i]; }
      ess = 1.0 / s2;
    }
    out(t, 2) = ess;
    out(t, 3) = N;
    const int Nn = (t < T - 1) ? sched_at10(Nv, t + 1) : N;
    if ((ess < ess_thr * N || Nn != N) && t < T - 1) {
      double u0 = R::unif_rand() / Nn, cum = Wn[0];
      int j = 0;
      for (int i = 0; i < Nn; ++i) {
        double ui = u0 + (double)i / Nn;
        while (ui > cum && j < N - 1) { ++j; cum += Wn[j]; }
        for (int k = 0; k < K; ++k) hp_new[i * K + k] = hp[j * K + k];
        gp_new[i] = gp[j]; ap_new[i] = ap[j]; Pp_new[i] = Pp[j];
      }
      hp.swap(hp_new); gp.swap(gp_new); ap.swap(ap_new); Pp.swap(Pp_new);
      N = Nn;
      std::fill(logW.begin(), logW.begin() + N, -std::log((double)N));
    }
  }
  return out;
}
