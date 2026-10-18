// Filters for m13 (stan/m13_serial_latent_sv.stan): m12 + stochastic volatility on the
// latent innovation.
//   mu_t = c + rho mu_{t-1} + exp(g_t/2) eta_t,  mu_1 ~ N(c/(1-rho), exp(g_1)/(1-rho^2))
//   g_1 ~ N(mu_g, omega_g^2/(1-phi_g^2)),  g_t = mu_g + phi_g (g_{t-1} - mu_g) + omega_g zeta_t
//   y_t = mu_t 1 + eps_t,  eps_t = diag(psi) eps_{t-1} + u_t,  u_t ~ N(0, R_t),
//   R_t = D_t C D_t,  D_t = diag(exp(h_t/2)),  C = [1 r 0; r 1 0; 0 0 1]
//   eps_1 ~ N(0, V_1),  V_1[j,k] = R_1[j,k] / (1 - psi_j psi_k)
//   h_{1,k} ~ N(mu_h_k, omega_k^2/(1-phi_k^2)),  h_{t,k} = mu_h_k + phi_k (h_{t-1,k} - mu_h_k) + omega_k xi
//
// m13_step is m12_step (R/lib/filter_m12.cpp, a transcription of one period of kf_m12) with
// the state innovation sd tau replaced by the period's tau_t = exp(g_t/2); it is a
// transcription of one period of kf_m13 in the Stan file.
//
// kf_m13_cpp   : conditional-on-(g, h) Kalman filter; T x 6 output identical to Stan's kf_m13
//                (argument tau_t = exp(g/2), as in the Stan call).
// rbpf_m13_cpp : Rao-Blackwellised particle filter, particles over (h_t, g_t), exact Kalman
//                over (mu_t, mu_{t-1}) per particle (carrying the filtered moments (m, P) of
//                mu_{t-1}; y_{t-1} is common to all particles).
//   h: bootstrap (AR(1) transition), capped-sd h_1 proposal with exact importance
//      correction, drawn in exactly the order of rbpf_m12_cpp (all particles' h first,
//      then the per-particle g draws and Kalman steps).
//   g: prop = 0 bootstrap (transition; g_1 with the same capped-sd proposal as h_1), or
//      prop = 1 the locally adapted proposal of R/lib/filter_m10.cpp: given the particle's
//      h_t and Kalman moments, the predictive covariance of the quasi-differenced y_t is
//        F(g) = F0 + B exp(g) 1 1'   (t >= 2: F0 = R_t + P (rho - psi)(rho - psi)', B = 1;
//                                     t  = 1: F0 = V_1, B = 1/(1 - rho^2)),
//      so log p(y_t | g, ...) = const - 0.5 [log(1 + P s) - P q^2 / (1 + P s)], P = B e^g,
//      with s = 1'F0^-1 1 and q = 1'F0^-1 v. The proposal is
//        q(g) = (1 - eps) t_4(g^, 1.2/sqrt(-f''(g^))) + eps N(m_t, s_t^2)
//      (g^ = mode of N(g; m_t, s_t^2) p(y_t | g), found by safeguarded Newton steps) and the
//      weight carries the exact correction N(g; m_t, s_t^2) / q(g); the defensive share
//      bounds it by 1/eps. Either proposal gives an unbiased estimate of p(y_t | y_{1:t-1}).
//      With omega_g = 0 no g is drawn (g_t = mu_g), so the random-number stream, and with
//      mu_g = log(tau^2) the whole output, is that of rbpf_m12_cpp for the same seed.
//   Adaptive systematic resampling at ESS < ess_thr * N; N is a scalar or a length-T
//   schedule fixed in advance (R/lib/lfo.R), handled as in rbpf_m12_cpp (resample to N_{t+1}
//   when it differs from N_t). Output columns: (1) log p-hat(y_t | y_{1:t-1}, theta),
//   (2) the same (no sub-block), (3) ESS of the weights at t before resampling, (4) N_t.

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

static const double LOG2PI_13 = std::log(2.0 * M_PI);

struct M13Par {
  double c, rho;
  double ep, em, log_det_C;          // 1 + r, 1 - r, log det C
  double psi[3], s1[3], l1m[3];      // psi, sqrt(1 - psi^2), log(1 - psi^2)
  double ep1, em1;                   // 1 + r_1, 1 - r_1 (correlation of V_1)
};

static M13Par make_m13(double c, double rho, double ep, double em, const NumericVector& psi) {
  if (psi.size() != 3) stop("psi must have length 3");
  M13Par p;
  p.c = c; p.rho = rho; p.ep = ep; p.em = em;
  p.log_det_C = std::log(ep) + std::log(em);
  for (int k = 0; k < 3; ++k) {
    p.psi[k] = psi[k];
    p.l1m[k] = std::log1p(-psi[k] * psi[k]);
    p.s1[k] = std::exp(0.5 * p.l1m[k]);
  }
  const double d1 = 1 - psi[0] * psi[1];
  const double k1 = std::sqrt((1 - psi[0] * psi[0]) * (1 - psi[1] * psi[1]));
  const double dk1 = (psi[0] - psi[1]) * (psi[0] - psi[1]) / (d1 + k1);
  p.ep1 = (dk1 + ep * k1) / d1;
  p.em1 = (dk1 + em * k1) / d1;
  return p;
}

static inline void wh3_13(const double* x, const double* sc, double ep, double em, double* o) {
  const double xE = x[0] * sc[0], xI = x[1] * sc[1];
  o[0] = (xE + xI) / std::sqrt(2 * ep);
  o[1] = (xE - xI) / std::sqrt(2 * em);
  o[2] = x[2] * sc[2];
}
static inline void cross3_13(const double* a, const double* b, double* o) {
  o[0] = a[1] * b[2] - a[2] * b[1];
  o[1] = a[2] * b[0] - a[0] * b[2];
  o[2] = a[0] * b[1] - a[1] * b[0];
}
static inline double dot3_13(const double* a, const double* b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// One period (m12_step with tau -> tau). first = (t == 1). On entry (m, P) are the filtered
// moments of mu_{t-1} (ignored if first); on exit those of mu_t. out6 (may be null) receives
// the 6 columns of Stan's kf_m13 row. Returns log p(y_t | y_{1:t-1}, g, h, theta); if it is
// not finite the moments are left unchanged (the particle then gets weight 0).
static inline double m13_step(const M13Par& p, const double* yt, const double* yprev,
                              const double* h, double tau, bool first, double& m, double& P,
                              double* out6) {
  double sc[3], w[3], b1[3], b2[3], x[3], a1, g11, g12, sP = 0, lH;
  for (int k = 0; k < 3; ++k) sc[k] = std::exp(-0.5 * h[k]);
  const double sumh = h[0] + h[1] + h[2];
  if (first) {
    double sc1[3];
    for (int k = 0; k < 3; ++k) sc1[k] = sc[k] * p.s1[k];
    a1 = p.c / (1 - p.rho);
    g11 = 0;
    g12 = tau / std::sqrt(1 - p.rho * p.rho);
    for (int k = 0; k < 3; ++k) x[k] = yt[k] - a1;
    wh3_13(x, sc1, p.ep1, p.em1, w);
    b1[0] = b1[1] = b1[2] = 0;
    for (int k = 0; k < 3; ++k) x[k] = g12;
    wh3_13(x, sc1, p.ep1, p.em1, b2);
    lH = sumh - (p.l1m[0] + p.l1m[1] + p.l1m[2]) + std::log(p.ep1) + std::log(p.em1);
  } else {
    sP = std::sqrt(P);
    a1 = p.c + p.rho * m;
    g11 = p.rho * sP;
    g12 = tau;
    for (int k = 0; k < 3; ++k) x[k] = (yt[k] - a1) - p.psi[k] * (yprev[k] - m);
    wh3_13(x, sc, p.ep, p.em, w);
    for (int k = 0; k < 3; ++k) x[k] = sP * (p.rho - p.psi[k]);
    wh3_13(x, sc, p.ep, p.em, b1);
    for (int k = 0; k < 3; ++k) x[k] = tau;
    wh3_13(x, sc, p.ep, p.em, b2);
    lH = sumh + p.log_det_C;
  }
  double x12[3], x1w[3], x2w[3], gb[3];
  cross3_13(b1, b2, x12); cross3_13(b1, w, x1w); cross3_13(b2, w, x2w);
  const double d3 = dot3_13(x12, w);
  const double den = 1 + dot3_13(b1, b1) + dot3_13(b2, b2) + dot3_13(x12, x12);
  const double num = dot3_13(w, w) + dot3_13(x1w, x1w) + dot3_13(x2w, x2w) + d3 * d3;
  const double zh1 = (dot3_13(b1, w) - dot3_13(x12, x2w)) / den;
  const double zh2 = (dot3_13(b2, w) + dot3_13(x12, x1w)) / den;
  for (int k = 0; k < 3; ++k) gb[k] = g11 * b2[k] - g12 * b1[k];
  const double g1aG1 = g11 * g11 + g12 * g12 + dot3_13(gb, gb);
  const double lp = -0.5 * (3 * LOG2PI_13 + lH + std::log(den) + num / den);
  if (!std::isfinite(lp)) return -INFINITY;
  const double m_new = a1 + g11 * zh1 + g12 * zh2;
  const double P_new = g1aG1 / den;
  if (out6) {
    out6[0] = lp; out6[1] = m_new; out6[2] = P_new;
    if (first) { out6[3] = 0; out6[4] = 0; out6[5] = 0; }
    else {
      const double cov12 = sP * (g11 + dot3_13(b2, gb)) / den;
      out6[3] = m + sP * zh1;
      out6[4] = cov12 / P_new;
      out6[5] = (tau * sP) * (tau * sP) / g1aG1;
    }
  }
  m = m_new; P = P_new;
  return lp;
}

// s = 1'F0^-1 1, q = 1'F0^-1 v and B for the g proposal (see the header). Stable forms:
// with o = W 1, w = W v, b = W sP (rho - psi), (I + b b')^-1 = I - b b' / (1 + |b|^2) and the
// Lagrange / Binet-Cauchy identities give s = (|o|^2 + |b x o|^2) / (1 + |b|^2) and
// q = (o.w + (b x o).(b x w)) / (1 + |b|^2).
static inline void m13_sq(const M13Par& p, const double* yt, const double* yprev,
                          const double* h, bool first, double m, double P,
                          double& s, double& q, double& B) {
  double sc[3], x[3], o[3], w[3];
  for (int k = 0; k < 3; ++k) sc[k] = std::exp(-0.5 * h[k]);
  const double ones[3] = {1, 1, 1};
  if (first) {
    double sc1[3];
    for (int k = 0; k < 3; ++k) sc1[k] = sc[k] * p.s1[k];
    const double a1 = p.c / (1 - p.rho);
    for (int k = 0; k < 3; ++k) x[k] = yt[k] - a1;
    wh3_13(x, sc1, p.ep1, p.em1, w);
    wh3_13(ones, sc1, p.ep1, p.em1, o);
    s = dot3_13(o, o);
    q = dot3_13(o, w);
    B = 1.0 / (1 - p.rho * p.rho);
  } else {
    double b[3], bo[3], bw[3];
    const double sP = std::sqrt(P), a1 = p.c + p.rho * m;
    for (int k = 0; k < 3; ++k) x[k] = (yt[k] - a1) - p.psi[k] * (yprev[k] - m);
    wh3_13(x, sc, p.ep, p.em, w);
    wh3_13(ones, sc, p.ep, p.em, o);
    for (int k = 0; k < 3; ++k) x[k] = sP * (p.rho - p.psi[k]);
    wh3_13(x, sc, p.ep, p.em, b);
    cross3_13(b, o, bo); cross3_13(b, w, bw);
    const double d0 = 1 + dot3_13(b, b);
    s = (dot3_13(o, o) + dot3_13(bo, bo)) / d0;
    q = (dot3_13(o, w) + dot3_13(bo, bw)) / d0;
    B = 1.0;
  }
}

// [[Rcpp::export]]
NumericMatrix kf_m13_cpp(NumericMatrix y, NumericMatrix h, NumericVector tau_t, double c,
                         double rho, double ep, double em, NumericVector psi) {
  if (y.ncol() != 3 || h.ncol() != 3 || h.nrow() != y.nrow()) stop("y and h must be T x 3");
  if (tau_t.size() != y.nrow()) stop("tau_t must have length nrow(y)");
  M13Par p = make_m13(c, rho, ep, em, psi);
  const int T = y.nrow();
  NumericMatrix out(T, 6);
  double m = 0, P = 0, yt[3], yp[3], ht[3], o6[6];
  for (int t = 0; t < T; ++t) {
    for (int k = 0; k < 3; ++k) { yt[k] = y(t, k); ht[k] = h(t, k); if (t > 0) yp[k] = y(t - 1, k); }
    double lp = m13_step(p, yt, yp, ht, tau_t[t], t == 0, m, P, o6);
    if (!std::isfinite(lp)) { for (int j = 0; j < 6; ++j) out(t, j) = NA_REAL; out(t, 0) = lp; continue; }
    for (int j = 0; j < 6; ++j) out(t, j) = o6[j];
  }
  return out;
}

static inline double lse_vec13(const std::vector<double>& x, int n) {
  double mx = -INFINITY;
  for (int i = 0; i < n; ++i) if (x[i] > mx) mx = x[i];
  if (!std::isfinite(mx)) return -INFINITY;
  double s = 0;
  for (int i = 0; i < n; ++i) s += std::exp(x[i] - mx);
  return mx + std::log(s);
}
static inline int sched_max13(const IntegerVector& Nv, int T) {
  if (Nv.size() != 1 && Nv.size() != T) stop("N must have length 1 or nrow(y)");
  int mx = 0;
  for (int i = 0; i < Nv.size(); ++i) { if (Nv[i] < 1) stop("N must be >= 1"); if (Nv[i] > mx) mx = Nv[i]; }
  return mx;
}
static inline int sched_at13(const IntegerVector& Nv, int t) { return Nv.size() == 1 ? Nv[0] : Nv[t]; }

static inline double ldt4_13(double x) {           // log density of Student t_4
  return std::log(0.375) - 2.5 * std::log1p(x * x / 4.0);
}
// Locally adapted proposal for g (R/lib/filter_m10.cpp draw_adapted, with s and q given):
// l(P) = -0.5 [log(1 + P s) - P q^2 / (1 + P s)], P = A + B e^g. Draws g and returns
// log N(g; m, sg^2) - log q(g).
static double draw_adapted13(double s, double q, double A, double B, double m, double sg,
                             double eps, double& g) {
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
  double lp = -0.5 * (g - m) * (g - m) * ip - std::log(sg) - 0.5 * std::log(2 * M_PI);
  double lt = ldt4_13((g - x) / scl) - std::log(scl);
  double mx = std::max(lt, lp);
  double lq = mx + std::log((1 - eps) * std::exp(lt - mx) + eps * std::exp(lp - mx));
  return lp - lq;
}

// [[Rcpp::export]]
NumericMatrix rbpf_m13_cpp(NumericMatrix y, double c, double rho, double ep, double em,
                           NumericVector psi, double mu_g, double phi_g, double omega_g,
                           NumericVector mu_h, NumericVector phi, NumericVector omega,
                           IntegerVector Nv, double ess_thr, double h1_sd_max,
                           int prop, double eps_def) {
  if (y.ncol() != 3) stop("y must be T x 3");
  if (prop != 0 && prop != 1) stop("prop must be 0 or 1");
  M13Par p = make_m13(c, rho, ep, em, psi);
  const int T = y.nrow(), K = 3;
  const int Nmax = sched_max13(Nv, T);
  int N = sched_at13(Nv, 0);
  NumericMatrix out(T, 4);
  std::vector<double> hp(Nmax * K), hp_new(Nmax * K), gp(Nmax, mu_g), gp_new(Nmax);
  std::vector<double> mp(Nmax, 0.0), Pp(Nmax, 0.0), mp_new(Nmax), Pp_new(Nmax);
  std::vector<double> logW(Nmax, -std::log((double)N)), lw(Nmax), Wn(Nmax);
  double yt[3], yp[3];
  const double sd0g = omega_g / std::sqrt(1 - phi_g * phi_g);
  // Initial h_1: identical to rbpf_m12_cpp (same proposal, same correction, same RNG order).
  for (int i = 0; i < N; ++i) {
    for (int k = 0; k < K; ++k) {
      double sd0 = omega[k] / std::sqrt(1 - phi[k] * phi[k]);
      double sp = std::min(sd0, h1_sd_max);
      double e = R::norm_rand();
      hp[i * K + k] = mu_h[k] + sp * e;
      if (sp < sd0) logW[i] += std::log(sp / sd0) + 0.5 * e * e * (1.0 - (sp * sp) / (sd0 * sd0));
    }
  }
  // Initial g_1 under the bootstrap proposal (capped sd, exact correction), after all h_1.
  if (prop == 0 && sd0g > 0) {
    for (int i = 0; i < N; ++i) {
      double sp = std::min(sd0g, h1_sd_max);
      double e = R::norm_rand();
      gp[i] = mu_g + sp * e;
      if (sp < sd0g) logW[i] += std::log(sp / sd0g) + 0.5 * e * e * (1.0 - (sp * sp) / (sd0g * sd0g));
    }
  }
  for (int t = 0; t < T; ++t) {
    if (t > 0) {
      for (int i = 0; i < N; ++i)
        for (int k = 0; k < K; ++k) {
          double& hk = hp[i * K + k];
          hk = mu_h[k] + phi[k] * (hk - mu_h[k]) + omega[k] * R::norm_rand();
        }
    }
    for (int k = 0; k < K; ++k) { yt[k] = y(t, k); if (t > 0) yp[k] = y(t - 1, k); }
    for (int i = 0; i < N; ++i) {
      // transition (prior) moments of g_t for this particle
      const double mg = (t == 0) ? mu_g : mu_g + phi_g * (gp[i] - mu_g);
      const double sg = (t == 0) ? sd0g : omega_g;
      if (!(sg > 0)) {
        gp[i] = mg;
      } else if (prop == 1) {
        double s, q, B;
        m13_sq(p, yt, yp, &hp[i * K], t == 0, mp[i], Pp[i], s, q, B);
        logW[i] += draw_adapted13(s, q, 0.0, B, mg, sg, eps_def, gp[i]);
      } else if (t > 0) {
        gp[i] = mg + sg * R::norm_rand();
      }                                            // prop == 0, t == 0: drawn above
      double lf = m13_step(p, yt, yp, &hp[i * K], std::exp(0.5 * gp[i]), t == 0, mp[i], Pp[i],
                           nullptr);
      lw[i] = logW[i] + lf;
    }
    out(t, 0) = lse_vec13(lw, N);
    out(t, 1) = out(t, 0);
    double ess;
    if (!std::isfinite(out(t, 0))) {
      std::fill(logW.begin(), logW.end(), -std::log((double)N));
      std::fill(Wn.begin(), Wn.end(), 1.0 / N);
      ess = 0;
    } else {
      double s2 = 0;
      for (int i = 0; i < N; ++i) { logW[i] = lw[i] - out(t, 0); Wn[i] = std::exp(logW[i]); s2 += Wn[i] * Wn[i]; }
      ess = 1.0 / s2;
    }
    out(t, 2) = ess;
    out(t, 3) = N;
    const int Nn = (t < T - 1) ? sched_at13(Nv, t + 1) : N;
    if ((ess < ess_thr * N || Nn != N) && t < T - 1) {   // systematic resampling (as rbpf_m12_cpp)
      double u0 = R::unif_rand() / Nn, cum = Wn[0];
      int j = 0;
      for (int i = 0; i < Nn; ++i) {
        double ui = u0 + (double)i / Nn;
        while (ui > cum && j < N - 1) { ++j; cum += Wn[j]; }
        for (int k = 0; k < K; ++k) hp_new[i * K + k] = hp[j * K + k];
        gp_new[i] = gp[j]; mp_new[i] = mp[j]; Pp_new[i] = Pp[j];
      }
      hp.swap(hp_new); gp.swap(gp_new); mp.swap(mp_new); Pp.swap(Pp_new);
      N = Nn;
      std::fill(logW.begin(), logW.begin() + N, -std::log((double)N));
    }
  }
  return out;
}
