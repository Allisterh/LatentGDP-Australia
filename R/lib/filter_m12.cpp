// Filters for m12 (stan/m12_struct_sv_serial.stan): m9 + AR(1) measurement errors.
//   mu_t = c + rho mu_{t-1} + tau eta_t,  mu_1 ~ N(c/(1-rho), tau^2/(1-rho^2))
//   y_t = mu_t 1 + eps_t,  eps_t = diag(psi) eps_{t-1} + u_t,  u_t ~ N(0, R_t),
//   R_t = D_t C D_t,  D_t = diag(exp(h_t/2)),  C = [1 r 0; r 1 0; 0 0 1]
//   eps_1 ~ N(0, V_1),  V_1[j,k] = R_1[j,k] / (1 - psi_j psi_k)
//   h_{1,k} ~ N(mu_h_k, omega_k^2/(1-phi_k^2)),  h_{t,k} = mu_h_k + phi_k (h_{t-1,k} - mu_h_k) + omega_k xi
//
// m12_step is a line-by-line transcription of one period of kf_m12() in the Stan file
// (same quasi-differenced observation, same 2-dim predicted state written as a + G z, same
// Cauchy-Binet form of the quadratic form and log-determinant); see the comments there.
//
// kf_m12_cpp   : conditional-on-h Kalman filter, T x 6 output identical to Stan's kf_m12.
// rbpf_m12_cpp : Rao-Blackwellised particle filter (particles over h, Kalman over mu given
//                h). Same particle scheme as rbpf_cpp in filters.cpp (bootstrap proposal,
//                capped-sd h_1 proposal with exact importance correction, adaptive
//                systematic resampling at ESS < ess_thr * N) and the same order of random
//                number draws, so at psi = 0 it reproduces rbpf_cpp for m9 for the same
//                seed. Output columns: (1) log p-hat(y_t | y_{1:t-1}, theta), (2) the same
//                (all three measures are GDP measures, there is no sub-block), (3) ESS of
//                the weights at t before resampling.
//                Per particle the filter carries (m, P) = moments of mu_{t-1} | y_{1:t-1}, h;
//                the observation at t also needs y_{t-1}, which is common to all particles.
//                N is a scalar or a length-T schedule fixed in advance, handled exactly as in
//                rbpf_cpp (resample to N_{t+1} when it differs from N_t); column (4) is N_t.

#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

static const double LOG2PI_12 = std::log(2.0 * M_PI);

struct M12Par {
  double c, rho, tau;
  double ep, em, log_det_C;          // 1 + r, 1 - r, log det C
  double psi[3], s1[3], l1m[3];      // psi, sqrt(1 - psi^2), log(1 - psi^2)
  double ep1, em1;                   // 1 + r_1, 1 - r_1 (correlation of V_1)
};

// ep = 1 + r and em = 1 - r are passed in (computed stably by the R wrapper).
static M12Par make_m12(double c, double rho, double tau, double ep, double em,
                       const NumericVector& psi) {
  if (psi.size() != 3) stop("psi must have length 3");
  M12Par p;
  p.c = c; p.rho = rho; p.tau = tau; p.ep = ep; p.em = em;
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

static inline void wh3(const double* x, const double* sc, double ep, double em, double* o) {
  const double xE = x[0] * sc[0], xI = x[1] * sc[1];
  o[0] = (xE + xI) / std::sqrt(2 * ep);
  o[1] = (xE - xI) / std::sqrt(2 * em);
  o[2] = x[2] * sc[2];
}
static inline void cross3(const double* a, const double* b, double* o) {
  o[0] = a[1] * b[2] - a[2] * b[1];
  o[1] = a[2] * b[0] - a[0] * b[2];
  o[2] = a[0] * b[1] - a[1] * b[0];
}
static inline double dot3(const double* a, const double* b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// One period. first = (t == 1). On entry (m, P) are the filtered moments of mu_{t-1}
// (ignored if first); on exit those of mu_t. out6 (may be null) receives the 6 columns of
// Stan's kf_m12 row. Returns log p(y_t | y_{1:t-1}, h, theta); if it is not finite the
// moments are left unchanged (the particle then gets weight 0).
static inline double m12_step(const M12Par& p, const double* yt, const double* yprev,
                              const double* h, bool first, double& m, double& P,
                              double* out6) {
  double sc[3], w[3], b1[3], b2[3], x[3], a1, g11, g12, sP = 0, lH;
  for (int k = 0; k < 3; ++k) sc[k] = std::exp(-0.5 * h[k]);
  const double sumh = h[0] + h[1] + h[2];
  if (first) {
    double sc1[3];
    for (int k = 0; k < 3; ++k) sc1[k] = sc[k] * p.s1[k];
    a1 = p.c / (1 - p.rho);
    g11 = 0;
    g12 = p.tau / std::sqrt(1 - p.rho * p.rho);
    for (int k = 0; k < 3; ++k) x[k] = yt[k] - a1;
    wh3(x, sc1, p.ep1, p.em1, w);
    b1[0] = b1[1] = b1[2] = 0;
    for (int k = 0; k < 3; ++k) x[k] = g12;
    wh3(x, sc1, p.ep1, p.em1, b2);
    lH = sumh - (p.l1m[0] + p.l1m[1] + p.l1m[2]) + std::log(p.ep1) + std::log(p.em1);
  } else {
    sP = std::sqrt(P);
    a1 = p.c + p.rho * m;
    g11 = p.rho * sP;
    g12 = p.tau;
    for (int k = 0; k < 3; ++k) x[k] = (yt[k] - a1) - p.psi[k] * (yprev[k] - m);
    wh3(x, sc, p.ep, p.em, w);
    for (int k = 0; k < 3; ++k) x[k] = sP * (p.rho - p.psi[k]);
    wh3(x, sc, p.ep, p.em, b1);
    for (int k = 0; k < 3; ++k) x[k] = p.tau;
    wh3(x, sc, p.ep, p.em, b2);
    lH = sumh + p.log_det_C;
  }
  double x12[3], x1w[3], x2w[3], gb[3];
  cross3(b1, b2, x12); cross3(b1, w, x1w); cross3(b2, w, x2w);
  const double d3 = dot3(x12, w);
  const double den = 1 + dot3(b1, b1) + dot3(b2, b2) + dot3(x12, x12);
  const double num = dot3(w, w) + dot3(x1w, x1w) + dot3(x2w, x2w) + d3 * d3;
  const double zh1 = (dot3(b1, w) - dot3(x12, x2w)) / den;
  const double zh2 = (dot3(b2, w) + dot3(x12, x1w)) / den;
  for (int k = 0; k < 3; ++k) gb[k] = g11 * b2[k] - g12 * b1[k];
  const double g1aG1 = g11 * g11 + g12 * g12 + dot3(gb, gb);
  const double lp = -0.5 * (3 * LOG2PI_12 + lH + std::log(den) + num / den);
  if (!std::isfinite(lp)) return -INFINITY;
  const double m_new = a1 + g11 * zh1 + g12 * zh2;
  const double P_new = g1aG1 / den;
  if (out6) {
    out6[0] = lp; out6[1] = m_new; out6[2] = P_new;
    if (first) { out6[3] = 0; out6[4] = 0; out6[5] = 0; }
    else {
      const double cov12 = sP * (g11 + dot3(b2, gb)) / den;
      out6[3] = m + sP * zh1;
      out6[4] = cov12 / P_new;
      out6[5] = (p.tau * sP) * (p.tau * sP) / g1aG1;
    }
  }
  m = m_new; P = P_new;
  return lp;
}

// [[Rcpp::export]]
NumericMatrix kf_m12_cpp(NumericMatrix y, NumericMatrix h, double c, double rho, double tau,
                         double ep, double em, NumericVector psi) {
  if (y.ncol() != 3 || h.ncol() != 3 || h.nrow() != y.nrow()) stop("y and h must be T x 3");
  M12Par p = make_m12(c, rho, tau, ep, em, psi);
  const int T = y.nrow();
  NumericMatrix out(T, 6);
  double m = 0, P = 0, yt[3], yp[3], ht[3], o6[6];
  for (int t = 0; t < T; ++t) {
    for (int k = 0; k < 3; ++k) { yt[k] = y(t, k); ht[k] = h(t, k); if (t > 0) yp[k] = y(t - 1, k); }
    double lp = m12_step(p, yt, yp, ht, t == 0, m, P, o6);
    if (!std::isfinite(lp)) { for (int j = 0; j < 6; ++j) out(t, j) = NA_REAL; out(t, 0) = lp; continue; }
    for (int j = 0; j < 6; ++j) out(t, j) = o6[j];
  }
  return out;
}

static inline double lse_vec12(const std::vector<double>& x, int n) {
  double mx = -INFINITY;
  for (int i = 0; i < n; ++i) if (x[i] > mx) mx = x[i];
  if (!std::isfinite(mx)) return -INFINITY;
  double s = 0;
  for (int i = 0; i < n; ++i) s += std::exp(x[i] - mx);
  return mx + std::log(s);
}
static inline int sched_max12(const IntegerVector& Nv, int T) {
  if (Nv.size() != 1 && Nv.size() != T) stop("N must have length 1 or nrow(y)");
  int mx = 0;
  for (int i = 0; i < Nv.size(); ++i) { if (Nv[i] < 1) stop("N must be >= 1"); if (Nv[i] > mx) mx = Nv[i]; }
  return mx;
}
static inline int sched_at12(const IntegerVector& Nv, int t) { return Nv.size() == 1 ? Nv[0] : Nv[t]; }

// [[Rcpp::export]]
NumericMatrix rbpf_m12_cpp(NumericMatrix y, double c, double rho, double tau,
                           double ep, double em, NumericVector psi,
                           NumericVector mu_h, NumericVector phi, NumericVector omega,
                           IntegerVector Nv, double ess_thr, double h1_sd_max) {
  if (y.ncol() != 3) stop("y must be T x 3");
  M12Par p = make_m12(c, rho, tau, ep, em, psi);
  const int T = y.nrow(), K = 3;
  const int Nmax = sched_max12(Nv, T);
  int N = sched_at12(Nv, 0);
  NumericMatrix out(T, 4);
  std::vector<double> hp(Nmax * K), hp_new(Nmax * K), mp(Nmax, 0.0), Pp(Nmax, 0.0), mp_new(Nmax), Pp_new(Nmax);
  std::vector<double> logW(Nmax, -std::log((double)N)), lw(Nmax), Wn(Nmax);
  double yt[3], yp[3];
  // Initial h_1: identical to rbpf_cpp (same proposal, same correction, same RNG order).
  for (int i = 0; i < N; ++i) {
    for (int k = 0; k < K; ++k) {
      double sd0 = omega[k] / std::sqrt(1 - phi[k] * phi[k]);
      double sp = std::min(sd0, h1_sd_max);
      double e = R::norm_rand();
      hp[i * K + k] = mu_h[k] + sp * e;
      if (sp < sd0) logW[i] += std::log(sp / sd0) + 0.5 * e * e * (1.0 - (sp * sp) / (sd0 * sd0));
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
      double lf = m12_step(p, yt, yp, &hp[i * K], t == 0, mp[i], Pp[i], nullptr);
      lw[i] = logW[i] + lf;
    }
    out(t, 0) = lse_vec12(lw, N);
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
    const int Nn = (t < T - 1) ? sched_at12(Nv, t + 1) : N;
    if ((ess < ess_thr * N || Nn != N) && t < T - 1) {   // systematic resampling (as rbpf_cpp)
      double u0 = R::unif_rand() / Nn, cum = Wn[0];
      int j = 0;
      for (int i = 0; i < Nn; ++i) {
        double ui = u0 + (double)i / Nn;
        while (ui > cum && j < N - 1) { ++j; cum += Wn[j]; }
        for (int k = 0; k < K; ++k) hp_new[i * K + k] = hp[j * K + k];
        mp_new[i] = mp[j]; Pp_new[i] = Pp[j];
      }
      hp.swap(hp_new); mp.swap(mp_new); Pp.swap(Pp_new);
      N = Nn;
      std::fill(logW.begin(), logW.begin() + N, -std::log((double)N));
    }
  }
  return out;
}
