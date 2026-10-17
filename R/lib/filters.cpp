// Filters for marginal-likelihood model comparison (latent states integrated out).
// Plain Rcpp (no Armadillo / LAPACK): all matrices are tiny (state dim <= 3,
// observation dim <= 4), so Cholesky and triangular solves are written out.
//
// kf_loglik_cpp : Kalman filter for a general linear-Gaussian state-space model
//                   alpha_1 ~ N(a1, P1)
//                   alpha_{t+1} = cst + Tm alpha_t + eta_t,   eta_t ~ N(0, Q)
//                   y_t = d_t + Z_t alpha_t + e_t,             e_t ~ N(0, H_t)
//                 Missing entries of y (NA) are skipped. Observations are processed
//                 in blocks blk = 1..nblk in that order within each period, so
//                 column b of the output is log p(y_{t,b} | y_{1:t-1}, y_{t,<b}).
//                 This is exact only if H_t is block-diagonal across blocks
//                 (checked by the R wrapper). The row sum is log p(y_t | y_{1:t-1}).
//                 Z: array k x m x (1 or T); H: array k x k x (1 or T); d: (1 or T) x k.
//
// rbpf_cpp      : Rao-Blackwellised particle filter for the SV specs (m4, m9, m10):
//                   mu_1 ~ N(c/(1-rho), tau^2/(1-rho^2)),  mu_t = c + rho mu_{t-1} + tau eta_t
//                   h_{1,k} ~ N(mu_h_k, omega_k^2/(1-phi_k^2)),
//                   h_{t,k} = mu_h_k + phi_k (h_{t-1,k} - mu_h_k) + omega_k xi_{t,k}
//                   y_t = d + Z mu_t + e_t,  e_t ~ N(0, D_t C D_t),  D_t = diag(exp(h_t/2))
//                 Particles carry h_t (bootstrap proposal = the AR(1) transition);
//                 mu is integrated out per particle by a scalar Kalman filter.
//                 Output columns: (1) log p-hat(y_t | y_{1:t-1}, theta);
//                 (2) log p-hat(y_{t,1:nsub} | y_{1:t-1}, theta), the marginal of the
//                 first nsub measures;
//                 (3) effective sample size of the weights at t (before resampling).
//                 (4) the particle count N_t used at t.
//                 The per-period estimates multiply to an unbiased estimate of
//                 p(y_{1:T} | theta) (adaptive systematic resampling at ESS < thr * N).
//                 N may be a scalar or a length-T schedule N_t fixed in advance (the LFO
//                 code sets it from a pilot run, raising N where the ESS collapses).
//                 When N_{t+1} != N_t the weighted particle set at t is resampled
//                 (systematic) to N_{t+1} particles before propagation; resampling to any
//                 deterministic count keeps the estimate unbiased. With a constant N the
//                 random-number sequence, and hence the output, is unchanged.
//
// kf_given_h_cpp: the conditional (given h) scalar Kalman filter used inside the
//                 RBPF, exposed so it can be checked against the Stan programs'
//                 log_lik_pred output.

#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

static const double LOG2PI = std::log(2.0 * M_PI);

// In-place lower Cholesky of an n x n row-major matrix A. Returns false if not PD.
static bool chol_lower(std::vector<double>& A, int n) {
  for (int j = 0; j < n; ++j) {
    double s = A[j * n + j];
    for (int k = 0; k < j; ++k) s -= A[j * n + k] * A[j * n + k];
    if (!(s > 0)) return false;
    double ljj = std::sqrt(s);
    A[j * n + j] = ljj;
    for (int i = j + 1; i < n; ++i) {
      double v = A[i * n + j];
      for (int k = 0; k < j; ++k) v -= A[i * n + k] * A[j * n + k];
      A[i * n + j] = v / ljj;
    }
    for (int k = j + 1; k < n; ++k) A[j * n + k] = 0.0;
  }
  return true;
}
// Solve L x = b in place (L lower, row-major n x n), b has `ncol` columns stored
// row-major as n x ncol.
static void fwd_solve(const std::vector<double>& L, int n, std::vector<double>& B, int ncol) {
  for (int c = 0; c < ncol; ++c)
    for (int i = 0; i < n; ++i) {
      double v = B[i * ncol + c];
      for (int k = 0; k < i; ++k) v -= L[i * n + k] * B[k * ncol + c];
      B[i * ncol + c] = v / L[i * n + i];
    }
}

// [[Rcpp::export]]
NumericMatrix kf_loglik_cpp(NumericMatrix y, NumericVector Z, NumericMatrix d,
                            NumericVector H, NumericMatrix Tm, NumericVector cst,
                            NumericMatrix Q, NumericVector a1, NumericMatrix P1,
                            IntegerVector blk, int nblk) {
  const int T = y.nrow(), k = y.ncol(), m = a1.size();
  IntegerVector zd = Z.attr("dim"), hd = H.attr("dim");
  const int nZ = zd[2], nH = hd[2];
  if (zd[0] != k || zd[1] != m || hd[0] != k || hd[1] != k) stop("dimension mismatch");
  NumericMatrix out(T, nblk);
  std::vector<double> a(m), P(m * m), an(m), Pn(m * m);
  for (int i = 0; i < m; ++i) { a[i] = a1[i]; for (int j = 0; j < m; ++j) P[i * m + j] = P1(i, j); }
  std::vector<int> idx; idx.reserve(k);
  std::vector<double> F, v, ZP, Zb;
  for (int t = 0; t < T; ++t) {
    const int zs = (nZ == 1 ? 0 : t), hs = (nH == 1 ? 0 : t), ds = (d.nrow() == 1 ? 0 : t);
    for (int b = 1; b <= nblk; ++b) {
      idx.clear();
      for (int j = 0; j < k; ++j) if (blk[j] == b && !ISNAN(y(t, j))) idx.push_back(j);
      const int n = idx.size();
      if (n == 0) continue;
      Zb.assign(n * m, 0.0); v.assign(n, 0.0); ZP.assign(n * m, 0.0); F.assign(n * n, 0.0);
      for (int r = 0; r < n; ++r) {
        int j = idx[r];
        double pred = d(ds, j);
        for (int c = 0; c < m; ++c) {
          Zb[r * m + c] = Z[j + k * c + k * m * zs];
          pred += Zb[r * m + c] * a[c];
        }
        v[r] = y(t, j) - pred;
      }
      for (int r = 0; r < n; ++r)
        for (int c = 0; c < m; ++c) {
          double s = 0;
          for (int l = 0; l < m; ++l) s += Zb[r * m + l] * P[l * m + c];
          ZP[r * m + c] = s;
        }
      for (int r = 0; r < n; ++r)
        for (int c = 0; c <= r; ++c) {
          double s = 0;
          for (int l = 0; l < m; ++l) s += ZP[r * m + l] * Zb[c * m + l];
          s += 0.5 * (H[idx[r] + k * idx[c] + k * k * hs] + H[idx[c] + k * idx[r] + k * k * hs]);
          F[r * n + c] = s; F[c * n + r] = s;
        }
      if (!chol_lower(F, n)) stop("prediction covariance not positive definite at t = %d", t + 1);
      fwd_solve(F, n, v, 1);       // u = L^{-1} v
      fwd_solve(F, n, ZP, m);      // W = L^{-1} Z P  (n x m)
      double ld = 0, qf = 0;
      for (int r = 0; r < n; ++r) { ld += std::log(F[r * n + r]); qf += v[r] * v[r]; }
      out(t, b - 1) = -0.5 * (n * LOG2PI + 2.0 * ld + qf);
      for (int c = 0; c < m; ++c) {
        double s = 0;
        for (int r = 0; r < n; ++r) s += ZP[r * m + c] * v[r];
        a[c] += s;
      }
      for (int i = 0; i < m; ++i)
        for (int j = 0; j <= i; ++j) {
          double s = 0;
          for (int r = 0; r < n; ++r) s += ZP[r * m + i] * ZP[r * m + j];
          double val = 0.5 * (P[i * m + j] + P[j * m + i]) - s;
          P[i * m + j] = val; P[j * m + i] = val;
        }
    }
    // predict
    for (int i = 0; i < m; ++i) {
      double s = cst[i];
      for (int j = 0; j < m; ++j) s += Tm(i, j) * a[j];
      an[i] = s;
    }
    for (int i = 0; i < m; ++i)
      for (int j = 0; j < m; ++j) {
        double s = 0;
        for (int l = 0; l < m; ++l) s += Tm(i, l) * P[l * m + j];
        Pn[i * m + j] = s;                       // (Tm P)
      }
    for (int i = 0; i < m; ++i)
      for (int j = 0; j <= i; ++j) {
        double s = 0;
        for (int l = 0; l < m; ++l) s += Pn[i * m + l] * Tm(j, l);
        s += 0.5 * (Q(i, j) + Q(j, i));
        P[i * m + j] = s; P[j * m + i] = s;
      }
    a = an;
  }
  return out;
}

// ---- SV pieces ------------------------------------------------------------------
// Numerically stable conditional step. With D_t = diag(exp(h_t/2)) and C = L L',
// whiten: zt = L^{-1} D^{-1} Z, wt = L^{-1} D^{-1} e (e = y - d - a Z). Then
//   s = zt'zt, q = zt'wt, Q = wt'wt, g = 1 + P s, and
//   e' F^{-1} e = Q - P q^2 / g = (Q + P * sum_{i<j} (zt_i wt_j - zt_j wt_i)^2) / g
// (Lagrange identity: every term is non-negative, so there is no cancellation when a
// particle has a tiny measurement variance). log|F| = sum(h) + log|C| + log(g).
// The first nsub whitened components depend only on the first nsub measures (L is
// lower-triangular and its leading block is the Cholesky factor of C's leading block),
// so the same formulas over i, j < nsub give the marginal density of those measures.
struct SVSys {
  int K, nsub;
  std::vector<double> Z, d, Linv;
  double logdetC, logdetC_sub;
};

static inline void sv_quad(const double* zt, const double* wt, int n, double P,
                           double& s, double& q, double& quad, double& g) {
  double Q = 0, cross = 0;
  s = 0; q = 0;
  for (int i = 0; i < n; ++i) {
    s += zt[i] * zt[i]; q += zt[i] * wt[i]; Q += wt[i] * wt[i];
    for (int j = i + 1; j < n; ++j) { double c = zt[i] * wt[j] - zt[j] * wt[i]; cross += c * c; }
  }
  g = 1.0 + P * s;
  quad = (Q + P * cross) / g;
}

// On entry (a, P) are the predicted moments; on exit the filtered moments (left
// unchanged if the particle's density is not finite, which then gets weight 0).
static inline void sv_step(const SVSys& S, const double* yt, const double* h,
                           double& a, double& P, double& lp_full, double& lp_sub,
                           double* z, double* w, double* zt, double* wt) {
  const int K = S.K, nsub = S.nsub;
  double sumh = 0.0, sumh_sub = 0.0;
  for (int k = 0; k < K; ++k) {
    double uk = std::exp(-0.5 * h[k]);
    z[k] = S.Z[k] * uk;
    w[k] = (yt[k] - S.d[k] - a * S.Z[k]) * uk;
    sumh += h[k];
    if (k < nsub) sumh_sub += h[k];
  }
  for (int i = 0; i < K; ++i) {
    double zi = 0, wi = 0;
    for (int j = 0; j <= i; ++j) { zi += S.Linv[i * K + j] * z[j]; wi += S.Linv[i * K + j] * w[j]; }
    zt[i] = zi; wt[i] = wi;
  }
  double s, q, quad, g;
  sv_quad(zt, wt, K, P, s, q, quad, g);
  lp_full = -0.5 * (K * LOG2PI + sumh + S.logdetC + std::log(g) + quad);
  if (nsub < K) {
    double ss, qs, quads, gs;
    sv_quad(zt, wt, nsub, P, ss, qs, quads, gs);
    lp_sub = -0.5 * (nsub * LOG2PI + sumh_sub + S.logdetC_sub + std::log(gs) + quads);
  } else {
    lp_sub = lp_full;
  }
  if (!std::isfinite(lp_full) || !std::isfinite(lp_sub)) {
    lp_full = -INFINITY; lp_sub = -INFINITY;
    return;
  }
  a = a + P * q / g;
  P = P / g;
}

static SVSys make_sys(NumericVector Z, NumericVector d, NumericMatrix Linv, double logdetC,
                      double logdetC_sub, int nsub) {
  SVSys S;
  S.K = Z.size(); S.nsub = nsub;
  S.Z.assign(Z.begin(), Z.end()); S.d.assign(d.begin(), d.end());
  S.Linv.resize(S.K * S.K);
  for (int i = 0; i < S.K; ++i) for (int j = 0; j < S.K; ++j) S.Linv[i * S.K + j] = Linv(i, j);
  S.logdetC = logdetC; S.logdetC_sub = logdetC_sub;
  return S;
}

// [[Rcpp::export]]
NumericMatrix kf_given_h_cpp(NumericMatrix y, NumericMatrix h, NumericVector Z, NumericVector d,
                             NumericMatrix Linv, double logdetC, double logdetC_sub, int nsub,
                             double c, double rho, double tau) {
  SVSys S = make_sys(Z, d, Linv, logdetC, logdetC_sub, nsub);
  const int T = y.nrow(), K = y.ncol();
  NumericMatrix out(T, 2);
  double a = c / (1 - rho), P = tau * tau / (1 - rho * rho);
  std::vector<double> z(K), w(K), zt(K), wt(K), yt(K), ht(K);
  for (int t = 0; t < T; ++t) {
    for (int k = 0; k < K; ++k) { yt[k] = y(t, k); ht[k] = h(t, k); }
    double lf, ls;
    sv_step(S, yt.data(), ht.data(), a, P, lf, ls, z.data(), w.data(), zt.data(), wt.data());
    out(t, 0) = lf; out(t, 1) = ls;
    a = c + rho * a;
    P = rho * rho * P + tau * tau;
  }
  return out;
}

// log sum_i exp(x_i); -Inf if every term is -Inf.
static inline double lse_vec(const std::vector<double>& x, int n) {
  double mx = -INFINITY;
  for (int i = 0; i < n; ++i) if (x[i] > mx) mx = x[i];
  if (!std::isfinite(mx)) return -INFINITY;
  double s = 0;
  for (int i = 0; i < n; ++i) s += std::exp(x[i] - mx);
  return mx + std::log(s);
}

// Particle-count schedule: scalar or length T, every entry >= 1.
static inline int sched_max(const IntegerVector& Nv, int T) {
  if (Nv.size() != 1 && Nv.size() != T) stop("N must have length 1 or nrow(y)");
  int mx = 0;
  for (int i = 0; i < Nv.size(); ++i) { if (Nv[i] < 1) stop("N must be >= 1"); if (Nv[i] > mx) mx = Nv[i]; }
  return mx;
}
static inline int sched_at(const IntegerVector& Nv, int t) { return Nv.size() == 1 ? Nv[0] : Nv[t]; }
// Systematic resampling of N weighted particles (normalised weights W) to M offspring:
// writes the ancestor index of each offspring to anc[0..M-1].
static inline void sys_resample(const std::vector<double>& W, int N, int M, std::vector<int>& anc) {
  double u0 = R::unif_rand() / M, cum = W[0];
  int j = 0;
  for (int i = 0; i < M; ++i) {
    double ui = u0 + (double)i / M;
    while (ui > cum && j < N - 1) { ++j; cum += W[j]; }
    anc[i] = j;
  }
}

// h_{1,k} is drawn from N(mu_h_k, sp_k^2), sp_k = min(sd0_k, h1_sd_max), with the exact
// importance correction N(h; mu_h, sd0^2) / N(h; mu_h, sp^2) in the initial weights.
// This equals the bootstrap filter when sd0 <= h1_sd_max; it matters only for draws
// with phi_k very close to 1, whose stationary sd0 = omega/sqrt(1 - phi^2) is huge and
// would otherwise scatter the initial particles over absurd log-variances.
// Caveat: the correction is exact (the period-1 likelihood estimate stays unbiased), but
// the weights have finite variance only if sp^2 > sd0^2 / 2, i.e. sd0 < sqrt(2) * h1_sd_max
// (4.24 at the default 3). For draws with larger sd0 the estimate is heavy-tailed and in
// practice biased low at any feasible N, and seed-to-seed spread does not reveal it
// (toy check: sd0 = 21 gives about -0.66 nats at N = 20000). Only period 1 (1959 Q4,
// outside the LFO window) is affected directly; recomputing the m9 LFO with the plain
// bootstrap initialisation changed its total elpd by 0.14, within Monte Carlo error.
// [[Rcpp::export]]
NumericMatrix rbpf_cpp(NumericMatrix y, NumericVector Z, NumericVector d,
                       NumericMatrix Linv, double logdetC, double logdetC_sub, int nsub,
                       double c, double rho, double tau,
                       NumericVector mu_h, NumericVector phi, NumericVector omega,
                       IntegerVector Nv, double ess_thr, double h1_sd_max) {
  SVSys S = make_sys(Z, d, Linv, logdetC, logdetC_sub, nsub);
  const int T = y.nrow(), K = y.ncol();
  const int Nmax = sched_max(Nv, T);
  int N = sched_at(Nv, 0);
  NumericMatrix out(T, 4);
  std::vector<double> hp(Nmax * K), hp_new(Nmax * K), ap(Nmax), Pp(Nmax), ap_new(Nmax), Pp_new(Nmax);
  std::vector<double> logW(Nmax, -std::log((double)N)), lw(Nmax), lws(Nmax), Wn(Nmax);
  std::vector<int> anc(Nmax);
  std::vector<double> z(K), w(K), zt(K), wt(K), yt(K);
  const double a0 = c / (1 - rho), P0 = tau * tau / (1 - rho * rho);
  for (int i = 0; i < N; ++i) {
    ap[i] = a0; Pp[i] = P0;
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
    for (int k = 0; k < K; ++k) yt[k] = y(t, k);
    for (int i = 0; i < N; ++i) {
      double lf, ls;
      sv_step(S, yt.data(), &hp[i * K], ap[i], Pp[i], lf, ls, z.data(), w.data(), zt.data(), wt.data());
      lw[i] = logW[i] + lf;
      lws[i] = logW[i] + ls;
    }
    out(t, 0) = lse_vec(lw, N);  // log sum_i W_{t-1}^i p(y_t | y_{1:t-1}, h^i_{1:t})
    out(t, 1) = lse_vec(lws, N);
    double ess;
    if (!std::isfinite(out(t, 0))) {
      // every particle has zero density: the estimate is 0 (log = -Inf) for this draw;
      // restart from equal weights so the remaining periods are still defined.
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
    for (int i = 0; i < N; ++i) { ap[i] = c + rho * ap[i]; Pp[i] = rho * rho * Pp[i] + tau * tau; }
    if (t < T - 1) {
      const int Nn = sched_at(Nv, t + 1);
      if (ess < ess_thr * N || Nn != N) {           // systematic resampling (to Nn)
        sys_resample(Wn, N, Nn, anc);
        for (int i = 0; i < Nn; ++i) {
          const int j = anc[i];
          for (int k = 0; k < K; ++k) hp_new[i * K + k] = hp[j * K + k];
          ap_new[i] = ap[j]; Pp_new[i] = Pp[j];
        }
        hp.swap(hp_new); ap.swap(ap_new); Pp.swap(Pp_new);
        N = Nn;
        std::fill(logW.begin(), logW.begin() + N, -std::log((double)N));
      }
    }
  }
  return out;
}
