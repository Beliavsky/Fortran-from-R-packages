// varmapack_sim — simulate burn-in-free AR, VAR, ARMA or VARMA time series
//
// See varmapack.h for parameter descriptions and storage conventions
//
// GENERAL VARMA TIME SERIES MODEL:
//   The following is called a VARMA(p,q) model:
//
//             x(t) - mu(t) = A1·(x(t-1) - mu(t-1)) + ... +
//                            Ap·(x(t-p) - mu(t-p)) + y(t)
//   where:
//             y(t) = eps(t) + B1·eps(t-1) + ... + Bq·eps(t-q),
//
//   x(t) is r-dimensional and eps(t) is r-variate normal with mean 0 and
//   covariance Sig. The eps(t) are called shocks and they are uncorrelated in
//   time (i.e. eps(t) and eps(s) are independent if t and s differ). VARMA
//   stands for vector-autoregressive-moving-average.
//
// SPECIAL CASES, ARMA, VAR and AR:
//   If r=1 the model above is an ARMA model, if q=0 and r>1 it is a VAR model
//   and if q=0 and r=1 it is an AR model. Thus VAR and AR models may be
//   written:
//            x(t) - mu(t) = A1·(x(t-1) - mu(t-1)) + ... +
//                           Ap·(x(t-p) - mu(t-p)) + eps(t)
//
//   where x(t) and eps(t) are scalar for an AR series and vectors for a VAR
//   series.
//
// ABOUT THE PROGRAM:
//   The program is based on the Matlab program varmapack_sim described in
//   references [1], [2] and [3] (the report [3] contains some details not found
//   in the published papers), but uses the C Randompack interface for random
//   number generation and may use either VYW or Lyapunov setup internally.
//
// RANDOM NUMBER GENERATION
//   The random number generator used for the shocks is specified by the rng
//   parameter. Random number generation is provided by the external Randompack
//   library. See GitHub.com/jonasson2/randompack, and [4] for details.
//
// MEANS
//   If mu is zero or nmu is zero, the time series has zero mean. Otherwise mu
//   gives r by nmu time series means: E[x(t)] = mu(:,t) for t < nmu, and
//   E[x(t)] = mu(:,nmu-1) for t >= nmu.
//
// STARTING VALUES:
//   When starting values are specified in X0, each generated series starts with
//   X0 and subsequent terms are simulated conditionally on X0. Without X0, the
//   first h x(t) and eps(t) values, h = max(p,q), are drawn from their
//   stationary joint distribution, so there is no initial burn-in segment to
//   throw away. With X0 and a stationary model, h = max(p,q,nX0), and the first
//   h eps(t) values are drawn from the conditional distribution of the startup
//   shocks given the supplied x values. This conditional path requires the
//   startup covariance of x to be positive definite. With X0 and a
//   nonstationary pure AR model, no latent startup shocks are needed; future
//   shocks are drawn and the recurrence is run forward from the supplied X0.
//   For a nonstationary model with MA terms, startup shocks are drawn
//   conditionally on the residual equations implied by X0. This requires
//   positive-definite Sig.
//
// MATRIX STORAGE SCHEME
//   All matrices are stored in column-major order, so that for example the
//   elements of A1, A2... are stored in the order A1(1,1),
//   A1(2,1),...,A1(r,r),..., Ap(r,r).
//
// STATIONARITY
//   Stationary models can be simulated without X0, including models with
//   singular positive-semidefinite startup covariance. Supplying X0 requires
//   this covariance to be positive definite. Nonstationary pure AR models can
//   be simulated with X0, because no latent startup shocks are needed for
//   forward recursion. Nonstationary models with MA terms require X0 and
//   positive-definite Sig.
//
// FAILURES
//   If the series is nonstationary without X0, or another recoverable error is
//   detected, the function returns false and sets varmapack_last_error().
//
// REFERENCES:
//   [1] K Jonasson and SE Ferrando 2008. Evaluating exact VARMA likelihood
//       and its gradient when data are incomplete. ACM Trans. Math. Softw.
//       35, No 1.
//   [2] K Jonasson 2008. Algorithm 878: Exact VARMA likelihood and its
//       gradient for complete and incomplete data with Matlab. ACM Trans.
//       Math. Softw. 35, No 1.
//   [3] K Jonasson and SE Ferrando 2006. Efficient likelihood evaluation for
//       VARMA processes with missing values. Report VHI-01-2006, Engineering
//       Research Institute, University of Iceland.
//   [4] K Jonasson 2026. Randompack: Cross-Platform Reproducible Random Number
//       Generation and Distribution Sampling, arXiv:2605.05099. Submitted to
//       ACM Trans. Math. Softw.

#include <math.h>
#include <stdbool.h>
#include "error.h"
#include "BlasGateway.h"
#include "VarmaUtilities.h"
#include "VarmaPackUtil.h"
#include "RandompackGateway.h"
#include "varmapack.h"

static void addMean(double X[], double mu[], int nmu, int r, int n, int M);
static void subtractMean(double X[], double mu[], int nmu, int r, int n, int M, int ldX);
static bool buildStartupCovar(double A[], double B[], double Sig[], int p, int q,
                              int r, int h, double C[], double SS[]);
static bool startUnconditional(double A[], double B[], double Sig[], double SS[],
                               double E[], int ldE, double X[], int p, int q,
                               int r, int h, int n, int M, randompack_rng *rng);
static bool startFromX0(double A[], double B[], double Sig[], double X0[], double mu[],
                        int nmu,
                        double C[], double SS[], double E[], int ldE, bool returnE,
                        bool stationary, double X[], int p, int q, int r, int h,
                        int n, int M, int MX0, randompack_rng *rng);
static void startupResidual(double A[], double X[], int p, int r, int h, int t0,
                            int rn, int M, double residual[]);
static bool tailSimulate(double Aflp[], double Bflp[], double Sig[], double E[], bool
                         rollingE, int ldE, double X[], int p, int q, int r,
                         int n, int M, int h, randompack_rng *rng);
static bool validSimSizes(int p, int q, int r, int n, int M, int nX0, bool returnE,
                          int *h, int *rn, int *rh, int *ldE, size_t *nM,
                          size_t *eCount, size_t *cCount, size_t *ssCount,
                          size_t *aCount, size_t *bCount);

bool varmapack_sim(double A[], double B[], double Sig[], double mu[], int nmu,
                   int p, int q, int r, int n, int M, double X0[], int nX0,
                   int MX0, double X[], double E[], randompack_rng *rng)
{
  double *C = 0, *SS = 0;
  double *Aflp = 0, *Bflp = 0;
  int h, rn, rh, ldE;
  size_t nM, eCount, cCount, ssCount, aCount, bCount;
  bool returnE = E != 0;
  bool X0Given = X0 != 0;
  bool rollingE = false;
  clear_error();
  if ((p > 0 && A == 0) || (q > 0 && B == 0) || Sig == 0 || X == 0 || rng == 0) {
    return fail_error("invalid argument");
  }
  if (p < 0 || q < 0 || r <= 0 || n <= 0 || M <= 0 || n < imax(p, q)) {
    return fail_error("invalid argument");
  }
  if (nmu < 0 || nmu > n || (mu == 0 && nmu > 0)) {
    return fail_error("invalid argument");
  }
  if ((!X0Given && nX0 != 0) ||
      (X0Given && (nX0 <= 0 || nX0 < imax(p, q) || nX0 > n)) ||
      (X0Given && MX0 != 1 && MX0 != M)) {
    return fail_error("invalid argument");
  }
  if (!validSimSizes(p, q, r, n, M, nX0, returnE, &h, &rn, &rh, &ldE, &nM,
                     &eCount, &cCount, &ssCount, &aCount, &bCount)) {
    return fail_error("problem size too large");
  }
  double rho = p == 0 ? 0 : varmapack_specrad(A, r, p);
  if (isnan(rho)) return false;
  bool stationary = rho < 1;
  if (!X0Given && !stationary)
    return fail_error("nonstationary model");
  if (!X0Given && p == 0 && q == 0) { // White noise can be drawn directly into X.
    if (!randompack_mvn("T", 0, Sig, r, nM, X, r, 0, rng)) goto random_fail;
    if (returnE) lacpy("All", rn, M, X, rn, E, rn); // and copy it to E
    addMean(X, mu, nmu, r, n, M);
    return true;
  }
  if (!returnE)
    if (!ALLOC(E, eCount)) goto alloc_fail;
  if (!X0Given || stationary) { // compute C and SS
    if (!ALLOC(C, cCount)) goto alloc_fail;
    if (rh > 0)
      if (!ALLOC(SS, ssCount)) goto alloc_fail;
    if (!buildStartupCovar(A, B, Sig, p, q, r, h, C, SS)) goto fail;
  }
  if (X0Given) {
    if (!startFromX0(A, B, Sig, X0, mu, nmu, C, SS, E, ldE, returnE, stationary, X,
                     p, q, r, h, n, M, MX0, rng)) goto fail;
  }
  else { // Unconditional exact startup: draw X[:h] given E[:h].
    if (!startUnconditional(A, B, Sig, SS, E, ldE, X, p, q, r, h, n, M, rng)) {
      goto fail;
    }
  }
  FREE(C); FREE(SS);
  if (p > 0 && !ALLOC(Aflp, aCount)) goto alloc_fail;
  if (q > 0 && !ALLOC(Bflp, bCount)) goto alloc_fail;
  flipmat(A, Aflp, r, p);
  flipmat(B, Bflp, r, q);
  rollingE = !returnE;
  if (!tailSimulate(Aflp, Bflp, Sig, E, rollingE, ldE, X, p, q, r, n, M, h, rng))
    goto fail;
  addMean(X, mu, nmu, r, n, M);
  FREE(Bflp); FREE(Aflp);
  if (!returnE) FREE(E);
  return true;
random_fail:
  fail_error(randompack_last_error(rng));
  goto fail;
alloc_fail:
  fail_error("allocation failed");
fail:
  FREE(Bflp); FREE(Aflp); FREE(SS); FREE(C);
  if (!returnE) FREE(E);
  return false;
}

static bool validSimSizes(int p, int q, int r, int n, int M, int nX0, bool returnE,
                          int *h, int *rn, int *rh, int *ldE, size_t *nM,
                          size_t *eCount, size_t *cCount, size_t *ssCount,
                          size_t *aCount, size_t *bCount) {
  int r2, cInt, aInt, bInt, square, rollingLdE, m, sColLd, sColCount;
  *h = imax(imax(p, q), nX0);
  if (*h == INT_MAX || !intProduct(r, n, rn) || !intProduct(r, *h, rh) ||
      !intProduct(r, *h + 1, &rollingLdE) || !intProduct(r, r, &r2) ||
      !intProduct(*rh, *rh, &square) || !intProduct(r2, q + 1, &cInt) ||
      !intProduct(r2, p, &aInt) || !intProduct(r2, q, &bInt)) return false;
  m = imax(*h, p + 1);
  if (!intProduct(r, m, &sColLd) || !intProduct(sColLd, r, &sColCount) ||
      !sizeProduct((size_t)n, (size_t)M, nM)) return false;
  *ldE = returnE ? *rn : rollingLdE;
  if (!sizeProduct((size_t)*ldE, (size_t)M, eCount)) return false;
  *cCount = cInt;
  *ssCount = square;
  *aCount = aInt;
  *bCount = bInt;
  return true;
}

static void addMean(double X[], double mu[], int nmu, int r, int n, int M) {
  size_t rn = (size_t)r*n;
  if (mu == 0 || nmu == 0) return;
  for (int j=0; j<M; j++) {
    for (int t=0; t<n; t++) {
      int imu = t < nmu ? t : nmu - 1;
      axpy(r, 1.0, mu + (size_t)imu*r, 1, X + (size_t)j*rn + (size_t)t*r, 1);
    }
  }
}

static void subtractMean(double X[], double mu[], int nmu, int r, int n, int M, int ldX) {
  if (mu == 0 || nmu == 0) return;
  for (int j=0; j<M; j++) {
    for (int t=0; t<n; t++) {
      int imu = t < nmu ? t : nmu - 1;
      axpy(r, -1.0, mu + (size_t)imu*r, 1,
           X + (size_t)j*ldX + (size_t)t*r, 1);
    }
  }
}

static bool buildStartupCovar(double A[], double B[], double Sig[], int p,
                              int q, int r, int h, double C[], double SS[]) {
  double *G = 0;
  double *S = 0;
  size_t rr, gCount, sCount;
  if (!sizeProduct((size_t)r, (size_t)r, &rr) ||
      !sizeProduct(rr, (size_t)(q + 1), &gCount) ||
      !sizeProduct(rr, (size_t)(p + 1), &sCount)) {
    return fail_error("problem size too large");
  }
  if (!ALLOC(G, gCount)) goto alloc_fail;
  if (!ALLOC(S, sCount)) goto alloc_fail;
  if (!FindS(A, B, Sig, p, q, r, S, C, G)) {
    FREE(S); FREE(G);
    return varmapack_last_error() ? false : fail_error("internal error");
  }
  if (!SBuild("Low", S, A, G, p, q, r, h, SS)) goto alloc_fail;
  FREE(S); FREE(G);
  return true;
alloc_fail:
  FREE(S); FREE(G);
  return fail_error("allocation failed");
}

static bool startUnconditional(double A[], double B[], double Sig[], double SS[],
                               double E[], int ldE, double X[], int p, int q, int r,
                               int h, int n, int M, randompack_rng *rng) {
  // Simulate M replicates X[:m] and E[:m] using joint distribution of (X,E)
  int rh = r*h;
  int rn = r*n;
  bool ok = true;
  double *Wrk = 0;
  double *Psi = 0;
  double *PsiHat = 0;
  double *R = 0;
  size_t wrkCount, squareCount;
  if (!sizeProduct((size_t)rh, (size_t)M, &wrkCount) ||
      !sizeProduct((size_t)rh, (size_t)rh, &squareCount)) {
    return fail_error("problem size too large");
  }
  if (!ALLOC(Wrk, wrkCount)) goto alloc_fail;
  if (!ALLOC(Psi, squareCount)) goto alloc_fail;
  if (!ALLOC(PsiHat, squareCount)) goto alloc_fail;
  if (!ALLOC(R, squareCount)) goto alloc_fail;
  for (int j=0; j<M; j++) {
    if (!randompack_mvn("T", 0, Sig, r, h, E + (size_t)j*ldE, r, 0, rng)) {
      ok = fail_error(randompack_last_error(rng));
      goto done;
    }
  }
  FindPsi(A, B, Psi, p, q, r);
  if (!FindPsiHat(Psi, PsiHat, Sig, r, h)) {
    ok = false;
    goto done;
  }
  lacpy("Low", rh, rh, SS, rh, R, rh);
  syrk("Low", "NoT", rh, rh, -1.0, PsiHat, rh, 1.0, R, rh);
  if (!randompack_mvn("T", 0, R, rh, M, Wrk, rh, 0, rng)) {
    ok = fail_error(randompack_last_error(rng));
    goto done;
  }
  lacpy("All", rh, M, Wrk, rh, X, rn);
  gemm("NoT", "NoT", rh, M, rh, 1.0, Psi, rh, E, ldE, 1.0, X, rn);
  goto done;
alloc_fail:
  ok = fail_error("allocation failed");
done:
  FREE(R); FREE(PsiHat); FREE(Psi); FREE(Wrk);
  return ok;
}

static bool startFromX0(double A[], double B[], double Sig[], double X0[], double mu[],
                        int nmu,
                        double C[], double SS[], double E[], int ldE, bool returnE,
                        bool stationary, double X[], int p, int q, int r, int h,
                        int n, int M, int MX0, randompack_rng *rng) {
  // Simulate M replicates from supplied X0.
  int info;
  int rh = r*h;
  int rn = r*n;
  bool ok = true;
  double *CC = 0;
  double *x0bar = 0;
  double *e = 0;
  double *wrk = 0;
  double *R = 0;
  double *Chat;
  size_t squareCount;
  if (!sizeProduct((size_t)rh, (size_t)rh, &squareCount)) {
    return fail_error("problem size too large");
  }
  for (int j=0; j<M; j++) {
    double *X0j = X0 + (MX0 == 1 ? 0 : (size_t)j*rh);
    copy(rh, X0j, 1, X + (size_t)j*rn, 1);
    subtractMean(X + (size_t)j*rn, mu, nmu, r, h, 1, rn);
  }
  if (!stationary && q == 0) {
    if (returnE) {
      for (int j=0; j<M; j++) setzero(rh, E + (size_t)j*ldE);
    }
    return true;
  }
  if (!stationary) {
    double *Eall = 0, *residual = 0;
    int firstActiveShock = p - q;
    int firstShockTime = imin(firstActiveShock, 0);
    int residualLength, ldEall;
    size_t residualCount, eallCount;
    if (firstShockTime < h - INT_MAX ||
        !intProduct(r, h - p, &residualLength) ||
        !intProduct(r, h - firstShockTime, &ldEall) ||
        !sizeProduct((size_t)residualLength, (size_t)M, &residualCount) ||
        !sizeProduct((size_t)ldEall, (size_t)M, &eallCount)) {
      return fail_error("problem size too large");
    }
    if (!requirePositiveDefinite(Sig, r)) return false;
    if (!ALLOC(Eall, eallCount)) goto nonstationary_alloc_fail;
    if (residualLength > 0 && !ALLOC(residual, residualCount)) {
      goto nonstationary_alloc_fail;
    }
    if (residualLength > 0)
      startupResidual(A, X, p, r, h, p, rn, M, residual);
    if (!drawConditionalStartupShocks(B, Sig, residual, q, r, h, p,
                                      firstShockTime, firstActiveShock, M,
                                      Eall, ldEall, rng)) goto nonstationary_fail;
    for (int j=0; j<M; j++) {
      copy(rh, Eall + (size_t)j*ldEall + (size_t)(-firstShockTime)*r, 1,
           E + (size_t)j*ldE, 1);
    }
    FREE(residual); FREE(Eall);
    return true;
nonstationary_alloc_fail:
    fail_error("allocation failed");
nonstationary_fail:
    FREE(residual); FREE(Eall);
    return false;
  }
  if (!ALLOC(CC, squareCount)) goto alloc_fail;
  if (!ALLOC(x0bar, rh)) goto alloc_fail;
  if (!ALLOC(e, rh)) goto alloc_fail;
  if (!ALLOC(wrk, rh)) goto alloc_fail;
  if (!ALLOC(R, squareCount)) goto alloc_fail;
  potrf("Low", rh, SS, rh, &info);
  if (info > 0) {
    ok = fail_error(
      "cannot condition on X0: startup covariance is not positive definite");
    goto fail;
  }
  if (info < 0) { ok = fail_error("internal error"); goto fail; }
  CCBuild(A, C, p, q, r, h, CC);
  Chat = CC;
  trsm("Left", "Low", "NT", "NotUD", rh, rh, 1.0, SS, rh, Chat, rh);
  for (int j=0; j<h; j++) {
    lacpy("Low", r, r, Sig, r, R + (size_t)j*r*(rh + 1), rh);
  }
  syrk("L", "T", rh, rh, -1.0, Chat, rh, 1.0, R, rh);
  for (int j=0; j<M; j++) {
    copy(rh, X + (size_t)j*rn, 1, x0bar, 1);
    copy(rh, x0bar, 1, wrk, 1);
    trsv("Lo", "NT", "NotUD", rh, SS, rh, wrk, 1);
    gemv("T", rh, rh, 1.0, Chat, rh, wrk, 1, 0.0, e, 1);
    if (!randompack_mvn("T", e, R, rh, 1, E + (size_t)j*ldE, ldE, 0, rng)) {
      ok = fail_error(randompack_last_error(rng));
      goto fail;
    }
  }
  FREE(R); FREE(wrk); FREE(e); FREE(x0bar); FREE(CC);
  return true;
alloc_fail:
  ok = fail_error("allocation failed");
fail:
  FREE(R); FREE(wrk); FREE(e); FREE(x0bar); FREE(CC);
  return ok;
}

static void startupResidual(double A[], double X[], int p, int r, int h, int t0,
                            int rn, int M, double residual[]) {
  int rm = r*(h - t0);
  for (int j=0; j<M; j++) {
    double *Xj = X + (size_t)j*rn;
    double *rj = residual + (size_t)j*rm;
    for (int t=t0; t<h; t++) {
      double *rt = rj + r*(t - t0);
      copy(r, Xj + r*t, 1, rt, 1);
      for (int i=1; i<=p; i++) {
        gemv("NoT", r, r, -1, A + (size_t)(i - 1)*r*r, r, Xj + r*(t - i),
             1, 1, rt, 1);
      }
    }
  }
}

static bool tailSimulate(double Aflp[], double Bflp[], double Sig[], double E[], bool
                         rollingE, int ldE, double X[], int p, int q, int r,
                         int n, int M, int h, randompack_rng *rng) {
  // Simulate X[h:] with forward computation from X[:h] and E[:h]
  int rn = r*n;
  double *LSig = 0;
  if (n > h && !ALLOC(LSig, (size_t)r*r)) {
    return fail_error("allocation failed");
  }
  for (int t=h; t<n; t++) {
    int iX = r*t;
    int iA = r*(t - p);
    int iB = r*(t - q);
    int iE = rollingE ? r*(t%(h + 1)) : iX;
    if (!randompack_mvn("T", 0, t == h ? Sig : 0, r, M, E + iE, ldE, LSig, rng)) {
      FREE(LSig);
      return fail_error(randompack_last_error(rng));
    }
    lacpy("All", r, M, E + iE, ldE, X + iX, rn);
    if (p > 0) // Apply AR
      gemm("NoT", "NoT", r, M, r*p, 1.0, Aflp, r, X + iA, rn, 1.0, X + iX, rn);
    if (q > 0 && rollingE) { // rolling MA to avoid allocating full E
      int firstSlot = (t - q)%(h + 1);
      int firstCount = h + 1 - firstSlot;
      if (firstCount > q) firstCount = q;
      gemm("NoT", "NoT", r, M, r*firstCount, 1.0, Bflp, r,
           E + r*firstSlot, ldE, 1.0, X + iX, rn);
      if (firstCount < q) {
        int secondCount = q - firstCount;
        gemm("NoT", "NoT", r, M, r*secondCount, 1.0,
             Bflp + (size_t)firstCount*r*r, r, E, ldE, 1.0, X + iX, rn);
      }
    }
    else if (q > 0) // Using full E
      gemm("NoT", "NoT", r, M, r*q, 1.0, Bflp, r, E + iB, rn, 1.0, X + iX, rn);
  }
  FREE(LSig);
  return true;
}
