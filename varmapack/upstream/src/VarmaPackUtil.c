#include "BlasGateway.h"
#include "Lyapunov.h"
#include "error.h"
#include "printX.h"
#include "varmapack.h"
#include "VarmaPackUtil.h"
#include "varmapack_config.h"
#include "VarmaUtilities.h"
#include "VYW.h"

static int slicotCutoff(int p, int q);
static void SExtend(double A[], double G[], double S[], double Scol[], int p, int q,
                    int r, int n);
static void buildStartupH(double B[], int q, int r, int h, int t0,
                          int firstShockTime, double H[]);

HIDDEN void FindC( // Calculate Ci = cov(x(t), eps(t-i))
  double A[],   // in   r×r×p, autoregressive parameter matrices
  double B[],   // in   r×r×q, moving average parameter matrices
  double Sig[], // in   r×r, covariance of the shock terms eps(t)
  int p,        // in   number of autoregressive terms
  int q,        // in   number of moving average terms
  int r,        // in   dimension of each x(t)
  double C[])   // out  r×r×(q+1) with C0, C1...Cq
{
  int i, j, rr = r*r;
  double *Bj, *Cj, *Ai, *Cjmi;
  copy(rr, Sig, 1, C, 1);
  for (j=1; j<=q; j++) {
    Cj = C + j*rr;
    Bj = B + (j-1)*rr;
    gemm("NoT", "NoT", r, r, r, 1, Bj, r, Sig, r, 0, Cj, r);
    for (i=1; i<=p && i<=j; i++) {
      Cjmi = C + (j-i)*rr;
      Ai = A + (i-1)*rr;
      gemm("NoT", "NoT", r, r, r, 1, Ai, r, Cjmi, r, 1, Cj, r);
    }
  }
}

HIDDEN void FindG( // Calculate Gi = cov(y(t), x(t-i))
  double B[],   // in   r×r×q, moving average parameter matrices
  double C[],   // in   r×r×(q+1), Ci = cov(x(t), eps(t-i))
  int q,        // in   number of moving average terms
  int r,        // in   dimension of each x(t)
  double G[])   // out  r×r×(q+1) with G0, G1...Gq
{
  int i, j, rr = r*r;
  double *Bi, *Cimj, *Gj;
  setzero(rr*(q+1), G);
  for (j=0; j<=q; j++) {
    Gj = G + j*rr;
    for (i=j; i<=q; i++) {
      Cimj = C + (i-j)*rr;
      if (i==0) {
        addmat("All", r, r, Cimj, r, Gj, r);
      }
      else {
        Bi = B + (i-1)*rr;
        gemm("NoT", "T", r, r, r, 1, Bi, r, Cimj, r, 1, Gj, r);
      }
    }
  }
}

HIDDEN bool FindS(double A[], double B[], double Sig[], int p, int q, int r,
                  double S[], double C[], double G[]) {
  int rc = slicotCutoff(p, q);
  if (rc > 0 && r >= rc) return LyapunovFactorizeSolve(A, B, Sig, p, q, r, S, C, G);
  return VYWFactorizeSolve(A, B, Sig, p, q, r, S, C, G);
}

HIDDEN bool SBuild(char *uplo, double S[], double A[], double G[], int p, int q, int r,
                   int n, double SS[]) {
  double *Scol = 0;
  double *SSj, *SSi;
  int j, m;
  if (n == 0) return true;
  m = imax(p+1, n);
  if (!ALLOC(Scol, (r*m)*r)) return false;
  SExtend(A, G, S, Scol, p, q, r, m);
  for (j = 0; j < n; j++) {
    SSj = SS + j*r*n*r + j*r;
    if (uplo[0] == 'A') {
      SSi = SSj + r*n*r;
      lacpy("All", r*(n-j), r, Scol, r*m, SSj, r*n);
      if (j < n-1) copytranspose(r*(n-j-1), r, Scol+r, r*m, SSi, r*n);
    }
    else {
      lacpy("Low", r*(n-j), r, Scol, r*m, SSj, r*n);
    }
  }
  FREE(Scol);
  return true;
}

HIDDEN void CCBuild( // Build covariance between terms and shocks
  double A[],  // in   r × r × p, A=[A1...Ap], autoregressive parameter matrices
  double C[],  // in   r × r·(q+1), [C0...Cq], Ci = cov(x(t),eps(t-i))
  int p,       // in   number of autoregressive terms
  int q,       // in   number of moving average terms
  int r,       // in   dimension of xt
  int n,       // in   length of series
  double CC[]) // out  r·n × r·n, cov(x1'...xn',eps1'...epsn')
{
  int i, j;
  double *CCj;
  setzero(r*n*r*n, CC);
  for (j=0; j<=q && j<n; j++) {
    lacpy("All", r, r, C+j*r*r, r, CC + j*r, r*n);
  }
  for (j=q+1; j<n; j++) {
    CCj = CC + j*r;
    for (i=1; i<=j && i<=p; i++) {
      gemm("N", "N", r, r, r, 1, A+(i-1)*r*r, r, CC+(j-i)*r, r*n, 1, CCj, r*n);
    }
  }
  for (j=1; j<n; j++) {
    CCj = CC + j*r*n*r + j*r;
    lacpy("All", r*(n-j), r, CC, r*n, CCj, r*n);
  }
}

static int slicotCutoff(int p, int q) {
  int pc, qc;
  // Average break-even r for choosing SLICOT over VYW; 0 means always use VYW.
  // The numbers are averages across Mac, XEON, Ubuntu; OpenBLAS, Accelerate, MKL.
  int rc[7][8] = {
    {16, 26, 0, 0, 0, 0, 0, 0}, {10, 13, 18, 21, 24, 24, 26, 28},
    {9, 13, 15, 17, 19, 21, 22, 23}, {11, 13, 14, 15, 16, 18, 18, 19},
    {11, 13, 14, 15, 15, 16, 17, 18}, {11, 13, 13, 14, 14, 15, 15, 16},
    {12, 13, 13, 13, 14, 14, 14, 15}
  };
  if (p <= 0 || q < 0) return 0;
  pc = p < 7 ? p : 7;
  qc = q < 7 ? q : 7;
  return rc[pc-1][qc];
}

static void SExtend(double A[], double G[], double S[], double Scol[], int p, int q,
                    int r, int n) {
  int iScol = n*r, i, j;
  double *Scolj, *Ai, *Scoli;
  for (j = 0; j < p+1; j++) {
    lacpy("All", r, r, S + j*r*r, r, Scol + j*r, iScol);
  }
  for (j = p+1; j < n; j++) {
    Scolj = Scol + j*r;
    if (j <= q) {
      lacpy("All", r, r, G + j*r*r, r, Scolj, iScol);
    }
    for (i = 0; i < p; i++) {
      Ai = A+i*r*r;
      Scoli = Scolj - (i+1)*r;
      gemm("N", "N", r, r, r, 1, Ai, r, Scoli, iScol, 1, Scolj, iScol);
    }
  }
}

HIDDEN void FindPsi(double *A, double *B, double *Psi, int p, int q, int r) {
  // Prepare
  int h = imax(p, q), rr = r*r, hr = h*r, i, j, k, l;
  double *Aflp = 0, *Psi_jj;
  if (h == 0) return;
  if (p > 0) xAssert(ALLOC(Aflp, p*rr));
  flipmat(A, Aflp, r, p);
  laset("All", hr, hr, 0.0, 0.0, Psi, hr);

  // Compute first block column
  setEye(r, Psi, hr);
  for (i=1; i<h; i++) {
    l = i*r;
    k = p*r - i*r;
    if (i <= p) {
      gemm("NoT", "NoT", r, r, l, 1.0, Aflp + k*r, r, Psi, hr, 1.0, Psi + i*r, hr);
    }
    if (i <= q) {
      addmat("All", r, r, B + (i - 1)*rr, r, Psi + i*r, hr);
    }
  }
  // Copy to remaining lower block columns
  for (j=1; j<h; j++) {
    Psi_jj = Psi + j*hr*r + j*r;
    lacpy("All", hr - j*r, r, Psi, hr, Psi_jj, hr);
  }
  FREE(Aflp);
}

HIDDEN bool FindPsiHat(double *Psi, double *Psi_hat, double *Sig, int r, int h) {
  double *LSig, *Psi_hat_kk, *W;
  int hr = h*r, rr = r*r, k, nrow;
  bool triangular;
  bool ok = false;
  LSig = 0;
  W = 0;
  if (!ALLOC(LSig, rr) || !ALLOC(W, hr*r)) {
    fail_error("allocation failed");
    goto done;
  }
  if (!psdFactor(Sig, r, LSig, &triangular)) goto done;
  lacpy("Low", hr, hr, Psi, hr, Psi_hat, hr);
  for (k=0; k<h; k++) {
    Psi_hat_kk = Psi_hat + k*h*rr + k*r;
    nrow = hr - k*r;
    if (triangular) {
      trmm("Right", "Low", "NoT", "NonUnit", nrow, r, 1, LSig, r, Psi_hat_kk, hr);
    }
    else {
      lacpy("All", nrow, r, Psi_hat_kk, hr, W, nrow);
      gemm("NoT", "NoT", nrow, r, r, 1, W, nrow, LSig, r, 0, Psi_hat_kk, hr);
    }
  }
  ok = true;
done:
  FREE(W);
  FREE(LSig);
  return ok;
}

HIDDEN bool psdFactor(double Sig[], int r, double L[], bool *triangular) {
  int info, rank, rr = r*r;
  int *piv = 0;
  double *C = 0;
  double *work = 0;
  double *recon = 0;
  double tol;
  bool ok = true;
  if (triangular) *triangular = false;
  for (int i=0; i<r; i++) {
    if (!isfinite(Sig[i + i*r])) return fail_error("matrix is not positive semidefinite");
  }
  if (!ALLOC(C, rr)) goto alloc_fail;
  if (!ALLOC(piv, r)) goto alloc_fail;
  if (!ALLOC(work, 2*r)) goto alloc_fail;
  if (!ALLOC(recon, rr)) goto alloc_fail;
  lacpy("All", r, r, Sig, r, C, r);
  potrf("Low", r, C, r, &info);
  if (info < 0) { ok = fail_error("internal error"); goto fail; }
  if (info == 0) {
    setzero(rr, L);
    lacpy("Low", r, r, C, r, L, r);
    if (triangular) *triangular = true;
    goto fail;
  }
  lacpy("All", r, r, Sig, r, C, r);
  pstrf("Low", r, C, r, piv, &rank, -1, work, &info);
  if (info < 0) { ok = fail_error("internal error"); goto fail; }
  setzero(rr, L);
  for (int k=0; k<rank; k++) {
    for (int i=k; i<r; i++) {
      L[piv[i] + k*r] = C[i + k*r];
    }
  }
  setzero(rr, recon);
  if (rank > 0) syrk("Low", "NoT", r, rank, 1, L, r, 0, recon, r);
  copylowertoupper(r, recon, r);
  tol = 100*r*lamch("E");
  if (relabsdiff(Sig, recon, rr) > tol) {
    ok = fail_error("matrix is not positive semidefinite");
  }
fail:
  FREE(recon);
  FREE(work);
  FREE(piv);
  FREE(C);
  return ok;
alloc_fail:
  ok = fail_error("allocation failed");
  goto fail;
}

HIDDEN bool drawConditionalStartupShocks(
  double B[], double Sig[], double residual[], int q, int r, int h, int t0,
  int firstShockTime, int firstActiveShock, int M, double E[], int ldE,
  randompack_rng *rng) {
  int info, m = h - t0, nPrefix = imax(firstActiveShock, 0);
  int activeOffset = r*(firstActiveShock - firstShockTime);
  int nActive = h - firstActiveShock, rm = r*m, re = r*nActive;
  bool ok = true;
  double *H = 0, *HD = 0, *Wlag = 0, *W = 0, *K = 0, *R = 0, *work = 0;
  double *Ehat = 0, *residualWork = 0;
  size_t hCount, wlagCount, wCount, rCount, workCount;
  if (m == 0) {
    int nStartup = h - firstShockTime;
    for (int j=0; j<M; j++) {
      if (nStartup > 0 && !randompack_mvn("T", 0, Sig, r, nStartup,
          E + (size_t)j*ldE, r, 0, rng)) return fail_error(randompack_last_error(rng));
    }
    return true;
  }
  for (int j=0; j<M; j++) {
    if (nPrefix > 0 && !randompack_mvn("T", 0, Sig, r, nPrefix,
        E + (size_t)j*ldE, r, 0, rng)) {
      ok = fail_error(randompack_last_error(rng));
      goto done;
    }
  }
  if (q == 0) {
    for (int j=0; j<M; j++) {
      copy(rm, residual + (size_t)j*rm, 1,
           E + (size_t)j*ldE + activeOffset, 1);
    }
    return true;
  }
  if (!sizeProduct((size_t)rm, (size_t)re, &hCount) ||
      !sizeProduct((size_t)r, (size_t)r, &wlagCount) ||
      !sizeProduct(wlagCount, (size_t)(q + 1), &wlagCount) ||
      !sizeProduct((size_t)rm, (size_t)rm, &wCount) ||
      !sizeProduct((size_t)re, (size_t)re, &rCount) ||
      !sizeProduct((size_t)imax(rm, re), (size_t)r, &workCount)) {
    return fail_error("problem size too large");
  }
  if (!ALLOC(H, hCount)) goto alloc_fail;
  if (!ALLOC(HD, hCount)) goto alloc_fail;
  if (!ALLOC(Wlag, wlagCount)) goto alloc_fail;
  if (!ALLOC(W, wCount)) goto alloc_fail;
  if (!ALLOC(K, hCount)) goto alloc_fail;
  if (!ALLOC(R, rCount)) goto alloc_fail;
  if (!ALLOC(work, workCount)) goto alloc_fail;
  if (!ALLOC(Ehat, re)) goto alloc_fail;
  if (!ALLOC(residualWork, rm)) goto alloc_fail;
  buildStartupH(B, q, r, h, t0, firstActiveShock, H);
  lacpy("All", rm, re, H, rm, HD, rm);
  postmultiplySigmaPrime(HD, rm, rm, nActive, Sig, r, work);
  if (!FindW(B, Sig, q, r, Wlag)) {
    ok = varmapack_last_error() ? false : fail_error("allocation failed");
    goto done;
  }
  WBuild(Wlag, q, r, m, W);
  potrf("Low", rm, W, rm, &info);
  if (info > 0) { ok = fail_error("singular matrix"); goto done; }
  if (info < 0) { ok = fail_error("internal error"); goto done; }
  lacpy("All", rm, re, HD, rm, K, rm);
  trsm("Left", "Low", "NoT", "NonUnit", rm, re, 1, W, rm, K, rm);
  trsm("Left", "Low", "Trans", "NonUnit", rm, re, 1, W, rm, K, rm);
  laset("All", re, re, 0, 0, R, re);
  for (int j=0; j<nActive; j++)
    lacpy("All", r, r, Sig, r, R + (size_t)j*r*((size_t)re + 1), re);
  gemm("Trans", "NoT", re, re, rm, -1, HD, rm, K, rm, 1, R, re);
  for (int j=0; j<M; j++) {
    double *rj = residual + (size_t)j*rm;
    copy(rm, rj, 1, residualWork, 1);
    trsv("Low", "NoT", "NonUnit", rm, W, rm, residualWork, 1);
    trsv("Low", "Trans", "NonUnit", rm, W, rm, residualWork, 1);
    gemv("Trans", rm, re, 1, HD, rm, residualWork, 1, 0, Ehat, 1);
    if (!randompack_mvn("T", Ehat, R, re, 1,
        E + (size_t)j*ldE + activeOffset, ldE, 0, rng)) {
      ok = fail_error(randompack_last_error(rng));
      goto done;
    }
    for (int t=t0; t<h; t++) {
      double *Et = E + (size_t)j*ldE + r*(t - firstShockTime);
      copy(r, rj + r*(t - t0), 1, Et, 1);
      for (int i=1; i<=q; i++) {
        gemv("NoT", r, r, -1, B + (size_t)(i - 1)*r*r, r,
             E + (size_t)j*ldE + r*(t - i - firstShockTime), 1, 1, Et, 1);
      }
    }
  }
  goto done;
alloc_fail:
  ok = fail_error("allocation failed");
done:
  FREE(residualWork); FREE(Ehat); FREE(work); FREE(R); FREE(K); FREE(W); FREE(Wlag);
  FREE(HD); FREE(H);
  return ok;
}

HIDDEN bool requirePositiveDefinite(double Sig[], int r) {
  double *L = 0;
  bool triangular;
  size_t count;
  if (!sizeProduct((size_t)r, (size_t)r, &count))
    return fail_error("problem size too large");
  if (!ALLOC(L, count)) return fail_error("allocation failed");
  if (!psdFactor(Sig, r, L, &triangular)) {
    FREE(L);
    return false;
  }
  FREE(L);
  if (!triangular) {
    return fail_error(
      "cannot condition on X0: shock covariance is not positive definite");
  }
  return true;
}

static void buildStartupH(double B[], int q, int r, int h, int t0,
                          int firstShockTime, double H[]) {
  int m = h - t0, nE0 = h - firstShockTime, rm = r*m;
  laset("All", rm, r*nE0, 0, 0, H, rm);
  for (int t=t0; t<h; t++) {
    int row = r*(t - t0);
    for (int ell=firstShockTime; ell<h; ell++) {
      int lag = t - ell, col = r*(ell - firstShockTime);
      double *Hij = H + (size_t)col*rm + row;
      if (lag == 0) {
        for (int k=0; k<r; k++) Hij[k + (size_t)k*rm] = 1;
      }
      else if (1 <= lag && lag <= q) {
        lacpy("All", r, r, B + (size_t)(lag - 1)*r*r, r, Hij, rm);
      }
    }
  }
}
