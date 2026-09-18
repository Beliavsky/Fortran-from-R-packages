// varmapack_simx — simulate VARMAX time series with fixed exogenous data
//
// See varmapack.h for parameter descriptions and storage conventions.

#include "BlasGateway.h"
#include "VarmaPackUtil.h"
#include "VarmaUtilities.h"
#include "error.h"
#include "RandompackGateway.h"
#include "varmapack.h"
#include "varmapack_config.h"

static void startupResidual(double A[], double C[], double z[], double X0[], int p,
                            int s, int d, int r, int h, int t0, double residual[]);
static bool forwardSimx(double Aflp[], double Bflp[], double C[], double Sig[],
                        double z[], int Mz, double Eall[], int ldE, double X[], int p,
                        int q, int s, int d, int r, int n, int M, int h, int rn,
                        int firstShockTime, randompack_rng *rng);
static bool validSimxSizes(int p, int q, int s, int d, int r, int n, int M, int h,
                           int *t0, int *firstActiveShock, int *firstShockTime,
                           int *rh, int *rn, int *ldE, size_t *eCount,
                           size_t *aCount, size_t *bCount);

bool varmapack_simx(double A[], double B[], double C[], double Sig[], double z[],
                    int Mz, int p, int q, int s, int d, int r, int n, int M, double X0[],
                    int h, int MX0, double X[], double E[], randompack_rng *rng)
{
  double *Eall = 0, *Aflp = 0, *Bflp = 0, *residual = 0;
  int t0, firstActiveShock, firstShockTime, residualLength, rh, ldE, rn;
  size_t eCount, aCount, bCount, residualCount;
  clear_error();
  if ((p > 0 && A == 0) || (q > 0 && B == 0) || (s > 0 && C == 0) ||
      (s > 0 && z == 0) || Sig == 0 || (h > 0 && X0 == 0) || X == 0 || rng == 0) {
    return fail_error("invalid argument");
  }
  if (p < 0 || q < 0 || s < 0 || d < 0 || (s == 0 && d != 0) ||
      (s > 0 && d == 0) || r <= 0 || n <= 0 || M <= 0 ||
      (Mz != 1 && Mz != M) || (MX0 != 1 && MX0 != M)) {
    return fail_error("invalid argument");
  }
  if (h < imax(imax(p, q), s > 0 ? s - 1 : 0) || n < h) {
    return fail_error("invalid argument");
  }
  if (!validSimxSizes(p, q, s, d, r, n, M, h, &t0, &firstActiveShock,
                      &firstShockTime, &rh, &rn, &ldE, &eCount, &aCount, &bCount)) {
    return fail_error("problem size too large");
  }
  if (!intProduct(r, h - t0, &residualLength) ||
      !sizeProduct((size_t)residualLength, (size_t)M, &residualCount)) {
    return fail_error("problem size too large");
  }
  if (!requirePositiveDefinite(Sig, r)) return false;
  if (!ALLOC(Eall, eCount)) goto alloc_fail;
  if (p > 0 && !ALLOC(Aflp, aCount)) goto alloc_fail;
  if (q > 0 && !ALLOC(Bflp, bCount)) goto alloc_fail;
  if (residualLength > 0 && !ALLOC(residual, residualCount)) goto alloc_fail;
  flipmat(A, Aflp, r, p);
  flipmat(B, Bflp, r, q);
  if (h > 0) {
    for (int j=0; j<M; j++) {
      double *X0j = X0 + (MX0 == 1 ? 0 : (size_t)j*rh);
      copy(rh, X0j, 1, X + (size_t)j*rn, 1);
    }
  }
  for (int j=0; j<M; j++) {
    double *zj = s > 0 ? z + (Mz == 1 ? 0 : (size_t)j*d*n) : 0;
    double *X0j = h > 0 ? X0 + (MX0 == 1 ? 0 : (size_t)j*rh) : 0;
    startupResidual(A, C, zj, X0j, p, s, d, r, h, t0,
                    residualLength > 0 ? residual + (size_t)j*residualLength : 0);
  }
  if (!drawConditionalStartupShocks(B, Sig, residual, q, r, h, t0,
                                    firstShockTime, firstActiveShock, M, Eall,
                                    ldE, rng)) {
    goto fail;
  }
  if (!forwardSimx(Aflp, Bflp, C, Sig, z, Mz, Eall, ldE, X, p, q, s, d, r, n, M, h,
                   rn, firstShockTime, rng)) {
    goto fail;
  }
  if (E != 0) {
    for (int j=0; j<M; j++) {
      lacpy("All", rn, 1, Eall + (size_t)j*ldE + (size_t)(-firstShockTime)*r, ldE,
            E + (size_t)j*rn, rn);
    }
  }
  FREE(residual); FREE(Bflp); FREE(Aflp); FREE(Eall);
  return true;
alloc_fail:
  fail_error("allocation failed");
fail:
  FREE(residual); FREE(Bflp); FREE(Aflp); FREE(Eall);
  return false;
}

static bool validSimxSizes(int p, int q, int s, int d, int r, int n, int M, int h,
                           int *t0, int *firstActiveShock, int *firstShockTime,
                           int *rh, int *rn, int *ldE, size_t *eCount,
                           size_t *aCount, size_t *bCount) {
  int r2, rp, rq, m, nActive, rm, re, nE;
  size_t hCount, wlagCount, wCount, rCount, workCount, cCount, zCount;
  *t0 = imax(p, s > 0 ? s - 1 : 0);
  if (q == INT_MAX) return false;
  *firstActiveShock = *t0 - q;
  *firstShockTime = imin(*firstActiveShock, 0);
  if (*firstShockTime < n - INT_MAX) return false;
  nE = n - *firstShockTime;
  m = h - *t0;
  nActive = h - *firstActiveShock;
  if (!intProduct(r, h, rh) || !intProduct(r, n, rn) || !intProduct(r, nE, ldE) ||
      !intProduct(r, r, &r2) || !intProduct(r, p, &rp) || !intProduct(r, q, &rq) ||
      !intProduct(r, m, &rm) || !intProduct(r, nActive, &re)) return false;
  if (!sizeProduct((size_t)*ldE, (size_t)M, eCount) ||
      !sizeProduct((size_t)r2, (size_t)p, aCount) ||
      !sizeProduct((size_t)r2, (size_t)q, bCount) ||
      !sizeProduct((size_t)r, (size_t)d, &cCount) ||
      !sizeProduct(cCount, (size_t)s, &cCount) ||
      !sizeProduct((size_t)d, (size_t)n, &zCount) ||
      !sizeProduct(zCount, (size_t)M, &zCount)) return false;
  if (q > 0 && m > 0 &&
      (!sizeProduct((size_t)rm, (size_t)re, &hCount) ||
       !sizeProduct((size_t)r2, (size_t)(q + 1), &wlagCount) ||
       !sizeProduct((size_t)rm, (size_t)rm, &wCount) ||
       !sizeProduct((size_t)re, (size_t)re, &rCount) ||
       !sizeProduct((size_t)imax(rm, re), (size_t)r, &workCount))) return false;
  return true;
}

static void startupResidual(double A[], double C[], double z[], double X0[], int p,
                            int s, int d, int r, int h, int t0, double residual[]) {
  for (int t=t0; t<h; t++) {
    double *rt = residual + r*(t - t0);
    copy(r, X0 + r*t, 1, rt, 1);
    for (int i=1; i<=p; i++) {
      gemv("NoT", r, r, -1, A + (size_t)(i-1)*r*r, r, X0 + r*(t-i), 1, 1, rt, 1);
    }
    for (int i=1; i<=s; i++) {
      gemv("NoT", r, d, -1, C + (size_t)(i-1)*r*d, r, z + (size_t)(t-i+1)*d, 1,
           1, rt, 1);
    }
  }
}

static bool forwardSimx(double Aflp[], double Bflp[], double C[], double Sig[],
                        double z[], int Mz, double Eall[], int ldE, double X[], int p,
                        int q, int s, int d, int r, int n, int M, int h, int rn,
                        int firstShockTime, randompack_rng *rng) {
  if (n > h) {
    for (int t=h; t<n; t++) {
      int iE = r*(t - firstShockTime);
      if (!randompack_mvn("T", 0, Sig, r, M, Eall + iE, ldE, 0, rng)) {
        return fail_error(randompack_last_error(rng));
      }
    }
  }
  for (int t=h; t<n; t++) {
    int iX = r*t;
    int iA = r*(t - p);
    int iB = r*(t - q - firstShockTime);
    lacpy("All", r, M, Eall + r*(t - firstShockTime), ldE, X + iX, rn);
    if (p > 0) {
      gemm("NoT", "NoT", r, M, r*p, 1, Aflp, r, X + iA, rn, 1, X + iX, rn);
    }
    if (q > 0) {
      gemm("NoT", "NoT", r, M, r*q, 1, Bflp, r, Eall + iB, ldE, 1, X + iX, rn);
    }
    for (int i=1; i<=s; i++) {
      for (int j=0; j<M; j++) {
        double *zj = z + (Mz == 1 ? 0 : (size_t)j*d*n);
        gemv("NoT", r, d, 1, C + (size_t)(i-1)*r*d, r, zj + (size_t)(t-i+1)*d, 1,
             1, X + (size_t)j*rn + iX, 1);
      }
    }
  }
  return true;
}
