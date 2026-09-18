// varmapack_psi - compute VARMA impulse-response coefficient matrices

#include "varmapack.h"
#include "BlasGateway.h"
#include "error.h"
#include "VarmaPackUtil.h"
#include "VarmaUtilities.h"

bool varmapack_psi(double A[], double B[], int p, int q, int r,
                   int maxlag, double Psi[]) {
  int rr, maxAi;
  double *Psij;
  clear_error();
  if (p < 0 || q < 0 || r <= 0 || maxlag < 0 || Psi == 0) {
    return fail_error("invalid argument");
  }
  if ((p > 0 && A == 0) || (q > 0 && B == 0)) {
    return fail_error("invalid argument");
  }
  rr = r*r;
  setzero((maxlag + 1)*rr, Psi);
  setEye(r, Psi, r);
  for (int j=1; j<=maxlag; j++) {
    Psij = Psi + j*rr;
    if (j <= q) {
      lacpy("All", r, r, B + (j - 1)*rr, r, Psij, r);
    }
    maxAi = imin(p, j);
    for (int i=1; i<=maxAi; i++) {
      gemm("NoT", "NoT", r, r, r, 1, A + (i - 1)*rr, r, Psi + (j - i)*rr, r, 1, Psij, r);
    }
  }
  return true;
}

bool varmapack_irf(double A[], double B[], double Sig[], int p,
                   int q, int r, int maxlag, double Theta[]) {
  int rr;
  double *L = 0;
  double *W = 0;
  bool triangular;
  clear_error();
  if (Sig == 0 || Theta == 0) return fail_error("invalid argument");
  if (!varmapack_psi(A, B, p, q, r, maxlag, Theta)) return false;
  rr = r*r;
  if (!ALLOC(L, rr)) goto alloc_fail;
  if (!psdFactor(Sig, r, L, &triangular)) goto fail;
  if (!triangular && !ALLOC(W, rr)) goto alloc_fail;
  for (int j=0; j<=maxlag; j++) {
    if (triangular) {
      trmm("Right", "Low", "NoT", "NonUnit", r, r, 1, L, r, Theta + j*rr, r);
    }
    else {
      lacpy("All", r, r, Theta + j*rr, r, W, r);
      gemm("NoT", "NoT", r, r, r, 1, W, r, L, r, 0, Theta + j*rr, r);
    }
  }
  FREE(W);
  FREE(L);
  return true;
alloc_fail:
  fail_error("allocation failed");
fail:
  FREE(W);
  FREE(L);
  return false;
}
