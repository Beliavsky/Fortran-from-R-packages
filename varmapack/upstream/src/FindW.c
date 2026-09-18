#include "BlasGateway.h"
#include "error.h"
#include "varmapack_config.h"
#include "VarmaPackUtil.h"

HIDDEN bool FindW(
  double B[],   // in   r×r×q, moving average parameter matrices
  double Sig[], // in   r×r, covariance of the shock terms eps(t)
  int q,        // in   number of moving average terms
  int r,        // in   dimension of each x(t)
  double W[])   // out  r×r×(q+1) with Wi = cov(y(t), y(t-i))
{
  int i, k, rr;
  size_t bSigCount;
  double *BSig = 0;
  double *Wk, *Bi, *BipkSig;
  if (r <= 0 || q < 0) return fail_error("invalid argument");
  if (q == INT_MAX || !intProduct(r, r, &rr) ||
      !sizeProduct((size_t)rr, (size_t)(q + 1), &bSigCount)) {
    return fail_error("problem size too large");
  }
  if (!ALLOC(BSig, bSigCount)) return fail_error("allocation failed");
  copy(rr, Sig, 1, BSig, 1);
  for (i=1; i<=q; i++) {
    Bi = B + (size_t)(i-1)*rr;
    gemm("NoT", "NoT", r, r, r, 1, Bi, r, Sig, r, 0, BSig + (size_t)i*rr, r);
  }
  for (k=0; k<=q; k++) {
    Wk = W + (size_t)k*rr;
    copy(rr, BSig + (size_t)k*rr, 1, Wk, 1);
    for (i=1; i<=q-k; i++) {
      Bi = B + (size_t)(i-1)*rr;
      BipkSig = BSig + (size_t)(i+k)*rr;
      gemm("NoT", "T", r, r, r, 1, BipkSig, r, Bi, r, 1, Wk, r);
    }
  }
  FREE(BSig);
  return true;
}
