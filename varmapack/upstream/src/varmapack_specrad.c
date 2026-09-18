// varmapack_specrad  Spectral radius of companion matrix of a VAR process
#include <math.h>
#include "BlasGateway.h"
#include "VarmaPackUtil.h"
#include "VarmaUtilities.h"
#include "error.h"

double varmapack_specrad(double *A, int r, int p) {
  clear_error();
  return CompanionSpecrad(A, r, p, 1);
}

HIDDEN double CompanionSpecrad(double P[], int r, int nblocks, double sign) {
  int i, j, lwork, info;
  int n;
  size_t pCount, acCount;
  double *Ac, *wr, *wi, *work, alwork, hyp, rho, dummy[1];
  if (nblocks < 0 || r <= 0 || (nblocks > 0 && P == 0)) {
    fail_error("invalid argument");
    return NAN;
  }
  if (nblocks == 0) return 0;
  if (!intProduct(r, nblocks, &n) ||
      !sizeProduct((size_t)r, (size_t)n, &pCount) ||
      !sizeProduct((size_t)n, (size_t)n, &acCount)) {
    fail_error("invalid dimension(s)");
    return NAN;
  }
  for (size_t k=0; k<pCount; k++) {
    if (!isfinite(P[k])) {
      fail_error("invalid argument");
      return NAN;
    }
  }
  Ac = wr = wi = work = 0;
  if (!ALLOC(Ac, acCount)) goto alloc_fail;
  if (!ALLOC(wr, n)) goto alloc_fail;
  if (!ALLOC(wi, n)) goto alloc_fail;
  lacpy("All", r, n, P, r, Ac, n);
  if (sign != 1) {
    for (j=0; j<n; j++) scal(r, sign, Ac + j*n, 1);
  }
  for (i = j = r; i < n; i++, j += n + 1) Ac[j] = 1;
  geev("N", "N", n, Ac, n, wr, wi, dummy, 1, dummy, 1, &alwork, -1, &info);
  if (info != 0) goto internal_fail;
  if (!isfinite(alwork) || alwork < 1 || alwork > INT_MAX) goto internal_fail;
  lwork = (int)ceil(alwork);
  if (!ALLOC(work, lwork)) goto alloc_fail;
  geev("N", "N", n, Ac, n, wr, wi, dummy, 1, dummy, 1, work, lwork, &info);
  if (info != 0) goto internal_fail;
  rho = 0;
  for (j=0; j<n; j++) {
    hyp = hypot(wr[j], wi[j]);
    rho = fmax(rho, hyp);
  }
  FREE(work);
  FREE(wi);
  FREE(wr);
  FREE(Ac);
  return rho;
alloc_fail:
  fail_error("allocation failed");
  goto fail;
internal_fail:
  fail_error("internal error");
fail:
  FREE(work);
  FREE(wi);
  FREE(wr);
  FREE(Ac);
  return NAN;
}
