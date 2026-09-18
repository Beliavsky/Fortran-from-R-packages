// varmapack_autocov — sample autocovariance of observed data

#include <stdbool.h>
#include "varmapack.h"
#include "error.h"
#include "BlasGateway.h"
#include "VarmaUtilities.h"

bool varmapack_autocov(const char *transp, const char *norm, int r,
                       int n, double X[], int maxlag, double C[]) {
  clear_error();
  if (transp == 0 || norm == 0 || X == 0 || C == 0) {
    return fail_error("invalid argument");
  }
  if (r <= 0 || n <= 0 || maxlag < 0 || maxlag > n-1) {
    return fail_error("invalid argument");
  }
  bool ML = (norm[0] == 'M' || norm[0] == 'm');
  bool CORRECTED = (norm[0] == 'C' || norm[0] == 'c');
  if (!ML && !CORRECTED) return fail_error("invalid argument");
  bool TRANSP = !(transp[0] == 'N' || transp[0] == 'n');
  if (!(!TRANSP || transp[0] == 'T' || transp[0] == 't')) {
    return fail_error("invalid argument");
  }
  double *Xc = 0;
  double *mu = 0;
  if (!ALLOC(Xc, r*n)) return fail_error("allocation failed");
  if (!ALLOC(mu, r)) {
    FREE(Xc);
    return fail_error("allocation failed");
  }
  copy(r*n, X, 1, Xc, 1);
  if (TRANSP) {
    meanmat("N", n, r, Xc, n, mu);
    for (int i=0; i<n; i++) {
      double *row = Xc + i;
      axpy(r, -1.0, mu, 1, row, n);
    }
  }
  else {
    meanmat("T", r, n, Xc, r, mu);
    for (int j=0; j<n; j++) {
      double *col = Xc + j*r;
      axpy(r, -1.0, mu, 1, col, 1);
    }
  }
  double *Xlag = Xc;
  for (int k=0; k<=maxlag; k++) {
    double fctr = ML ? 1.0/n : 1.0/(n - k);
    double *Ck = C + r*r*k;
    if (TRANSP) {
      double *Xnow = Xc + k;
      gemm("T", "N", r, r, n-k, fctr, Xnow, n, Xlag, n, 0.0, Ck, r);
    }
    else {
      double *Xnow = Xc + r*k;
      gemm("N", "T", r, r, n-k, fctr, Xnow, r, Xlag, r, 0.0, Ck, r);
    }
  }
  FREE(mu);
  FREE(Xc);
  return true;
}
