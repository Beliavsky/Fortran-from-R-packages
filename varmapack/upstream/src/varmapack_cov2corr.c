// Convert covariance matrices at successive lags to correlation matrices.

#include <math.h>
#include <stdbool.h>
#include "varmapack.h"
#include "error.h"

bool varmapack_cov2corr(double cov[], int r, int maxlag, double corr[]) {
  double *sd = 0;
  clear_error();
  if (cov == 0 || corr == 0 || r <= 0 || maxlag < 0) {
    return fail_error("invalid argument");
  }
  int r2, n;
  if (maxlag == INT_MAX || !intProduct(r, r, &r2) ||
      !intProduct(r2, maxlag + 1, &n)) {
    return fail_error("problem size too large");
  }
  for (int i=0; i<n; i++) {
    if (!isfinite(cov[i])) return fail_error("nonfinite covariance");
  }
  for (int i=0; i<r; i++) {
    if (cov[i + i*r] <= 0) return fail_error("nonpositive variance");
  }
  if (!ALLOC(sd, r)) return fail_error("allocation failed");
  for (int i=0; i<r; i++) sd[i] = sqrt(cov[i + i*r]);
  for (int k=0; k<=maxlag; k++) {
    for (int j=0; j<r; j++) {
      for (int i=0; i<r; i++) {
        int index = i + j*r + k*r2;
        corr[index] = cov[index]/sd[i]/sd[j];
      }
    }
  }
  // Set only the lag-zero self-correlations exactly to one.
  for (int i=0; i<r; i++) corr[i + i*r] = 1;
  FREE(sd);
  return true;
}
