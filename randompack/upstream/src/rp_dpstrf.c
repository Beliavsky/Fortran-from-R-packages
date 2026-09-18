#include <math.h>
#include <stdbool.h>
#include "randompack_config.h"
#include "BlasGateway.h"

// isnan portable stand-in (matches Fortran DISNAN)
#define disnan(x) isnan(x)

// Column-major element access: A(i,j) in Fortran -> A[i + lda*j] in C (0-indexed)
#define A(i,j) a[(i) + lda*(j)]
#define WORK(i) work[(i)]

// ----------------------------------------------------------------------------
// pstf2 -- unblocked pivoted Cholesky factorization
//
// Computes the Cholesky factorization with complete pivoting of a real
// symmetric positive semidefinite matrix A:
//   P^T * A * P = U^T * U  (uplo='U')
//   P^T * A * P = L * L^T  (uplo='L')
//
// Arguments (all 0-indexed on entry/exit):
//   uplo  - 'U': upper triangle stored; 'L': lower triangle stored
//   n     - order of A
//   a     - matrix of dimension (lda, n); overwritten with factor
//   lda   - leading dimension of a (>= n)
//   piv   - pivot indices (1-indexed): column/row j was interchanged with piv[j]-1
//   rank  - estimated rank of A on output
//   tol   - user tolerance; if < 0 use n*eps*maxi(diag(A))
//   work  - workspace of length 2*n
//   info  - 0: success; 1: matrix not positive definite (rank < n)
// ----------------------------------------------------------------------------
static void rp_dpstf2(char *uplo, int n, double *a, int lda, int *piv,
  int *rank, double tol, double *work, int *info)
{
  *info = 0;
  bool upper = (uplo[0] == 'U' || uplo[0] == 'u');
  if (n == 0)
    return;
  for (int i = 0; i < n; i++)
    piv[i] = i + 1;
  int pvt = 0;
  double ajj = A(0, 0);
  for (int i = 1; i < n; i++) {
    if (A(i, i) > ajj) {
      pvt = i;
      ajj = A(pvt, pvt);
    }
  }
  if (ajj <= 0.0 || disnan(ajj)) {
    *rank = 0;
    *info = 1;
    return;
  }
  double dstop = (tol < 0.0) ? n * lamch("E") * ajj : tol;
  for (int i = 0; i < n; i++)
    WORK(i) = 0.0;
  if (upper) {
    for (int j = 0; j < n; j++) {
      for (int i = j; i < n; i++) {
        if (j > 0)
          WORK(i) += A(j-1, i) * A(j-1, i);
        WORK(n+i) = A(i, i) - WORK(i);
      }
      if (j > 0) {
        pvt = j;
        ajj = WORK(n+j);
        for (int i = j+1; i < n; i++) {
          if (WORK(n+i) > ajj) {
            pvt = i;
            ajj = WORK(n+i);
          }
        }
        if (ajj <= dstop || disnan(ajj)) {
          A(j, j) = ajj;
          *rank = j;
          *info = 1;
          return;
        }
      }
      if (j != pvt) {
        A(pvt, pvt) = A(j, j);
        swap(j,       &A(0, j),     1,   &A(0, pvt),     1);
        if (pvt < n-1)
          swap(n-pvt-1, &A(j, pvt+1), lda, &A(pvt, pvt+1), lda);
        swap(pvt-j-1, &A(j, j+1),   lda, &A(j+1, pvt),   1);
        double dtemp = WORK(j);
        WORK(j)   = WORK(pvt);
        WORK(pvt) = dtemp;
        int itemp = piv[pvt];
        piv[pvt]  = piv[j];
        piv[j]    = itemp;
      }
      ajj = sqrt(ajj);
      A(j, j) = ajj;
      if (j < n-1) {
        gemv("Trans", j, n-j-1, -1.0, &A(0, j+1), lda,
             &A(0, j), 1, 1.0, &A(j, j+1), lda);
        scal(n-j-1, 1.0/ajj, &A(j, j+1), lda);
      }
    }
  }
  else {
    for (int j = 0; j < n; j++) {
      for (int i = j; i < n; i++) {
        if (j > 0)
          WORK(i) += A(i, j-1) * A(i, j-1);
        WORK(n+i) = A(i, i) - WORK(i);
      }
      if (j > 0) {
        pvt = j;
        ajj = WORK(n+j);
        for (int i = j+1; i < n; i++) {
          if (WORK(n+i) > ajj) {
            pvt = i;
            ajj = WORK(n+i);
          }
        }
        if (ajj <= dstop || disnan(ajj)) {
          A(j, j) = ajj;
          *rank = j;
          *info = 1;
          return;
        }
      }
      if (j != pvt) {
        A(pvt, pvt) = A(j, j);
        swap(j,       &A(j,     0), lda, &A(pvt,   0),   lda);
        if (pvt < n-1)
          swap(n-pvt-1, &A(pvt+1, j), 1,   &A(pvt+1, pvt), 1);
        swap(pvt-j-1, &A(j+1,   j), 1,   &A(pvt,   j+1), lda);
        double dtemp = WORK(j);
        WORK(j)   = WORK(pvt);
        WORK(pvt) = dtemp;
        int itemp = piv[pvt];
        piv[pvt]  = piv[j];
        piv[j]    = itemp;
      }
      ajj = sqrt(ajj);
      A(j, j) = ajj;
      if (j < n-1) {
        gemv("No Trans", n-j-1, j, -1.0, &A(j+1, 0), lda,
             &A(j, 0), lda, 1.0, &A(j+1, j), 1);
        scal(n-j-1, 1.0/ajj, &A(j+1, j), 1);
      }
    }
  }
  *rank = n;
}

// ----------------------------------------------------------------------------
// pstrf -- blocked pivoted Cholesky factorization
//
// Same interface as pstf2 but uses level-3 BLAS (dsyrk) for the trailing
// submatrix update, falling back to pstf2 when the block size NB >= N.
//
// Additional workspace: work must be length 2*n.
// NB is obtained from ilaenv (same query as DPOTRF uses).
// ----------------------------------------------------------------------------
HIDDEN void rp_dpstrf(char *uplo, int n, double *a, int lda, int *piv,
  int *rank, double tol, double *work, int *info)
{
  *info = 0;
  bool upper = (uplo[0] == 'U' || uplo[0] == 'u');
  if (n == 0)
    return;
  int nb = ilaenv(1, "DPOTRF", uplo, n, -1, -1, -1);
  if (nb <= 1 || nb >= n) {
    rp_dpstf2(uplo, n, a, lda, piv, rank, tol, work, info);
    return;
  }
  for (int i = 0; i < n; i++)
    piv[i] = i + 1;
  int pvt = 0;
  double ajj = A(0, 0);
  for (int i = 1; i < n; i++) {
    if (A(i, i) > ajj) {
      pvt = i;
      ajj = A(pvt, pvt);
    }
  }
  if (ajj <= 0.0 || disnan(ajj)) {
    *rank = 0;
    *info = 1;
    return;
  }
  double dstop = (tol < 0.0) ? n * lamch("E") * ajj : tol;
  if (upper) {
    for (int k = 0; k < n; k += nb) {
      int jb = (nb < n-k) ? nb : n-k;
      for (int i = k; i < n; i++)
        WORK(i) = 0.0;
      int j;
      for (j = k; j < k+jb; j++) {
        for (int i = j; i < n; i++) {
          if (j > k)
            WORK(i) += A(j-1, i) * A(j-1, i);
          WORK(n+i) = A(i, i) - WORK(i);
        }
        if (j > 0) {
          pvt = j;
          ajj = WORK(n+j);
          for (int i = j+1; i < n; i++) {
            if (WORK(n+i) > ajj) {
              pvt = i;
              ajj = WORK(n+i);
            }
          }
          if (ajj <= dstop || disnan(ajj)) {
            A(j, j) = ajj;
            *rank = j;
            *info = 1;
            return;
          }
        }
        if (j != pvt) {
          A(pvt, pvt) = A(j, j);
          swap(j,       &A(0, j),     1,   &A(0, pvt),     1);
          if (pvt < n-1)
            swap(n-pvt-1, &A(j, pvt+1), lda, &A(pvt, pvt+1), lda);
          swap(pvt-j-1, &A(j, j+1),   lda, &A(j+1, pvt),   1);
          double dtemp = WORK(j);
          WORK(j)   = WORK(pvt);
          WORK(pvt) = dtemp;
          int itemp = piv[pvt];
          piv[pvt]  = piv[j];
          piv[j]    = itemp;
        }
        ajj = sqrt(ajj);
        A(j, j) = ajj;
        if (j < n-1) {
          gemv("Trans", j-k, n-j-1, -1.0, &A(k, j+1), lda,
               &A(k, j), 1, 1.0, &A(j, j+1), lda);
          scal(n-j-1, 1.0/ajj, &A(j, j+1), lda);
        }
      }
      // Trailing submatrix update: A(j:n, j:n) -= A(k:j, j:n)^T * A(k:j, j:n)
      if (k+jb <= n-1) {
        int ntrail = n - j;
        syrk("Upper", "Trans", ntrail, jb, -1.0, &A(k, j), lda,
             1.0, &A(j, j), lda);
      }
    }
  }
  else {
    for (int k = 0; k < n; k += nb) {
      int jb = (nb < n-k) ? nb : n-k;
      for (int i = k; i < n; i++)
        WORK(i) = 0.0;
      int j;
      for (j = k; j < k+jb; j++) {
        for (int i = j; i < n; i++) {
          if (j > k)
            WORK(i) += A(i, j-1) * A(i, j-1);
          WORK(n+i) = A(i, i) - WORK(i);
        }
        if (j > 0) {
          pvt = j;
          ajj = WORK(n+j);
          for (int i = j+1; i < n; i++) {
            if (WORK(n+i) > ajj) {
              pvt = i;
              ajj = WORK(n+i);
            }
          }
          if (ajj <= dstop || disnan(ajj)) {
            A(j, j) = ajj;
            *rank = j;
            *info = 1;
            return;
          }
        }
        if (j != pvt) {
          A(pvt, pvt) = A(j, j);
          swap(j,       &A(j,     0), lda, &A(pvt,   0),   lda);
          if (pvt < n-1)
            swap(n-pvt-1, &A(pvt+1, j), 1,   &A(pvt+1, pvt), 1);
          swap(pvt-j-1, &A(j+1,   j), 1,   &A(pvt,   j+1), lda);
          double dtemp = WORK(j);
          WORK(j)   = WORK(pvt);
          WORK(pvt) = dtemp;
          int itemp = piv[pvt];
          piv[pvt]  = piv[j];
          piv[j]    = itemp;
        }
        ajj = sqrt(ajj);
        A(j, j) = ajj;
        if (j < n-1) {
          gemv("No Trans", n-j-1, j-k, -1.0, &A(j+1, k), lda,
               &A(j, k), lda, 1.0, &A(j+1, j), 1);
          scal(n-j-1, 1.0/ajj, &A(j+1, j), 1);
        }
      }
      if (k+jb <= n-1) {
        int ntrail = n - j;
        syrk("Lower", "No Trans", ntrail, jb, -1.0, &A(j, k), lda,
             1.0, &A(j, j), lda);
      }
    }
  }
  *rank = n;
}

#undef A
#undef WORK
