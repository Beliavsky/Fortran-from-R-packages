#ifndef BLASGATEWAY_H
#define BLASGATEWAY_H
// Gateway to Fortran Blas functions used by the VarmaLoglik / VarLoglik package
#include "blasref.h"
#include <string.h>
//#include "printX.h"

static inline void axpy(int n, double alpha, double x[], int incx, double y[], int incy) {
  daxpy_(&n, &alpha, x, &incx, y, &incy);
}

static inline void copy(int n, double x[], int incx, double y[], int incy) {
  dcopy_(&n, x, &incx, y, &incy);
}

static inline double dot(int n, double x[], int incx, double y[], int incy) {
  return ddot_(&n, x, &incx, y, &incy);
}

static inline double nrm2(int n, double x[], int incx) {
  return dnrm2_(&n, x, &incx);
}

static inline void gemm(char *transa, char *transb, int m, int n, int k, double alpha,
          double a[], int lda, double b[], int ldb, double beta, double c[], int ldc) {
  dgemm_(transa, transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc, 1, 1);
}

static inline void gemv(char *transa, int m, int n, double alpha, double a[], int lda,
          double x[], int incx, double beta, double y[], int incy) {
  dgemv_(transa, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy, 1);
}

static inline void ger(int m, int n, double alpha, double x[], int incx, double y[], int incy,
         double a[], int lda) {
  dger_(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
}

static inline void getrf(int m, int n, double a[], int lda, int ipiv[], int *info) {
  dgetrf_(&m, &n, a, &lda, ipiv, info);
}

static inline void getrs(char *transa, int n, int nrhs, double a[], int lda, int ipiv[],
           double b[], int ldb, int *info) {
  dgetrs_(transa, &n, &nrhs, a, &lda, ipiv, b, &ldb, info, 1);
}

static inline void geev(char *jobvl, char *jobvr, int n, double a[], int lda, double wr[],
          double wi[], double vl[], int ldvl, double vr[], int ldvr, double work[], int
          lwork, int *info) {
  dgeev_(jobvl, jobvr, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork, info, 1, 1);
}

static inline int iamax(int n, double dx[], int incx) { // NOTE: Returns 0 for 1st elem. etc.
  return idamax_(&n, dx, &incx) - 1;
}

static inline void lacpy(char *uplo, int m, int n, double a[], int lda, double b[], int ldb) {
  dlacpy_(uplo, &m, &n, a, &lda, b, &ldb, 1);
}

static inline void laset(char *uplo, int m, int n, double alpha, double beta, double a[],
                         int lda) {
  dlaset_(uplo, &m, &n, &alpha, &beta, a, &lda, 1);
}

static inline void potrf(char *uplo, int n, double a[], int lda, int *info) {
  dpotrf_(uplo, &n, a, &lda, info, 1);
}

static inline void posv(char *uplo, int n, int nrhs, double a[], int lda, double b[],
                        int ldb, int *info) {
  dposv_(uplo, &n, &nrhs, a, &lda, b, &ldb, info, 1);
}

static inline void pstrf(char *uplo, int n, double a[], int lda, int piv[], int *rank,
								 double tol, double work[], int *info) {
#if defined(LOCAL_DPSTRF)
  rp_dpstrf(uplo, n, a, lda, piv, rank, tol, work, info);
#else
  dpstrf_(uplo, &n, a, &lda, piv, rank, &tol, work, info, 1);
#endif  
  for (int i=0; i<n; i++) piv[i]--;
}

static inline void scal(int m, double alpha, double *x, int incx) {
  dscal_(&m, &alpha, x, &incx);
}

static inline void swap(int n, double *x, int incx, double *y, int incy) {
  dswap_(&n, x, &incx, y, &incy);
}

static inline void syev(char *jobz, char *uplo, int n, double a[], int lda, double w[],
          double work[], int lwork, int *info) {
  dsyev_(jobz, uplo, &n, a, &lda, w, work, &lwork, info, 1, 1);
}

static inline void symm(char *side, char *uplo, int m, int n, double alpha, double a[],
          int lda, double b[], int ldb, double beta, double c[], int ldc) {
  dsymm_(side, uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc, 1, 1);
}

static inline void symv(char *uplo, int n, double alpha, double a[], int lda, double x[],
          int incx, double beta, double y[], int incy) {
  dsymv_(uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy, 1);
}

static inline void syr(char *uplo, int n, double alpha, double x[], int incx, double a[],
                       int lda) {
  dsyr_(uplo, &n, &alpha, x, &incx, a, &lda, 1);
}

static inline void syr2k(char *uplo, char *trans, int n, int k, double alpha, double a[],
           int lda, double b[], int ldb, double beta, double c[], int ldc) {
  dsyr2k_(uplo, trans, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc, 1, 1);
}

static inline void syrk(char *uplo, char *trans, int m, int n, double alpha, double a[], 
          int lda, double beta, double c[], int ldc) {
  dsyrk_(uplo, trans, &m, &n, &alpha, a, &lda, &beta, c, &ldc, 1, 1);
}

static inline void trmm(char *side, char *uplo, char *transa, char *diag, int m, int n,
          double alpha, double a[], int lda, double b[], int ldb) {
  dtrmm_(side, uplo, transa, diag, &m, &n, &alpha, a, &lda, b, &ldb, 1, 1, 1, 1);
}

static inline void trmv(char *uplo, char *transa, char *diag, int n, double a[], int lda, 
          double x[], int incx) {
  dtrmv_(uplo, transa, diag, &n, a, &lda, x, &incx, 1, 1, 1);
}

static inline void trsm(char *side, char *uplo, char *transa, char *diag, int m, int n,
          double alpha, double a[], int lda, double b[], int ldb) {
  dtrsm_(side, uplo, transa, diag, &m, &n, &alpha, a, &lda, b, &ldb, 1, 1, 1, 1);
}

static inline void trsv(char *uplo, char *transa, char *diag, int n, double a[], int lda, 
          double x[], int incx) {
  dtrsv_(uplo, transa, diag, &n, a, &lda, x, &incx, 1, 1, 1);
}

static inline int ilaenv(int ispec, char *name, char *opts, int n1, int n2, int n3, int n4) {
  return ilaenv_(&ispec, name, opts, &n1, &n2, &n3, &n4, (fstrlen)strlen(name), 1);
}

static inline double lamch(char *cmach) {
  return dlamch_(cmach, 1);
}

#endif
