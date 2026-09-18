// Include file declaring the reference blas functions
//
#ifndef BLASREF_H
#define BLASREF_H

#if defined(__GNUC__)
typedef long int fstrlen;
#else
typedef int fstrlen;
#endif

void daxpy_(int *n, double *alpha, double x[], int *incx, double y[], int *incy);

void dcopy_(int *n, double x[], int *incx, double y[], int *incy);

void dswap_(int *n, double x[], int *incx, double y[], int *incy);

double ddot_(int *n, double x[], int *incx, double y[], int *incy);

void dgeev_(char *jobvl, char *jobvr, int *n, double a[], int *lda, double wr[], double
	    wi[], double vl[], int *ldvl, double vr[], int *ldvr, double work[], int
	    *lwork, int *info, fstrlen lenjobvl, fstrlen lenjobvr);

void dgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha,
            double a[], int *lda, double b[], int *ldb, double *beta, double c[], 
            int *ldc, fstrlen lentransa, fstrlen lentransb);

void dgemv_(char *transa, int *m, int *n, double *alpha, double a[], int *lda,
            double x[], int *incx, double *beta, double y[], int *incy, fstrlen lentransa);

void dger_(int *m, int *n, double *alpha, double x[], int *incx, double y[], 
           int *incy, double a[], int *lda);

void dgetrf_(int *m, int *n, double a[], int *lda, int ipiv[], int *info);

void dgetrs_(char *transa, int *n, int *nrhs, double a[], int *lda, int ipiv[],
             double b[], int *ldb, int *info, fstrlen lentransa);

void dlacpy_(char *uplo, int *m, int *n, double a[], int *lda, double b[], int *ldb,
	     fstrlen lenuplo);

void dlaset_(char *uplo, int *m, int *n, double *alpha, double *beta, double a[], int
	     *lda, fstrlen lenuplo);

double dnrm2_(int *n, double x[], int *incx);

void dpotrf_(char *uplo, int *n, double a[], int *lda, int *info, fstrlen lenuplo);

double dlamch_(char *cmach, fstrlen lencmach);

void dpstrf_(char *uplo, int *n, double a[], int *lda, int piv[], int *rank, double *tol,
	     double work[], int *info, fstrlen lenuplo);

void dposv_(char *uplo, int *n, int *nrhs, double a[], int *lda, double b[], int *ldb,
            int *info, fstrlen lenuplo);

void dscal_(int *m, double *alpha, double *x, int *incx);

void dspr_(char *uplo, int *n, double *alpha, double x[], int *incx, double ap[],
           fstrlen lenuplo);

void dspr2_(char *uplo, int *n, double *alpha, double x[], int *incx, double y[],
            int *incy, double ap[], fstrlen lenuplo);

void dsyev_(char *jobz, char *uplo, int *n, double a[], int *lda, double w[],
            double work[], int *lwork, int *info, fstrlen lenjobz, fstrlen lenuplo);

void dsymm_(char *side, char *uplo, int *m, int *n, double *alpha, double a[], 
            int *lda, double b[], int *ldb, double *beta, double c[], int *ldc,
            fstrlen lenside, fstrlen lenuplo);

void dsymv_(char *uplo, int *n, double *alpha, double a[], int *lda, double x[],
            int *incx, double *beta, double y[], int *incy, fstrlen lenuplo);
void dsyr_(char *uplo, int *n, double *alpha, double x[], int *incx, double a[],
          int *lda, fstrlen lenuplo);

void dsyr2k_(char *uplo, char *trans, int *m, int *n, double *alpha, double a[], 
             int *lda, double b[], int *ldb, double *beta, double c[], int *ldc, 
             fstrlen lenuplo, fstrlen lentrans);

void dsyrk_(char *uplo, char *trans, int *m, int *n, double *alpha, double a[], int *lda,
            double *beta, double c[], int *ldc, fstrlen lenuplo, fstrlen lentrans);

void dtrmm_(char *side, char *uplo, char *transa, char *diag, int *m, int *n,
            double *alpha, double a[], int *lda, double b[], int *ldb, fstrlen lenside, 
            fstrlen lenuplo, fstrlen lentransa, fstrlen lendiag);

void dtrmv_(char *uplo, char *trans, char *diag, int *n, double A[], int *lda, double x[],
            int *incx, fstrlen lenuplo, fstrlen lentrans, fstrlen lendiag);

void dtrsm_(char *side, char *uplo, char *transa, char *diag, int *m, int *n, double
            *alpha, double a[], int *lda, double b[], int *ldb, fstrlen lenside, fstrlen
            lenuplo, fstrlen lentransa, fstrlen lendiag);

void dtrsv_(char *uplo, char *transa, char *diag, int *n, double a[], int *lda, double
            x[], int *incx, fstrlen lenuplo, fstrlen lentransa, fstrlen lendiag);

int idamax_(int *n, double dx[], int *incx);

int ilaenv_(int *ispec, char *name, char *opts, int *n1, int *n2, int *n3, int *n4,
            fstrlen lenname, fstrlen lenopts);

void rp_dpstrf(char *uplo, int n, double a[], int lda, int piv[], int *rank, double tol,
               double work[], int *info);

#endif // BLASREF_H
