/*
 *   socp.h
 *
 *   Second-Order Cone Programming
 *   Header file for socp.c and socp_mex.c
 *
 *   mlobo@isl.stanford.edu -- 96/97
 */



/* basic macros */

#define SQR(x)    ((x)*(x))
#define MAX(x,y)  ((x) > (y) ? (x) : (y))
#define MIN(x,y)  ((x) < (y) ? (x) : (y))


/* constants for socp algorithm */

#define MAX_ITER_PLANE 20
#define MAX_LAMBDA2 1e-2  /* plane search stopping crit. */
#define DIV_ALPHA 2
#define MIN_ALPHA 1e-6    /* max. of about 20 line search iterations */



#ifdef nounderscores
#define ddot_ ddot
#define dcopy_ dcopy
#define daxpy_ daxpy
#define dscal_ dscal
#define dgemv_ dgemv
#define dsyr_ dsyr
#define dsyrk_ dsyrk
#define dposvx_ dposvx
#define dgelss_ dgelss
#endif


/* BLAS 1 */
extern double ddot_(int *n, double *x, int *incx, double *y, int *incy);
extern void dcopy_(int *n, double *x, int *incx, double *y, int *incy);
extern void daxpy_(int *n, double *alpha, double *x, int *incx, double *y, int *incy);
extern void dscal_(int *n, double *alpha, double *x, int *incx);

/* BLAS 2 */
extern void dgemv_(char *transa, int *m, int *n, double *alpha, double *a, int *lda, 
            double *x, int *incx, double *beta, double *y, int *incy);
extern void dsyr_(char *uplo, int *n, double *alpha, double *x, int *incx, double *a, 
           int *lda);
/* BLAS 3 */
extern void dsyrk_(char *uplo, char *trans, int *n, int *k, double *alpha, double *a, 
            int *lda, double *beta, double *c, int *ldc);
/* LAPACK */
extern void dposvx_(char *fact, char *uplo, int *n, int *nrhs, double *a, int *lda, 
                    double *af, int *ldaf, int *equed, double *s, double *b, 
                    int *ldb, double *x, int *ldx, double *rcond, double *ferr, 
                    double *berr, double *work, int *iwork, int *info);
extern void dgelss_(int *m, int *n, int *nrhs, double *a, int *lda, double *b, 
                    int *ldb, double *s, double *rcond, int *rank, double *work, 
                    int *lwork, int *info);



/* socp.c */

void socp_getwork(
	/* input args.:  problem dimensions and max. num. of iterations */
 int *L,
 int *N,
 int *n,
 int *max_iter,
 int *out_mode,
        /* output args.:  dimensions of history output matrix, */
        /*   number of doubles, pointers and ints required for workspace */
 int *mhist,
 int *nhist,
 int *ndbl,
 int *nint
);


int socp(
   /* problem dimensions */
 int *L,
 int *N,
 int *n,
   /* problem data */
 double *f,
 double *A,
 double *b,
   /* in:  initial primal and dual strictly feasible points */
   /* out: final points */
 double *x,
 double *z,
   /* stopping criteria, */
   /*   *iter on entry is max. number of iterations, */
   /*   on exit actual number performed */
 double *abs_tol,
 double *rel_tol,
 double *target,
 int *iter,
   /* algorithm parameter */
 double *Nu,
   /* reason for exit, output matrix with extra info */
   /* out_mode specifies what will be stored in *hist */
 int *info,
 int *out_mode,
 double *hist,
   /* workspace, use socp_getwork() to determine required sizes */
 double *dblwork,
 int *intwork
);
