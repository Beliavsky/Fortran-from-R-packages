/*
 *   Copyright (c) 1998-2000 Lucent Technologies.
 *   See README file for details.
 *
 *
 *   Headers for math utility functions.
 */

#ifndef I_MUT_H
#define I_MUT_H

#include <math.h>

typedef struct {
  double *Z;   /* jacobian matrix, length p*p          */
  double *Q;   /* eigenvalue matrix, length p*p        */
  double *wk;  /* work vector in eig_solve, length p   */
  double *dg;  /* diag vector in eigd, length p        */
  int p;       /* dimension */
  int st;      /* status    */
  int sm;      /* requested decomposition */
} jacobian;

/* m_jacob.c */
extern int jac_reqd(int p);
extern double *jac_alloc(jacobian *J, int p, double *wk);
extern void jacob_dec(jacobian *J, int meth);
extern int jacob_solve(jacobian *J, double *v);
extern int jacob_hsolve(jacobian *J, double *v);
extern double jacob_qf(jacobian *J, double *v);
extern void chol_dec(double *A, int n, int p);
extern int chol_solve(double *A, double *v, int n, int p);
extern int chol_hsolve(double *A, double *v, int n, int p);
extern double chol_qf(double *A, double *v, int n, int p);
extern void eig_dec(double *X, double *P, int d);
extern int eig_solve(jacobian *J, double *x);
extern int eig_hsolve(jacobian *J, double *v);
extern double eig_qf(jacobian *J, double *v);

/* m_max.c */
extern double max_grid(double (*f)(), double xlo, double xhi, int n, char flag);
extern double max_golden(double (*f)(), double xlo, double xhi, int n, double tol,
                  int *err, char flag);
extern double max_quad(double (*f)(), double xlo, double xhi, int n, double tol,
                int *err, char flag);
extern double max_nr(int (*F)(), double *coef, double *old_coef, double *f1, double *delta, 
              jacobian *J, int p, int maxit, double tol, int *err);

/* m_qr.c */
extern void qr(double *X, int n, int p, double *w);
extern void qrinvx(double *R, double *x, int n, int p);
extern void qrtinvx(double *R, double *x, int n, int p);
extern void qrsolv(double *R, double *x, int n, int p);

/* m_svd.c */
extern void svd(double *x, double *p, double *q, int d, int mxit);
extern void hsvdsolve(double *x, double *w, double *P, double *D, double *Q, int d,
               double tol);
extern int svdsolve(double *x, double *w, double *P, double *D, double *Q, int d, 
             double tol);

/* m_solve.c */
extern double solve_secant(double (*f)(), double c, double xlo, double xhi, double tol,
                    int bd_flag, int *err);
extern double solve_nr(double (*f)(), double (*f1)(), double c, double x0, double tol,
                int *err);
extern double solve_fp(double (*f)(), double x0, double tol, int maxit);

/* m_vector.c */
extern void setzero(double *v, int p);
extern void unitvec(double *x, int k, int p);
extern void addouter(double *A, double *v1, double *v2, int p, double c);
extern void multmatscal(double *A, double z, int n);
extern void transpose(double *x, int m, int n);
extern double innerprod(double *v1, double *v2, int p);
extern double m_trace(double *x, int n);

#define BDF_NONE  0
#define BDF_EXPLEFT  1
#define BDF_EXPRIGHT 2

/* return codes for functions optimized by max_nr */
#define NR_OK 0
#define NR_INVALID 1
#define NR_BREAK   2
#define NR_REDUCE  3
#define NR_NCON  10
#define NR_NDIV  11


/* jacobian status definitions */
#define JAC_RAW 0
#define JAC_CHOL 1
#define JAC_EIG  2
#define JAC_EIGD 3

/*  Numerical Integration Stuff
 */
#define MXRESULT 5
#define MXIDIM  10  /* max. dimension */
extern void simpsonm(), simpson4(), integ_disc(), integ_circ();
extern void integ_sphere(), monte(), rn3();
extern double simpson(), sptarea();

/*  Density, distribution stuff
 */

#ifndef PI
#define PI  3.141592653589793238462643
#endif
#define PIx2 6.283185307179586476925286        /* 2*pi */
#define HF_LG_PIx2  0.918938533204672741780329736406    /* 0.5*log(2*pi) */
#define SQRT2 1.4142135623730950488

#define LOG_ZERO -1e100
#define D_0 ((give_log) ? LOG_ZERO : 0.0)
#define D_1 ((give_log) ? 0.0 : 1.0)
#define DEXP(x)   ((give_log) ? (x) : exp(x))
#define FEXP(f,x) ((give_log) ? -0.5*log(f)+(x) : exp(x)/sqrt(f))

#define INVALID_PARAMS 0.0

extern double stirlerr(double), bd0(double, double);
extern double dbinom_raw(double, double, double, double, int), 
              dpois_raw(double, double, int);
extern double dbinom(int, int, double, int), 
              dpois(int, double, int), 
              dnbinom(int, double, double, int), 
              dbeta(double, double, double, int), 
              dgamma(double, double, double, int), 
              dt(double, double, int), 
              df(double, double, double, int), 
              dhyper(int, int, int, int, int);
extern double dchisq(double, double, int);

extern double igamma(double, double), ibeta(double, double, double);
extern double pf(double, double, double), pchisq(), mut_pnorm(double, double, double);
#define pchisq(x,df) igamma((x)/2.0,(df)/2.0)

#endif  /* define I_MUT_H */
