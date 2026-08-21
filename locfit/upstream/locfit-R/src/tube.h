/*
 *   Copyright (c) 1998-2001 Catherine Loader, Jiayang Sun
 *   See README file for details.
 *
 *
 *   Headers for the tube library.
 */

#ifndef I_TUBE_H
#define I_TUBE_H

/*
 * public functions needed by routines calling the tube library.
 */

/* from scb_crit.c */
extern double area(int d);
extern double tailp_uniform(double c, double *k0, int m, int d, int s, double n);
extern double tailp_gaussian(double c, double *k0, int m, int d, int s, double n);
extern double tailp_tprocess(double c, double *k0, int m, int d, int s, double n);
extern double taild_uniform(double c, double *k0, int m, int d, int s, double n);
extern double taild_gaussian(double c, double *k0, int m, int d, int s, double n);
extern double taild_tprocess(double c, double *k0, int m, int d, int s, double n);
extern double tailp(double c, double *k0, int m, int d, int s, double nu, int process);
extern double taild(double c, double *k0, int m, int d, int s, double nu, int process);
extern double critval(double alpha, double *k0, int m, int d, int s, double nu, int process);

/* from scb_cons.c */
int k0_reqd(int d, int n, int uc);
void assignk0(double *z, int d, int n);
void rproject(double *y, double *A, double *R, int n, int p);
double k2c(double *lij, double *A, int m, int dd, int d);
double k2x(double *lij, double *A, int m, int d, int dd);
void d2c(double *ll, double *nn, double *li, double *ni, double *lij, double *nij, double *M, int m, int dd, int d);
void d2x(double *li, double *lij, double *nij, double *M, int m, int dd, int d);
int k0x(double *x, int d, double *kap, double *M);
void d1c(double *li, double *ni, int m, int d, double *M);
void d1x(double *li, double *ni, int m, int d, double *M);
int l1x(double *x, int d, double *lap, double *M);
int m0x(double *x, int d, double *m0, double *M);
int n0x(double *x, int d, double *n0, double *M);
int kodf(double *ll, double *ur, int *mg, double *kap, double *lap);
int tube_constants(int (*f)(), int d, int m, int ev, int *mg, double *fl, double *kap, double *wk, int terms, int uc);

/*
 * stuff used internally.
 */

#include "mutil.h"

#define TUBE_MXDIM 10

/*
 * definitions for integration methods.
 * these match locfit evaluation structures where applicable.
 */

#define ISIMPSON  4    /* grid */
#define ISPHERIC 11    /* circle or sphere */
#define IDERFREE 25    /* derivative free */
#define IMONTE   30    /* monte carlo */

#ifndef PI
#define PI    3.141592653589793238462643

#endif

#define ONE_SIDED 1
#define TWO_SIDED 2

#define UNIF    400
#define GAUSS   401
#define TPROC   402
#endif  /* define I_TUBE_H */
