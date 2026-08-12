
#ifndef R_PKG_EXPM_H
#define R_PKG_EXPM_H

#define USE_FC_LEN_T
#include <Rconfig.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include "R_NLS_locale.h"

typedef enum {Ward_2, Ward_1, Ward_buggy_octave} precond_type;

SEXP R_dgebal(SEXP x, SEXP type);
SEXP R_zgebal(SEXP x, SEXP type);

SEXP do_expm(SEXP x, SEXP kind);
void expm(double *x, int n, double *z, precond_type precond_kind);

SEXP R_matexp_MH09(SEXP x, SEXP p);
void matexp_MH09(double *x, int n, const int p, double *ret);

// The legacy code: -----------------------------

// matexp.f: matexpRBS  << is what I'd want
int F77_NAME(matexprbs)(int *ideg, int *m, double *t,
			double *H, int *iflag);

// matrexp.f:
int F77_NAME(matrexp)(double* a, int* n, int* ntaylor, int* npade,
		      double* accuracy);
// matrexpO.f: matrexpO << is what I'd want
int F77_NAME(matrexpo)(double* a, int* n, int* ntaylor, int* npade,
		       double* accuracy);

#endif /* R_PKG_EXPM_H */
