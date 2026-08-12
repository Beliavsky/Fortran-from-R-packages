#include "expm.h"

/* The C API :*/
void matpow  (double *x, int n, int k, double *z);
void matpow_z(Rcomplex *x, int n, int k, Rcomplex *z);


/* as .Call()ed from R */
SEXP R_matpow(SEXP x, SEXP k);

