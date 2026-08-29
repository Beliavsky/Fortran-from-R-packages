#include <R.h>
#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */
extern void F77_NAME(gausq2)(void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"gausq2", (DL_FUNC) &F77_NAME(gausq2), 5},
    {NULL, NULL, 0}
};

/* .C calls */
extern void mpoisdev (double *, double *, double *, int *);
extern void mbinomdev (double *, int *, double *, double *, int *, int *, int *);
extern void mnbinomdev (double *, double *, double *, double *, int *, int *);

static const R_CMethodDef CEntries[] = {
  {"mpoisdev", (DL_FUNC) &mpoisdev, 4},
  {"mbinomdev", (DL_FUNC) &mbinomdev, 7},
  {"mnbinomdev", (DL_FUNC) &mnbinomdev, 6},
  {NULL, NULL, 0}
};

void R_init_statmod(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
