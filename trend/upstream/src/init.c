#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* Created by
 * tools::package_native_routine_registration_skeleton
 * and manually modified by 
 * Thorsten Pohlert 2020-01-11
*/

/* .Fortran calls */
extern void F77_NAME(mcbr)(double *r, int *n, int *m, double *pval);
extern void F77_NAME(mcbu)(double *r, int *n, int *m, double *pval);
extern void F77_NAME(mcsnht)(double *r, int *n, int *m, double *pval);

static const R_FortranMethodDef FortranEntries[] = {
    {"mcbr",   (DL_FUNC) &F77_NAME(mcbr),   4},
    {"mcbu",   (DL_FUNC) &F77_NAME(mcbu),   4},
    {"mcsnht", (DL_FUNC) &F77_NAME(mcsnht), 4},
    {NULL, NULL, 0}
};

void R_init_trend(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
