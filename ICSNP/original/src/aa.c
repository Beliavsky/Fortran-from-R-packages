#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void norming(void *, void *, void *);
extern void pairdiff(void *, void *, void *);
extern void pairprod(void *, void *, void *);
extern void pairsum(void *, void *, void *);
extern void sum_of_diff_sign_outers(void *, void *, void *);
extern void sum_of_sign_outers(void *, void *, void *);
extern void symm_huber(void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"norming",                 (DL_FUNC) &norming,                 3},
    {"pairdiff",                (DL_FUNC) &pairdiff,                3},
    {"pairprod",                (DL_FUNC) &pairprod,                3},
    {"pairsum",                 (DL_FUNC) &pairsum,                 3},
    {"sum_of_diff_sign_outers", (DL_FUNC) &sum_of_diff_sign_outers, 3},
    {"sum_of_sign_outers",      (DL_FUNC) &sum_of_sign_outers,      3},
    {"symm_huber",              (DL_FUNC) &symm_huber,              6},
    {NULL, NULL, 0}
};

void R_init_ICSNP(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
