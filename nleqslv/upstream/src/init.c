#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void deactivatenleq(void);

/* .Call calls */
extern SEXP nleqslv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"deactivatenleq", (DL_FUNC) &deactivatenleq, 0},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"nleqslv", (DL_FUNC) &nleqslv, 9},
    {NULL, NULL, 0}
};

void R_init_nleqslv(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
