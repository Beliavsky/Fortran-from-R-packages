#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP rewens_c(SEXP, SEXP);

extern SEXP gcrp_c(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"rewens_c", (DL_FUNC) &rewens_c, 2},
    {"gcrp_c",   (DL_FUNC) &gcrp_c,   3},
    {NULL, NULL, 0}
};

void R_init_ewens(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
