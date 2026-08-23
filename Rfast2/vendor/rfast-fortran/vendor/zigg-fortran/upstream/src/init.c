
#include <R.h>
#include <Rinternals.h>         /* for SEXP */
#include <R_ext/Rdynload.h>

SEXP zrnorm(SEXP ns);
SEXP zrexp(SEXP ns);
SEXP zrunif(SEXP ns);
SEXP zsetseed(SEXP ns);

static const R_CallMethodDef CallEntries[] = {
    {"zrnorm",   (DL_FUNC) &zrnorm,   1},
    {"zrexp",    (DL_FUNC) &zrexp,    1},
    {"zrunif",   (DL_FUNC) &zrunif,   1},
    {"zsetseed", (DL_FUNC) &zsetseed, 1},
    {NULL, NULL, 0}
};

void R_init_zigg(DllInfo *dll) {
    /* used by external packages linking to internal serialization code from C */
    R_RegisterCCallable("zigg", "zrnorm",   (DL_FUNC) &zrnorm);
    R_RegisterCCallable("zigg", "zrexp",    (DL_FUNC) &zrexp);
    R_RegisterCCallable("zigg", "zrunif",   (DL_FUNC) &zrunif);
    R_RegisterCCallable("zigg", "zsetseed", (DL_FUNC) &zsetseed);

    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
