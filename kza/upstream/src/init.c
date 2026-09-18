#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/RS.h>
#include <R_ext/Applic.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>

#include "kz.h"

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

/* kz.c */
extern SEXP kz1d(SEXP x, SEXP window, SEXP iterations);
extern SEXP kz2d(SEXP x, SEXP window, SEXP iterations);
/* extern void copyArray(SEXP destination, SEXP source); */
extern SEXP kz3d(SEXP x, SEXP window, SEXP iterations);
extern SEXP kz(SEXP x, SEXP window, SEXP iterations);
extern SEXP kzs(SEXP ans, SEXP y, SEXP x, SEXP dx, SEXP q, SEXP iterations, SEXP minx, SEXP maxx);
/* kza.c */
extern SEXP kza1d(SEXP v, SEXP y, SEXP window, SEXP iterations, SEXP minimum_window_length, SEXP tolerance, SEXP normalize_prob);
extern SEXP kza2d(SEXP v, SEXP kz, SEXP window, SEXP iterations, SEXP minimum_window_length, SEXP tolerance, SEXP normalize_prob);
extern SEXP kza3d(SEXP v, SEXP kz, SEXP window, SEXP iterations, SEXP minimum_window_length, SEXP tolerance, SEXP normalize_prob);
extern SEXP kza(SEXP x, SEXP kz, SEXP window, SEXP iterations, SEXP min_window, SEXP tol, SEXP normalize_prob);
/* kzf.c */
/* extern double R_maximum(SEXP v); */
/* extern void R_differenced(SEXP y, SEXP d, SEXP dprime, int q); */
/* extern void differenced(double y[], double d[], double dprime[], long n, int q); */
/* kzsv.c */
extern SEXP kzsv(SEXP kza_data, SEXP kz_data, SEXP window, SEXP minimum_window_length, SEXP iterations, SEXP tolerance);
extern SEXP R_kzsv(SEXP kza_data, SEXP kz_data, SEXP window, SEXP minimum_window_length, SEXP tolerance);

static const R_CallMethodDef CallEntries[] = {
    {"kz", (DL_FUNC) &kz, 3},
    {"kz1d", (DL_FUNC) &kz1d, 3},
    {"kz2d", (DL_FUNC) &kz2d, 3},
    /* {"copyArray", (DL_FUNC) &copyArray, 2}, */
    {"kz3d", (DL_FUNC) &kz3d, 3},
    {"kza", (DL_FUNC) &kza, 7},
    {"kza1d", (DL_FUNC) &kza1d, 7},
    {"kza2d", (DL_FUNC) &kza2d, 7},
    {"kza3d", (DL_FUNC) &kza3d, 7},
    /* {"R_maximum", (DL_FUNC) &R_maximum, 1}, */
    /* {"R_differenced", (DL_FUNC) &R_differenced, 4}, */
    /* {"differenced", (DL_FUNC) &differenced, 5}, */
    {"kzs", (DL_FUNC) &kzs, 8},
    {"kzsv", (DL_FUNC) &kzsv, 6},
    {"R_kzsv", (DL_FUNC) &R_kzsv, 5},
    {NULL, NULL, 0}
};

void R_init_kza(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

