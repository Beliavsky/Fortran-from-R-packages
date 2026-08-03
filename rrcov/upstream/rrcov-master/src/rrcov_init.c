//  VT::21.05.2018 - this file was added to fix the warning
//      "Found no calls to: 'R_registerRoutines', 'R_useDynamicSymbols'"
//
//  About registration of native symbols see for example: https://www.r-bloggers.com/1-easy-package-registration/
//      also here http://r.789695.n4.nabble.com/Registration-of-native-routines-td4728874.html
//
//	There is a new function in 'tools' package: 
//		package_native_routine_registration_skeleton()
//
//      - about Windows - take the 64 bit version of mingw!
//


#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void r_fast_mve(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void sest(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP covOPW(SEXP, SEXP, SEXP, SEXP);

/* .Fortran calls */
extern void F77_NAME(fsada)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(rlds)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"r_fast_mve", (DL_FUNC) &r_fast_mve, 12},
    {"sest",       (DL_FUNC) &sest,       14},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"covOPW", (DL_FUNC) &covOPW, 4},
    {NULL, NULL, 0}
};

static const R_FortranMethodDef FortranEntries[] = {
    {"fsada", (DL_FUNC) &F77_NAME(fsada), 15},
    {"rlds",  (DL_FUNC) &F77_NAME(rlds),  14},
    {NULL, NULL, 0}
};

void R_init_rrcov(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
