/*
 *  Native routines registration
 */

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "expm-eigen.h"
#include "expm.h"
#include "logm-eigen.h"
#include "matpow.h"

static const R_CallMethodDef CallEntries[] = {
    {"do_expm", (DL_FUNC) &do_expm, 2},
    {"R_matpow", (DL_FUNC) &R_matpow, 2},
    {"R_dgebal", (DL_FUNC) &R_dgebal, 2},
    {"R_zgebal", (DL_FUNC) &R_zgebal, 2},
    {"do_expm_eigen", (DL_FUNC) &do_expm_eigen, 2},
    {"do_logm_eigen", (DL_FUNC) &do_logm_eigen, 2},
    {"R_matexp_MH09", (DL_FUNC) &R_matexp_MH09, 2},
    {NULL, NULL, 0}
};

static const R_FortranMethodDef FortEntries[] = {
    {"matexpRBS", (DL_FUNC) &F77_SUB(matexprbs), 5}, // ./matexp.f
    {"matrexp",   (DL_FUNC) &F77_SUB(matrexp),   5}, // ./matrexp.f
    {"matrexpO",  (DL_FUNC) &F77_SUB(matrexpo),  5}, // ./matrexpO.f
    {NULL, NULL, 0}
};


void R_init_expm(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, FortEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
    /* callable C code from other packages C code :*/
    R_RegisterCCallable("expm",        "expm",        (DL_FUNC) expm);
    R_RegisterCCallable("matpow",      "matpow",      (DL_FUNC) matpow);
    R_RegisterCCallable("expm_eigen",  "expm_eigen",  (DL_FUNC) expm_eigen);
    R_RegisterCCallable("logm_eigen",  "logm_eigen",  (DL_FUNC) logm_eigen);
    R_RegisterCCallable("matexp_MH09", "matexp_MH09", (DL_FUNC) matexp_MH09);
}
