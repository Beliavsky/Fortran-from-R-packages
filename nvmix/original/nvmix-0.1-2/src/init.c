/* Register routines with R ***************************************************/

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>

#include "eval_nvmix_integral.h"
#include "eval_densmix_integrand.h"
#include "eval_gdensmix_integrand.h"
#include "logsumexp2.h"
#include "precond.h"

/* For .Call methods */
static const R_CallMethodDef callMethods[] = {
    {"eval_nvmix_integral", (DL_FUNC) &eval_nvmix_integral, 15},
    {"eval_densmix_integrand", (DL_FUNC) &eval_densmix_integrand, 7},
    {"eval_gdensmix_integrand_returnall", (DL_FUNC) &eval_gdensmix_integrand_returnall, 8},
        {"eval_gdensmix_integrand_LSE", (DL_FUNC) &eval_gdensmix_integrand_LSE, 8},
    {"logsumexp2", (DL_FUNC) &logsumexp2, 3},
    {NULL, NULL, 0}
};

/* Argument types for the .C method 'precond' */
static R_NativePrimitiveArgType precond_t[] = {
    REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP
};

/* For .C methods */
static const R_CMethodDef cMethods[] = {
    {"precond", (DL_FUNC) &precond, 9, precond_t},
    {NULL, NULL, 0}
};

void R_init_nvmix(DllInfo *dll)
{
    R_registerRoutines(dll, cMethods, callMethods, NULL, NULL); /* s. WRE (2015, Section 5.4) */
    R_useDynamicSymbols(dll, FALSE);
}
