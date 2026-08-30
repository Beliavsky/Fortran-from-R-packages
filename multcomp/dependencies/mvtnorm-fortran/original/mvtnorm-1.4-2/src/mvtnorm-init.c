
#include "mvtnorm.h"
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>


void C_mvtdst(int *n, int *nu, double *lower, double *upper,
              int *infin, double *corr, double *delta,
              int *maxpts, double *abseps, double *releps,
              double *error, double *value, int *inform, int *rnd)
{

    if (rnd[0]) GetRNGstate();

    /* call FORTRAN subroutine */
    F77_CALL(mvtdst)(n, nu, lower, upper, 
                     infin, corr, delta,
                     maxpts, abseps, releps, 
                     error, value, inform);

    if (rnd[0]) PutRNGstate();

}

// TVPACK n=3
void C_tvtlr(int *NU, double *H, double *R, double *EPSI, double *TVTL) {

    F77_CALL(tvtlrcall)(NU, H, R, EPSI, TVTL); 
}

// TVPACK n=2
void C_bvtlr(int *NU, double *DH, double *DK, double *R, double *BVTL) {

    F77_CALL(bvtlrcall)(NU, DH, DK, R, BVTL );
}


static const R_CMethodDef cMethods[] = {
    {"C_mvtdst", (DL_FUNC) &C_mvtdst, 14, (R_NativePrimitiveArgType[14]){INTSXP, INTSXP, REALSXP, REALSXP, 
                                           INTSXP, REALSXP, REALSXP, 
                                           INTSXP, REALSXP, REALSXP, 
                                           REALSXP, REALSXP, INTSXP, INTSXP}}, 
    {"C_tvtlr", (DL_FUNC) &C_tvtlr, 5, (R_NativePrimitiveArgType[13]){INTSXP, REALSXP, REALSXP, REALSXP, REALSXP}},
    {"C_bvtlr", (DL_FUNC) &C_bvtlr, 5, (R_NativePrimitiveArgType[13]){INTSXP, REALSXP, REALSXP, REALSXP, REALSXP}},
    {NULL, NULL, 0, NULL}
};

static const R_CallMethodDef callMethods[] = {
    {"R_miwa", (DL_FUNC) &R_miwa, 5},
    {"R_ltMatrices_solve", (DL_FUNC) &R_ltMatrices_solve, 6},
    {"R_ltMatrices_solve_C", (DL_FUNC) &R_ltMatrices_solve_C, 5},
    {"R_ltMatrices_logdet", (DL_FUNC) &R_ltMatrices_logdet, 5},
    {"R_ltMatrices_Mult", (DL_FUNC) &R_ltMatrices_Mult, 6},
    {"R_ltMatrices_colSumsdnorm", (DL_FUNC) &R_ltMatrices_colSumsdnorm, 3},
    {"R_lpmvnorm", (DL_FUNC) &R_lpmvnorm, 11},
    {"R_slpmvnorm", (DL_FUNC) &R_slpmvnorm, 10},
    {"R_syMatrices_chol", (DL_FUNC) &R_syMatrices_chol, 3},
    {"R_syMatrices_invchol", (DL_FUNC) &R_syMatrices_invchol, 3},
    {NULL, NULL, 0}
};

void attribute_visible R_init_mvtnorm(DllInfo *dll)
{
    R_registerRoutines(dll, cMethods, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
    R_RegisterCCallable("mvtnorm", "C_mvtdst", (DL_FUNC) &C_mvtdst);
    R_RegisterCCallable("mvtnorm", "R_miwa", (DL_FUNC) &R_miwa);
    R_RegisterCCallable("mvtnorm", "R_ltMatrices_solve", (DL_FUNC) &R_ltMatrices_solve);
    R_RegisterCCallable("mvtnorm", "R_ltMatrices_solve_C", (DL_FUNC) &R_ltMatrices_solve_C);
    R_RegisterCCallable("mvtnorm", "R_ltMatrices_logdet", (DL_FUNC) &R_ltMatrices_logdet);
    R_RegisterCCallable("mvtnorm", "R_ltMatrices_Mult", (DL_FUNC) &R_ltMatrices_Mult);
    R_RegisterCCallable("mvtnorm", "R_ltMatrices_colSumsdnorm", (DL_FUNC) &R_ltMatrices_colSumsdnorm);
    R_RegisterCCallable("mvtnorm", "R_lpmvnorm", (DL_FUNC) &R_lpmvnorm);
    R_RegisterCCallable("mvtnorm", "R_slpmvnorm", (DL_FUNC) &R_slpmvnorm);
    R_RegisterCCallable("mvtnorm", "R_syMatrices_chol", (DL_FUNC) &R_syMatrices_chol);
    R_RegisterCCallable("mvtnorm", "R_syMatrices_invchol", (DL_FUNC) &R_syMatrices_invchol);
}
