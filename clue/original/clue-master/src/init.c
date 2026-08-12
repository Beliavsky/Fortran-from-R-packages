#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "clue.h"

static R_NativePrimitiveArgType solve_LSAP_t[3] = { REALSXP, INTSXP, INTSXP };
static R_NativePrimitiveArgType clue_dissimilarity_count_inversions_t[4] = { REALSXP, REALSXP, INTSXP, REALSXP };

static R_NativePrimitiveArgType deviation_from_ultrametricity_t[4] = { REALSXP, INTSXP, REALSXP, LGLSXP }; 
static R_NativePrimitiveArgType deviation_from_ultrametricity_gradient_t[3] = { REALSXP, INTSXP, REALSXP };
static R_NativePrimitiveArgType deviation_from_additivity_t[4] = { REALSXP, INTSXP, REALSXP, LGLSXP };
static R_NativePrimitiveArgType deviation_from_additivity_gradient_t[3] = { REALSXP, INTSXP, REALSXP };

static R_NativePrimitiveArgType ls_fit_ultrametric_by_iterative_reduction_t[7] = { REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, LGLSXP };
static R_NativePrimitiveArgType ls_fit_ultrametric_by_iterative_projection_t[7] = { REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, LGLSXP };
static R_NativePrimitiveArgType ls_fit_addtree_by_iterative_reduction_t[7] = { REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, LGLSXP };
static R_NativePrimitiveArgType ls_fit_addtree_by_iterative_projection_t[7] = { REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, LGLSXP };

#define CDEF(name)  {#name, (DL_FUNC) &name, sizeof(name ## _t)/sizeof(name ## _t[0]), name ##_t}

static const R_CMethodDef cMethods[] = {
    CDEF(solve_LSAP),
    CDEF(clue_dissimilarity_count_inversions),
    CDEF(deviation_from_ultrametricity),
    CDEF(deviation_from_ultrametricity_gradient),
    CDEF(deviation_from_additivity),
    CDEF(deviation_from_additivity_gradient),
    CDEF(ls_fit_ultrametric_by_iterative_reduction),
    CDEF(ls_fit_ultrametric_by_iterative_projection),
    CDEF(ls_fit_addtree_by_iterative_reduction),
    CDEF(ls_fit_addtree_by_iterative_projection),
    {NULL, NULL, 0}
};

void R_init_clue(DllInfo *dll)
{
    R_registerRoutines(dll, cMethods, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
