#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* Layer 1: One-shot solver */
SEXP R_scip_solve(SEXP obj, SEXP Ai, SEXP Ap, SEXP Ax,
                  SEXP m, SEXP n,
                  SEXP lhs, SEXP rhs,
                  SEXP vtype, SEXP lb, SEXP ub,
                  SEXP ctrl);

/* Layer 2: Model-building API */
SEXP R_scip_model_create(SEXP name);
SEXP R_scip_model_add_var(SEXP ext, SEXP obj, SEXP lb, SEXP ub,
                           SEXP vtype, SEXP name);
SEXP R_scip_model_add_vars(SEXP ext, SEXP obj, SEXP lb, SEXP ub,
                            SEXP vtype, SEXP names);
SEXP R_scip_model_add_linear_cons(SEXP ext, SEXP vars, SEXP coefs,
                                   SEXP lhs, SEXP rhs, SEXP name);
SEXP R_scip_model_add_quadratic_cons(SEXP ext,
                                      SEXP linvars, SEXP lincoefs,
                                      SEXP quadvars1, SEXP quadvars2,
                                      SEXP quadcoefs,
                                      SEXP lhs, SEXP rhs, SEXP name);
SEXP R_scip_model_add_sos1_cons(SEXP ext, SEXP vars, SEXP weights,
                                 SEXP name);
SEXP R_scip_model_add_sos2_cons(SEXP ext, SEXP vars, SEXP weights,
                                 SEXP name);
SEXP R_scip_model_add_indicator_cons(SEXP ext, SEXP binvar,
                                      SEXP vars, SEXP coefs,
                                      SEXP rhs, SEXP name);
SEXP R_scip_model_set_param(SEXP ext, SEXP name, SEXP value);
SEXP R_scip_model_set_objective_sense(SEXP ext, SEXP sense);
SEXP R_scip_model_optimize(SEXP ext);
SEXP R_scip_model_get_status(SEXP ext);
SEXP R_scip_model_get_solution(SEXP ext);
SEXP R_scip_model_get_nsols(SEXP ext);
SEXP R_scip_model_get_sol(SEXP ext, SEXP k);
SEXP R_scip_model_get_info(SEXP ext);
SEXP R_scip_model_get_nvars(SEXP ext);
SEXP R_scip_model_get_nconss(SEXP ext);
SEXP R_scip_model_free(SEXP ext);

static const R_CallMethodDef CallEntries[] = {
    /* Layer 1 */
    {"R_scip_solve",                      (DL_FUNC) &R_scip_solve,                      12},
    /* Layer 2 */
    {"R_scip_model_create",               (DL_FUNC) &R_scip_model_create,                1},
    {"R_scip_model_add_var",              (DL_FUNC) &R_scip_model_add_var,               6},
    {"R_scip_model_add_vars",             (DL_FUNC) &R_scip_model_add_vars,              6},
    {"R_scip_model_add_linear_cons",      (DL_FUNC) &R_scip_model_add_linear_cons,       6},
    {"R_scip_model_add_quadratic_cons",   (DL_FUNC) &R_scip_model_add_quadratic_cons,    9},
    {"R_scip_model_add_sos1_cons",        (DL_FUNC) &R_scip_model_add_sos1_cons,         4},
    {"R_scip_model_add_sos2_cons",        (DL_FUNC) &R_scip_model_add_sos2_cons,         4},
    {"R_scip_model_add_indicator_cons",   (DL_FUNC) &R_scip_model_add_indicator_cons,    6},
    {"R_scip_model_set_param",            (DL_FUNC) &R_scip_model_set_param,             3},
    {"R_scip_model_set_objective_sense",  (DL_FUNC) &R_scip_model_set_objective_sense,   2},
    {"R_scip_model_optimize",             (DL_FUNC) &R_scip_model_optimize,              1},
    {"R_scip_model_get_status",           (DL_FUNC) &R_scip_model_get_status,            1},
    {"R_scip_model_get_solution",         (DL_FUNC) &R_scip_model_get_solution,          1},
    {"R_scip_model_get_nsols",            (DL_FUNC) &R_scip_model_get_nsols,             1},
    {"R_scip_model_get_sol",              (DL_FUNC) &R_scip_model_get_sol,               2},
    {"R_scip_model_get_info",             (DL_FUNC) &R_scip_model_get_info,              1},
    {"R_scip_model_get_nvars",            (DL_FUNC) &R_scip_model_get_nvars,             1},
    {"R_scip_model_get_nconss",           (DL_FUNC) &R_scip_model_get_nconss,            1},
    {"R_scip_model_free",                 (DL_FUNC) &R_scip_model_free,                  1},
    {NULL, NULL, 0}
};

void R_init_scip(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
