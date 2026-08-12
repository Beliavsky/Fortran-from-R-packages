#include <R.h>
#include <Rinternals.h>
#include <scip/scip.h>
#include <scip/scipdefplugins.h>

/* =====================================================================
 * SCIP message handler routed through R's I/O
 *
 * SCIP has a built-in message handler system for embedding. We use it
 * to route all SCIP output through Rprintf/REprintf instead of
 * stdout/stderr. This is the CRAN-compliant approach — no preprocessor
 * hacks, no fake FILE* streams, no -include flags.
 * ===================================================================== */

static SCIP_DECL_MESSAGEWARNING(r_message_warning) {
    if (msg != NULL) REprintf("%s", msg);
}

static SCIP_DECL_MESSAGEDIALOG(r_message_dialog) {
    if (msg != NULL) Rprintf("%s", msg);
}

static SCIP_DECL_MESSAGEINFO(r_message_info) {
    if (msg != NULL) Rprintf("%s", msg);
}

static SCIP_DECL_ERRORPRINTING(r_error_printer) {
    if (msg != NULL) REprintf("%s", msg);
}

/* Install our message handler on a SCIP instance.
 * Call right after SCIPcreate(), before SCIPincludeDefaultPlugins(). */
static void install_r_message_handler(SCIP *scip) {
    SCIP_MESSAGEHDLR *hdlr = NULL;

    /* Route error messages through REprintf (global, not per-instance) */
    SCIPmessageSetErrorPrinting(r_error_printer, NULL);

    /* Create and install a message handler for info/warning/dialog */
    SCIPmessagehdlrCreate(&hdlr, TRUE, NULL, FALSE,
                          r_message_warning,
                          r_message_dialog,
                          r_message_info,
                          NULL, NULL);
    SCIPsetMessagehdlr(scip, hdlr);
    /* Release our reference; SCIP now owns it */
    SCIPmessagehdlrRelease(&hdlr);
}

/* Status string from SCIP status code */
static const char* scip_status_string(SCIP_STATUS status) {
    switch (status) {
        case SCIP_STATUS_OPTIMAL:        return "optimal";
        case SCIP_STATUS_INFEASIBLE:     return "infeasible";
        case SCIP_STATUS_UNBOUNDED:      return "unbounded";
        case SCIP_STATUS_INFORUNBD:      return "infeasible_or_unbounded";
        case SCIP_STATUS_TIMELIMIT:      return "timelimit";
        case SCIP_STATUS_NODELIMIT:      return "nodelimit";
        case SCIP_STATUS_SOLLIMIT:       return "sollimit";
        case SCIP_STATUS_GAPLIMIT:       return "gaplimit";
        case SCIP_STATUS_MEMLIMIT:       return "memlimit";
        case SCIP_STATUS_STALLNODELIMIT: return "stallnodelimit";
        case SCIP_STATUS_USERINTERRUPT:  return "userinterrupt";
        default:                         return "unknown";
    }
}

/* Macro for SCIP call with error handling - use our own name to avoid clash with SCIP's macro */
#undef SCIP_CALL_R
#define SCIP_CALL_R(x) do { \
    SCIP_RETCODE _retcode = (x); \
    if (_retcode != SCIP_OKAY) { \
        error("SCIP error %d in %s at line %d", _retcode, __FILE__, __LINE__); \
    } \
} while(0)

/* Get a named element from an R list (NULL if not found) */
static SEXP getListElement(SEXP list, const char *name) {
    SEXP names = getAttrib(list, R_NamesSymbol);
    if (names == R_NilValue) return R_NilValue;
    for (int i = 0; i < length(list); i++) {
        if (strcmp(CHAR(STRING_ELT(names, i)), name) == 0)
            return VECTOR_ELT(list, i);
    }
    return R_NilValue;
}

/* =====================================================================
 * Model-building API (Layer 2)
 *
 * Uses R external pointers to hold SCIP model state across calls.
 * ===================================================================== */

/* Model state held in external pointer */
typedef struct {
    SCIP       *scip;
    SCIP_VAR  **vars;
    SCIP_CONS **conss;
    int         nvars;
    int         nconss;
    int         vars_cap;
    int         conss_cap;
} ScipModel;

#define MODEL_INIT_CAP 64

static void scip_model_finalizer(SEXP ext) {
    ScipModel *model = (ScipModel *)R_ExternalPtrAddr(ext);
    if (model == NULL) return;

    if (model->scip != NULL) {
        for (int i = 0; i < model->nconss; i++) {
            SCIPreleaseCons(model->scip, &model->conss[i]);
        }
        for (int j = 0; j < model->nvars; j++) {
            SCIPreleaseVar(model->scip, &model->vars[j]);
        }
        SCIPfree(&model->scip);
    }
    if (model->vars) free(model->vars);
    if (model->conss) free(model->conss);
    free(model);
    R_ClearExternalPtr(ext);
}

static ScipModel* get_model(SEXP ext) {
    ScipModel *model = (ScipModel *)R_ExternalPtrAddr(ext);
    if (model == NULL)
        error("SCIP model has been freed");
    return model;
}

static void ensure_var_cap(ScipModel *model) {
    if (model->nvars >= model->vars_cap) {
        model->vars_cap *= 2;
        model->vars = (SCIP_VAR **)realloc(model->vars,
            model->vars_cap * sizeof(SCIP_VAR *));
        if (!model->vars) error("Failed to allocate variable array");
    }
}

static void ensure_cons_cap(ScipModel *model) {
    if (model->nconss >= model->conss_cap) {
        model->conss_cap *= 2;
        model->conss = (SCIP_CONS **)realloc(model->conss,
            model->conss_cap * sizeof(SCIP_CONS *));
        if (!model->conss) error("Failed to allocate constraint array");
    }
}

/* Create a new SCIP model */
SEXP R_scip_model_create(SEXP s_name) {
    const char *name = CHAR(STRING_ELT(s_name, 0));

    ScipModel *model = (ScipModel *)calloc(1, sizeof(ScipModel));
    if (!model) error("Failed to allocate ScipModel");

    SCIP_CALL_R(SCIPcreate(&model->scip));
    install_r_message_handler(model->scip);
    SCIP_CALL_R(SCIPincludeDefaultPlugins(model->scip));
    SCIP_CALL_R(SCIPcreateProbBasic(model->scip, name));
    /* Silence output by default; user can override via set_param */
    SCIP_CALL_R(SCIPsetIntParam(model->scip, "display/verblevel", 0));

    model->vars_cap = MODEL_INIT_CAP;
    model->conss_cap = MODEL_INIT_CAP;
    model->vars = (SCIP_VAR **)calloc(model->vars_cap, sizeof(SCIP_VAR *));
    model->conss = (SCIP_CONS **)calloc(model->conss_cap, sizeof(SCIP_CONS *));
    if (!model->vars || !model->conss) {
        if (model->vars) free(model->vars);
        if (model->conss) free(model->conss);
        SCIPfree(&model->scip);
        free(model);
        error("Failed to allocate arrays");
    }

    SEXP ext = PROTECT(R_MakeExternalPtr(model, R_NilValue, R_NilValue));
    R_RegisterCFinalizerEx(ext, scip_model_finalizer, TRUE);
    UNPROTECT(1);
    return ext;
}

/* Add a variable. Returns 1-based variable index. */
SEXP R_scip_model_add_var(SEXP ext, SEXP s_obj, SEXP s_lb, SEXP s_ub,
                           SEXP s_vtype, SEXP s_name) {
    ScipModel *model = get_model(ext);
    double obj_coef = asReal(s_obj);
    double lb_val = asReal(s_lb);
    double ub_val = asReal(s_ub);
    const char *vt = CHAR(STRING_ELT(s_vtype, 0));
    const char *name = CHAR(STRING_ELT(s_name, 0));

    SCIP_VARTYPE vartype;
    if (vt[0] == 'B')      vartype = SCIP_VARTYPE_BINARY;
    else if (vt[0] == 'I') vartype = SCIP_VARTYPE_INTEGER;
    else                    vartype = SCIP_VARTYPE_CONTINUOUS;

    ensure_var_cap(model);
    int idx = model->nvars;
    SCIP_CALL_R(SCIPcreateVarBasic(model->scip, &model->vars[idx],
                                    name, lb_val, ub_val, obj_coef, vartype));
    SCIP_CALL_R(SCIPaddVar(model->scip, model->vars[idx]));
    model->nvars++;

    return ScalarInteger(idx + 1);  /* 1-based */
}

/* Add multiple variables at once. Returns 1-based index of first var added. */
SEXP R_scip_model_add_vars(SEXP ext, SEXP s_obj, SEXP s_lb, SEXP s_ub,
                            SEXP s_vtype, SEXP s_names) {
    ScipModel *model = get_model(ext);
    int nv = length(s_obj);
    double *obj = REAL(s_obj);
    double *lb = REAL(s_lb);
    double *ub = REAL(s_ub);

    int first_idx = model->nvars + 1;  /* 1-based */

    for (int j = 0; j < nv; j++) {
        const char *vt = CHAR(STRING_ELT(s_vtype, j));
        const char *nm = CHAR(STRING_ELT(s_names, j));

        SCIP_VARTYPE vartype;
        if (vt[0] == 'B')      vartype = SCIP_VARTYPE_BINARY;
        else if (vt[0] == 'I') vartype = SCIP_VARTYPE_INTEGER;
        else                    vartype = SCIP_VARTYPE_CONTINUOUS;

        ensure_var_cap(model);
        int idx = model->nvars;
        SCIP_CALL_R(SCIPcreateVarBasic(model->scip, &model->vars[idx],
                                        nm, lb[j], ub[j], obj[j], vartype));
        SCIP_CALL_R(SCIPaddVar(model->scip, model->vars[idx]));
        model->nvars++;
    }

    return ScalarInteger(first_idx);
}

/* Add a linear constraint: lhs <= sum(coefs * vars) <= rhs
 * s_vars: integer vector of 1-based variable indices
 * s_coefs: double vector of coefficients */
SEXP R_scip_model_add_linear_cons(SEXP ext, SEXP s_vars, SEXP s_coefs,
                                   SEXP s_lhs, SEXP s_rhs, SEXP s_name) {
    ScipModel *model = get_model(ext);
    int nv = length(s_vars);
    int *var_idx = INTEGER(s_vars);
    double *coefs = REAL(s_coefs);
    double lhs_val = asReal(s_lhs);
    double rhs_val = asReal(s_rhs);
    const char *name = CHAR(STRING_ELT(s_name, 0));

    SCIP_VAR **cvars = (SCIP_VAR **)R_alloc(nv, sizeof(SCIP_VAR *));
    for (int k = 0; k < nv; k++) {
        int idx = var_idx[k] - 1;  /* convert to 0-based */
        if (idx < 0 || idx >= model->nvars)
            error("Variable index %d out of range [1, %d]", var_idx[k], model->nvars);
        cvars[k] = model->vars[idx];
    }

    ensure_cons_cap(model);
    int cidx = model->nconss;
    SCIP_CALL_R(SCIPcreateConsBasicLinear(model->scip, &model->conss[cidx],
                                          name, nv, cvars, coefs,
                                          lhs_val, rhs_val));
    SCIP_CALL_R(SCIPaddCons(model->scip, model->conss[cidx]));
    model->nconss++;

    return ScalarInteger(cidx + 1);  /* 1-based */
}

/* Add a quadratic constraint: lhs <= bilinexpr + linexpr <= rhs
 * Linear part: s_linvars (1-based), s_lincoefs
 * Quadratic part: s_quadvars1, s_quadvars2 (1-based), s_quadcoefs
 *   represents sum of quadcoefs[k] * quadvars1[k] * quadvars2[k] */
SEXP R_scip_model_add_quadratic_cons(SEXP ext,
                                      SEXP s_linvars, SEXP s_lincoefs,
                                      SEXP s_quadvars1, SEXP s_quadvars2,
                                      SEXP s_quadcoefs,
                                      SEXP s_lhs, SEXP s_rhs, SEXP s_name) {
    ScipModel *model = get_model(ext);
    int nlin = length(s_linvars);
    int nquad = length(s_quadvars1);
    double lhs_val = asReal(s_lhs);
    double rhs_val = asReal(s_rhs);
    const char *name = CHAR(STRING_ELT(s_name, 0));

    /* Build the expression: sum of quadratic terms + linear terms */
    /* We use SCIPcreateConsBasicQuadraticNonlinear for SCIP >= 8 */

    /* Collect linear variables */
    SCIP_VAR **linvars = NULL;
    double *lincoefs = NULL;
    if (nlin > 0) {
        linvars = (SCIP_VAR **)R_alloc(nlin, sizeof(SCIP_VAR *));
        lincoefs = REAL(s_lincoefs);
        int *lv = INTEGER(s_linvars);
        for (int k = 0; k < nlin; k++) {
            int idx = lv[k] - 1;
            if (idx < 0 || idx >= model->nvars)
                error("Variable index %d out of range", lv[k]);
            linvars[k] = model->vars[idx];
        }
    }

    /* Collect quadratic variables */
    SCIP_VAR **qvars1 = NULL, **qvars2 = NULL;
    double *qcoefs = NULL;
    if (nquad > 0) {
        qvars1 = (SCIP_VAR **)R_alloc(nquad, sizeof(SCIP_VAR *));
        qvars2 = (SCIP_VAR **)R_alloc(nquad, sizeof(SCIP_VAR *));
        qcoefs = REAL(s_quadcoefs);
        int *qv1 = INTEGER(s_quadvars1);
        int *qv2 = INTEGER(s_quadvars2);
        for (int k = 0; k < nquad; k++) {
            int i1 = qv1[k] - 1, i2 = qv2[k] - 1;
            if (i1 < 0 || i1 >= model->nvars || i2 < 0 || i2 >= model->nvars)
                error("Quadratic variable index out of range");
            qvars1[k] = model->vars[i1];
            qvars2[k] = model->vars[i2];
        }
    }

    ensure_cons_cap(model);
    int cidx = model->nconss;
    SCIP_CALL_R(SCIPcreateConsBasicQuadraticNonlinear(model->scip,
        &model->conss[cidx], name,
        nlin, linvars, lincoefs,
        nquad, qvars1, qvars2, qcoefs,
        lhs_val, rhs_val));
    SCIP_CALL_R(SCIPaddCons(model->scip, model->conss[cidx]));
    model->nconss++;

    return ScalarInteger(cidx + 1);
}

/* Add SOS1 constraint: at most one variable in set is nonzero */
SEXP R_scip_model_add_sos1_cons(SEXP ext, SEXP s_vars, SEXP s_weights,
                                 SEXP s_name) {
    ScipModel *model = get_model(ext);
    int nv = length(s_vars);
    int *var_idx = INTEGER(s_vars);
    double *weights = REAL(s_weights);
    const char *name = CHAR(STRING_ELT(s_name, 0));

    SCIP_VAR **cvars = (SCIP_VAR **)R_alloc(nv, sizeof(SCIP_VAR *));
    for (int k = 0; k < nv; k++) {
        int idx = var_idx[k] - 1;
        if (idx < 0 || idx >= model->nvars)
            error("Variable index %d out of range", var_idx[k]);
        cvars[k] = model->vars[idx];
    }

    ensure_cons_cap(model);
    int cidx = model->nconss;
    SCIP_CALL_R(SCIPcreateConsBasicSOS1(model->scip, &model->conss[cidx],
                                         name, nv, cvars, weights));
    SCIP_CALL_R(SCIPaddCons(model->scip, model->conss[cidx]));
    model->nconss++;

    return ScalarInteger(cidx + 1);
}

/* Add SOS2 constraint: at most two adjacent variables in set are nonzero */
SEXP R_scip_model_add_sos2_cons(SEXP ext, SEXP s_vars, SEXP s_weights,
                                 SEXP s_name) {
    ScipModel *model = get_model(ext);
    int nv = length(s_vars);
    int *var_idx = INTEGER(s_vars);
    double *weights = REAL(s_weights);
    const char *name = CHAR(STRING_ELT(s_name, 0));

    SCIP_VAR **cvars = (SCIP_VAR **)R_alloc(nv, sizeof(SCIP_VAR *));
    for (int k = 0; k < nv; k++) {
        int idx = var_idx[k] - 1;
        if (idx < 0 || idx >= model->nvars)
            error("Variable index %d out of range", var_idx[k]);
        cvars[k] = model->vars[idx];
    }

    ensure_cons_cap(model);
    int cidx = model->nconss;
    SCIP_CALL_R(SCIPcreateConsBasicSOS2(model->scip, &model->conss[cidx],
                                         name, nv, cvars, weights));
    SCIP_CALL_R(SCIPaddCons(model->scip, model->conss[cidx]));
    model->nconss++;

    return ScalarInteger(cidx + 1);
}

/* Add indicator constraint: binvar = 1 => sum(coefs * vars) <= rhs */
SEXP R_scip_model_add_indicator_cons(SEXP ext, SEXP s_binvar,
                                      SEXP s_vars, SEXP s_coefs,
                                      SEXP s_rhs, SEXP s_name) {
    ScipModel *model = get_model(ext);
    int bin_idx = asInteger(s_binvar) - 1;
    int nv = length(s_vars);
    int *var_idx = INTEGER(s_vars);
    double *coefs = REAL(s_coefs);
    double rhs_val = asReal(s_rhs);
    const char *name = CHAR(STRING_ELT(s_name, 0));

    if (bin_idx < 0 || bin_idx >= model->nvars)
        error("Binary variable index %d out of range", asInteger(s_binvar));

    SCIP_VAR **cvars = (SCIP_VAR **)R_alloc(nv, sizeof(SCIP_VAR *));
    for (int k = 0; k < nv; k++) {
        int idx = var_idx[k] - 1;
        if (idx < 0 || idx >= model->nvars)
            error("Variable index %d out of range", var_idx[k]);
        cvars[k] = model->vars[idx];
    }

    ensure_cons_cap(model);
    int cidx = model->nconss;
    SCIP_CALL_R(SCIPcreateConsBasicIndicator(model->scip, &model->conss[cidx],
                                              name, model->vars[bin_idx],
                                              nv, cvars, coefs, rhs_val));
    SCIP_CALL_R(SCIPaddCons(model->scip, model->conss[cidx]));
    model->nconss++;

    return ScalarInteger(cidx + 1);
}

/* Set a SCIP parameter */
SEXP R_scip_model_set_param(SEXP ext, SEXP s_name, SEXP s_value) {
    ScipModel *model = get_model(ext);
    const char *pname = CHAR(STRING_ELT(s_name, 0));

    SCIP_PARAM *param = SCIPgetParam(model->scip, pname);
    if (param == NULL)
        error("Unknown SCIP parameter: %s", pname);

    switch (SCIPparamGetType(param)) {
        case SCIP_PARAMTYPE_BOOL:
            SCIP_CALL_R(SCIPsetBoolParam(model->scip, pname,
                (SCIP_Bool)asLogical(s_value)));
            break;
        case SCIP_PARAMTYPE_INT:
            SCIP_CALL_R(SCIPsetIntParam(model->scip, pname,
                asInteger(s_value)));
            break;
        case SCIP_PARAMTYPE_LONGINT:
            SCIP_CALL_R(SCIPsetLongintParam(model->scip, pname,
                (SCIP_Longint)asReal(s_value)));
            break;
        case SCIP_PARAMTYPE_REAL:
            SCIP_CALL_R(SCIPsetRealParam(model->scip, pname,
                asReal(s_value)));
            break;
        case SCIP_PARAMTYPE_CHAR:
            SCIP_CALL_R(SCIPsetCharParam(model->scip, pname,
                CHAR(asChar(s_value))[0]));
            break;
        case SCIP_PARAMTYPE_STRING:
            SCIP_CALL_R(SCIPsetStringParam(model->scip, pname,
                CHAR(asChar(s_value))));
            break;
    }
    return R_NilValue;
}

/* Set objective sense: "minimize" or "maximize" */
SEXP R_scip_model_set_objective_sense(SEXP ext, SEXP s_sense) {
    ScipModel *model = get_model(ext);
    const char *sense = CHAR(STRING_ELT(s_sense, 0));

    if (strcmp(sense, "maximize") == 0 || strcmp(sense, "max") == 0) {
        SCIP_CALL_R(SCIPsetObjsense(model->scip, SCIP_OBJSENSE_MAXIMIZE));
    } else {
        SCIP_CALL_R(SCIPsetObjsense(model->scip, SCIP_OBJSENSE_MINIMIZE));
    }
    return R_NilValue;
}

/* Solve the model */
SEXP R_scip_model_optimize(SEXP ext) {
    ScipModel *model = get_model(ext);
    SCIP_CALL_R(SCIPsolve(model->scip));
    return R_NilValue;
}

/* Get solver status */
SEXP R_scip_model_get_status(SEXP ext) {
    ScipModel *model = get_model(ext);
    SCIP_STATUS status = SCIPgetStatus(model->scip);
    return mkString(scip_status_string(status));
}

/* Get primal solution and objective value */
SEXP R_scip_model_get_solution(SEXP ext) {
    ScipModel *model = get_model(ext);
    int nsols = SCIPgetNSols(model->scip);

    SEXP result = PROTECT(allocVector(VECSXP, 2));
    SEXP result_names = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(result_names, 0, mkChar("objval"));
    SET_STRING_ELT(result_names, 1, mkChar("x"));
    setAttrib(result, R_NamesSymbol, result_names);

    if (nsols > 0) {
        SCIP_SOL *sol = SCIPgetBestSol(model->scip);
        SET_VECTOR_ELT(result, 0,
            ScalarReal(SCIPgetSolOrigObj(model->scip, sol)));

        SEXP r_x = PROTECT(allocVector(REALSXP, model->nvars));
        double *px = REAL(r_x);
        for (int j = 0; j < model->nvars; j++) {
            px[j] = SCIPgetSolVal(model->scip, sol, model->vars[j]);
        }
        SET_VECTOR_ELT(result, 1, r_x);
        UNPROTECT(1);
    } else {
        SET_VECTOR_ELT(result, 0, ScalarReal(NA_REAL));
        SET_VECTOR_ELT(result, 1, R_NilValue);
    }

    UNPROTECT(2);
    return result;
}

/* Get number of solutions */
SEXP R_scip_model_get_nsols(SEXP ext) {
    ScipModel *model = get_model(ext);
    return ScalarInteger(SCIPgetNSols(model->scip));
}

/* Get k-th solution (1-based) */
SEXP R_scip_model_get_sol(SEXP ext, SEXP s_k) {
    ScipModel *model = get_model(ext);
    int k = asInteger(s_k) - 1;  /* to 0-based */
    int nsols = SCIPgetNSols(model->scip);

    if (k < 0 || k >= nsols)
        error("Solution index %d out of range [1, %d]", k + 1, nsols);

    SCIP_SOL **sols = SCIPgetSols(model->scip);
    SCIP_SOL *sol = sols[k];

    SEXP result = PROTECT(allocVector(VECSXP, 2));
    SEXP result_names = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(result_names, 0, mkChar("objval"));
    SET_STRING_ELT(result_names, 1, mkChar("x"));
    setAttrib(result, R_NamesSymbol, result_names);

    SET_VECTOR_ELT(result, 0,
        ScalarReal(SCIPgetSolOrigObj(model->scip, sol)));

    SEXP r_x = PROTECT(allocVector(REALSXP, model->nvars));
    double *px = REAL(r_x);
    for (int j = 0; j < model->nvars; j++) {
        px[j] = SCIPgetSolVal(model->scip, sol, model->vars[j]);
    }
    SET_VECTOR_ELT(result, 1, r_x);

    UNPROTECT(3);
    return result;
}

/* Get solving statistics */
SEXP R_scip_model_get_info(SEXP ext) {
    ScipModel *model = get_model(ext);

    SEXP info = PROTECT(allocVector(VECSXP, 5));
    SEXP info_names = PROTECT(allocVector(STRSXP, 5));
    SET_STRING_ELT(info_names, 0, mkChar("solve_time"));
    SET_STRING_ELT(info_names, 1, mkChar("nodes"));
    SET_STRING_ELT(info_names, 2, mkChar("iterations"));
    SET_STRING_ELT(info_names, 3, mkChar("gap"));
    SET_STRING_ELT(info_names, 4, mkChar("sol_count"));
    setAttrib(info, R_NamesSymbol, info_names);

    SET_VECTOR_ELT(info, 0, ScalarReal(SCIPgetSolvingTime(model->scip)));
    SET_VECTOR_ELT(info, 1, ScalarReal((double)SCIPgetNNodes(model->scip)));
    SET_VECTOR_ELT(info, 2, ScalarReal((double)SCIPgetNLPIterations(model->scip)));

    int nsols = SCIPgetNSols(model->scip);
    SCIP_STATUS status = SCIPgetStatus(model->scip);
    double gap = (nsols > 0 && status != SCIP_STATUS_OPTIMAL) ?
                  SCIPgetGap(model->scip) : 0.0;
    SET_VECTOR_ELT(info, 3, ScalarReal(gap));
    SET_VECTOR_ELT(info, 4, ScalarInteger(nsols));

    UNPROTECT(2);
    return info;
}

/* Get number of variables */
SEXP R_scip_model_get_nvars(SEXP ext) {
    ScipModel *model = get_model(ext);
    return ScalarInteger(model->nvars);
}

/* Get number of constraints */
SEXP R_scip_model_get_nconss(SEXP ext) {
    ScipModel *model = get_model(ext);
    return ScalarInteger(model->nconss);
}

/* Free model explicitly */
SEXP R_scip_model_free(SEXP ext) {
    scip_model_finalizer(ext);
    return R_NilValue;
}

/* =====================================================================
 * One-shot solver (Layer 1)
 * ===================================================================== */

/*
 * One-shot solver: R_scip_solve
 *
 * Arguments:
 *   obj   - double[n]: objective coefficients
 *   Ai    - int[nnz]: row indices (0-based) of CSC constraint matrix
 *   Ap    - int[n+1]: column pointers of CSC constraint matrix
 *   Ax    - double[nnz]: values of CSC constraint matrix
 *   m     - int: number of constraints
 *   n     - int: number of variables
 *   lhs   - double[m]: left-hand side of constraints
 *   rhs   - double[m]: right-hand side of constraints
 *   vtype - character[n]: variable types ("C", "B", "I")
 *   lb    - double[n]: variable lower bounds
 *   ub    - double[n]: variable upper bounds
 *   ctrl  - list: control parameters
 */
SEXP R_scip_solve(SEXP obj, SEXP Ai, SEXP Ap, SEXP Ax,
                  SEXP s_m, SEXP s_n,
                  SEXP lhs, SEXP rhs,
                  SEXP vtype, SEXP lb, SEXP ub,
                  SEXP ctrl) {
    int n = INTEGER(s_n)[0];
    int m = INTEGER(s_m)[0];
    double *c_obj = REAL(obj);
    int *c_Ai = INTEGER(Ai);
    int *c_Ap = INTEGER(Ap);
    double *c_Ax = REAL(Ax);
    double *c_lhs = REAL(lhs);
    double *c_rhs = REAL(rhs);
    double *c_lb = REAL(lb);
    double *c_ub = REAL(ub);

    SCIP *scip = NULL;
    SCIP_VAR **vars = NULL;
    SCIP_CONS **conss = NULL;

    /* Create SCIP instance */
    SCIP_CALL_R(SCIPcreate(&scip));
    install_r_message_handler(scip);
    SCIP_CALL_R(SCIPincludeDefaultPlugins(scip));

    /* Apply control parameters from the scip_control structure.
     * The R layer produces: list(verbose=T/F, scip_params=list(...),
     *   heuristics_emphasis=...) */

    /* Handle heuristics emphasis (must be set before individual params) */
    SEXP hem = getListElement(ctrl, "heuristics_emphasis");
    if (hem != R_NilValue) {
        const char *emphasis = CHAR(STRING_ELT(hem, 0));
        if (strcmp(emphasis, "aggressive") == 0) {
            SCIP_CALL_R(SCIPsetHeuristics(scip, SCIP_PARAMSETTING_AGGRESSIVE, TRUE));
        } else if (strcmp(emphasis, "fast") == 0) {
            SCIP_CALL_R(SCIPsetHeuristics(scip, SCIP_PARAMSETTING_FAST, TRUE));
        } else if (strcmp(emphasis, "off") == 0) {
            SCIP_CALL_R(SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE));
        }
    }

    /* Apply all native SCIP parameters from scip_params list */
    SEXP params = getListElement(ctrl, "scip_params");
    if (params != R_NilValue && length(params) > 0) {
        SEXP param_names = getAttrib(params, R_NamesSymbol);
        for (int j = 0; j < length(params); j++) {
            const char *pname = CHAR(STRING_ELT(param_names, j));
            SEXP pval = VECTOR_ELT(params, j);
            SCIP_PARAM *param = SCIPgetParam(scip, pname);
            if (param == NULL) {
                warning("Unknown SCIP parameter: %s", pname);
                continue;
            }
            switch (SCIPparamGetType(param)) {
                case SCIP_PARAMTYPE_BOOL:
                    SCIP_CALL_R(SCIPsetBoolParam(scip, pname,
                        (SCIP_Bool)asLogical(pval)));
                    break;
                case SCIP_PARAMTYPE_INT:
                    SCIP_CALL_R(SCIPsetIntParam(scip, pname,
                        asInteger(pval)));
                    break;
                case SCIP_PARAMTYPE_LONGINT:
                    SCIP_CALL_R(SCIPsetLongintParam(scip, pname,
                        (SCIP_Longint)asReal(pval)));
                    break;
                case SCIP_PARAMTYPE_REAL:
                    SCIP_CALL_R(SCIPsetRealParam(scip, pname,
                        asReal(pval)));
                    break;
                case SCIP_PARAMTYPE_CHAR:
                    SCIP_CALL_R(SCIPsetCharParam(scip, pname,
                        CHAR(asChar(pval))[0]));
                    break;
                case SCIP_PARAMTYPE_STRING:
                    SCIP_CALL_R(SCIPsetStringParam(scip, pname,
                        CHAR(asChar(pval))));
                    break;
            }
        }
    }

    /* Create problem */
    SCIP_CALL_R(SCIPcreateProbBasic(scip, "scip_r"));

    /* Allocate variable and constraint arrays */
    vars = (SCIP_VAR **)R_alloc(n, sizeof(SCIP_VAR *));
    conss = (SCIP_CONS **)R_alloc(m, sizeof(SCIP_CONS *));

    /* Create variables */
    for (int j = 0; j < n; j++) {
        SCIP_VARTYPE vartype;
        const char *vt = CHAR(STRING_ELT(vtype, j));
        if (vt[0] == 'B') {
            vartype = SCIP_VARTYPE_BINARY;
        } else if (vt[0] == 'I') {
            vartype = SCIP_VARTYPE_INTEGER;
        } else {
            vartype = SCIP_VARTYPE_CONTINUOUS;
        }
        char varname[32];
        snprintf(varname, sizeof(varname), "x%d", j);
        SCIP_CALL_R(SCIPcreateVarBasic(scip, &vars[j], varname,
                                            c_lb[j], c_ub[j], c_obj[j],
                                            vartype));
        SCIP_CALL_R(SCIPaddVar(scip, vars[j]));
    }

    /* Create constraints from CSC matrix */
    /* First, build row-wise representation */
    /* Count entries per row */
    int *row_count = (int *)R_alloc(m, sizeof(int));
    memset(row_count, 0, m * sizeof(int));
    int nnz = c_Ap[n];
    for (int k = 0; k < nnz; k++) {
        row_count[c_Ai[k]]++;
    }

    /* Row pointers and storage */
    int *row_ptr = (int *)R_alloc(m + 1, sizeof(int));
    row_ptr[0] = 0;
    for (int i = 0; i < m; i++) {
        row_ptr[i + 1] = row_ptr[i] + row_count[i];
    }

    int *row_col = (int *)R_alloc(nnz, sizeof(int));
    double *row_val = (double *)R_alloc(nnz, sizeof(double));
    int *row_pos = (int *)R_alloc(m, sizeof(int));
    memcpy(row_pos, row_ptr, m * sizeof(int));

    /* Fill row-wise data from CSC */
    for (int j = 0; j < n; j++) {
        for (int k = c_Ap[j]; k < c_Ap[j + 1]; k++) {
            int i = c_Ai[k];
            int pos = row_pos[i]++;
            row_col[pos] = j;
            row_val[pos] = c_Ax[k];
        }
    }

    /* Create linear constraints */
    for (int i = 0; i < m; i++) {
        int nrow_vars = row_count[i];
        SCIP_VAR **con_vars = (SCIP_VAR **)R_alloc(nrow_vars, sizeof(SCIP_VAR *));
        double *con_vals = (double *)R_alloc(nrow_vars, sizeof(double));
        for (int k = 0; k < nrow_vars; k++) {
            con_vars[k] = vars[row_col[row_ptr[i] + k]];
            con_vals[k] = row_val[row_ptr[i] + k];
        }
        char consname[32];
        snprintf(consname, sizeof(consname), "c%d", i);
        SCIP_CALL_R(SCIPcreateConsBasicLinear(scip, &conss[i], consname,
                                                   nrow_vars, con_vars, con_vals,
                                                   c_lhs[i], c_rhs[i]));
        SCIP_CALL_R(SCIPaddCons(scip, conss[i]));
    }

    /* Solve */
    SCIP_CALL_R(SCIPsolve(scip));

    /* Extract results */
    SCIP_STATUS status = SCIPgetStatus(scip);
    const char *status_str = scip_status_string(status);
    int nsols = SCIPgetNSols(scip);

    /* Build result list */
    SEXP result = PROTECT(allocVector(VECSXP, 6));
    SEXP result_names = PROTECT(allocVector(STRSXP, 6));
    SET_STRING_ELT(result_names, 0, mkChar("status"));
    SET_STRING_ELT(result_names, 1, mkChar("objval"));
    SET_STRING_ELT(result_names, 2, mkChar("x"));
    SET_STRING_ELT(result_names, 3, mkChar("sol_count"));
    SET_STRING_ELT(result_names, 4, mkChar("gap"));
    SET_STRING_ELT(result_names, 5, mkChar("info"));
    setAttrib(result, R_NamesSymbol, result_names);

    /* Status */
    SET_VECTOR_ELT(result, 0, mkString(status_str));

    /* Solution values */
    if (nsols > 0) {
        SCIP_SOL *sol = SCIPgetBestSol(scip);

        /* Objective value */
        SEXP r_objval = PROTECT(ScalarReal(SCIPgetSolOrigObj(scip, sol)));
        SET_VECTOR_ELT(result, 1, r_objval);
        UNPROTECT(1);

        /* Primal solution */
        SEXP r_x = PROTECT(allocVector(REALSXP, n));
        double *px = REAL(r_x);
        for (int j = 0; j < n; j++) {
            px[j] = SCIPgetSolVal(scip, sol, vars[j]);
        }
        SET_VECTOR_ELT(result, 2, r_x);
        UNPROTECT(1);
    } else {
        SET_VECTOR_ELT(result, 1, ScalarReal(NA_REAL));
        SET_VECTOR_ELT(result, 2, R_NilValue);
    }

    /* Solution count */
    SET_VECTOR_ELT(result, 3, ScalarInteger(nsols));

    /* Gap */
    double gap = (nsols > 0 && status != SCIP_STATUS_OPTIMAL) ?
                  SCIPgetGap(scip) : 0.0;
    SET_VECTOR_ELT(result, 4, ScalarReal(gap));

    /* Info sub-list */
    SEXP info = PROTECT(allocVector(VECSXP, 3));
    SEXP info_names = PROTECT(allocVector(STRSXP, 3));
    SET_STRING_ELT(info_names, 0, mkChar("solve_time"));
    SET_STRING_ELT(info_names, 1, mkChar("nodes"));
    SET_STRING_ELT(info_names, 2, mkChar("iterations"));
    setAttrib(info, R_NamesSymbol, info_names);
    SET_VECTOR_ELT(info, 0, ScalarReal(SCIPgetSolvingTime(scip)));
    SET_VECTOR_ELT(info, 1, ScalarReal((double)SCIPgetNNodes(scip)));
    SET_VECTOR_ELT(info, 2, ScalarReal((double)SCIPgetNLPIterations(scip)));
    SET_VECTOR_ELT(result, 5, info);
    UNPROTECT(2);

    /* Release constraints and variables */
    for (int i = 0; i < m; i++) {
        SCIP_CALL_R(SCIPreleaseCons(scip, &conss[i]));
    }
    for (int j = 0; j < n; j++) {
        SCIP_CALL_R(SCIPreleaseVar(scip, &vars[j]));
    }

    /* Free SCIP */
    SCIP_CALL_R(SCIPfree(&scip));

    UNPROTECT(2); /* result, result_names */
    return result;
}
