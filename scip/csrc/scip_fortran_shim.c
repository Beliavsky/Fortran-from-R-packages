/*
 * Apache-2.0
 * Plain-C ABI shim for the Fortran translation of the R package scip.
 * The solver implementation remains the vendored SCIP Optimization Suite.
 */
#include <scip/scip.h>
#include <scip/scipdefplugins.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

typedef struct {
    SCIP *scip;
    SCIP_VAR **vars;
    SCIP_CONS **conss;
    int nvars;
    int nconss;
    int vars_cap;
    int conss_cap;
} scipf_model;

enum { SCIPF_ALLOC_ERROR = -1001, SCIPF_BAD_INDEX = -1002,
       SCIPF_NO_SOLUTION = -1003, SCIPF_BAD_ARGUMENT = -1004 };

static int ensure_var_cap(scipf_model *m)
{
    SCIP_VAR **tmp;
    int cap;
    if (m->nvars < m->vars_cap) return 0;
    cap = m->vars_cap > 0 ? 2 * m->vars_cap : 64;
    tmp = (SCIP_VAR **)realloc(m->vars, (size_t)cap * sizeof(*tmp));
    if (tmp == NULL) return SCIPF_ALLOC_ERROR;
    m->vars = tmp;
    m->vars_cap = cap;
    return 0;
}

static int ensure_cons_cap(scipf_model *m)
{
    SCIP_CONS **tmp;
    int cap;
    if (m->nconss < m->conss_cap) return 0;
    cap = m->conss_cap > 0 ? 2 * m->conss_cap : 64;
    tmp = (SCIP_CONS **)realloc(m->conss, (size_t)cap * sizeof(*tmp));
    if (tmp == NULL) return SCIPF_ALLOC_ERROR;
    m->conss = tmp;
    m->conss_cap = cap;
    return 0;
}

static SCIP_VARTYPE vartype_from_char(char vtype)
{
    if (vtype == 'B' || vtype == 'b') return SCIP_VARTYPE_BINARY;
    if (vtype == 'I' || vtype == 'i') return SCIP_VARTYPE_INTEGER;
    return SCIP_VARTYPE_CONTINUOUS;
}

static const char *status_string(SCIP_STATUS status)
{
    switch (status) {
    case SCIP_STATUS_OPTIMAL: return "optimal";
    case SCIP_STATUS_INFEASIBLE: return "infeasible";
    case SCIP_STATUS_UNBOUNDED: return "unbounded";
    case SCIP_STATUS_INFORUNBD: return "infeasible_or_unbounded";
    case SCIP_STATUS_TIMELIMIT: return "timelimit";
    case SCIP_STATUS_NODELIMIT: return "nodelimit";
    case SCIP_STATUS_SOLLIMIT: return "sollimit";
    case SCIP_STATUS_GAPLIMIT: return "gaplimit";
    case SCIP_STATUS_MEMLIMIT: return "memlimit";
    case SCIP_STATUS_STALLNODELIMIT: return "stallnodelimit";
    case SCIP_STATUS_USERINTERRUPT: return "userinterrupt";
    default: return "unknown";
    }
}

int scipf_model_create(const char *name, void **out)
{
    scipf_model *m;
    SCIP_RETCODE rc;
    if (out == NULL) return SCIPF_BAD_ARGUMENT;
    *out = NULL;
    m = (scipf_model *)calloc(1, sizeof(*m));
    if (m == NULL) return SCIPF_ALLOC_ERROR;
    rc = SCIPcreate(&m->scip);
    if (rc != SCIP_OKAY) { free(m); return (int)rc; }
    rc = SCIPincludeDefaultPlugins(m->scip);
    if (rc != SCIP_OKAY) { SCIPfree(&m->scip); free(m); return (int)rc; }
    rc = SCIPcreateProbBasic(m->scip, name != NULL ? name : "scip_model");
    if (rc != SCIP_OKAY) { SCIPfree(&m->scip); free(m); return (int)rc; }
    (void)SCIPsetIntParam(m->scip, "display/verblevel", 0);
    m->vars_cap = 64;
    m->conss_cap = 64;
    m->vars = (SCIP_VAR **)calloc((size_t)m->vars_cap, sizeof(*m->vars));
    m->conss = (SCIP_CONS **)calloc((size_t)m->conss_cap, sizeof(*m->conss));
    if (m->vars == NULL || m->conss == NULL) {
        free(m->vars); free(m->conss); SCIPfree(&m->scip); free(m);
        return SCIPF_ALLOC_ERROR;
    }
    *out = m;
    return 0;
}

int scipf_model_free(void *handle)
{
    scipf_model *m = (scipf_model *)handle;
    int i;
    if (m == NULL) return 0;
    if (m->scip != NULL) {
        for (i = 0; i < m->nconss; ++i)
            if (m->conss[i] != NULL) (void)SCIPreleaseCons(m->scip, &m->conss[i]);
        for (i = 0; i < m->nvars; ++i)
            if (m->vars[i] != NULL) (void)SCIPreleaseVar(m->scip, &m->vars[i]);
        (void)SCIPfree(&m->scip);
    }
    free(m->vars);
    free(m->conss);
    free(m);
    return 0;
}

int scipf_add_var(void *handle, double obj, double lb, double ub,
                  int vtype, const char *name, int *index1)
{
    scipf_model *m = (scipf_model *)handle;
    SCIP_RETCODE rc;
    int idx, er;
    if (m == NULL || index1 == NULL) return SCIPF_BAD_ARGUMENT;
    er = ensure_var_cap(m); if (er != 0) return er;
    idx = m->nvars;
    rc = SCIPcreateVarBasic(m->scip, &m->vars[idx], name, lb, ub, obj,
                            vartype_from_char((char)vtype));
    if (rc != SCIP_OKAY) return (int)rc;
    rc = SCIPaddVar(m->scip, m->vars[idx]);
    if (rc != SCIP_OKAY) { SCIPreleaseVar(m->scip, &m->vars[idx]); return (int)rc; }
    m->nvars++;
    *index1 = idx + 1;
    return 0;
}

static int collect_vars(scipf_model *m, int nv, const int *idx1, SCIP_VAR ***out)
{
    SCIP_VAR **v;
    int k, idx;
    if (nv == 0) { *out = NULL; return 0; }
    v = (SCIP_VAR **)malloc((size_t)nv * sizeof(*v));
    if (v == NULL) return SCIPF_ALLOC_ERROR;
    for (k = 0; k < nv; ++k) {
        idx = idx1[k] - 1;
        if (idx < 0 || idx >= m->nvars) { free(v); return SCIPF_BAD_INDEX; }
        v[k] = m->vars[idx];
    }
    *out = v;
    return 0;
}

int scipf_add_linear_cons(void *handle, int nv, const int *vars1,
                          const double *coefs, double lhs, double rhs,
                          const char *name, int *index1)
{
    scipf_model *m = (scipf_model *)handle;
    SCIP_VAR **vars = NULL;
    SCIP_RETCODE rc;
    int er, idx;
    if (m == NULL || nv < 0 || index1 == NULL) return SCIPF_BAD_ARGUMENT;
    er = collect_vars(m, nv, vars1, &vars); if (er != 0) return er;
    er = ensure_cons_cap(m); if (er != 0) { free(vars); return er; }
    idx = m->nconss;
    rc = SCIPcreateConsBasicLinear(m->scip, &m->conss[idx], name, nv,
                                   vars, (SCIP_Real *)coefs, lhs, rhs);
    free(vars);
    if (rc != SCIP_OKAY) return (int)rc;
    rc = SCIPaddCons(m->scip, m->conss[idx]);
    if (rc != SCIP_OKAY) { SCIPreleaseCons(m->scip, &m->conss[idx]); return (int)rc; }
    m->nconss++;
    *index1 = idx + 1;
    return 0;
}

int scipf_add_quadratic_cons(void *handle,
                             int nlin, const int *linvars1, const double *lincoefs,
                             int nquad, const int *qvars1, const int *qvars2,
                             const double *qcoefs, double lhs, double rhs,
                             const char *name, int *index1)
{
    scipf_model *m = (scipf_model *)handle;
    SCIP_VAR **linvars = NULL, **v1 = NULL, **v2 = NULL;
    SCIP_RETCODE rc;
    int er, idx;
    if (m == NULL || nlin < 0 || nquad < 0 || index1 == NULL) return SCIPF_BAD_ARGUMENT;
    er = collect_vars(m, nlin, linvars1, &linvars); if (er != 0) return er;
    er = collect_vars(m, nquad, qvars1, &v1); if (er != 0) { free(linvars); return er; }
    er = collect_vars(m, nquad, qvars2, &v2);
    if (er != 0) { free(linvars); free(v1); return er; }
    er = ensure_cons_cap(m);
    if (er != 0) { free(linvars); free(v1); free(v2); return er; }
    idx = m->nconss;
    rc = SCIPcreateConsBasicQuadraticNonlinear(m->scip, &m->conss[idx], name,
         nlin, linvars, (SCIP_Real *)lincoefs,
         nquad, v1, v2, (SCIP_Real *)qcoefs, lhs, rhs);
    free(linvars); free(v1); free(v2);
    if (rc != SCIP_OKAY) return (int)rc;
    rc = SCIPaddCons(m->scip, m->conss[idx]);
    if (rc != SCIP_OKAY) { SCIPreleaseCons(m->scip, &m->conss[idx]); return (int)rc; }
    m->nconss++;
    *index1 = idx + 1;
    return 0;
}

static int add_sos(void *handle, int sos2, int nv, const int *vars1,
                   const double *weights, const char *name, int *index1)
{
    scipf_model *m = (scipf_model *)handle;
    SCIP_VAR **vars = NULL;
    SCIP_RETCODE rc;
    int er, idx;
    if (m == NULL || nv < 0 || index1 == NULL) return SCIPF_BAD_ARGUMENT;
    er = collect_vars(m, nv, vars1, &vars); if (er != 0) return er;
    er = ensure_cons_cap(m); if (er != 0) { free(vars); return er; }
    idx = m->nconss;
    if (sos2)
        rc = SCIPcreateConsBasicSOS2(m->scip, &m->conss[idx], name, nv, vars,
                                     (SCIP_Real *)weights);
    else
        rc = SCIPcreateConsBasicSOS1(m->scip, &m->conss[idx], name, nv, vars,
                                     (SCIP_Real *)weights);
    free(vars);
    if (rc != SCIP_OKAY) return (int)rc;
    rc = SCIPaddCons(m->scip, m->conss[idx]);
    if (rc != SCIP_OKAY) { SCIPreleaseCons(m->scip, &m->conss[idx]); return (int)rc; }
    m->nconss++;
    *index1 = idx + 1;
    return 0;
}

int scipf_add_sos1_cons(void *h, int nv, const int *v, const double *w,
                        const char *name, int *idx)
{ return add_sos(h, 0, nv, v, w, name, idx); }

int scipf_add_sos2_cons(void *h, int nv, const int *v, const double *w,
                        const char *name, int *idx)
{ return add_sos(h, 1, nv, v, w, name, idx); }

int scipf_add_indicator_cons(void *handle, int binvar1, int nv,
                             const int *vars1, const double *coefs,
                             double rhs, const char *name, int *index1)
{
    scipf_model *m = (scipf_model *)handle;
    SCIP_VAR **vars = NULL;
    SCIP_RETCODE rc;
    int er, idx, b = binvar1 - 1;
    if (m == NULL || b < 0 || b >= m->nvars || index1 == NULL) return SCIPF_BAD_INDEX;
    er = collect_vars(m, nv, vars1, &vars); if (er != 0) return er;
    er = ensure_cons_cap(m); if (er != 0) { free(vars); return er; }
    idx = m->nconss;
    rc = SCIPcreateConsBasicIndicator(m->scip, &m->conss[idx], name,
                                      m->vars[b], nv, vars,
                                      (SCIP_Real *)coefs, rhs);
    free(vars);
    if (rc != SCIP_OKAY) return (int)rc;
    rc = SCIPaddCons(m->scip, m->conss[idx]);
    if (rc != SCIP_OKAY) { SCIPreleaseCons(m->scip, &m->conss[idx]); return (int)rc; }
    m->nconss++;
    *index1 = idx + 1;
    return 0;
}

static int normalize_rc(SCIP_RETCODE rc) { return rc == SCIP_OKAY ? 0 : (int)rc; }

int scipf_set_param_bool(void *h, const char *name, int value)
{ return normalize_rc(SCIPsetBoolParam(((scipf_model *)h)->scip, name, value ? TRUE : FALSE)); }
int scipf_set_param_int(void *h, const char *name, int value)
{ return normalize_rc(SCIPsetIntParam(((scipf_model *)h)->scip, name, value)); }
int scipf_set_param_long(void *h, const char *name, long long value)
{ return normalize_rc(SCIPsetLongintParam(((scipf_model *)h)->scip, name, (SCIP_Longint)value)); }
int scipf_set_param_real(void *h, const char *name, double value)
{ return normalize_rc(SCIPsetRealParam(((scipf_model *)h)->scip, name, value)); }
int scipf_set_param_char(void *h, const char *name, char value)
{ return normalize_rc(SCIPsetCharParam(((scipf_model *)h)->scip, name, value)); }
int scipf_set_param_string(void *h, const char *name, const char *value)
{ return normalize_rc(SCIPsetStringParam(((scipf_model *)h)->scip, name, value)); }

int scipf_set_heuristics(void *h, int setting)
{
    SCIP_PARAMSETTING s = SCIP_PARAMSETTING_DEFAULT;
    if (setting == 1) s = SCIP_PARAMSETTING_AGGRESSIVE;
    else if (setting == 2) s = SCIP_PARAMSETTING_FAST;
    else if (setting == 3) s = SCIP_PARAMSETTING_OFF;
    return normalize_rc(SCIPsetHeuristics(((scipf_model *)h)->scip, s, TRUE));
}

int scipf_set_objective_sense(void *h, int maximize)
{
    return normalize_rc(SCIPsetObjsense(((scipf_model *)h)->scip,
        maximize ? SCIP_OBJSENSE_MAXIMIZE : SCIP_OBJSENSE_MINIMIZE));
}

int scipf_optimize(void *h)
{ return normalize_rc(SCIPsolve(((scipf_model *)h)->scip)); }

int scipf_get_status(void *h, char *buf, int nbuf)
{
    const char *s;
    if (h == NULL || buf == NULL || nbuf <= 0) return SCIPF_BAD_ARGUMENT;
    s = status_string(SCIPgetStatus(((scipf_model *)h)->scip));
    snprintf(buf, (size_t)nbuf, "%s", s);
    return 0;
}

int scipf_get_nsols(void *h)
{ return SCIPgetNSols(((scipf_model *)h)->scip); }
int scipf_get_nvars(void *h)
{ return ((scipf_model *)h)->nvars; }
int scipf_get_nconss(void *h)
{ return ((scipf_model *)h)->nconss; }

static int get_solution(scipf_model *m, SCIP_SOL *sol, double *obj, double *x, int n)
{
    int j;
    if (sol == NULL) return SCIPF_NO_SOLUTION;
    if (n < m->nvars) return SCIPF_BAD_ARGUMENT;
    if (obj != NULL) *obj = SCIPgetSolOrigObj(m->scip, sol);
    if (x != NULL)
        for (j = 0; j < m->nvars; ++j) x[j] = SCIPgetSolVal(m->scip, sol, m->vars[j]);
    return 0;
}

int scipf_get_best_solution(void *h, double *obj, double *x, int n)
{
    scipf_model *m = (scipf_model *)h;
    if (m == NULL || SCIPgetNSols(m->scip) <= 0) return SCIPF_NO_SOLUTION;
    return get_solution(m, SCIPgetBestSol(m->scip), obj, x, n);
}

int scipf_get_solution_k(void *h, int k1, double *obj, double *x, int n)
{
    scipf_model *m = (scipf_model *)h;
    SCIP_SOL **sols;
    int ns;
    if (m == NULL) return SCIPF_BAD_ARGUMENT;
    ns = SCIPgetNSols(m->scip);
    if (k1 < 1 || k1 > ns) return SCIPF_BAD_INDEX;
    sols = SCIPgetSols(m->scip);
    return get_solution(m, sols[k1 - 1], obj, x, n);
}

int scipf_get_info(void *h, double *solve_time, long long *nodes,
                   long long *iterations, double *gap, int *sol_count)
{
    scipf_model *m = (scipf_model *)h;
    SCIP_STATUS status;
    int ns;
    if (m == NULL) return SCIPF_BAD_ARGUMENT;
    ns = SCIPgetNSols(m->scip);
    status = SCIPgetStatus(m->scip);
    if (solve_time) *solve_time = SCIPgetSolvingTime(m->scip);
    if (nodes) *nodes = (long long)SCIPgetNNodes(m->scip);
    if (iterations) *iterations = (long long)SCIPgetNLPIterations(m->scip);
    if (gap) *gap = (ns > 0 && status != SCIP_STATUS_OPTIMAL) ? SCIPgetGap(m->scip) : 0.0;
    if (sol_count) *sol_count = ns;
    return 0;
}

/* Type-detecting setters corresponding to the R wrapper's behavior. */
int scipf_set_param_integer_auto(void *h, const char *name, long long value)
{
    scipf_model *m = (scipf_model *)h;
    SCIP_PARAM *p;
    if (m == NULL) return SCIPF_BAD_ARGUMENT;
    p = SCIPgetParam(m->scip, name);
    if (p == NULL) return SCIPF_BAD_ARGUMENT;
    switch (SCIPparamGetType(p)) {
    case SCIP_PARAMTYPE_BOOL: return normalize_rc(SCIPsetBoolParam(m->scip, name, value != 0));
    case SCIP_PARAMTYPE_INT: return normalize_rc(SCIPsetIntParam(m->scip, name, (int)value));
    case SCIP_PARAMTYPE_LONGINT:
        return normalize_rc(SCIPsetLongintParam(m->scip, name, (SCIP_Longint)value));
    case SCIP_PARAMTYPE_REAL: return normalize_rc(SCIPsetRealParam(m->scip, name, (double)value));
    case SCIP_PARAMTYPE_CHAR: return normalize_rc(SCIPsetCharParam(m->scip, name, (char)value));
    default: return SCIPF_BAD_ARGUMENT;
    }
}

int scipf_set_param_real_auto(void *h, const char *name, double value)
{
    scipf_model *m = (scipf_model *)h;
    SCIP_PARAM *p;
    if (m == NULL) return SCIPF_BAD_ARGUMENT;
    p = SCIPgetParam(m->scip, name);
    if (p == NULL) return SCIPF_BAD_ARGUMENT;
    switch (SCIPparamGetType(p)) {
    case SCIP_PARAMTYPE_BOOL: return normalize_rc(SCIPsetBoolParam(m->scip, name, value != 0.0));
    case SCIP_PARAMTYPE_INT: return normalize_rc(SCIPsetIntParam(m->scip, name, (int)value));
    case SCIP_PARAMTYPE_LONGINT:
        return normalize_rc(SCIPsetLongintParam(m->scip, name, (SCIP_Longint)value));
    case SCIP_PARAMTYPE_REAL: return normalize_rc(SCIPsetRealParam(m->scip, name, value));
    default: return SCIPF_BAD_ARGUMENT;
    }
}

int scipf_set_param_text_auto(void *h, const char *name, const char *value)
{
    scipf_model *m = (scipf_model *)h;
    SCIP_PARAM *p;
    if (m == NULL || value == NULL) return SCIPF_BAD_ARGUMENT;
    p = SCIPgetParam(m->scip, name);
    if (p == NULL) return SCIPF_BAD_ARGUMENT;
    if (SCIPparamGetType(p) == SCIP_PARAMTYPE_CHAR)
        return normalize_rc(SCIPsetCharParam(m->scip, name, value[0]));
    if (SCIPparamGetType(p) == SCIP_PARAMTYPE_STRING)
        return normalize_rc(SCIPsetStringParam(m->scip, name, value));
    return SCIPF_BAD_ARGUMENT;
}
