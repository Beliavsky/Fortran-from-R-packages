/* SPDX-License-Identifier: Apache-2.0 */
#include "osqp.h"
#include "osqp_fortran_abi.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef _WIN32
#define OF_EXPORT __declspec(dllexport)
#else
#define OF_EXPORT __attribute__((visibility("default")))
#endif

typedef struct of_handle {
  OSQPSolver* solver;
  OSQPInt n;
  OSQPInt m;
  OSQPInt p_nnz;
  OSQPInt a_nnz;
  OSQPFloat* px;
  OSQPInt* pi;
  OSQPInt* pp;
  OSQPFloat* ax;
  OSQPInt* ai;
  OSQPInt* ap;
  OSQPFloat* q;
  OSQPFloat* l;
  OSQPFloat* u;
  OSQPSettings settings;
} of_handle;

static void copy_settings_to_osqp(OSQPSettings* out, const of_settings_c* in) {
  osqp_set_default_settings(out);
  if (!in) return;
  out->device = in->device;
  out->linsys_solver = (enum osqp_linsys_solver_type)in->linsys_solver;
  out->allocate_solution = in->allocate_solution;
  out->verbose = in->verbose;
  out->profiler_level = in->profiler_level;
  out->warm_starting = in->warm_starting;
  out->scaling = in->scaling;
  out->polishing = in->polishing;
  out->rho = in->rho;
  out->rho_is_vec = in->rho_is_vec;
  out->sigma = in->sigma;
  out->alpha = in->alpha;
  out->cg_max_iter = in->cg_max_iter;
  out->cg_tol_reduction = in->cg_tol_reduction;
  out->cg_tol_fraction = in->cg_tol_fraction;
  out->cg_precond = (osqp_precond_type)in->cg_precond;
  out->adaptive_rho = in->adaptive_rho;
  out->adaptive_rho_interval = in->adaptive_rho_interval;
  out->adaptive_rho_fraction = in->adaptive_rho_fraction;
  out->adaptive_rho_tolerance = in->adaptive_rho_tolerance;
  out->max_iter = in->max_iter;
  out->eps_abs = in->eps_abs;
  out->eps_rel = in->eps_rel;
  out->eps_prim_inf = in->eps_prim_inf;
  out->eps_dual_inf = in->eps_dual_inf;
  out->scaled_termination = in->scaled_termination;
  out->check_termination = in->check_termination;
  out->check_dualgap = in->check_dualgap;
  out->time_limit = in->time_limit;
  out->delta = in->delta;
  out->polish_refine_iter = in->polish_refine_iter;
}

static void copy_settings_from_osqp(of_settings_c* out, const OSQPSettings* in) {
  if (!out || !in) return;
  out->device = in->device;
  out->linsys_solver = (int)in->linsys_solver;
  out->allocate_solution = in->allocate_solution;
  out->verbose = in->verbose;
  out->profiler_level = in->profiler_level;
  out->warm_starting = in->warm_starting;
  out->scaling = in->scaling;
  out->polishing = in->polishing;
  out->rho = in->rho;
  out->rho_is_vec = in->rho_is_vec;
  out->sigma = in->sigma;
  out->alpha = in->alpha;
  out->cg_max_iter = in->cg_max_iter;
  out->cg_tol_reduction = in->cg_tol_reduction;
  out->cg_tol_fraction = in->cg_tol_fraction;
  out->cg_precond = (int)in->cg_precond;
  out->adaptive_rho = in->adaptive_rho;
  out->adaptive_rho_interval = in->adaptive_rho_interval;
  out->adaptive_rho_fraction = in->adaptive_rho_fraction;
  out->adaptive_rho_tolerance = in->adaptive_rho_tolerance;
  out->max_iter = in->max_iter;
  out->eps_abs = in->eps_abs;
  out->eps_rel = in->eps_rel;
  out->eps_prim_inf = in->eps_prim_inf;
  out->eps_dual_inf = in->eps_dual_inf;
  out->scaled_termination = in->scaled_termination;
  out->check_termination = in->check_termination;
  out->check_dualgap = in->check_dualgap;
  out->time_limit = in->time_limit;
  out->delta = in->delta;
  out->polish_refine_iter = in->polish_refine_iter;
}

static void free_handle(of_handle* h) {
  if (!h) return;
  if (h->solver) osqp_cleanup(h->solver);
  free(h->px); free(h->pi); free(h->pp);
  free(h->ax); free(h->ai); free(h->ap);
  free(h->q); free(h->l); free(h->u);
  free(h);
}

static int copy_int_array(OSQPInt** dst, const int* src, int n) {
  int k;
  if (n <= 0) { *dst = NULL; return 1; }
  *dst = (OSQPInt*)malloc((size_t)n * sizeof(OSQPInt));
  if (!*dst) return 0;
  for (k = 0; k < n; ++k) (*dst)[k] = (OSQPInt)src[k];
  return 1;
}

static int copy_float_array(OSQPFloat** dst, const double* src, int n) {
  int k;
  if (n <= 0) { *dst = NULL; return 1; }
  *dst = (OSQPFloat*)malloc((size_t)n * sizeof(OSQPFloat));
  if (!*dst) return 0;
  for (k = 0; k < n; ++k) (*dst)[k] = (OSQPFloat)src[k];
  return 1;
}

static void clamp_bounds(OSQPFloat* x, int n) {
  int k;
  for (k = 0; k < n; ++k) {
    if (isinf((double)x[k]) || x[k] > OSQP_INFTY) x[k] = OSQP_INFTY;
    else if (x[k] < -OSQP_INFTY) x[k] = -OSQP_INFTY;
  }
}

OF_EXPORT void* of_backend_create(
    int n, int m,
    int p_nnz, const int* p_colptr, const int* p_rowind, const double* p_value,
    int a_nnz, const int* a_colptr, const int* a_rowind, const double* a_value,
    const double* q, const double* l, const double* u,
    const of_settings_c* settings, int* error_code) {
  of_handle* h = NULL;
  OSQPCscMatrix P;
  OSQPCscMatrix A;
  OSQPInt rc;
  if (error_code) *error_code = OSQP_MEM_ALLOC_ERROR;
  if (n < 1 || m < 0 || p_nnz < 0 || a_nnz < 0 || !p_colptr || !a_colptr || !q) {
    if (error_code) *error_code = OSQP_DATA_VALIDATION_ERROR;
    return NULL;
  }
  h = (of_handle*)calloc(1, sizeof(of_handle));
  if (!h) return NULL;
  h->n = (OSQPInt)n;
  h->m = (OSQPInt)m;
  h->p_nnz = (OSQPInt)p_nnz;
  h->a_nnz = (OSQPInt)a_nnz;
  if (!copy_int_array(&h->pp, p_colptr, n + 1) ||
      !copy_int_array(&h->pi, p_rowind, p_nnz) ||
      !copy_float_array(&h->px, p_value, p_nnz) ||
      !copy_int_array(&h->ap, a_colptr, n + 1) ||
      !copy_int_array(&h->ai, a_rowind, a_nnz) ||
      !copy_float_array(&h->ax, a_value, a_nnz) ||
      !copy_float_array(&h->q, q, n) ||
      !copy_float_array(&h->l, l, m) ||
      !copy_float_array(&h->u, u, m)) {
    free_handle(h);
    return NULL;
  }
  clamp_bounds(h->l, m);
  clamp_bounds(h->u, m);
  copy_settings_to_osqp(&h->settings, settings);
  OSQPCscMatrix_set_data(&P, n, n, p_nnz, h->px, h->pi, h->pp);
  OSQPCscMatrix_set_data(&A, m, n, a_nnz, h->ax, h->ai, h->ap);
  rc = osqp_setup(&h->solver, &P, h->q, &A, h->l, h->u, m, n, &h->settings);
  if (error_code) *error_code = (int)rc;
  if (rc != OSQP_NO_ERROR || !h->solver) {
    free_handle(h);
    return NULL;
  }
  return h;
}

OF_EXPORT void of_backend_destroy(void* p) { free_handle((of_handle*)p); }

OF_EXPORT int of_backend_solve(void* p) {
  of_handle* h = (of_handle*)p;
  return h && h->solver ? (int)osqp_solve(h->solver) : OSQP_WORKSPACE_NOT_INIT_ERROR;
}

OF_EXPORT int of_backend_get_dimensions(void* p, int* n, int* m) {
  of_handle* h = (of_handle*)p;
  if (!h) return OSQP_WORKSPACE_NOT_INIT_ERROR;
  if (n) *n = (int)h->n;
  if (m) *m = (int)h->m;
  return 0;
}

OF_EXPORT int of_backend_get_solution(void* p, double* x, double* y,
                                      double* prim_inf_cert, double* dual_inf_cert,
                                      of_info_c* out) {
  int k;
  of_handle* h = (of_handle*)p;
  OSQPSolver* s;
  if (!h || !(s = h->solver) || !s->info || !s->solution) return OSQP_WORKSPACE_NOT_INIT_ERROR;
  if (x) for (k = 0; k < h->n; ++k) x[k] = s->solution->x ? s->solution->x[k] : 0.0;
  if (y) for (k = 0; k < h->m; ++k) y[k] = s->solution->y ? s->solution->y[k] : 0.0;
  if (prim_inf_cert) for (k = 0; k < h->m; ++k) prim_inf_cert[k] = s->solution->prim_inf_cert ? s->solution->prim_inf_cert[k] : 0.0;
  if (dual_inf_cert) for (k = 0; k < h->n; ++k) dual_inf_cert[k] = s->solution->dual_inf_cert ? s->solution->dual_inf_cert[k] : 0.0;
  if (out) {
    out->status_val = (int)s->info->status_val;
    out->status_polish = (int)s->info->status_polish;
    out->obj_val = s->info->obj_val;
    out->dual_obj_val = s->info->dual_obj_val;
    out->prim_res = s->info->prim_res;
    out->dual_res = s->info->dual_res;
    out->duality_gap = s->info->duality_gap;
    out->iter = (int)s->info->iter;
    out->rho_updates = (int)s->info->rho_updates;
    out->rho_estimate = s->info->rho_estimate;
    out->setup_time = s->info->setup_time;
    out->solve_time = s->info->solve_time;
    out->update_time = s->info->update_time;
    out->polish_time = s->info->polish_time;
    out->run_time = s->info->run_time;
    out->primdual_int = s->info->primdual_int;
    out->rel_kkt_error = s->info->rel_kkt_error;
  }
  return 0;
}

static int copy_string(const char* s, char* buffer, int n) {
  int len = s ? (int)strlen(s) : 0;
  if (buffer && n > 0) {
    int m = len < n - 1 ? len : n - 1;
    if (m > 0) memcpy(buffer, s, (size_t)m);
    buffer[m] = '\0';
  }
  return len;
}

OF_EXPORT int of_backend_status_message(void* p, char* buffer, int n) {
  of_handle* h = (of_handle*)p;
  if (!h || !h->solver || !h->solver->info) return copy_string("backend solver is not initialized", buffer, n);
  return copy_string(h->solver->info->status, buffer, n);
}

OF_EXPORT int of_backend_version(char* buffer, int n) { return copy_string(osqp_version(), buffer, n); }
OF_EXPORT int of_backend_capabilities(void) { return (int)osqp_capabilities(); }
OF_EXPORT double of_backend_infinity(void) { return (double)OSQP_INFTY; }
OF_EXPORT int of_backend_error_message(int code, char* buffer, int n) { return copy_string(osqp_error_message((OSQPInt)code), buffer, n); }

OF_EXPORT int of_backend_get_settings(void* p, of_settings_c* out) {
  of_handle* h = (of_handle*)p;
  if (!h || !h->solver || !h->solver->settings || !out) return OSQP_WORKSPACE_NOT_INIT_ERROR;
  copy_settings_from_osqp(out, h->solver->settings);
  return 0;
}

OF_EXPORT int of_backend_update_data_vec(void* p,
                                         const double* q_new, int update_q,
                                         const double* l_new, int update_l,
                                         const double* u_new, int update_u) {
  int k;
  OSQPInt rc;
  OSQPFloat *qv = NULL, *lv = NULL, *uv = NULL;
  of_handle* h = (of_handle*)p;
  if (!h || !h->solver) return OSQP_WORKSPACE_NOT_INIT_ERROR;
  if (update_q) {
    if (!q_new || !copy_float_array(&qv, q_new, h->n)) return OSQP_MEM_ALLOC_ERROR;
  }
  if (update_l) {
    if (!l_new || !copy_float_array(&lv, l_new, h->m)) { free(qv); return OSQP_MEM_ALLOC_ERROR; }
    clamp_bounds(lv, h->m);
  }
  if (update_u) {
    if (!u_new || !copy_float_array(&uv, u_new, h->m)) { free(qv); free(lv); return OSQP_MEM_ALLOC_ERROR; }
    clamp_bounds(uv, h->m);
  }
  rc = osqp_update_data_vec(h->solver, qv, lv, uv);
  if (rc == 0) {
    if (update_q) for (k = 0; k < h->n; ++k) h->q[k] = qv[k];
    if (update_l) for (k = 0; k < h->m; ++k) h->l[k] = lv[k];
    if (update_u) for (k = 0; k < h->m; ++k) h->u[k] = uv[k];
  }
  free(qv); free(lv); free(uv);
  return (int)rc;
}

OF_EXPORT int of_backend_update_data_mat(void* p,
                                         const double* px_new, const int* px_idx, int p_new_n, int update_p,
                                         const double* ax_new, const int* ax_idx, int a_new_n, int update_a) {
  int k;
  OSQPInt rc;
  OSQPFloat *pv = NULL, *av = NULL;
  OSQPInt *pi = NULL, *ai = NULL;
  of_handle* h = (of_handle*)p;
  if (!h || !h->solver) return OSQP_WORKSPACE_NOT_INIT_ERROR;
  if (update_p) {
    if (!px_new || p_new_n < 0 || !copy_float_array(&pv, px_new, p_new_n)) return OSQP_MEM_ALLOC_ERROR;
    if (px_idx && !copy_int_array(&pi, px_idx, p_new_n)) { free(pv); return OSQP_MEM_ALLOC_ERROR; }
  }
  if (update_a) {
    if (!ax_new || a_new_n < 0 || !copy_float_array(&av, ax_new, a_new_n)) { free(pv); free(pi); return OSQP_MEM_ALLOC_ERROR; }
    if (ax_idx && !copy_int_array(&ai, ax_idx, a_new_n)) { free(pv); free(pi); free(av); return OSQP_MEM_ALLOC_ERROR; }
  }
  rc = osqp_update_data_mat(h->solver, pv, pi, update_p ? p_new_n : 0,
                            av, ai, update_a ? a_new_n : 0);
  if (rc == 0) {
    if (update_p) {
      if (pi) for (k = 0; k < p_new_n; ++k) h->px[pi[k]] = pv[k];
      else for (k = 0; k < p_new_n; ++k) h->px[k] = pv[k];
    }
    if (update_a) {
      if (ai) for (k = 0; k < a_new_n; ++k) h->ax[ai[k]] = av[k];
      else for (k = 0; k < a_new_n; ++k) h->ax[k] = av[k];
    }
  }
  free(pv); free(pi); free(av); free(ai);
  return (int)rc;
}

OF_EXPORT int of_backend_warm_start(void* p, const double* x, int has_x,
                                    const double* y, int has_y) {
  OSQPFloat *xv = NULL, *yv = NULL;
  OSQPInt rc;
  of_handle* h = (of_handle*)p;
  if (!h || !h->solver) return OSQP_WORKSPACE_NOT_INIT_ERROR;
  if (has_x && (!x || !copy_float_array(&xv, x, h->n))) return OSQP_MEM_ALLOC_ERROR;
  if (has_y && (!y || !copy_float_array(&yv, y, h->m))) { free(xv); return OSQP_MEM_ALLOC_ERROR; }
  rc = osqp_warm_start(h->solver, xv, yv);
  free(xv); free(yv);
  return (int)rc;
}

OF_EXPORT int of_backend_cold_start(void* p) {
  of_handle* h = (of_handle*)p;
  if (!h || !h->solver) return OSQP_WORKSPACE_NOT_INIT_ERROR;
  osqp_cold_start(h->solver);
  return 0;
}

OF_EXPORT int of_backend_update_settings(void* p, const of_settings_c* settings) {
  OSQPInt rc;
  OSQPSettings s;
  of_handle* h = (of_handle*)p;
  if (!h || !h->solver || !settings) return OSQP_WORKSPACE_NOT_INIT_ERROR;
  copy_settings_to_osqp(&s, settings);
  rc = osqp_update_rho(h->solver, s.rho);
  if (rc != 0) return (int)rc;
  rc = osqp_update_settings(h->solver, &s);
  return (int)rc;
}

