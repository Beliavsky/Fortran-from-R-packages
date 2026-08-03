/* SPDX-License-Identifier: Apache-2.0 */
#ifndef CLARABEL_FORTRAN_BRIDGE_H
#define CLARABEL_FORTRAN_BRIDGE_H
#include <stddef.h>
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    size_t nrows, ncols, nnz;
    const size_t *colptr;
    const size_t *rowind;
    const double *values;
} clarabel_csc_c;

typedef struct {
    uint8_t tag;
    size_t dim;
    double parameter;
    const double *alpha;
    size_t alpha_len;
} clarabel_cone_c;

typedef struct {
    int32_t max_iter;
    double time_limit;
    int32_t verbose;
    double max_step_fraction;
    double tol_gap_abs, tol_gap_rel, tol_feas;
    double tol_infeas_abs, tol_infeas_rel, tol_ktratio;
    double reduced_tol_gap_abs, reduced_tol_gap_rel, reduced_tol_feas;
    double reduced_tol_infeas_abs, reduced_tol_infeas_rel, reduced_tol_ktratio;
    int32_t equilibrate_enable, equilibrate_max_iter;
    double equilibrate_min_scaling, equilibrate_max_scaling;
    double linesearch_backtrack_step, min_switch_step_length, min_terminate_step_length;
    int32_t max_threads, direct_kkt_solver, direct_solve_method;
    int32_t static_regularization_enable;
    double static_regularization_constant, static_regularization_proportional;
    int32_t dynamic_regularization_enable;
    double dynamic_regularization_eps, dynamic_regularization_delta;
    int32_t iterative_refinement_enable;
    double iterative_refinement_reltol, iterative_refinement_abstol;
    int32_t iterative_refinement_max_iter;
    double iterative_refinement_stop_ratio;
    int32_t presolve_enable, input_sparse_dropzeros;
    int32_t chordal_decomposition_enable, chordal_decomposition_merge_method;
    int32_t chordal_decomposition_compact, chordal_decomposition_complete_dual;
} clarabel_settings_c;

typedef struct {
    int32_t status, iterations;
    double obj_val, obj_val_dual, solve_time, r_prim, r_dual;
    double mu, sigma, step_length;
    double cost_primal, cost_dual, res_primal, res_dual;
    double res_primal_inf, res_dual_inf, gap_abs, gap_rel, ktratio;
    size_t linear_solver_threads, linear_solver_nnz_a, linear_solver_nnz_l;
} clarabel_result_c;

void clarabel_settings_default(clarabel_settings_c *settings);
int32_t clarabel_solver_create(const clarabel_csc_c *P, const double *q, size_t q_len,
                               const clarabel_csc_c *A, const double *b, size_t b_len,
                               const clarabel_cone_c *cones, size_t ncones,
                               const clarabel_settings_c *settings, void **out,
                               char *error, size_t error_capacity);
int32_t clarabel_solver_solve(void *solver, double *x, size_t x_len,
                              double *z, size_t z_len, double *s, size_t s_len,
                              clarabel_result_c *result, char *error, size_t error_capacity);
int32_t clarabel_solver_update(void *solver,
                               const double *p_values, size_t p_len,
                               const double *a_values, size_t a_len,
                               const double *q, size_t q_len,
                               const double *b, size_t b_len,
                               char *error, size_t error_capacity);
int32_t clarabel_solver_is_update_allowed(const void *solver);
void clarabel_solver_free(void *solver);
#ifdef __cplusplus
}
#endif
#endif
