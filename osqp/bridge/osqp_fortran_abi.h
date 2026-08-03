/* SPDX-License-Identifier: Apache-2.0 */
#ifndef OSQP_FORTRAN_ABI_H
#define OSQP_FORTRAN_ABI_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct of_settings_c {
  int device;
  int linsys_solver;
  int allocate_solution;
  int verbose;
  int profiler_level;
  int warm_starting;
  int scaling;
  int polishing;
  double rho;
  int rho_is_vec;
  double sigma;
  double alpha;
  int cg_max_iter;
  int cg_tol_reduction;
  double cg_tol_fraction;
  int cg_precond;
  int adaptive_rho;
  int adaptive_rho_interval;
  double adaptive_rho_fraction;
  double adaptive_rho_tolerance;
  int max_iter;
  double eps_abs;
  double eps_rel;
  double eps_prim_inf;
  double eps_dual_inf;
  int scaled_termination;
  int check_termination;
  int check_dualgap;
  double time_limit;
  double delta;
  int polish_refine_iter;
} of_settings_c;

typedef struct of_info_c {
  int status_val;
  int status_polish;
  double obj_val;
  double dual_obj_val;
  double prim_res;
  double dual_res;
  double duality_gap;
  int iter;
  int rho_updates;
  double rho_estimate;
  double setup_time;
  double solve_time;
  double update_time;
  double polish_time;
  double run_time;
  double primdual_int;
  double rel_kkt_error;
} of_info_c;

#ifdef __cplusplus
}
#endif
#endif
