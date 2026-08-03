/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef HIGHS_FORTRAN_ABI_H
#define HIGHS_FORTRAN_ABI_H
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif

typedef struct hf_info_c {
  int valid;
  int64_t mip_node_count;
  int simplex_iteration_count;
  int ipm_iteration_count;
  int crossover_iteration_count;
  int qp_iteration_count;
  int primal_solution_status;
  int dual_solution_status;
  int basis_validity;
  double objective_function_value;
  double mip_dual_bound;
  double mip_gap;
  double max_integrality_violation;
  int num_primal_infeasibilities;
  double max_primal_infeasibility;
  double sum_primal_infeasibilities;
  int num_dual_infeasibilities;
  double max_dual_infeasibility;
  double sum_dual_infeasibilities;
  double run_time;
} hf_info_c;

#ifdef __cplusplus
}
#endif
#endif
