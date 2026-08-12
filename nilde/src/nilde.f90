! SPDX-License-Identifier: GPL-2.0-or-later
module nilde
   use nilde_kinds, only : i8, dp
   use nilde_types, only : integer_solutions_t, knapsack_result_t, bin_packing_result_t, tsp_result_t
   use nilde_diophantine, only : nlde, get_subsetsum
   use nilde_partitions, only : get_partitions
   use nilde_knapsack, only : get_knapsack
   use nilde_binpacking, only : bin_packing
   use nilde_tsp, only : tsp_solver, assignment_lower_bound
   implicit none
   public
end module nilde
