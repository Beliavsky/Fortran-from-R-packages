! SPDX-License-Identifier: BSD-2-Clause
module piqp
   use piqp_kinds, only : dp
   use piqp_types
   use piqp_solver, only : piqp_model_type, solve_piqp_dense
   use piqp_matrix_adapter, only : solve_piqp_sparse, csc_to_dense_piqp
   implicit none
   public
   interface solve_piqp
      module procedure solve_piqp_dense
      module procedure solve_piqp_sparse
   end interface solve_piqp
end module piqp
