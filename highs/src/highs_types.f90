! SPDX-License-Identifier: GPL-2.0-or-later
module highs_types
   use, intrinsic :: iso_c_binding, only : c_ptr, c_null_ptr
   use highs_kinds, only : dp, highs_int, highs_int64
   use highs_constants, only : highs_minimize, highs_var_continuous, highs_model_notset
   use highs_sparse, only : highs_sparse_matrix, highs_empty_matrix
   implicit none
   private

   type, public :: highs_control
      integer(highs_int) :: threads = 1
      real(dp) :: time_limit = -1.0_dp
      logical :: log_to_console = .false.
      character(len=32) :: solver = "choose"
      character(len=16) :: presolve = "choose"
      logical :: parallel = .false.
      real(dp) :: mip_rel_gap = -1.0_dp
      real(dp) :: mip_abs_gap = -1.0_dp
      integer(highs_int) :: random_seed = -1
   end type highs_control

   type, public :: highs_model
      integer(highs_int) :: num_col = 0
      integer(highs_int) :: num_row = 0
      integer(highs_int) :: sense = highs_minimize
      real(dp) :: offset = 0.0_dp
      real(dp), allocatable :: col_cost(:)
      real(dp), allocatable :: col_lower(:)
      real(dp), allocatable :: col_upper(:)
      real(dp), allocatable :: row_lower(:)
      real(dp), allocatable :: row_upper(:)
      integer(highs_int), allocatable :: integrality(:)
      type(highs_sparse_matrix) :: a
      type(highs_sparse_matrix) :: q
      logical :: has_hessian = .false.
   contains
      procedure :: valid => model_valid
   end type highs_model

   type, public :: highs_info
      logical :: valid = .false.
      integer(highs_int64) :: mip_node_count = 0
      integer(highs_int) :: simplex_iteration_count = 0
      integer(highs_int) :: ipm_iteration_count = 0
      integer(highs_int) :: crossover_iteration_count = 0
      integer(highs_int) :: qp_iteration_count = 0
      integer(highs_int) :: primal_solution_status = 0
      integer(highs_int) :: dual_solution_status = 0
      integer(highs_int) :: basis_validity = 0
      real(dp) :: objective_function_value = 0.0_dp
      real(dp) :: mip_dual_bound = 0.0_dp
      real(dp) :: mip_gap = 0.0_dp
      real(dp) :: max_integrality_violation = 0.0_dp
      integer(highs_int) :: num_primal_infeasibilities = 0
      real(dp) :: max_primal_infeasibility = 0.0_dp
      real(dp) :: sum_primal_infeasibilities = 0.0_dp
      integer(highs_int) :: num_dual_infeasibilities = 0
      real(dp) :: max_dual_infeasibility = 0.0_dp
      real(dp) :: sum_dual_infeasibilities = 0.0_dp
      real(dp) :: run_time = 0.0_dp
   end type highs_info

   type, public :: highs_solution
      integer(highs_int) :: call_status = 0
      integer(highs_int) :: model_status = highs_model_notset
      character(len=128) :: status_message = "not set"
      logical :: value_valid = .false.
      logical :: dual_valid = .false.
      real(dp) :: objective_value = 0.0_dp
      real(dp), allocatable :: col_value(:)
      real(dp), allocatable :: col_dual(:)
      real(dp), allocatable :: row_value(:)
      real(dp), allocatable :: row_dual(:)
      type(highs_info) :: info
   end type highs_solution

   type, public :: highs_basis
      logical :: valid = .false.
      integer(highs_int), allocatable :: col_status(:)
      integer(highs_int), allocatable :: row_status(:)
   end type highs_basis

   type, public :: highs_solver
      type(c_ptr) :: handle = c_null_ptr
      integer(highs_int) :: num_col = 0
      integer(highs_int) :: num_row = 0
   contains
      final :: finalize_solver_stub
   end type highs_solver

contains

   pure logical function model_valid(self) result(ok)
      class(highs_model), intent(in) :: self
      ok = self%num_col >= 0 .and. self%num_row >= 0
      if (.not. ok) return
      ok = allocated(self%col_cost) .and. size(self%col_cost) == self%num_col
      ok = ok .and. allocated(self%col_lower) .and. size(self%col_lower) == self%num_col
      ok = ok .and. allocated(self%col_upper) .and. size(self%col_upper) == self%num_col
      ok = ok .and. allocated(self%row_lower) .and. size(self%row_lower) == self%num_row
      ok = ok .and. allocated(self%row_upper) .and. size(self%row_upper) == self%num_row
      ok = ok .and. allocated(self%integrality) .and. size(self%integrality) == self%num_col
      ok = ok .and. self%a%nrow == self%num_row .and. self%a%ncol == self%num_col
      ok = ok .and. self%a%valid()
      if (self%has_hessian) ok = ok .and. self%q%nrow == self%num_col .and. self%q%valid()
   end function model_valid

   subroutine finalize_solver_stub(self)
      use, intrinsic :: iso_c_binding, only : c_associated
      type(highs_solver), intent(inout) :: self
      interface
         subroutine hf_api_destroy_local(handle) bind(c, name="hf_api_destroy")
            use, intrinsic :: iso_c_binding, only : c_ptr
            type(c_ptr), value :: handle
         end subroutine hf_api_destroy_local
      end interface
      if (c_associated(self%handle)) call hf_api_destroy_local(self%handle)
      self%handle = c_null_ptr
      self%num_col = 0
      self%num_row = 0
   end subroutine finalize_solver_stub

end module highs_types
