! SPDX-License-Identifier: Apache-2.0
module osqp_types
   use, intrinsic :: iso_c_binding, only : c_ptr, c_null_ptr, c_associated
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use osqp_kinds, only : dp, osqp_int
   use osqp_constants, only : osqp_direct_solver, osqp_diagonal_preconditioner, &
      osqp_adaptive_rho_iterations, osqp_unsolved
   use osqp_sparse, only : osqp_sparse_matrix
   implicit none
   private

   type, public :: osqp_settings
      integer(osqp_int) :: device = 0
      integer(osqp_int) :: linsys_solver = osqp_direct_solver
      logical :: allocate_solution = .true.
      logical :: verbose = .true.
      integer(osqp_int) :: profiler_level = 0
      logical :: warm_starting = .true.
      integer(osqp_int) :: scaling = 10
      logical :: polishing = .false.
      real(dp) :: rho = 0.1_dp
      logical :: rho_is_vec = .true.
      real(dp) :: sigma = 1.0e-6_dp
      real(dp) :: alpha = 1.6_dp
      integer(osqp_int) :: cg_max_iter = 20
      integer(osqp_int) :: cg_tol_reduction = 10
      real(dp) :: cg_tol_fraction = 0.15_dp
      integer(osqp_int) :: cg_precond = osqp_diagonal_preconditioner
      integer(osqp_int) :: adaptive_rho = osqp_adaptive_rho_iterations
      integer(osqp_int) :: adaptive_rho_interval = 50
      real(dp) :: adaptive_rho_fraction = 0.4_dp
      real(dp) :: adaptive_rho_tolerance = 5.0_dp
      integer(osqp_int) :: max_iter = 4000
      real(dp) :: eps_abs = 1.0e-3_dp
      real(dp) :: eps_rel = 1.0e-3_dp
      real(dp) :: eps_prim_inf = 1.0e-4_dp
      real(dp) :: eps_dual_inf = 1.0e-4_dp
      logical :: scaled_termination = .false.
      integer(osqp_int) :: check_termination = 25
      logical :: check_dualgap = .true.
      real(dp) :: time_limit = 1.0e10_dp
      real(dp) :: delta = 1.0e-6_dp
      integer(osqp_int) :: polish_refine_iter = 3
   contains
      procedure :: valid => settings_valid
   end type osqp_settings

   type, public :: osqp_model
      integer(osqp_int) :: n = 0
      integer(osqp_int) :: m = 0
      type(osqp_sparse_matrix) :: p
      type(osqp_sparse_matrix) :: a
      real(dp), allocatable :: q(:)
      real(dp), allocatable :: l(:)
      real(dp), allocatable :: u(:)
   contains
      procedure :: valid => model_valid
   end type osqp_model

   type, public :: osqp_info
      integer(osqp_int) :: status_val = osqp_unsolved
      integer(osqp_int) :: status_polish = 0
      real(dp) :: obj_val = 0.0_dp
      real(dp) :: dual_obj_val = 0.0_dp
      real(dp) :: prim_res = 0.0_dp
      real(dp) :: dual_res = 0.0_dp
      real(dp) :: duality_gap = 0.0_dp
      integer(osqp_int) :: iter = 0
      integer(osqp_int) :: rho_updates = 0
      real(dp) :: rho_estimate = 0.0_dp
      real(dp) :: setup_time = 0.0_dp
      real(dp) :: solve_time = 0.0_dp
      real(dp) :: update_time = 0.0_dp
      real(dp) :: polish_time = 0.0_dp
      real(dp) :: run_time = 0.0_dp
      real(dp) :: primdual_int = 0.0_dp
      real(dp) :: rel_kkt_error = 0.0_dp
   end type osqp_info

   type, public :: osqp_solution
      integer(osqp_int) :: call_status = 0
      character(len=64) :: status = "unsolved"
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: prim_inf_cert(:)
      real(dp), allocatable :: dual_inf_cert(:)
      type(osqp_info) :: info
   contains
      procedure :: solved => solution_solved
   end type osqp_solution

   type, public :: osqp_solver
      type(c_ptr) :: handle = c_null_ptr
      logical :: initialized = .false.
      type(osqp_model) :: model
      type(osqp_settings) :: settings
   contains
      final :: finalize_solver
   end type osqp_solver

contains

   pure logical function settings_valid(self) result(ok)
      class(osqp_settings), intent(in) :: self
      ok = self%rho > 0.0_dp .and. self%sigma > 0.0_dp
      ok = ok .and. self%alpha > 0.0_dp .and. self%alpha < 2.0_dp
      ok = ok .and. self%max_iter > 0 .and. self%check_termination >= 0
      ok = ok .and. self%eps_abs >= 0.0_dp .and. self%eps_rel >= 0.0_dp
      ok = ok .and. self%eps_prim_inf > 0.0_dp .and. self%eps_dual_inf > 0.0_dp
      ok = ok .and. self%delta > 0.0_dp .and. self%polish_refine_iter >= 0
      ok = ok .and. self%scaling >= 0 .and. self%time_limit > 0.0_dp
      ok = ok .and. self%adaptive_rho >= 0 .and. self%adaptive_rho <= 3
      ok = ok .and. self%adaptive_rho_interval >= 0
      ok = ok .and. self%adaptive_rho_fraction > 0.0_dp
      ok = ok .and. self%adaptive_rho_tolerance >= 1.0_dp
      ok = ok .and. self%cg_max_iter > 0 .and. self%cg_tol_reduction > 0
      ok = ok .and. self%cg_tol_fraction > 0.0_dp
   end function settings_valid

   pure logical function model_valid(self) result(ok)
      class(osqp_model), intent(in) :: self
      integer :: j, k
      ok = self%n >= 1 .and. self%m >= 0
      if (.not. ok) return
      ok = allocated(self%q) .and. allocated(self%l) .and. allocated(self%u)
      if (.not. ok) return
      ok = size(self%q) == self%n .and. size(self%l) == self%m .and. size(self%u) == self%m
      if (.not. ok) return
      ok = self%p%valid() .and. self%a%valid()
      ok = ok .and. self%p%nrow == self%n .and. self%p%ncol == self%n
      ok = ok .and. self%a%nrow == self%m .and. self%a%ncol == self%n
      if (.not. ok) return
      if (any(.not. ieee_is_finite(self%q)) .or. any(.not. ieee_is_finite(self%p%value)) .or. &
          any(.not. ieee_is_finite(self%a%value))) then
         ok = .false.
         return
      end if
      if (any(self%l > self%u)) then
         ok = .false.
         return
      end if
      do j = 1, self%n
         do k = self%p%col_ptr(j), self%p%col_ptr(j+1) - 1
            if (self%p%row_index(k) > j) then
               ok = .false.
               return
            end if
         end do
      end do
   end function model_valid

   pure logical function solution_solved(self) result(ok)
      class(osqp_solution), intent(in) :: self
      ok = self%info%status_val == 1 .or. self%info%status_val == 2
   end function solution_solved

   subroutine finalize_solver(self)
      type(osqp_solver), intent(inout) :: self
      interface
         subroutine c_destroy(handle) bind(c, name="of_api_destroy")
            use, intrinsic :: iso_c_binding, only : c_ptr
            type(c_ptr), value :: handle
         end subroutine c_destroy
      end interface
      if (c_associated(self%handle)) call c_destroy(self%handle)
      self%handle = c_null_ptr
      self%initialized = .false.
   end subroutine finalize_solver

end module osqp_types
