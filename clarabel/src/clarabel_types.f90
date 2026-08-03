! SPDX-License-Identifier: Apache-2.0
module clarabel_types
   use, intrinsic :: iso_c_binding, only : c_int8_t, c_int32_t, c_size_t
   use clarabel_kinds, only : dp
   implicit none
   private

   integer(c_int8_t), parameter, public :: cone_zero = 0_c_int8_t
   integer(c_int8_t), parameter, public :: cone_nonnegative = 1_c_int8_t
   integer(c_int8_t), parameter, public :: cone_second_order = 2_c_int8_t
   integer(c_int8_t), parameter, public :: cone_exponential = 3_c_int8_t
   integer(c_int8_t), parameter, public :: cone_power = 4_c_int8_t
   integer(c_int8_t), parameter, public :: cone_generalized_power = 5_c_int8_t
   integer(c_int8_t), parameter, public :: cone_psd_triangle = 6_c_int8_t

   integer(c_int32_t), parameter, public :: status_unsolved = 0_c_int32_t
   integer(c_int32_t), parameter, public :: status_solved = 1_c_int32_t
   integer(c_int32_t), parameter, public :: status_primal_infeasible = 2_c_int32_t
   integer(c_int32_t), parameter, public :: status_dual_infeasible = 3_c_int32_t
   integer(c_int32_t), parameter, public :: status_almost_solved = 4_c_int32_t
   integer(c_int32_t), parameter, public :: status_almost_primal_infeasible = 5_c_int32_t
   integer(c_int32_t), parameter, public :: status_almost_dual_infeasible = 6_c_int32_t
   integer(c_int32_t), parameter, public :: status_max_iterations = 7_c_int32_t
   integer(c_int32_t), parameter, public :: status_max_time = 8_c_int32_t
   integer(c_int32_t), parameter, public :: status_numerical_error = 9_c_int32_t
   integer(c_int32_t), parameter, public :: status_insufficient_progress = 10_c_int32_t
   integer(c_int32_t), parameter, public :: status_callback_terminated = 11_c_int32_t

   integer(c_int32_t), parameter, public :: direct_solver_auto = 0_c_int32_t
   integer(c_int32_t), parameter, public :: direct_solver_qdldl = 1_c_int32_t
   integer(c_int32_t), parameter, public :: direct_solver_faer = 2_c_int32_t
   integer(c_int32_t), parameter, public :: direct_solver_mkl = 3_c_int32_t
   integer(c_int32_t), parameter, public :: direct_solver_panua = 4_c_int32_t

   integer(c_int32_t), parameter, public :: chordal_merge_none = 0_c_int32_t
   integer(c_int32_t), parameter, public :: chordal_merge_parent_child = 1_c_int32_t
   integer(c_int32_t), parameter, public :: chordal_merge_clique_graph = 2_c_int32_t

   type, public :: clarabel_cone
      integer(c_int8_t) :: kind = cone_zero
      integer(c_size_t) :: dim = 0_c_size_t
      real(dp) :: parameter = 0.0_dp
      real(dp), allocatable :: alpha(:)
   contains
      procedure :: total_dimension => cone_total_dimension
      procedure :: validate => cone_validate
   end type clarabel_cone

   type, public :: clarabel_settings
      integer(c_int32_t) :: max_iter = 200_c_int32_t
      real(dp) :: time_limit = huge(1.0_dp)
      logical :: verbose = .true.
      real(dp) :: max_step_fraction = 0.99_dp
      real(dp) :: tol_gap_abs = 1.0e-8_dp
      real(dp) :: tol_gap_rel = 1.0e-8_dp
      real(dp) :: tol_feas = 1.0e-8_dp
      real(dp) :: tol_infeas_abs = 1.0e-8_dp
      real(dp) :: tol_infeas_rel = 1.0e-8_dp
      real(dp) :: tol_ktratio = 1.0e-6_dp
      real(dp) :: reduced_tol_gap_abs = 5.0e-5_dp
      real(dp) :: reduced_tol_gap_rel = 5.0e-5_dp
      real(dp) :: reduced_tol_feas = 1.0e-4_dp
      real(dp) :: reduced_tol_infeas_abs = 5.0e-5_dp
      real(dp) :: reduced_tol_infeas_rel = 5.0e-5_dp
      real(dp) :: reduced_tol_ktratio = 1.0e-4_dp
      logical :: equilibrate_enable = .true.
      integer(c_int32_t) :: equilibrate_max_iter = 10_c_int32_t
      real(dp) :: equilibrate_min_scaling = 1.0e-4_dp
      real(dp) :: equilibrate_max_scaling = 1.0e4_dp
      real(dp) :: linesearch_backtrack_step = 0.8_dp
      real(dp) :: min_switch_step_length = 1.0e-1_dp
      real(dp) :: min_terminate_step_length = 1.0e-4_dp
      integer(c_int32_t) :: max_threads = 0_c_int32_t
      logical :: direct_kkt_solver = .true.
      integer(c_int32_t) :: direct_solve_method = direct_solver_qdldl
      logical :: static_regularization_enable = .true.
      real(dp) :: static_regularization_constant = 1.0e-8_dp
      real(dp) :: static_regularization_proportional = epsilon(1.0_dp)**2
      logical :: dynamic_regularization_enable = .true.
      real(dp) :: dynamic_regularization_eps = 1.0e-13_dp
      real(dp) :: dynamic_regularization_delta = 2.0e-7_dp
      logical :: iterative_refinement_enable = .true.
      real(dp) :: iterative_refinement_reltol = 1.0e-13_dp
      real(dp) :: iterative_refinement_abstol = 1.0e-12_dp
      integer(c_int32_t) :: iterative_refinement_max_iter = 10_c_int32_t
      real(dp) :: iterative_refinement_stop_ratio = 5.0_dp
      logical :: presolve_enable = .true.
      logical :: input_sparse_dropzeros = .false.
      logical :: chordal_decomposition_enable = .false.
      integer(c_int32_t) :: chordal_decomposition_merge_method = chordal_merge_none
      logical :: chordal_decomposition_compact = .false.
      logical :: chordal_decomposition_complete_dual = .false.
   end type clarabel_settings

   type, public :: clarabel_info
      real(dp) :: mu = 0.0_dp
      real(dp) :: sigma = 0.0_dp
      real(dp) :: step_length = 0.0_dp
      integer(c_int32_t) :: iterations = 0_c_int32_t
      real(dp) :: cost_primal = 0.0_dp
      real(dp) :: cost_dual = 0.0_dp
      real(dp) :: res_primal = 0.0_dp
      real(dp) :: res_dual = 0.0_dp
      real(dp) :: res_primal_inf = 0.0_dp
      real(dp) :: res_dual_inf = 0.0_dp
      real(dp) :: gap_abs = 0.0_dp
      real(dp) :: gap_rel = 0.0_dp
      real(dp) :: ktratio = 0.0_dp
      real(dp) :: solve_time = 0.0_dp
      integer(c_int32_t) :: status = status_unsolved
      integer(c_size_t) :: linear_solver_threads = 0_c_size_t
      integer(c_size_t) :: linear_solver_nnz_a = 0_c_size_t
      integer(c_size_t) :: linear_solver_nnz_l = 0_c_size_t
   end type clarabel_info

   type, public :: clarabel_solution
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: z(:)
      real(dp), allocatable :: s(:)
      real(dp) :: obj_val = 0.0_dp
      real(dp) :: obj_val_dual = 0.0_dp
      real(dp) :: solve_time = 0.0_dp
      integer(c_int32_t) :: iterations = 0_c_int32_t
      real(dp) :: r_prim = 0.0_dp
      real(dp) :: r_dual = 0.0_dp
      integer(c_int32_t) :: status = status_unsolved
      type(clarabel_info) :: info
   contains
      procedure :: solved => solution_solved
   end type clarabel_solution

   public :: zero_cone, nonnegative_cone, second_order_cone
   public :: exponential_cone, power_cone, generalized_power_cone
   public :: psd_triangle_cone, cones_total_dimension
   public :: status_name, default_clarabel_settings

contains

   pure function default_clarabel_settings() result(settings)
      type(clarabel_settings) :: settings
      settings = clarabel_settings()
   end function default_clarabel_settings

   function zero_cone(dim) result(cone)
      integer, intent(in) :: dim
      type(clarabel_cone) :: cone
      cone%kind = cone_zero
      cone%dim = int(dim, c_size_t)
   end function zero_cone

   function nonnegative_cone(dim) result(cone)
      integer, intent(in) :: dim
      type(clarabel_cone) :: cone
      cone%kind = cone_nonnegative
      cone%dim = int(dim, c_size_t)
   end function nonnegative_cone

   function second_order_cone(dim) result(cone)
      integer, intent(in) :: dim
      type(clarabel_cone) :: cone
      cone%kind = cone_second_order
      cone%dim = int(dim, c_size_t)
   end function second_order_cone

   function exponential_cone() result(cone)
      type(clarabel_cone) :: cone
      cone%kind = cone_exponential
      cone%dim = 3_c_size_t
   end function exponential_cone

   function power_cone(alpha) result(cone)
      real(dp), intent(in) :: alpha
      type(clarabel_cone) :: cone
      cone%kind = cone_power
      cone%dim = 3_c_size_t
      cone%parameter = alpha
   end function power_cone

   function generalized_power_cone(alpha, dim) result(cone)
      real(dp), intent(in) :: alpha(:)
      integer, intent(in) :: dim
      type(clarabel_cone) :: cone
      cone%kind = cone_generalized_power
      cone%dim = int(dim, c_size_t)
      allocate(cone%alpha(size(alpha)))
      cone%alpha = alpha
   end function generalized_power_cone

   function psd_triangle_cone(matrix_order) result(cone)
      integer, intent(in) :: matrix_order
      type(clarabel_cone) :: cone
      cone%kind = cone_psd_triangle
      cone%dim = int(matrix_order, c_size_t)
   end function psd_triangle_cone

   pure integer(c_size_t) function cone_total_dimension(self) result(n)
      class(clarabel_cone), intent(in) :: self
      select case (self%kind)
      case (cone_zero, cone_nonnegative, cone_second_order)
         n = self%dim
      case (cone_exponential, cone_power)
         n = 3_c_size_t
      case (cone_generalized_power)
         if (allocated(self%alpha)) then
            n = int(size(self%alpha), c_size_t) + self%dim
         else
            n = self%dim
         end if
      case (cone_psd_triangle)
         n = self%dim * (self%dim + 1_c_size_t) / 2_c_size_t
      case default
         n = 0_c_size_t
      end select
   end function cone_total_dimension

   pure integer(c_size_t) function cones_total_dimension(cones) result(n)
      type(clarabel_cone), intent(in) :: cones(:)
      integer :: i
      n = 0_c_size_t
      do i = 1, size(cones)
         n = n + cones(i)%total_dimension()
      end do
   end function cones_total_dimension

   subroutine cone_validate(self, ok, message)
      class(clarabel_cone), intent(in) :: self
      logical, intent(out) :: ok
      character(len=:), allocatable, intent(out) :: message
      real(dp) :: s

      ok = .false.
      message = ""
      select case (self%kind)
      case (cone_zero, cone_nonnegative)
         if (self%dim < 0_c_size_t) then
            message = "cone dimension cannot be negative"
            return
         end if
      case (cone_second_order)
         if (self%dim < 2_c_size_t) then
            message = "a second-order cone must have dimension at least 2"
            return
         end if
      case (cone_exponential)
         continue
      case (cone_power)
         if (self%parameter <= 0.0_dp .or. self%parameter >= 1.0_dp) then
            message = "power-cone alpha must lie strictly between 0 and 1"
            return
         end if
      case (cone_generalized_power)
         if (.not. allocated(self%alpha) .or. size(self%alpha) < 2) then
            message = "generalized-power cone requires at least two exponents"
            return
         end if
         if (any(self%alpha <= 0.0_dp)) then
            message = "generalized-power exponents must be positive"
            return
         end if
         s = sum(self%alpha)
         if (abs(s - 1.0_dp) > 100.0_dp * epsilon(1.0_dp)) then
            message = "generalized-power exponents must sum to one"
            return
         end if
         if (self%dim < 1_c_size_t) then
            message = "generalized-power trailing dimension must be positive"
            return
         end if
      case (cone_psd_triangle)
         if (self%dim < 1_c_size_t) then
            message = "PSD matrix order must be positive"
            return
         end if
      case default
         message = "unknown cone kind"
         return
      end select
      ok = .true.
   end subroutine cone_validate

   pure logical function solution_solved(self) result(ok)
      class(clarabel_solution), intent(in) :: self
      ok = self%status == status_solved .or. self%status == status_almost_solved
   end function solution_solved

   pure function status_name(status) result(name)
      integer(c_int32_t), intent(in) :: status
      character(len=:), allocatable :: name
      select case (status)
      case (status_unsolved); name = "Unsolved"
      case (status_solved); name = "Solved"
      case (status_primal_infeasible); name = "PrimalInfeasible"
      case (status_dual_infeasible); name = "DualInfeasible"
      case (status_almost_solved); name = "AlmostSolved"
      case (status_almost_primal_infeasible); name = "AlmostPrimalInfeasible"
      case (status_almost_dual_infeasible); name = "AlmostDualInfeasible"
      case (status_max_iterations); name = "MaxIterations"
      case (status_max_time); name = "MaxTime"
      case (status_numerical_error); name = "NumericalError"
      case (status_insufficient_progress); name = "InsufficientProgress"
      case (status_callback_terminated); name = "CallbackTerminated"
      case default; name = "Unknown"
      end select
   end function status_name

end module clarabel_types
