! SPDX-License-Identifier: BSD-2-Clause
module piqp_types
   use piqp_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: PIQP_SOLVED = 1
   integer, parameter, public :: PIQP_MAX_ITER_REACHED = -1
   integer, parameter, public :: PIQP_PRIMAL_INFEASIBLE = -2
   integer, parameter, public :: PIQP_DUAL_INFEASIBLE = -3
   integer, parameter, public :: PIQP_NUMERICS = -8
   integer, parameter, public :: PIQP_UNSOLVED = -9
   integer, parameter, public :: PIQP_INVALID_SETTINGS = -10

   type, public :: piqp_settings_type
      real(dp) :: rho_init = 1.0e-6_dp
      real(dp) :: delta_init = 1.0e-4_dp
      real(dp) :: eps_abs = 1.0e-8_dp
      real(dp) :: eps_rel = 1.0e-9_dp
      logical :: check_duality_gap = .true.
      real(dp) :: eps_duality_gap_abs = 1.0e-8_dp
      real(dp) :: eps_duality_gap_rel = 1.0e-9_dp
      real(dp) :: infeasibility_threshold = 0.9_dp
      real(dp) :: reg_lower_limit = 1.0e-10_dp
      real(dp) :: reg_finetune_lower_limit = 1.0e-13_dp
      integer :: reg_finetune_primal_update_threshold = 7
      integer :: reg_finetune_dual_update_threshold = 7
      integer :: max_iter = 250
      integer :: max_factor_retires = 10
      logical :: preconditioner_scale_cost = .false.
      logical :: preconditioner_reuse_on_update = .false.
      integer :: preconditioner_iter = 10
      real(dp) :: tau = 0.99_dp
      logical :: iterative_refinement_always_enabled = .false.
      real(dp) :: iterative_refinement_eps_abs = 1.0e-12_dp
      real(dp) :: iterative_refinement_eps_rel = 1.0e-12_dp
      integer :: iterative_refinement_max_iter = 10
      real(dp) :: iterative_refinement_min_improvement_rate = 5.0_dp
      real(dp) :: iterative_refinement_static_regularization_eps = 1.0e-8_dp
      real(dp) :: iterative_refinement_static_regularization_rel = epsilon(1.0_dp)**2
      logical :: verbose = .false.
      logical :: compute_timings = .false.
   contains
      procedure :: valid => piqp_settings_valid
   end type piqp_settings_type

   type, public :: piqp_info_type
      integer :: status = PIQP_UNSOLVED
      integer :: iter = 0
      real(dp) :: rho = 0.0_dp
      real(dp) :: delta = 0.0_dp
      real(dp) :: mu = 0.0_dp
      real(dp) :: sigma = 0.0_dp
      real(dp) :: primal_step = 0.0_dp
      real(dp) :: dual_step = 0.0_dp
      real(dp) :: primal_res = huge(1.0_dp)
      real(dp) :: primal_res_rel = huge(1.0_dp)
      real(dp) :: dual_res = huge(1.0_dp)
      real(dp) :: dual_res_rel = huge(1.0_dp)
      real(dp) :: primal_res_reg = huge(1.0_dp)
      real(dp) :: primal_res_reg_rel = huge(1.0_dp)
      real(dp) :: dual_res_reg = huge(1.0_dp)
      real(dp) :: dual_res_reg_rel = huge(1.0_dp)
      real(dp) :: primal_prox_inf = 0.0_dp
      real(dp) :: dual_prox_inf = 0.0_dp
      real(dp) :: prev_primal_res = huge(1.0_dp)
      real(dp) :: prev_dual_res = huge(1.0_dp)
      real(dp) :: primal_obj = 0.0_dp
      real(dp) :: dual_obj = 0.0_dp
      real(dp) :: duality_gap = huge(1.0_dp)
      real(dp) :: duality_gap_rel = huge(1.0_dp)
      integer :: factor_retires = 0
      real(dp) :: reg_limit = 0.0_dp
      integer :: no_primal_update = 0
      integer :: no_dual_update = 0
      real(dp) :: setup_time = 0.0_dp
      real(dp) :: update_time = 0.0_dp
      real(dp) :: solve_time = 0.0_dp
      real(dp) :: kkt_factor_time = 0.0_dp
      real(dp) :: kkt_solve_time = 0.0_dp
      real(dp) :: run_time = 0.0_dp
   end type piqp_info_type

   type, public :: piqp_result_type
      real(dp), allocatable :: x(:), y(:)
      real(dp), allocatable :: z_l(:), z_u(:), z_bl(:), z_bu(:)
      real(dp), allocatable :: s_l(:), s_u(:), s_bl(:), s_bu(:)
      type(piqp_info_type) :: info
   end type piqp_result_type

   public :: status_description, status_to_string

contains

   logical function piqp_settings_valid(self) result(ok)
      class(piqp_settings_type), intent(in) :: self
      ok = self%rho_init > 0.0_dp .and. self%delta_init > 0.0_dp .and. &
           self%eps_abs > 0.0_dp .and. self%eps_rel >= 0.0_dp .and. &
           self%eps_duality_gap_abs > 0.0_dp .and. self%eps_duality_gap_rel >= 0.0_dp .and. &
           self%infeasibility_threshold >= 0.0_dp .and. self%reg_lower_limit > 0.0_dp .and. &
           self%reg_finetune_lower_limit > 0.0_dp .and. &
           self%reg_finetune_primal_update_threshold >= 0 .and. &
           self%reg_finetune_dual_update_threshold >= 0 .and. self%max_iter > 0 .and. &
           self%max_factor_retires > 0 .and. self%preconditioner_iter >= 0 .and. &
           self%tau > 0.0_dp .and. self%tau <= 1.0_dp .and. &
           self%iterative_refinement_eps_abs > 0.0_dp .and. &
           self%iterative_refinement_eps_rel >= 0.0_dp .and. &
           self%iterative_refinement_max_iter >= 0 .and. &
           self%iterative_refinement_min_improvement_rate >= 1.0_dp .and. &
           self%iterative_refinement_static_regularization_eps > 0.0_dp .and. &
           self%iterative_refinement_static_regularization_rel >= 0.0_dp
   end function piqp_settings_valid

   function status_to_string(code) result(text)
      integer, intent(in) :: code
      character(len=:), allocatable :: text
      select case (code)
      case (PIQP_SOLVED); text = 'solved'
      case (PIQP_MAX_ITER_REACHED); text = 'max iterations reached'
      case (PIQP_PRIMAL_INFEASIBLE); text = 'primal infeasible'
      case (PIQP_DUAL_INFEASIBLE); text = 'dual infeasible'
      case (PIQP_NUMERICS); text = 'numerics issue'
      case (PIQP_UNSOLVED); text = 'unsolved'
      case (PIQP_INVALID_SETTINGS); text = 'invalid settings'
      case default; text = 'unknown'
      end select
   end function status_to_string

   function status_description(code) result(text)
      integer, intent(in) :: code
      character(len=:), allocatable :: text
      select case (code)
      case (PIQP_SOLVED); text = 'Solver solved problem up to given tolerance.'
      case (PIQP_MAX_ITER_REACHED); text = 'Iteration limit was reached.'
      case (PIQP_PRIMAL_INFEASIBLE); text = 'The problem is primal infeasible.'
      case (PIQP_DUAL_INFEASIBLE); text = 'The problem is dual infeasible.'
      case (PIQP_NUMERICS); text = 'Numerical error occurred during solving.'
      case (PIQP_UNSOLVED); text = 'The problem is unsolved, i.e., solve was never called.'
      case (PIQP_INVALID_SETTINGS); text = 'Invalid settings were provided to the solver.'
      case default; text = 'Unknown solver status.'
      end select
   end function status_description
end module piqp_types
