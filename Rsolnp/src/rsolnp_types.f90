! SPDX-License-Identifier: GPL-2.0-only
module rsolnp_types
   use rsolnp_kinds, only : dp
   use rsolnp_callbacks, only : objective_callback, gradient_callback, vector_callback, jacobian_callback
   implicit none
   private

   integer, parameter, public :: solnp_success = 0
   integer, parameter, public :: solnp_max_iterations = 1
   integer, parameter, public :: solnp_invalid_problem = 2
   integer, parameter, public :: solnp_numerical_failure = 3
   integer, parameter, public :: solnp_infeasible = 4

   type, public :: solnp_control
      real(dp) :: rho = 1.0_dp
      integer :: max_iter = 400
      integer :: min_iter = 400
      real(dp) :: delta = 1.0e-6_dp
      real(dp) :: tol = 1.0e-8_dp
      integer :: trace = 0
      real(dp) :: penalty_growth = 5.0_dp
      real(dp) :: max_rho = 1.0e12_dp
      integer :: line_search_max = 32
      real(dp) :: armijo = 1.0e-4_dp
      real(dp) :: min_step = 1.0e-12_dp
      integer :: restoration_iter = 80
   end type solnp_control

   type, public :: solnp_problem
      character(len=64) :: name = 'unnamed'
      integer :: n = 0
      integer :: n_eq = 0
      integer :: n_ineq = 0
      integer :: raw_n_eq = 0
      integer :: raw_n_ineq = 0
      logical :: standard_form = .false.
      procedure(objective_callback), pointer, nopass :: fn => null()
      procedure(gradient_callback), pointer, nopass :: gr => null()
      procedure(vector_callback), pointer, nopass :: eq_fn => null()
      procedure(jacobian_callback), pointer, nopass :: eq_jac => null()
      procedure(vector_callback), pointer, nopass :: ineq_fn => null()
      procedure(jacobian_callback), pointer, nopass :: ineq_jac => null()
      class(*), allocatable :: data
      real(dp), allocatable :: start(:)
      real(dp), allocatable :: lower(:)
      real(dp), allocatable :: upper(:)
      real(dp), allocatable :: eq_b(:)
      real(dp), allocatable :: ineq_lower(:)
      real(dp), allocatable :: ineq_upper(:)
      real(dp), allocatable :: standard_eq_shift(:)
      real(dp), allocatable :: standard_ineq_lower(:)
      real(dp), allocatable :: standard_ineq_upper(:)
      real(dp) :: best_fn = huge(1.0_dp)
      real(dp), allocatable :: best_par(:)
   end type solnp_problem

   type, public :: kkt_diagnostics
      real(dp) :: stationarity = huge(1.0_dp)
      real(dp) :: eq_violation = 0.0_dp
      real(dp) :: ineq_violation = 0.0_dp
      real(dp) :: bound_violation = 0.0_dp
      real(dp) :: dual_feas_violation = 0.0_dp
      real(dp) :: complementarity = 0.0_dp
      logical :: primal_feasible = .false.
      logical :: first_order = .false.
   end type kkt_diagnostics

   type, public :: solnp_result
      real(dp), allocatable :: pars(:)
      real(dp) :: objective = huge(1.0_dp)
      real(dp), allocatable :: objective_history(:)
      real(dp), allocatable :: constraint_history(:)
      real(dp), allocatable :: step_history(:)
      real(dp), allocatable :: lagrange(:)
      real(dp), allocatable :: hessian(:, :)
      real(dp), allocatable :: ineq_slack(:)
      integer :: out_iterations = 0
      integer :: inner_iterations = 0
      integer :: convergence = solnp_invalid_problem
      integer :: n_eval = 0
      real(dp) :: elapsed = 0.0_dp
      character(len=160) :: message = 'problem not solved'
      type(kkt_diagnostics) :: kkt
   end type solnp_result

   type, public :: multistart_result
      type(solnp_result) :: best
      type(solnp_result), allocatable :: results(:)
      real(dp), allocatable :: starts(:, :)
      real(dp), allocatable :: objectives(:)
      integer, allocatable :: convergence(:)
      integer :: best_index = 0
   end type multistart_result

   type, public :: problem_table_entry
      character(len=24) :: suite = ''
      character(len=32) :: name = ''
      integer :: number = 0
      logical :: implemented = .false.
   end type problem_table_entry

end module rsolnp_types
