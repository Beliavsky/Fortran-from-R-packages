! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2012-2014 Alexios Galanos and Bernhard Pfaff.
! Modern Fortran translation Copyright (C) 2026 OpenAI.
! This derivative work is distributed under GPL-3.0-or-later.
! Derived from parma 1.7, Copyright (C) 2012-2014 Alexios Galanos and Bernhard Pfaff.
module parma_types
   use parma_kinds, only: dp
   implicit none
   private

   integer, parameter, public :: risk_mad = 1
   integer, parameter, public :: risk_ev = 2
   integer, parameter, public :: risk_minimax = 3
   integer, parameter, public :: risk_cvar = 4
   integer, parameter, public :: risk_cdar = 5
   integer, parameter, public :: risk_lpm = 6
   integer, parameter, public :: risk_upm = 7
   integer, parameter, public :: risk_rachev = 8

   integer, parameter, public :: solve_min_risk = 1
   integer, parameter, public :: solve_max_reward = 2
   integer, parameter, public :: solve_max_ratio = 3
   integer, parameter, public :: solve_utility = 4

   type, public :: parma_options
      integer :: max_iter = 2000
      integer :: population = 0
      integer :: seed = 12345
      integer :: print_level = 0
      real(dp) :: tol = 1.0e-8_dp
      real(dp) :: sigma0 = 0.20_dp
      real(dp) :: penalty = 1.0e5_dp
      real(dp) :: step0 = 0.10_dp
      logical :: use_cmaes = .true.
      logical :: polish = .true.
   end type parma_options

   type, public :: parma_spec
      integer :: risk = risk_ev
      integer :: objective = solve_min_risk
      real(dp) :: budget = 1.0_dp
      real(dp) :: leverage = -1.0_dp
      real(dp) :: target = 0.0_dp
      real(dp) :: alpha = 0.05_dp
      real(dp) :: moment = 1.0_dp
      real(dp) :: threshold = 0.0_dp
      real(dp) :: risk_aversion = 1.0_dp
      real(dp) :: ratio_rf = 0.0_dp
      real(dp) :: turnover_limit = huge(1.0_dp)
      real(dp) :: buy_turnover_limit = huge(1.0_dp)
      real(dp) :: sell_turnover_limit = huge(1.0_dp)
      real(dp) :: variance_limit = huge(1.0_dp)
      logical :: target_is_equality = .false.
      logical :: lpm_legacy = .false.
      integer :: max_positions = 0
      real(dp), allocatable :: data(:,:)
      real(dp), allocatable :: mu(:)
      real(dp), allocatable :: cov(:,:)
      real(dp), allocatable :: lb(:)
      real(dp), allocatable :: ub(:)
      real(dp), allocatable :: benchmark(:)
      real(dp), allocatable :: benchmark_cov(:)
      real(dp), allocatable :: initial(:)
      real(dp), allocatable :: eq_a(:,:)
      real(dp), allocatable :: eq_b(:)
      real(dp), allocatable :: ineq_a(:,:)
      real(dp), allocatable :: ineq_lb(:)
      real(dp), allocatable :: ineq_ub(:)
   end type parma_spec

   type, public :: parma_port
      integer :: status = 1
      integer :: iterations = 0
      real(dp) :: reward = 0.0_dp
      real(dp) :: risk = huge(1.0_dp)
      real(dp) :: objective = huge(1.0_dp)
      real(dp) :: multiplier = 1.0_dp
      real(dp) :: var_level = 0.0_dp
      real(dp) :: dar_level = 0.0_dp
      character(len=160) :: message = 'not solved'
      real(dp), allocatable :: weights(:)
   end type parma_port



   type, public :: qp_result
      integer :: status = 1
      integer :: iterations = 0
      real(dp) :: objective = huge(1.0_dp)
      character(len=160) :: message = 'not solved'
      real(dp), allocatable :: x(:)
   end type qp_result

   type, public :: lp_result
      integer :: status = 1
      integer :: iterations = 0
      real(dp) :: objective = huge(1.0_dp)
      character(len=160) :: message = 'not solved'
      real(dp), allocatable :: x(:)
   end type lp_result

   type, public :: milp_result
      integer :: status = 1
      integer :: evaluations = 0
      real(dp) :: objective = huge(1.0_dp)
      character(len=160) :: message = 'not solved'
      integer, allocatable :: x(:)
   end type milp_result

   type, public :: socp_result
      integer :: status = 1
      integer :: iterations = 0
      real(dp) :: objective = huge(1.0_dp)
      real(dp) :: max_violation = huge(1.0_dp)
      character(len=160) :: message = 'not solved'
      real(dp), allocatable :: x(:)
   end type socp_result

   type, public :: cmaes_result
      integer :: status = 1
      integer :: iterations = 0
      integer :: evaluations = 0
      real(dp) :: value = huge(1.0_dp)
      real(dp), allocatable :: x(:)
   end type cmaes_result

   abstract interface
      function objective_callback(x) result(value)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp) :: value
      end function objective_callback
   end interface
   public :: objective_callback

end module parma_types
