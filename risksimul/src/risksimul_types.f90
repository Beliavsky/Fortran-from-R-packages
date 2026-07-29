! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Derived from riskSimul 0.1.2 by Wolfgang Hormann and Ismail Basoglu.
module risksimul_types
   use ghyp_kinds, only : dp, i8
   use ghyp_model, only : ghyp_model_type
   implicit none
   private

   integer, parameter, public :: marginal_t = 1
   integer, parameter, public :: marginal_gh = 2

   integer, parameter, public :: objective_mse = 0
   integer, parameter, public :: objective_msre = -1
   integer, parameter, public :: objective_max_error = -2
   integer, parameter, public :: objective_max_relative = -3

   type, public :: inverse_table
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: cdf(:)
      logical :: ready = .false.
   end type inverse_table

   type, public :: portfolio_model
      integer :: marginal_family = marginal_t
      real(dp) :: copula_df = 0.0_dp
      real(dp), allocatable :: correlation(:,:)
      real(dp), allocatable :: cholesky(:,:)
      real(dp), allocatable :: marginal_parameters(:,:)
      real(dp), allocatable :: scale(:)
      real(dp), allocatable :: weight(:)
      type(ghyp_model_type), allocatable :: gh_models(:)
      type(inverse_table), allocatable :: gh_tables(:)
      logical :: ok = .false.
      character(len=192) :: message = ''
   contains
      procedure :: dimension => portfolio_dimension
   end type portfolio_model

   type, public :: estimate_row
      real(dp) :: estimate = 0.0_dp
      real(dp) :: halfwidth = 0.0_dp
      real(dp) :: ci_lower = 0.0_dp
      real(dp) :: ci_upper = 0.0_dp
      real(dp) :: variance = 0.0_dp
      real(dp) :: relative_error_percent = 0.0_dp
   end type estimate_row

   type, public :: simulation_result
      type(estimate_row), allocatable :: tail_probability(:)
      type(estimate_row), allocatable :: unconditional_excess(:)
      type(estimate_row), allocatable :: conditional_excess(:)
      real(dp), allocatable :: thresholds(:)
      integer :: samples_requested = 0
      integer :: samples_used = 0
      logical :: ok = .false.
      character(len=192) :: message = ''
   end type simulation_result

   type, public :: sis_control
      integer, allocatable :: allocations(:)
      integer :: normal_strata = 22
      integer :: gamma_strata = 22
      integer :: minimum_per_stratum = 10
      logical :: optimize_conditional_excess = .false.
      real(dp) :: intermediate_weight = 0.75_dp
      integer :: multi_objective = objective_msre
      real(dp) :: allocation_tolerance = 1.0e-6_dp
      integer :: direction_iterations = 160
      integer(i8) :: seed = 123456789_i8
      logical :: upstream_allocation_compatibility = .false.
   end type sis_control

   type, public :: importance_parameters
      real(dp), allocatable :: shift(:)
      real(dp) :: gamma_mean = 0.0_dp
      real(dp) :: objective = 0.0_dp
      logical :: ok = .false.
      character(len=192) :: message = ''
   end type importance_parameters

   type, public :: allocation_result
      real(dp), allocatable :: fractions(:)
      real(dp), allocatable :: objectives(:)
      real(dp) :: worst_objective = 0.0_dp
      logical :: ok = .false.
      character(len=192) :: message = ''
   end type allocation_result

contains

   pure function portfolio_dimension(self) result(d)
      class(portfolio_model), intent(in) :: self
      integer :: d
      if (allocated(self%weight)) then
         d = size(self%weight)
      else
         d = 0
      end if
   end function portfolio_dimension

end module risksimul_types
