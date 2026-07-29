! SPDX-License-Identifier: GPL-3.0-only
! Derived from riskParityPortfolio 0.2.2.9000, Copyright Ze Vinicius and Daniel P. Palomar.
module rpp_types
   use rpp_kinds, only: dp
   implicit none
   private

   integer, parameter, public :: RPP_OK = 0
   integer, parameter, public :: RPP_INVALID_INPUT = 1
   integer, parameter, public :: RPP_INFEASIBLE = 2
   integer, parameter, public :: RPP_LINEAR_SOLVE_FAILED = 3
   integer, parameter, public :: RPP_MAX_ITER = 4

   type, public :: risk_parity_result
      real(dp), allocatable :: weights(:)
      real(dp), allocatable :: relative_risk_contribution(:)
      real(dp), allocatable :: objective_history(:)
      real(dp) :: risk_concentration = 0.0_dp
      real(dp) :: mean_return = 0.0_dp
      real(dp) :: variance = 0.0_dp
      real(dp) :: theta = 0.0_dp
      integer :: iterations = 0
      integer :: status = RPP_OK
      logical :: converged = .false.
      logical :: feasible = .false.
      character(len=32) :: method = ''
      character(len=40) :: formulation = ''
   end type risk_parity_result
end module rpp_types
