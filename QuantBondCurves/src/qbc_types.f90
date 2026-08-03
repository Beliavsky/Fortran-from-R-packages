! SPDX-License-Identifier: GPL-3.0-or-later
module qbc_types
   use qbc_kinds, only : dp
   use qbc_dates, only : qbc_date
   implicit none
   private

   integer, parameter, public :: qbc_asset_tes = 1
   integer, parameter, public :: qbc_asset_fixed = 2
   integer, parameter, public :: qbc_asset_ibr = 3
   integer, parameter, public :: qbc_asset_libor = 4
   integer, parameter, public :: qbc_asset_irs = 5
   integer, parameter, public :: qbc_asset_ccs = 6

   integer, parameter, public :: qbc_rate_continuous = 0
   integer, parameter, public :: qbc_rate_discrete = 1

   integer, parameter, public :: qbc_schedule_short_first = 1
   integer, parameter, public :: qbc_schedule_long_first = 2
   integer, parameter, public :: qbc_schedule_short_last = 3
   integer, parameter, public :: qbc_schedule_long_last = 4

   integer, parameter, public :: qbc_leg_fixed_fixed = 1
   integer, parameter, public :: qbc_leg_fixed_float = 2
   integer, parameter, public :: qbc_leg_float_fixed = 3
   integer, parameter, public :: qbc_leg_float_float = 4

   type, public :: qbc_coupon_schedule
      type(qbc_date), allocatable :: dates(:)
      type(qbc_date), allocatable :: effective_dates(:)
   end type qbc_coupon_schedule

   type, public :: qbc_curve
      real(dp), allocatable :: terms(:)
      real(dp), allocatable :: rates(:)
      integer :: approximation = 1 ! 1 constant, 2 linear
      integer :: rate_type = qbc_rate_continuous
      integer :: frequency = 1
   end type qbc_curve

   type, public :: qbc_bond
      type(qbc_date) :: maturity
      type(qbc_date) :: analysis_date
      type(qbc_date) :: trade_date
      logical :: has_trade_date = .false.
      real(dp) :: coupon_rate = 0.0_dp
      real(dp) :: principal = 1.0_dp
      real(dp) :: spread = 0.0_dp
      integer :: asset_type = qbc_asset_tes
      integer :: frequency = 1
      integer :: rate_type = qbc_rate_discrete
      integer :: coupon_schedule = qbc_schedule_short_first
      character(len=16) :: daycount = 'ACT/365'
      character(len=2) :: business_convention = 'F'
   end type qbc_bond

   type, public :: qbc_swap
      type(qbc_date) :: maturity
      type(qbc_date) :: analysis_date
      integer :: frequency = 4
      integer :: rate_type = qbc_rate_discrete
      integer :: legs = qbc_leg_fixed_fixed
      real(dp) :: coupon_rate1 = 0.0_dp
      real(dp) :: coupon_rate2 = 0.0_dp
      real(dp) :: spread1 = 0.0_dp
      real(dp) :: spread2 = 0.0_dp
      real(dp) :: float_rate1 = 0.0_dp
      real(dp) :: float_rate2 = 0.0_dp
      real(dp) :: principal1 = 1.0_dp
      real(dp) :: principal2 = 1.0_dp
      real(dp) :: exchange_rate = 1.0_dp
      character(len=16) :: daycount = 'ACT/365'
   end type qbc_swap

   type, public :: qbc_bond_sensitivity
      real(dp) :: modified_duration = 0.0_dp
      real(dp) :: convexity = 0.0_dp
      real(dp) :: dv01 = 0.0_dp
      real(dp) :: price = 0.0_dp
   end type qbc_bond_sensitivity

   type, public :: qbc_calibration_result
      type(qbc_curve) :: curve
      real(dp) :: objective = huge(1.0_dp)
      integer :: convergence = 1
      integer :: iterations = 0
      character(len=160) :: message = 'not calibrated'
   end type qbc_calibration_result

end module qbc_types
