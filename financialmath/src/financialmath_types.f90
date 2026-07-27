! SPDX-License-Identifier: GPL-2.0-only
! Based on FinancialMath 0.1.1, Copyright (C) 2016 Kameron Penn and Jack Schmidt.
module financialmath_types
   use financialmath_kinds, only : dp
   implicit none
   private

   type, public :: status_t
      logical :: ok = .true.
      character(len=160) :: message = ''
   end type status_t

   type, public :: rate_conversion_t
      real(dp) :: effective_interest = 0.0_dp
      real(dp) :: effective_discount = 0.0_dp
      real(dp) :: force = 0.0_dp
      real(dp) :: input_nominal_interest = 0.0_dp
      real(dp) :: input_nominal_discount = 0.0_dp
      real(dp) :: target_nominal_interest = 0.0_dp
      real(dp) :: target_nominal_discount = 0.0_dp
   end type rate_conversion_t

   type, public :: tvm_result_t
      real(dp) :: present_value = 0.0_dp
      real(dp) :: future_value = 0.0_dp
      real(dp) :: periods = 0.0_dp
      real(dp) :: effective_rate = 0.0_dp
      real(dp) :: nominal_rate = 0.0_dp
      type(status_t) :: status
   end type tvm_result_t

   type, public :: cashflow_analysis_t
      real(dp) :: present_value = 0.0_dp
      real(dp) :: macaulay_duration = 0.0_dp
      real(dp) :: modified_duration = 0.0_dp
      real(dp) :: macaulay_convexity = 0.0_dp
      real(dp) :: modified_convexity = 0.0_dp
      type(status_t) :: status
   end type cashflow_analysis_t

   type, public :: annuity_result_t
      real(dp) :: present_value = 0.0_dp
      real(dp) :: future_value = 0.0_dp
      real(dp) :: first_payment = 0.0_dp
      real(dp) :: increment = 0.0_dp
      real(dp) :: growth_rate = 0.0_dp
      real(dp) :: periods = 0.0_dp
      real(dp) :: effective_rate = 0.0_dp
      real(dp) :: nominal_rate = 0.0_dp
      real(dp) :: payment_frequency_rate = 0.0_dp
      type(status_t) :: status
   end type annuity_result_t

   type, public :: amortization_result_t
      real(dp) :: loan = 0.0_dp
      real(dp) :: payment = 0.0_dp
      real(dp) :: periods = 0.0_dp
      real(dp) :: total_paid = 0.0_dp
      real(dp) :: total_interest = 0.0_dp
      real(dp) :: balloon_payment = 0.0_dp
      real(dp) :: drop_payment = 0.0_dp
      real(dp), allocatable :: time(:)
      real(dp), allocatable :: payments(:)
      real(dp), allocatable :: interest(:)
      real(dp), allocatable :: principal(:)
      real(dp), allocatable :: balance(:)
      type(status_t) :: status
   end type amortization_result_t

   type, public :: amort_period_result_t
      real(dp) :: loan = 0.0_dp
      real(dp) :: payment = 0.0_dp
      real(dp) :: periods = 0.0_dp
      real(dp) :: interest_paid = 0.0_dp
      real(dp) :: principal_paid = 0.0_dp
      real(dp) :: balance = 0.0_dp
      real(dp) :: effective_rate = 0.0_dp
      type(status_t) :: status
   end type amort_period_result_t

   type, public :: bond_result_t
      real(dp) :: price = 0.0_dp
      real(dp) :: premium = 0.0_dp
      real(dp) :: discount = 0.0_dp
      real(dp) :: coupon = 0.0_dp
      real(dp) :: effective_rate = 0.0_dp
      real(dp) :: macaulay_duration = 0.0_dp
      real(dp) :: modified_duration = 0.0_dp
      real(dp) :: macaulay_convexity = 0.0_dp
      real(dp) :: modified_convexity = 0.0_dp
      real(dp) :: price_at_time = 0.0_dp
      real(dp) :: full_price = 0.0_dp
      real(dp) :: clean_price = 0.0_dp
      real(dp) :: write_change = 0.0_dp
      type(status_t) :: status
   end type bond_result_t

   type, public :: option_order1_t
      real(dp) :: call_price = 0.0_dp
      real(dp) :: put_price = 0.0_dp
      real(dp) :: call_delta = 0.0_dp
      real(dp) :: put_delta = 0.0_dp
      real(dp) :: call_theta = 0.0_dp
      real(dp) :: put_theta = 0.0_dp
      real(dp) :: vega = 0.0_dp
      type(status_t) :: status
   end type option_order1_t

   type, public :: payoff_table_t
      real(dp), allocatable :: stock(:)
      real(dp), allocatable :: payoff(:)
      real(dp), allocatable :: profit(:)
      real(dp), allocatable :: premiums(:)
      type(status_t) :: status
   end type payoff_table_t

   type, public :: forward_result_t
      real(dp) :: delivery_price = 0.0_dp
      real(dp) :: prepaid_price = 0.0_dp
      real(dp), allocatable :: stock(:)
      real(dp), allocatable :: payoff(:)
      type(status_t) :: status
   end type forward_result_t

end module financialmath_types
