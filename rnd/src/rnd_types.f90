! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from RND 1.2, Copyright (C) 2017 Kam Hamidieh.
module rnd_types
   use rnd_kinds, only : dp
   implicit none
   private

   type, public :: option_prices
      real(dp), allocatable :: call(:)
      real(dp), allocatable :: put(:)
   end type option_prices

   type, public :: optimizer_result
      real(dp), allocatable :: par(:)
      real(dp) :: value = huge(1.0_dp)
      integer :: iterations = 0
      integer :: convergence = 1
      real(dp), allocatable :: hessian(:, :)
   end type optimizer_result

   type, public :: bsm_fit
      real(dp) :: mu = 0.0_dp
      real(dp) :: zeta = 0.0_dp
      type(optimizer_result) :: optimizer
   end type bsm_fit

   type, public :: gb_fit
      real(dp) :: a = 0.0_dp
      real(dp) :: b = 0.0_dp
      real(dp) :: v = 0.0_dp
      real(dp) :: w = 0.0_dp
      type(optimizer_result) :: optimizer
   end type gb_fit

   type, public :: mln_fit
      real(dp) :: alpha1 = 0.0_dp
      real(dp) :: meanlog1 = 0.0_dp
      real(dp) :: meanlog2 = 0.0_dp
      real(dp) :: sdlog1 = 0.0_dp
      real(dp) :: sdlog2 = 0.0_dp
      type(optimizer_result) :: optimizer
   end type mln_fit

   type, public :: ew_fit
      real(dp) :: sigma = 0.0_dp
      real(dp) :: skew = 0.0_dp
      real(dp) :: kurt = 0.0_dp
      type(optimizer_result) :: optimizer
   end type ew_fit

   type, public :: am_fit
      real(dp) :: w1 = 0.0_dp
      real(dp) :: w2 = 0.0_dp
      real(dp) :: u1 = 0.0_dp
      real(dp) :: u2 = 0.0_dp
      real(dp) :: u3 = 0.0_dp
      real(dp) :: sigma1 = 0.0_dp
      real(dp) :: sigma2 = 0.0_dp
      real(dp) :: sigma3 = 0.0_dp
      real(dp) :: p1 = 0.0_dp
      real(dp) :: p2 = 0.0_dp
      type(optimizer_result) :: optimizer
   end type am_fit

   type, public :: shimko_fit
      real(dp) :: a0 = 0.0_dp
      real(dp) :: a1 = 0.0_dp
      real(dp) :: a2 = 0.0_dp
      real(dp), allocatable :: implied_volatility(:)
      real(dp), allocatable :: density(:)
   end type shimko_fit

   type, public :: rate_result
      real(dp) :: risk_free_rate = 0.0_dp
      real(dp) :: dividend_yield = 0.0_dp
   end type rate_result

   type, public :: am_price_result
      real(dp) :: call = 0.0_dp
      real(dp) :: put = 0.0_dp
      real(dp) :: expected_f0 = 0.0_dp
      real(dp) :: prob_above = 0.0_dp
      real(dp) :: prob_below = 0.0_dp
      real(dp) :: conditional_above = 0.0_dp
      real(dp) :: conditional_below = 0.0_dp
   end type am_price_result

   type, public :: moe_result
      type(bsm_fit) :: bsm
      type(gb_fit) :: gb
      type(mln_fit) :: mln
      type(ew_fit) :: ew
      type(shimko_fit) :: shimko
      real(dp), allocatable :: point_density(:)
   end type moe_result
end module rnd_types
