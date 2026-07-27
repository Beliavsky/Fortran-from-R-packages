! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2013 Bruno Remillard
! Modern Fortran translation copyright (C) 2026 OpenAI
module opthedging_types
   use opthedging_interpolation, only : interpolation1d
   use opthedging_kinds, only : dp
   implicit none
   private

   type, public :: hedging_result
      real(dp), allocatable :: s(:)
      real(dp), allocatable :: c(:,:)
      real(dp), allocatable :: a(:,:)
      real(dp), allocatable :: phi1(:)
      real(dp) :: rho = 0.0_dp
      real(dp) :: discounted_strike = 0.0_dp
      logical :: ok = .false.
      character(len=:), allocatable :: message
   contains
      procedure :: auxiliary_at
      procedure :: initial_hedge_at
      procedure :: option_value_at
      procedure :: shares_at
   end type hedging_result

contains

   function option_value_at(self, period, spot) result(value)
      class(hedging_result), intent(in) :: self
      integer, intent(in) :: period
      real(dp), intent(in) :: spot
      real(dp) :: value

      if (.not. self%ok) then
         value = 0.0_dp
         return
      end if
      if (period < 1 .or. period > size(self%c, 1)) then
         value = 0.0_dp
         return
      end if
      value = interpolation1d(spot, self%c(period, :), self%s(1), self%s(size(self%s)))
   end function option_value_at

   function auxiliary_at(self, period, spot) result(value)
      class(hedging_result), intent(in) :: self
      integer, intent(in) :: period
      real(dp), intent(in) :: spot
      real(dp) :: value

      if (.not. self%ok) then
         value = 0.0_dp
         return
      end if
      if (period < 1 .or. period > size(self%a, 1)) then
         value = 0.0_dp
         return
      end if
      value = interpolation1d(spot, self%a(period, :), self%s(1), self%s(size(self%s)))
   end function auxiliary_at

   function shares_at(self, period, spot, portfolio_value) result(shares)
      class(hedging_result), intent(in) :: self
      integer, intent(in) :: period
      real(dp), intent(in) :: spot
      real(dp), intent(in) :: portfolio_value
      real(dp) :: shares

      if (.not. self%ok .or. abs(spot) <= tiny(1.0_dp)) then
         shares = 0.0_dp
         return
      end if
      shares = (self%auxiliary_at(period, spot) - portfolio_value * self%rho) / spot
   end function shares_at

   function initial_hedge_at(self, spot) result(shares)
      class(hedging_result), intent(in) :: self
      real(dp), intent(in) :: spot
      real(dp) :: shares
      real(dp) :: value

      value = self%option_value_at(1, spot)
      shares = self%shares_at(1, spot, value)
   end function initial_hedge_at

end module opthedging_types
