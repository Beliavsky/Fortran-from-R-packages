! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2013 Bruno Remillard
! Modern Fortran translation copyright (C) 2026 OpenAI
module opthedging_iid
   use opthedging_interpolation, only : interpolation1d
   use opthedging_kinds, only : dp
   use opthedging_statistics, only : mean_square, mean_value
   use opthedging_types, only : hedging_result
   implicit none
   private

   public :: call_payoff
   public :: hedging_iid
   public :: put_payoff

contains

   pure function call_payoff(spot, strike) result(value)
      real(dp), intent(in) :: spot
      real(dp), intent(in) :: strike
      real(dp) :: value

      value = max(spot - strike, 0.0_dp)
   end function call_payoff

   pure function put_payoff(spot, strike) result(value)
      real(dp), intent(in) :: spot
      real(dp), intent(in) :: strike
      real(dp) :: value

      value = max(strike - spot, 0.0_dp)
   end function put_payoff

   function hedging_iid(log_excess_returns, maturity, strike, rate, is_put, &
      n_periods, n_grid, min_s, max_s) result(out)
      real(dp), intent(in) :: log_excess_returns(:)
      real(dp), intent(in) :: maturity
      real(dp), intent(in) :: strike
      real(dp), intent(in) :: rate
      logical, intent(in) :: is_put
      integer, intent(in) :: n_periods
      integer, intent(in) :: n_grid
      real(dp), intent(in) :: min_s
      real(dp), intent(in) :: max_s
      type(hedging_result) :: out

      integer :: i
      integer :: j
      integer :: k
      integer :: n_returns
      real(dp) :: c2
      real(dp) :: dx
      real(dp) :: m1
      real(dp) :: m2
      real(dp) :: next_spot
      real(dp) :: payoff
      real(dp) :: sum_a
      real(dp) :: sum_c
      real(dp) :: z
      real(dp), allocatable :: xi(:)

      call initialize_failure(out, "invalid input")
      n_returns = size(log_excess_returns)
      if (n_returns < 2) then
         out%message = "at least two simulated returns are required"
         return
      end if
      if (n_periods < 1) then
         out%message = "n_periods must be positive"
         return
      end if
      if (n_grid < 2) then
         out%message = "n_grid must be at least two"
         return
      end if
      if (maturity < 0.0_dp) then
         out%message = "maturity must be nonnegative"
         return
      end if
      if (strike < 0.0_dp) then
         out%message = "strike must be nonnegative"
         return
      end if
      if (min_s <= 0.0_dp .or. max_s <= min_s) then
         out%message = "the price grid must satisfy 0 < min_s < max_s"
         return
      end if

      allocate(xi(n_returns))
      xi = exp(log_excess_returns) - 1.0_dp
      m1 = mean_value(xi)
      m2 = mean_square(xi)
      if (m2 <= epsilon(1.0_dp) * max(1.0_dp, abs(m1))) then
         out%message = "the return second moment is numerically zero"
         return
      end if

      out%rho = m1 / m2
      c2 = 1.0_dp - m1 * out%rho
      if (abs(c2) <= sqrt(epsilon(1.0_dp))) then
         out%message = "the return sample produces a singular hedging operator"
         return
      end if

      allocate(out%s(n_grid), out%c(n_periods, n_grid))
      allocate(out%a(n_periods, n_grid), out%phi1(n_grid))
      out%c = 0.0_dp
      out%a = 0.0_dp
      out%discounted_strike = strike * exp(-rate * maturity)
      dx = (max_s - min_s) / real(n_grid - 1, dp)
      do i = 1, n_grid
         out%s(i) = min_s + real(i - 1, dp) * dx
      end do

      k = n_periods
      do i = 1, n_grid
         sum_a = 0.0_dp
         sum_c = 0.0_dp
         do j = 1, n_returns
            next_spot = out%s(i) * (1.0_dp + xi(j))
            if (is_put) then
               payoff = put_payoff(next_spot, out%discounted_strike)
            else
               payoff = call_payoff(next_spot, out%discounted_strike)
            end if
            z = (1.0_dp - out%rho * xi(j)) / c2
            sum_a = sum_a + xi(j) * payoff
            sum_c = sum_c + z * payoff
         end do
         out%a(k, i) = sum_a / real(n_returns, dp) / m2
         out%c(k, i) = sum_c / real(n_returns, dp)
      end do

      do k = n_periods - 1, 1, -1
         do i = 1, n_grid
            sum_a = 0.0_dp
            sum_c = 0.0_dp
            do j = 1, n_returns
               next_spot = out%s(i) * (1.0_dp + xi(j))
               payoff = interpolation1d(next_spot, out%c(k + 1, :), min_s, max_s)
               z = (1.0_dp - out%rho * xi(j)) / c2
               sum_a = sum_a + xi(j) * payoff
               sum_c = sum_c + z * payoff
            end do
            out%a(k, i) = sum_a / real(n_returns, dp) / m2
            out%c(k, i) = sum_c / real(n_returns, dp)
         end do
      end do

      out%phi1 = (out%a(1, :) - out%c(1, :) * out%rho) / out%s
      out%ok = .true.
      out%message = "success"
   end function hedging_iid

   subroutine initialize_failure(out, message)
      type(hedging_result), intent(out) :: out
      character(len=*), intent(in) :: message

      out%ok = .false.
      out%rho = 0.0_dp
      out%discounted_strike = 0.0_dp
      out%message = message
   end subroutine initialize_failure

end module opthedging_iid
