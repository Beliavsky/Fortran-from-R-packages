! SPDX-License-Identifier: Artistic-2.0
module ldhmm_series
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use ldhmm_kinds, only : dp
   use ldhmm_math, only : quiet_nan, sample_sd
   implicit none
   private

   public :: ldhmm_prices_to_log_returns

contains

   function ldhmm_prices_to_log_returns(prices, randomize_zeros, zero_scale) result(returns)
      real(dp), intent(in) :: prices(:)
      logical, intent(in), optional :: randomize_zeros
      real(dp), intent(in), optional :: zero_scale
      real(dp), allocatable :: returns(:)
      logical :: add_noise
      real(dp) :: scale, sigma, u1, u2, noise
      integer :: i

      if (size(prices) < 2) then
         allocate(returns(0))
         return
      end if
      add_noise = .false.
      scale = 0.01_dp
      if (present(randomize_zeros)) add_noise = randomize_zeros
      if (present(zero_scale)) scale = zero_scale
      allocate(returns(size(prices)-1))
      do i = 2, size(prices)
         if (prices(i) > 0.0_dp .and. prices(i-1) > 0.0_dp .and. &
             ieee_is_finite(prices(i)) .and. ieee_is_finite(prices(i-1))) then
            returns(i-1) = log(prices(i)/prices(i-1))
         else
            returns(i-1) = quiet_nan()
         end if
      end do
      if (.not. add_noise) return
      sigma = sample_sd(returns)
      if (.not. ieee_is_finite(sigma)) return
      do i = 1, size(returns)
         if (abs(returns(i)) > tiny(1.0_dp)) cycle
         call random_number(u1)
         call random_number(u2)
         u1 = max(u1, 1.0e-300_dp)
         noise = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*acos(-1.0_dp)*u2)
         returns(i) = scale*sigma*noise
      end do
   end function ldhmm_prices_to_log_returns

end module ldhmm_series
