! SPDX-License-Identifier: GPL-2.0-or-later
program test_time_series
   use fints
   implicit none

   real(dp) :: x(5), arch_data(240), matrix_data(8, 2)
   type(acf_result) :: corr, covar, partial
   type(cross_acf_result) :: matrix_acf
   type(test_result) :: box, arch
   type(arma_acf_result) :: arma
   integer :: i

   x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
   call acf(x, corr, lag_max=3, acf_type='correlation')
   call check(corr%status == fints_ok, 'acf status')
   call check_close(corr%value(1), 1.0_dp, 1.0e-13_dp, 'acf lag zero')
   call check_close(corr%value(2), 0.4_dp, 1.0e-13_dp, 'acf lag one')
   call check_close(corr%value(3), -0.1_dp, 1.0e-13_dp, 'acf lag two')

   call acf(x, covar, lag_max=1, acf_type='covariance')
   call check_close(covar%value(1), 2.0_dp, 1.0e-13_dp, 'covariance lag zero')
   call check_close(covar%value(2), 0.8_dp, 1.0e-13_dp, 'covariance lag one')

   call acf(x, partial, lag_max=2, acf_type='partial')
   call check_close(partial%value(1), 0.4_dp, 1.0e-13_dp, 'pacf lag one')
   call check_close(partial%value(2), -0.3095238095238095_dp, 1.0e-12_dp, 'pacf lag two')

   call AutocorTest(x, box, lag=2, test_type='Box-Pierce', degrees_freedom=2.0_dp)
   call check_close(box%statistic, 0.85_dp, 1.0e-13_dp, 'box pierce statistic')
   call check(box%p_value > 0.0_dp .and. box%p_value < 1.0_dp, 'box p value')

   arch_data(1) = 0.2_dp
   do i = 2, size(arch_data)
      arch_data(i) = (0.25_dp + 0.70_dp * arch_data(i - 1) ** 2) * &
         (0.6_dp * sin(0.47_dp * real(i, dp)) + 0.8_dp * cos(0.19_dp * real(i, dp)))
   end do
   call ArchTest(arch_data, arch, lags=4, demean=.true.)
   call check(arch%status == fints_ok, 'arch status')
   call check(arch%statistic >= 0.0_dp, 'arch statistic')
   call check(arch%p_value >= 0.0_dp .and. arch%p_value <= 1.0_dp, 'arch p value')

   call arma_true_acf([0.8_dp], [real(dp) ::], 6, arma)
   call check(arma%status == fints_ok .and. arma%stationary, 'ar1 status')
   do i = 0, 6
      call check_close(arma%value(i + 1), 0.8_dp ** i, 2.0e-11_dp, 'ar1 acf')
   end do
   call check_close(real(arma%roots(1), dp), 0.8_dp, 1.0e-11_dp, 'ar1 root')

   call arma_true_acf([real(dp) ::], [0.5_dp], 4, arma)
   call check_close(arma%value(2), 0.4_dp, 1.0e-12_dp, 'ma1 acf one')
   call check(maxval(abs(arma%value(3:))) < 1.0e-12_dp, 'ma1 acf tail')

   call arma_true_acf([0.6_dp, -0.4_dp], [real(dp) ::], 5, arma)
   call check(arma%stationary .and. size(arma%period) == 1, 'complex roots periodicity')
   call check(arma%period(1) > 2.0_dp, 'positive periodicity')

   call arma_true_acf([1.1_dp], [real(dp) ::], 4, arma)
   call check(arma%status == fints_nonstationary, 'nonstationary detection')

   do i = 1, 8
      matrix_data(i, 1) = real(i, dp)
      matrix_data(i, 2) = 2.0_dp * real(i, dp) + 1.0_dp
   end do
   call cross_acf(matrix_data, matrix_acf, lag_max=2)
   call check(matrix_acf%status == fints_ok, 'cross acf status')
   call check_close(matrix_acf%value(1, 1, 2), 1.0_dp, 1.0e-13_dp, &
      'cross correlation at zero')

   print '(a)', 'test_time_series: PASS'

contains

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*, '(a,a)') 'FAIL: ', trim(label)
         error stop 1
      end if
   end subroutine check

   subroutine check_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: label
      call check(abs(actual - expected) <= tolerance * max(1.0_dp, abs(expected)), label)
   end subroutine check_close

end program test_time_series
