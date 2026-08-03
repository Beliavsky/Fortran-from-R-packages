! SPDX-License-Identifier: GPL-2.0-or-later
program test_quantile
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
   use gnorm
   implicit none
   real(dp), parameter :: probabilities(7) = [1.0e-8_dp, 0.01_dp, 0.1_dp, &
      0.5_dp, 0.9_dp, 0.99_dp, 1.0_dp - 1.0e-8_dp]
   real(dp), parameter :: shapes(4) = [0.5_dp, 1.0_dp, 2.0_dp, 5.0_dp]
   real(dp) :: quantile, back
   integer :: i, j

   call check_close(qgnorm(0.975_dp, alpha=sqrt(2.0_dp), beta=2.0_dp), &
      1.959963984540054_dp, 2.0e-12_dp, 'normal quantile')
   call check_close(qgnorm(0.75_dp, mu=3.0_dp, alpha=2.0_dp, beta=1.0_dp), &
      3.0_dp + 2.0_dp*log(2.0_dp), 2.0e-13_dp, 'laplace quantile')
   call check_close(qgnorm(log(0.025_dp), alpha=sqrt(2.0_dp), beta=2.0_dp, &
      lower_tail=.false., log_probability=.true.), 1.959963984540054_dp, &
      2.0e-12_dp, 'log upper-tail quantile')


   call check_close(qgnorm(1.0e-6_dp, mu=-0.2_dp, alpha=1.3_dp, beta=0.5_dp), &
      -331.0392170434414_dp, 2.0e-10_dp, 'shape 0.5 tail quantile')
   call check_close(qgnorm(0.123_dp, mu=0.4_dp, alpha=0.8_dp, beta=1.3_dp), &
      -0.4463685927720905_dp, 3.0e-12_dp, 'shape 1.3 quantile')
   call check_close(qgnorm(0.999_dp, mu=1.2_dp, alpha=2.1_dp, beta=5.0_dp), &
      3.8998896227755333_dp, 3.0e-12_dp, 'shape 5 quantile')

   do j = 1, size(shapes)
      do i = 1, size(probabilities)
         quantile = qgnorm(probabilities(i), mu=0.4_dp, alpha=1.7_dp, &
            beta=shapes(j))
         back = pgnorm(quantile, mu=0.4_dp, alpha=1.7_dp, beta=shapes(j))
         call check_close(back, probabilities(i), 2.0e-10_dp, 'round trip')
      end do
   end do

   if (ieee_is_finite(qgnorm(0.0_dp))) error stop 'FAIL: lower endpoint'
   if (ieee_is_finite(qgnorm(1.0_dp))) error stop 'FAIL: upper endpoint'
   if (.not. ieee_is_nan(qgnorm(-0.1_dp))) error stop 'FAIL: invalid probability'
   if (.not. ieee_is_nan(qgnorm(0.5_dp, alpha=0.0_dp))) error stop 'FAIL: invalid scale'

   print '(a)', 'test_quantile: PASS'

contains

   subroutine check_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual, expected, tolerance
      character(*), intent(in) :: label
      if (abs(actual - expected) > tolerance * max(1.0_dp, abs(expected))) then
         write(*, '(a,2es24.15)') 'FAIL: '//label//' ', actual, expected
         error stop 1
      end if
   end subroutine check_close

end program test_quantile
