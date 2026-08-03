! SPDX-License-Identifier: GPL-2.0-or-later
program test_special
   use gnorm
   implicit none
   real(dp) :: p, x
   integer :: i

   do i = 1, 9
      p = 0.1_dp * real(i, dp)
      x = inverse_regularized_gamma_p(1.0_dp, p)
      call check_close(x, -log(1.0_dp - p), 3.0e-13_dp, 'exponential quantile')
      call check_close(regularized_gamma_p(1.0_dp, x), p, 3.0e-13_dp, &
         'gamma inverse round trip')
   end do
   call check_close(regularized_gamma_p(0.5_dp, 1.0_dp), &
      0.8427007929497149_dp, 3.0e-13_dp, 'half gamma cdf')
   call check_close(regularized_gamma_p(4.0_dp, 2.0_dp) + &
      regularized_gamma_q(4.0_dp, 2.0_dp), 1.0_dp, 3.0e-14_dp, &
      'gamma complement')
   call check_close(gnorm_variance(alpha=sqrt(2.0_dp), beta=2.0_dp), &
      1.0_dp, 3.0e-14_dp, 'normal variance')
   call check_close(gnorm_variance(alpha=2.0_dp, beta=1.0_dp), &
      8.0_dp, 3.0e-14_dp, 'laplace variance')

   print '(a)', 'test_special: PASS'

contains

   subroutine check_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual, expected, tolerance
      character(*), intent(in) :: label
      if (abs(actual - expected) > tolerance * max(1.0_dp, abs(expected))) then
         write(*, '(a,2es24.15)') 'FAIL: '//label//' ', actual, expected
         error stop 1
      end if
   end subroutine check_close

end program test_special
