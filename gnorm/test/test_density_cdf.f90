! SPDX-License-Identifier: GPL-2.0-or-later
program test_density_cdf
   use gnorm
   implicit none
   real(dp), parameter :: tol = 3.0e-13_dp
   real(dp) :: x(3), y(3)

   call check_close(dgnorm(0.0_dp, alpha=sqrt(2.0_dp), beta=2.0_dp), &
      0.3989422804014327_dp, tol, 'normal density')
   call check_close(pgnorm(1.0_dp, alpha=sqrt(2.0_dp), beta=2.0_dp), &
      0.8413447460685429_dp, tol, 'normal cdf')
   call check_close(pgnorm(-1.0_dp, alpha=sqrt(2.0_dp), beta=2.0_dp), &
      0.1586552539314571_dp, tol, 'normal lower tail')
   call check_close(pgnorm(1.0_dp, alpha=sqrt(2.0_dp), beta=2.0_dp, &
      lower_tail=.false.), 0.1586552539314571_dp, tol, 'normal upper tail')

   call check_close(dgnorm(3.0_dp, mu=3.0_dp, alpha=2.0_dp, beta=1.0_dp), &
      0.25_dp, tol, 'laplace density')
   call check_close(pgnorm(3.0_dp, mu=3.0_dp, alpha=2.0_dp, beta=1.0_dp), &
      0.5_dp, tol, 'laplace center')
   call check_close(pgnorm(3.0_dp + 2.0_dp*log(2.0_dp), mu=3.0_dp, &
      alpha=2.0_dp, beta=1.0_dp), 0.75_dp, tol, 'laplace cdf')

   x = [-1.0_dp, 0.0_dp, 1.0_dp]
   y = pgnorm(x, alpha=sqrt(2.0_dp), beta=2.0_dp)
   call check_close(y(1) + y(3), 1.0_dp, tol, 'elemental symmetry')
   call check_close(y(2), 0.5_dp, tol, 'elemental center')


   call check_close(dgnorm(0.7_dp, mu=-0.2_dp, alpha=1.3_dp, beta=0.5_dp), &
      0.08368387938930931_dp, tol, 'shape 0.5 density')
   call check_close(pgnorm(1.1_dp, mu=0.4_dp, alpha=0.8_dp, beta=1.3_dp), &
      0.8388803569867330_dp, 5.0e-13_dp, 'shape 1.3 cdf')
   call check_close(pgnorm(-0.3_dp, mu=1.2_dp, alpha=2.1_dp, beta=5.0_dp), &
      0.1224948082441346_dp, 5.0e-13_dp, 'shape 5 cdf')

   call check_close(exp(dgnorm(0.7_dp, mu=0.2_dp, alpha=1.3_dp, beta=1.7_dp, &
      log_density=.true.)), dgnorm(0.7_dp, mu=0.2_dp, alpha=1.3_dp, beta=1.7_dp), &
      tol, 'log density')
   call check_close(exp(pgnorm(-4.0_dp, alpha=1.2_dp, beta=0.8_dp, &
      log_probability=.true.)), pgnorm(-4.0_dp, alpha=1.2_dp, beta=0.8_dp), &
      2.0e-13_dp, 'log cdf')

   print '(a)', 'test_density_cdf: PASS'

contains

   subroutine check_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual, expected, tolerance
      character(*), intent(in) :: label
      if (abs(actual - expected) > tolerance * max(1.0_dp, abs(expected))) then
         write(*, '(a,2es24.15)') 'FAIL: '//label//' ', actual, expected
         error stop 1
      end if
   end subroutine check_close

end program test_density_cdf
