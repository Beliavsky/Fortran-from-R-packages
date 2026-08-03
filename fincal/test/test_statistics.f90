! SPDX-License-Identifier: GPL-2.0-or-later
program test_statistics
   use fincal
   implicit none
   real(dp), parameter :: tol = 2.0e-12_dp
   integer :: status

   call assert_close(geometric_mean([-0.0934_dp, 0.2345_dp, 0.0892_dp]), 0.068246504876591_dp, 1.0e-10_dp, 'geometric mean')
   call assert_close(harmonic_mean([8.0_dp, 9.0_dp, 10.0_dp]), 8.92561983471_dp, 1.0e-10_dp, 'harmonic mean')
   call assert_close(sampling_error(0.45_dp, 0.5_dp), -0.05_dp, tol, 'sampling error')
   call assert_close(weighted_portfolio_return([0.12_dp, 0.07_dp, 0.03_dp], [0.5_dp, 0.4_dp, 0.1_dp], status), &
      0.091_dp, tol, 'weighted return')
   call assert_equal_int(status, fincal_ok, 'weighted return status')
   call assert_close(wpr([0.12_dp, 0.07_dp, 0.03_dp], [0.5_dp, 0.4_dp, 0.1_dp]), &
      0.091_dp, tol, 'wpr alias')
   call assert_close(twrr([120.0_dp, 260.0_dp], [100.0_dp, 240.0_dp], [2.0_dp, 4.0_dp]), &
      0.158447236605966_dp, 1.0e-10_dp, 'twrr')
   print '(a)', 'test_statistics: PASS'
contains
   subroutine assert_close(actual, expected, tolerance, message)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: message
      if (abs(actual - expected) > tolerance) then
         print *, trim(message), actual, expected, abs(actual - expected)
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_equal_int(actual, expected, message)
      integer, intent(in) :: actual, expected
      character(len=*), intent(in) :: message
      if (actual /= expected) then
         print *, trim(message), actual, expected
         error stop 1
      end if
   end subroutine assert_equal_int
end program test_statistics
