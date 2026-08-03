! SPDX-License-Identifier: GPL-2.0-or-later
program test_drawdowns
   use jfe
   implicit none

   real(dp), parameter :: r(5) = [0.02_dp, -0.01_dp, 0.03_dp, -0.02_dp, 0.01_dp]
   real(dp), parameter :: p(6) = [2.0_dp, -1.0_dp, -2.0_dp, 3.0_dp, -4.0_dp, 5.0_dp]
   real(dp), parameter :: expected_dd(5) = [0.0_dp, -0.01_dp, 0.0_dp, -0.02_dp, -0.0102_dp]
   real(dp), parameter :: expected_peak(6) = [0.0_dp, -1.0_dp, -2.98_dp, -0.0694_dp, &
      -4.066624_dp, 0.0_dp]
   real(dp), allocatable :: dd(:), peak(:)
   real(dp) :: ordinary, modified

   dd = drawdowns(r)
   call check_vector(dd, expected_dd, 1.0e-13_dp, 'standard drawdowns')
   call check_close(max_drawdown(r), 0.02_dp, 1.0e-13_dp, 'maximum drawdown')
   call check_close(calmar_ratio(r, 12.0_dp), 3.611517662489302_dp, 1.0e-12_dp, &
      'calmar ratio')
   call check_close(sterling_ratio(r, 12.0_dp), 0.902879415622326_dp, 1.0e-12_dp, &
      'sterling ratio')

   peak = drawdown_peak(p)
   call check_vector(peak, expected_peak, 1.0e-10_dp, 'peak drawdowns')
   call check_close(ulcer_index(p), 2.098517219267611_dp, 1.0e-12_dp, 'ulcer index')
   call check_close(pain_index(p), 1.3526706666666666_dp, 1.0e-12_dp, 'pain index')

   ordinary = burke_ratio(r, scale=12.0_dp)
   modified = burke_ratio(r, scale=12.0_dp, modified=.true.)
   call check_close(modified, ordinary*sqrt(real(size(r), dp)), 1.0e-12_dp, &
      'modified Burke scaling')

   print '(a)', 'test_drawdowns: PASS'

contains

   subroutine check_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual, expected, tolerance
      character(*), intent(in) :: label
      if (abs(actual - expected) > tolerance) then
         print '(a,2(1x,es24.16))', 'FAIL: '//trim(label), actual, expected
         error stop 1
      end if
   end subroutine check_close

   subroutine check_vector(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual(:), expected(:), tolerance
      character(*), intent(in) :: label
      if (size(actual) /= size(expected)) error stop 'FAIL: vector size'
      if (maxval(abs(actual - expected)) > tolerance) then
         print '(a)', 'FAIL: '//trim(label)
         print '(a,*(1x,es14.6))', 'actual:  ', actual
         print '(a,*(1x,es14.6))', 'expected:', expected
         error stop 1
      end if
   end subroutine check_vector

end program test_drawdowns
