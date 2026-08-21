program test_options
   use bridgedist, only : dp, dbridge, pbridge, qbridge
   implicit none
   real(dp) :: x, p, phi
   integer :: fails

   fails = 0
   x = 1.7_dp
   p = 0.23_dp
   phi = 0.72_dp

   call check(dbridge(x, phi, .true.), log(dbridge(x, phi)), 2.0e-14_dp, 'log density')
   call check(pbridge(x, phi, .false.), 1.0_dp - pbridge(x, phi), 2.0e-14_dp, 'upper tail')
   call check(pbridge(x, phi, .true., .true.), log(pbridge(x, phi)), 2.0e-14_dp, 'log cdf')
   call check(pbridge(x, phi, .false., .true.), log(pbridge(x, phi, .false.)), 2.0e-14_dp, 'log upper')
   call check(qbridge(log(p), phi, .true., .true.), qbridge(p, phi), 2.0e-14_dp, 'log quantile')
   call check(qbridge(p, phi, .false.), -qbridge(p, phi), 2.0e-14_dp, 'upper quantile')
   call check(qbridge(log(p), phi, .false., .true.), -qbridge(p, phi), 2.0e-14_dp, 'log upper quantile')

   if (fails /= 0) error stop 1
   print '(a)', 'test_options: PASS'

contains

   subroutine check(got, expected, rtol, label)
      real(dp), intent(in) :: got, expected, rtol
      character(len=*), intent(in) :: label
      if (abs(got - expected) > rtol * max(1.0_dp, abs(expected))) then
         print '(a,2(1x,es24.16))', trim(label)//' FAIL', got, expected
         fails = fails + 1
      end if
   end subroutine check

end program test_options
