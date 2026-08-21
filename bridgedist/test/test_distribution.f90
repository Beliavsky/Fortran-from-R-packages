program test_distribution
   use bridgedist, only : dp, dbridge, pbridge, qbridge, bridge_mean, bridge_variance
   implicit none
   real(dp), parameter :: tol = 5.0e-13_dp
   integer :: fails

   fails = 0
   call check(dbridge(-2.0_dp, 0.25_dp), 0.061338382792202408115_dp, tol, 'd1')
   call check(pbridge(-2.0_dp, 0.25_dp), 0.37127200025513721534_dp, tol, 'p1')
   call check(dbridge(0.0_dp, 0.5_dp), 0.15915494309189533577_dp, tol, 'd2')
   call check(pbridge(0.0_dp, 0.5_dp), 0.5_dp, tol, 'p2')
   call check(dbridge(1.3_dp, 0.7_dp), 0.15048307190253310466_dp, tol, 'd3')
   call check(pbridge(1.3_dp, 0.7_dp), 0.81665139320514155212_dp, tol, 'p3')
   call check(qbridge(0.01_dp, 0.2_dp), -22.64892824522124318_dp, 2.0e-12_dp, 'q1')
   call check(qbridge(0.37_dp, 0.5_dp), -0.84051630812518878126_dp, tol, 'q2')
   call check(qbridge(0.9_dp, 0.8_dp), 1.4135626933331296441_dp, tol, 'q3')
   call check(bridge_mean(0.6_dp), 0.0_dp, 0.0_dp, 'mean')
   call check(bridge_variance(0.5_dp), acos(-1.0_dp)**2, tol, 'variance')

   if (fails /= 0) error stop 1
   print '(a)', 'test_distribution: PASS'

contains

   subroutine check(got, expected, atol, label)
      real(dp), intent(in) :: got, expected, atol
      character(len=*), intent(in) :: label
      if (abs(got - expected) > atol * max(1.0_dp, abs(expected))) then
         print '(a,2(1x,es24.16))', trim(label)//' FAIL', got, expected
         fails = fails + 1
      end if
   end subroutine check

end program test_distribution
