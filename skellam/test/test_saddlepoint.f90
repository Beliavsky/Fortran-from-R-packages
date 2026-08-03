program test_saddlepoint
   use skellam, only : dp, i8, dskellam, pskellam, dskellam_sp, pskellam_sp
   implicit none

   integer(i8) :: k
   real(dp) :: exact, approximate

   do k = -8_i8, 12_i8, 4_i8
      exact = dskellam(k, 12.0_dp, 8.0_dp)
      approximate = dskellam_sp(k, 12.0_dp, 8.0_dp)
      if (abs(approximate - exact)/max(exact, tiny(1.0_dp)) > 0.012_dp) then
         print *, 'saddlepoint PMF failed:', k, exact, approximate
         error stop 1
      end if
      exact = pskellam(real(k, dp), 12.0_dp, 8.0_dp)
      approximate = pskellam_sp(real(k, dp), 12.0_dp, 8.0_dp)
      if (abs(approximate - exact) > 0.006_dp) then
         print *, 'saddlepoint CDF failed:', k, exact, approximate
         error stop 1
      end if
   end do
   call assert_close(log(dskellam_sp(2_i8, 12.0_dp, 8.0_dp)), &
      dskellam_sp(2_i8, 12.0_dp, 8.0_dp, log_p=.true.), 1.0e-13_dp)
   print '(a)', 'test_saddlepoint: PASS'

contains

   subroutine assert_close(actual, expected, tolerance)
      real(dp), intent(in) :: actual, expected, tolerance
      if (abs(actual - expected) > tolerance) error stop 'log saddlepoint failed'
   end subroutine assert_close

end program test_saddlepoint
