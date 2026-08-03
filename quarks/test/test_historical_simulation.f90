program test_historical_simulation
   use quarks
   implicit none
   real(dp), parameter :: tol = 2.0e-12_dp
   real(dp) :: x(8)
   real(dp), allocatable :: variance(:)
   type(risk_result) :: plain, age

   x = [0.01_dp, -0.02_dp, 0.015_dp, -0.03_dp, 0.005_dp, -0.01_dp, &
      0.02_dp, -0.025_dp]
   variance = ewma(x, 0.94_dp)
   call assert_close(variance(1), 0.0003745535714285714_dp, tol, 'EWMA initial')
   call assert_close(variance(8), 0.0003509995610802_dp, tol, 'EWMA final')

   plain = hs(x, 0.75_dp, method_plain)
   call assert_close(plain%var, 0.02125_dp, tol, 'plain VaR')
   call assert_close(plain%es, 0.0275_dp, tol, 'plain ES')

   age = hs(x, 0.75_dp, method_age, 0.98_dp)
   call assert_close(age%var, 0.0202845297113616_dp, tol, 'age VaR')
   call assert_close(age%es, 0.0273990413990211_dp, tol, 'age ES')
   print *, 'test_historical_simulation: PASS'

contains

   subroutine assert_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: label
      if (abs(actual - expected) > tolerance) then
         print *, trim(label), actual, expected
         error stop 1
      end if
   end subroutine assert_close

end program test_historical_simulation
