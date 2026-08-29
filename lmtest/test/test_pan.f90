program test_pan
   use lmtest, only : dp, pan_probability
   implicit none
   real(dp) :: a3(3), a5(5)

   a3 = [1.0_dp, 3.0_dp, 6.0_dp]
   call check(pan_probability(a3, 0.0_dp, 1.0_dp, 15), 0.0542_dp, 3.0e-4_dp)
   call check(pan_probability(a3, 0.0_dp, 7.0_dp, 15), 0.4936_dp, 3.0e-4_dp)
   call check(pan_probability(a3, 0.0_dp, 20.0_dp, 15), 0.8760_dp, 3.0e-4_dp)

   a5 = [1.0_dp, 3.0_dp, 5.0_dp, 7.0_dp, 9.0_dp]
   call check(pan_probability(a5, 0.0_dp, 5.0_dp, 15), 0.0544_dp, 3.0e-4_dp)
   call check(pan_probability(a5, 0.0_dp, 20.0_dp, 15), 0.4853_dp, 3.0e-4_dp)
   call check(pan_probability(a5, 0.0_dp, 50.0_dp, 15), 0.9069_dp, 3.0e-4_dp)
contains
   subroutine check(actual, expected, tol)
      real(dp), intent(in) :: actual, expected, tol
      if (abs(actual - expected) > tol) error stop 1
   end subroutine check
end program test_pan
