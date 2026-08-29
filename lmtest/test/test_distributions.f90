program test_distributions
   use lmtest, only : dp, normal_cdf, normal_quantile, chi_square_sf, &
      f_sf, student_t_cdf, student_t_quantile
   implicit none

   call check(normal_cdf(1.96_dp), 0.9750021048517795_dp, 2.0e-14_dp)
   call check(normal_quantile(0.975_dp), 1.959963984540054_dp, 3.0e-9_dp)
   call check(chi_square_sf(7.271207977024421_dp, 2.0_dp), &
      0.02636800360563612_dp, 2.0e-13_dp)
   call check(f_sf(1.641760171679299_dp, 37.0_dp, 37.0_dp), &
      0.06809147841060445_dp, 3.0e-13_dp)
   call check(student_t_cdf(1.4066348820050234_dp, 76.0_dp), &
      0.9181941395155292_dp, 2.0e-10_dp)
   call check(student_t_quantile(0.975_dp, 77.0_dp), &
      1.9912543953907817_dp, 2.0e-10_dp)
contains
   subroutine check(actual, expected, tol)
      real(dp), intent(in) :: actual, expected, tol
      if (abs(actual - expected) > tol) then
         print *, actual, expected
         error stop 1
      end if
   end subroutine check
end program test_distributions
