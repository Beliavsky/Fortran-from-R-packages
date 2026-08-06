program basic_fit
   use yieldcurve, only : dp, yc_success, nelson_siegel_fit, ns_rates
   implicit none

   real(dp), parameter :: maturity(8) = [0.25_dp, 0.5_dp, 1.0_dp, 2.0_dp, &
      3.0_dp, 5.0_dp, 7.0_dp, 10.0_dp]
   real(dp), parameter :: observed(8) = [2.72_dp, 2.92_dp, 3.23_dp, 3.64_dp, &
      3.86_dp, 4.04_dp, 4.08_dp, 4.08_dp]
   real(dp) :: coeff(4), fitted(8)
   integer :: i, stat
   character(len=160) :: message

   call nelson_siegel_fit(observed, maturity, coeff, stat, message)
   if (stat /= yc_success) error stop trim(message)
   call ns_rates(coeff, maturity, fitted, stat, message)
   if (stat /= yc_success) error stop trim(message)

   write(*, '(a, 4(1x, f10.6))') 'beta0 beta1 beta2 lambda:', coeff
   write(*, '(a)') ' maturity   observed     fitted'
   do i = 1, size(maturity)
      write(*, '(3f12.6)') maturity(i), observed(i), fitted(i)
   end do
end program basic_fit
