program basic
   use expint_fortran
   implicit none
   integer :: n
   real(dp) :: x

   x = 1.275_dp
   print '(a,f8.3)', 'x = ', x
   print '(a)', ' n          E_n(x)'
   do n = 1, 10
      print '(i2,2x,es16.8)', n, expint(x, n)
   end do

   print '(/,a,es16.8)', 'Ei(0.5)       = ', expint_ei(0.5_dp)
   print '(a,es16.8)', 'Gamma(-1.2,2.5) = ', gammainc(-1.2_dp, 2.5_dp)
end program basic
