program basic
   use rmutil
   implicit none
   real(dp) :: ans, p
   real(dp), allocatable :: nodes(:), weights(:), conc(:)
   integer :: status

   p = pinvgauss(1.4_dp,2.0_dp,0.3_dp)
   ans = integrate_romberg(fsin,0.0_dp,pi,1.0e-10_dp)
   call gauss_hermite(10,nodes,weights,status=status)
   conc = mu1_1o1c([log(10.0_dp),log(2.0_dp),log(0.4_dp)], &
      [0.1_dp,0.5_dp,1.0_dp],2.0_dp)

   print '(a,f12.8)', 'P(IG <= 1.4): ', p
   print '(a,f12.8)', 'Integral sin(x), 0..pi: ', ans
   print '(a,f12.8)', 'Gauss-Hermite weight sum: ', sum(weights)
   print '(a,3f12.8)', 'PK concentrations: ', conc
contains
   function fsin(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y
      y = sin(x)
   end function fsin
end program basic
