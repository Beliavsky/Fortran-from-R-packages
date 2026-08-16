program basic_cmp
   use compoissonreg
   implicit none
   integer, parameter :: n = 10
   integer :: y(n)
   real(dp) :: x(n,1), s(n,1)
   type(cmp_fit_t) :: fit
   type(cmp_init_t) :: init
   type(cmp_fixed_t) :: fixed

   print '(a,f12.8)', 'P(CMP(3,0.8) <= 4) = ', pcmp(4,3.0_dp,0.8_dp)
   print '(a,f12.8)', 'E[CMP(3,0.8)]       = ', ecmp(3.0_dp,0.8_dp)

   y = [0,1,2,1,3,2,4,1,0,2]
   x = 1.0_dp
   s = 1.0_dp
   init = default_init(1,1)
   fixed = default_fixed(1,1)
   fixed%gamma = .true.  ! nu = exp(0) = 1: Poisson submodel

   call fit_cmp_raw(y,x,s,fit,init,fixed)
   print '(a,f12.8)', 'Poisson-submodel lambda = ', exp(fit%beta(1))
   print '(a,l1)', 'Converged = ', fit%converged
end program basic_cmp
