module directional_special
   use directional_kinds, only : dp, pi
   implicit none
   private
   public :: log_i0, log_bessel_i, normal_pdf, normal_cdf, log_beta
contains
   pure real(dp) function log_i0(x) result(v)
      real(dp), intent(in) :: x
      real(dp) :: ax, y, p
      ax = abs(x)
      if (ax < 3.75_dp) then
         y = (x / 3.75_dp)**2
         p = 1.0_dp + y*(3.5156229_dp + y*(3.0899424_dp + y*(1.2067492_dp + &
             y*(0.2659732_dp + y*(0.0360768_dp + y*0.0045813_dp)))))
         v = log(p)
      else
         y = 3.75_dp / ax
         p = 0.39894228_dp + y*(0.01328592_dp + y*(0.00225319_dp + y*(-0.00157565_dp + &
             y*(0.00916281_dp + y*(-0.02057706_dp + y*(0.02635537_dp + &
             y*(-0.01647633_dp + y*0.00392377_dp)))))))
         v = ax - 0.5_dp*log(ax) + log(p)
      end if
   end function log_i0

   pure real(dp) function log_bessel_i(nu, x) result(v)
      real(dp), intent(in) :: nu, x
      real(dp) :: term, s, xx
      integer :: k
      if (x < 0.0_dp) then
         v = -huge(1.0_dp); return
      else if (x == 0.0_dp) then
         if (nu == 0.0_dp) then; v = 0.0_dp; else; v = -huge(1.0_dp); end if
         return
      else if (nu == 0.0_dp) then
         v = log_i0(x); return
      end if
      if (x > 40.0_dp) then
         v = x - 0.5_dp*log(2.0_dp*pi*x) - (4.0_dp*nu*nu - 1.0_dp)/(8.0_dp*x)
         return
      end if
      xx = 0.25_dp*x*x
      term = exp(nu*log(0.5_dp*x) - log_gamma(nu + 1.0_dp))
      s = term
      do k = 1, 500
         term = term*xx/(real(k,dp)*(nu + real(k,dp)))
         s = s + term
         if (abs(term) <= 1.0e-15_dp*abs(s)) exit
      end do
      v = log(s)
   end function log_bessel_i

   pure real(dp) function normal_pdf(x) result(v)
      real(dp), intent(in) :: x
      v = exp(-0.5_dp*x*x)/sqrt(2.0_dp*pi)
   end function normal_pdf

   pure real(dp) function normal_cdf(x) result(v)
      real(dp), intent(in) :: x
      v = 0.5_dp*erfc(-x/sqrt(2.0_dp))
   end function normal_cdf

   pure real(dp) function log_beta(a,b) result(v)
      real(dp), intent(in) :: a,b
      v = log_gamma(a)+log_gamma(b)-log_gamma(a+b)
   end function log_beta
end module directional_special
