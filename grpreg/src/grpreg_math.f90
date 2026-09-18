module grpreg_math
   use grpreg_kinds, only : dp
   implicit none
   private
   public :: soft_threshold, firm_threshold, scad_threshold
   public :: mcp_penalty, dmcp, logistic, log1pexp, normal_cdf
   public :: regularized_gamma_p, chi_square_cdf

contains

   pure elemental real(dp) function soft_threshold(z, lambda) result(value)
      real(dp), intent(in) :: z !! Scalar coefficient or group norm to threshold.
      real(dp), intent(in) :: lambda !! Nonnegative threshold magnitude.
      if (z > lambda) then
         value = z - lambda
      else if (z < -lambda) then
         value = z + lambda
      else
         value = 0.0_dp
      end if
   end function soft_threshold

   pure elemental real(dp) function firm_threshold(z, lambda1, lambda2, gamma) result(value)
      real(dp), intent(in) :: z !! Scalar argument to the MCP firm-threshold map.
      real(dp), intent(in) :: lambda1 !! MCP sparsity penalty, nonnegative.
      real(dp), intent(in) :: lambda2 !! Ridge component of the penalty, nonnegative.
      real(dp), intent(in) :: gamma !! MCP concavity parameter, required greater than one.
      real(dp) :: s
      s = sign(1.0_dp, z)
      if (abs(z) <= lambda1) then
         value = 0.0_dp
      else if (abs(z) <= gamma*lambda1*(1.0_dp + lambda2)) then
         value = s*(abs(z) - lambda1)/(1.0_dp + lambda2 - 1.0_dp/gamma)
      else
         value = z/(1.0_dp + lambda2)
      end if
   end function firm_threshold

   pure elemental real(dp) function scad_threshold(z, lambda1, lambda2, gamma) result(value)
      real(dp), intent(in) :: z !! Scalar argument to the SCAD threshold map.
      real(dp), intent(in) :: lambda1 !! SCAD sparsity penalty, nonnegative.
      real(dp), intent(in) :: lambda2 !! Ridge component of the penalty, nonnegative.
      real(dp), intent(in) :: gamma !! SCAD concavity parameter, required greater than two.
      real(dp) :: s
      s = sign(1.0_dp, z)
      if (abs(z) <= lambda1) then
         value = 0.0_dp
      else if (abs(z) <= lambda1*(2.0_dp + lambda2)) then
         value = s*(abs(z) - lambda1)/(1.0_dp + lambda2)
      else if (abs(z) <= gamma*lambda1*(1.0_dp + lambda2)) then
         value = s*(abs(z) - gamma*lambda1/(gamma - 1.0_dp))/ &
            (1.0_dp - 1.0_dp/(gamma - 1.0_dp) + lambda2)
      else
         value = z/(1.0_dp + lambda2)
      end if
   end function scad_threshold

   pure elemental real(dp) function mcp_penalty(theta, lambda, gamma) result(value)
      real(dp), intent(in) :: theta !! Coefficient magnitude argument to the MCP penalty.
      real(dp), intent(in) :: lambda !! MCP tuning parameter, nonnegative.
      real(dp), intent(in) :: gamma !! MCP concavity parameter, greater than one.
      real(dp) :: a
      a = abs(theta)
      if (a <= gamma*lambda) then
         value = lambda*a - a*a/(2.0_dp*gamma)
      else
         value = gamma*lambda*lambda/2.0_dp
      end if
   end function mcp_penalty

   pure elemental real(dp) function dmcp(theta, lambda, gamma) result(value)
      real(dp), intent(in) :: theta !! Coefficient magnitude argument to the MCP derivative.
      real(dp), intent(in) :: lambda !! MCP tuning parameter, nonnegative.
      real(dp), intent(in) :: gamma !! MCP concavity parameter, greater than one.
      if (abs(theta) < gamma*lambda) then
         value = lambda - abs(theta)/gamma
      else
         value = 0.0_dp
      end if
   end function dmcp

   pure elemental real(dp) function logistic(eta) result(value)
      real(dp), intent(in) :: eta !! Linear predictor for the logistic inverse link.
      if (eta >= 0.0_dp) then
         value = 1.0_dp/(1.0_dp + exp(-min(eta, 700.0_dp)))
      else
         value = exp(max(eta, -700.0_dp))/(1.0_dp + exp(max(eta, -700.0_dp)))
      end if
   end function logistic

   pure elemental real(dp) function log1pexp(x) result(value)
      real(dp), intent(in) :: x !! Real argument for stable log(1+exp(x)).
      if (x > 0.0_dp) then
         value = x + log(1.0_dp + exp(-x))
      else
         value = log(1.0_dp + exp(x))
      end if
   end function log1pexp

   pure elemental real(dp) function normal_cdf(x) result(value)
      real(dp), intent(in) :: x !! Standard-normal quantile at which to evaluate the CDF.
      value = 0.5_dp*erfc(-x/sqrt(2.0_dp))
   end function normal_cdf

   pure real(dp) function regularized_gamma_p(a, x) result(value)
      real(dp), intent(in) :: a !! Positive shape parameter of the lower incomplete gamma ratio.
      real(dp), intent(in) :: x !! Nonnegative evaluation point of the lower incomplete gamma ratio.
      integer, parameter :: itmax = 500
      real(dp), parameter :: eps = 1.0e-14_dp
      real(dp), parameter :: fpmin = tiny(1.0_dp)/eps
      integer :: n
      real(dp) :: ap, del, summ, b, c, d, h, an, gln
      if (a <= 0.0_dp .or. x < 0.0_dp) then
         value = 0.0_dp
         return
      end if
      if (x <= tiny(1.0_dp)) then
         value = 0.0_dp
         return
      end if
      gln = log_gamma(a)
      if (x < a + 1.0_dp) then
         ap = a
         summ = 1.0_dp/a
         del = summ
         do n = 1, itmax
            ap = ap + 1.0_dp
            del = del*x/ap
            summ = summ + del
            if (abs(del) <= abs(summ)*eps) exit
         end do
         value = summ*exp(-x + a*log(x) - gln)
      else
         b = x + 1.0_dp - a
         c = 1.0_dp/fpmin
         d = 1.0_dp/b
         h = d
         do n = 1, itmax
            an = -real(n, dp)*(real(n, dp) - a)
            b = b + 2.0_dp
            d = an*d + b
            if (abs(d) < fpmin) d = fpmin
            c = b + an/c
            if (abs(c) < fpmin) c = fpmin
            d = 1.0_dp/d
            del = d*c
            h = h*del
            if (abs(del - 1.0_dp) <= eps) exit
         end do
         value = 1.0_dp - exp(-x + a*log(x) - gln)*h
      end if
      value = min(max(value, 0.0_dp), 1.0_dp)
   end function regularized_gamma_p

   pure real(dp) function chi_square_cdf(x, df) result(value)
      real(dp), intent(in) :: x !! Nonnegative chi-square variate.
      real(dp), intent(in) :: df !! Positive chi-square degrees of freedom.
      if (x <= 0.0_dp) then
         value = 0.0_dp
      else
         value = regularized_gamma_p(0.5_dp*df, 0.5_dp*x)
      end if
   end function chi_square_cdf

end module grpreg_math
