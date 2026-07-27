! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from RND 1.2, Copyright (C) 2017 Kam Hamidieh.
module rnd_densities
   use rnd_kinds, only : dp
   use rnd_special, only : beta_pdf, regularized_beta, lognormal_pdf, normal_pdf
   implicit none
   private
   public :: approximate_max, dgb, pgb, dmln, dmln_am, dew, dshimko

contains

   elemental real(dp) function approximate_max(x, y, sharpness) result(value)
      real(dp), intent(in) :: x, y
      real(dp), intent(in), optional :: sharpness
      real(dp) :: k, z, weight
      k = 5.0_dp
      if (present(sharpness)) k = sharpness
      z = max(-700.0_dp, min(700.0_dp, -k*(x-y)))
      weight = 1.0_dp/(1.0_dp+exp(z))
      value = weight*x + (1.0_dp-weight)*y
   end function approximate_max

   elemental real(dp) function pgb(x, a, b, v, w) result(value)
      real(dp), intent(in) :: x, a, b, v, w
      real(dp) :: xa, xnew
      if (x <= 0.0_dp .or. a <= 0.0_dp .or. b <= 0.0_dp) then
         value = 0.0_dp
      else
         xa = (x/b)**a
         xnew = xa/(1.0_dp+xa)
         value = regularized_beta(xnew, v, w)
      end if
   end function pgb

   elemental real(dp) function dgb(x, a, b, v, w) result(value)
      real(dp), intent(in) :: x, a, b, v, w
      real(dp) :: xa, xnew, jacobian
      if (x <= 0.0_dp .or. a <= 0.0_dp .or. b <= 0.0_dp) then
         value = 0.0_dp
      else
         xa = (x/b)**a
         xnew = xa/(1.0_dp+xa)
         jacobian = a*b**a*x**(a-1.0_dp)/(x**a+b**a)**2
         value = beta_pdf(xnew, v, w)*jacobian
      end if
   end function dgb

   elemental real(dp) function dmln(x, alpha1, meanlog1, meanlog2, sdlog1, sdlog2) result(value)
      real(dp), intent(in) :: x, alpha1, meanlog1, meanlog2, sdlog1, sdlog2
      value = alpha1*lognormal_pdf(x, meanlog1, sdlog1) &
         + (1.0_dp-alpha1)*lognormal_pdf(x, meanlog2, sdlog2)
   end function dmln

   elemental real(dp) function dmln_am(x, u1, u2, u3, sigma1, sigma2, sigma3, p1, p2) result(value)
      real(dp), intent(in) :: x, u1, u2, u3, sigma1, sigma2, sigma3, p1, p2
      value = p1*lognormal_pdf(x, u1, sigma1) &
         + p2*lognormal_pdf(x, u2, sigma2) &
         + (1.0_dp-p1-p2)*lognormal_pdf(x, u3, sigma3)
   end function dmln_am

   elemental real(dp) function dew(x, r, dividend_yield, te, s0, sigma, skew, kurt) result(value)
      real(dp), intent(in) :: x, r, dividend_yield, te, s0, sigma, skew, kurt
      real(dp) :: v, m, skew_lognorm, kurt_lognorm, cumul_lognorm
      real(dp) :: density, first_derivative, second_derivative
      real(dp) :: third_derivative, fourth_derivative

      if (x <= 0.0_dp .or. sigma <= 0.0_dp .or. te <= 0.0_dp .or. s0 <= 0.0_dp) then
         value = 0.0_dp
         return
      end if
      v = sqrt(exp(sigma*sigma*te)-1.0_dp)
      m = log(s0)+(r-dividend_yield-0.5_dp*sigma*sigma)*te
      skew_lognorm = 3.0_dp*v+v**3
      kurt_lognorm = 16.0_dp*v**2+15.0_dp*v**4+6.0_dp*v**6+v**8
      cumul_lognorm = (s0*exp((r-dividend_yield)*te)*v)**2
      density = lognormal_pdf(x, m, sigma*sqrt(te))
      first_derivative = -(1.0_dp+(log(x)-m)/(te*sigma*sigma))*density/x
      second_derivative = -(2.0_dp+(log(x)-m)/(te*sigma*sigma))*first_derivative/x &
         - density/(x*x*sigma*sigma)
      third_derivative = -(3.0_dp+(log(x)-m)/(te*sigma*sigma))*second_derivative/x &
         - 2.0_dp*first_derivative/(x*x*sigma*sigma) &
         + density/(x**3*sigma*sigma)
      fourth_derivative = -(4.0_dp+(log(x)-m)/(te*sigma*sigma))*third_derivative/x &
         - 3.0_dp*second_derivative/(x*x*sigma*sigma) &
         + 3.0_dp*first_derivative/(x**3*sigma*sigma) &
         - 2.0_dp*density/(x**4*sigma*sigma)
      value = density - (skew-skew_lognorm)*cumul_lognorm**1.5_dp*third_derivative/6.0_dp &
         + (kurt-kurt_lognorm)*cumul_lognorm**2*fourth_derivative/24.0_dp
   end function dew

   elemental real(dp) function dshimko(r, te, s0, strike, dividend_yield, a0, a1, a2) result(value)
      real(dp), intent(in) :: r, te, s0, strike, dividend_yield, a0, a1, a2
      real(dp) :: sigma, vol_time, d1, d2, d1x, d2x, sigma_derivative
      if (strike <= 0.0_dp .or. te <= 0.0_dp .or. s0 <= 0.0_dp) then
         value = 0.0_dp
         return
      end if
      sigma = a0+a1*strike+a2*strike*strike
      if (sigma <= 0.0_dp) then
         value = 0.0_dp
         return
      end if
      vol_time = sigma*sqrt(te)
      d1 = (log(s0/strike)+(r-dividend_yield+0.5_dp*sigma*sigma)*te)/vol_time
      d2 = d1-vol_time
      sigma_derivative = a1+2.0_dp*a2*strike
      d1x = -1.0_dp/(strike*vol_time)+(1.0_dp-d1/vol_time)*sigma_derivative
      d2x = d1x-sigma_derivative
      value = -normal_pdf(d2)*(d2x-sigma_derivative*(1.0_dp-d2*d2x)-2.0_dp*a2*strike)
   end function dshimko

end module rnd_densities
