! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module vgam_links
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use vgam_kinds, only : dp, pi
   use vgam_special, only : logistic, logit, normal_cdf, normal_pdf, &
      normal_quantile, log1p_v, expm1_v
   implicit none
   private
   integer, parameter, public :: link_identity = 1
   integer, parameter, public :: link_log = 2
   integer, parameter, public :: link_logit = 3
   integer, parameter, public :: link_probit = 4
   integer, parameter, public :: link_cloglog = 5
   integer, parameter, public :: link_cauchit = 6
   integer, parameter, public :: link_reciprocal = 7
   integer, parameter, public :: link_sqrt = 8
   integer, parameter, public :: link_fisherz = 9
   integer, parameter, public :: link_clog = 10
   integer, parameter, public :: link_loglog = 11
   integer, parameter, public :: link_negidentity = 12
   integer, parameter, public :: link_neglog = 13
   integer, parameter, public :: link_negreciprocal = 14
   public :: link_value, link_inverse, link_derivative, inverse_link_derivative
   public :: loglink, logitlink, probitlink, clogloglink, cauchitlink
   public :: fisherzlink

contains

   elemental real(dp) function qnan() result(x)
      x = ieee_value(0.0_dp, ieee_quiet_nan)
   end function qnan

   elemental real(dp) function link_value(theta, link_id, power) result(eta)
      real(dp), intent(in) :: theta
      integer, intent(in) :: link_id
      real(dp), intent(in), optional :: power
      real(dp) :: p
      select case (link_id)
      case (link_identity)
         eta = theta
      case (link_negidentity)
         eta = -theta
      case (link_log)
         eta = merge(log(theta), qnan(), theta > 0.0_dp)
      case (link_neglog)
         eta = merge(-log(theta), qnan(), theta > 0.0_dp)
      case (link_logit)
         eta = logit(theta)
      case (link_probit)
         eta = normal_quantile(theta)
      case (link_cloglog)
         eta = merge(log(-log1p_v(-theta)), qnan(), &
                     theta > 0.0_dp .and. theta < 1.0_dp)
      case (link_clog)
         eta = merge(-log1p_v(-theta), qnan(), &
                     theta >= 0.0_dp .and. theta < 1.0_dp)
      case (link_cauchit)
         eta = tan(pi*(theta - 0.5_dp))
      case (link_reciprocal)
         eta = 1.0_dp/theta
      case (link_negreciprocal)
         eta = -1.0_dp/theta
      case (link_sqrt)
         eta = merge(sqrt(theta), qnan(), theta >= 0.0_dp)
      case (link_fisherz)
         eta = 0.5_dp*log((1.0_dp + theta)/(1.0_dp - theta))
      case (link_loglog)
         eta = merge(log(log(theta)), qnan(), theta > 1.0_dp)
      case default
         if (.not. present(power)) then
            eta = qnan()
         else
            p = power
            if (abs(p) < epsilon(1.0_dp)) then
               eta = log(theta)
            else
               eta = theta**p
            end if
         end if
      end select
   end function link_value

   elemental real(dp) function link_inverse(eta, link_id, power) result(theta)
      real(dp), intent(in) :: eta
      integer, intent(in) :: link_id
      real(dp), intent(in), optional :: power
      real(dp) :: p
      select case (link_id)
      case (link_identity)
         theta = eta
      case (link_negidentity)
         theta = -eta
      case (link_log)
         theta = exp(eta)
      case (link_neglog)
         theta = exp(-eta)
      case (link_logit)
         theta = logistic(eta)
      case (link_probit)
         theta = normal_cdf(eta)
      case (link_cloglog)
         theta = -expm1_v(-exp(eta))
      case (link_clog)
         theta = -expm1_v(-eta)
      case (link_cauchit)
         theta = 0.5_dp + atan(eta)/pi
      case (link_reciprocal)
         theta = 1.0_dp/eta
      case (link_negreciprocal)
         theta = -1.0_dp/eta
      case (link_sqrt)
         theta = eta*eta
      case (link_fisherz)
         theta = tanh(eta)
      case (link_loglog)
         theta = exp(exp(eta))
      case default
         if (.not. present(power)) then
            theta = qnan()
         else
            p = power
            if (abs(p) < epsilon(1.0_dp)) then
               theta = exp(eta)
            else
               theta = eta**(1.0_dp/p)
            end if
         end if
      end select
   end function link_inverse

   elemental real(dp) function link_derivative(theta, link_id, order, &
                                                power) result(d)
      real(dp), intent(in) :: theta
      integer, intent(in) :: link_id
      integer, intent(in), optional :: order
      real(dp), intent(in), optional :: power
      integer :: k
      real(dp) :: z, p, phi
      k = 1
      if (present(order)) k = order
      if (k < 1 .or. k > 3) then
         d = qnan()
         return
      end if
      select case (link_id)
      case (link_identity)
         d = merge(1.0_dp, 0.0_dp, k == 1)
      case (link_negidentity)
         d = merge(-1.0_dp, 0.0_dp, k == 1)
      case (link_log)
         select case (k)
         case (1); d = 1.0_dp/theta
         case (2); d = -1.0_dp/theta**2
         case (3); d = 2.0_dp/theta**3
         end select
      case (link_neglog)
         select case (k)
         case (1); d = -1.0_dp/theta
         case (2); d = 1.0_dp/theta**2
         case (3); d = -2.0_dp/theta**3
         end select
      case (link_logit)
         select case (k)
         case (1)
            d = 1.0_dp/(theta*(1.0_dp - theta))
         case (2)
            d = (2.0_dp*theta - 1.0_dp) / &
                (theta*(1.0_dp - theta))**2
         case (3)
            d = 2.0_dp*(1.0_dp - 3.0_dp*theta*(1.0_dp - theta)) / &
                (theta*(1.0_dp - theta))**3
         end select
      case (link_probit)
         z = normal_quantile(theta)
         phi = normal_pdf(z)
         select case (k)
         case (1); d = 1.0_dp/phi
         case (2); d = z/(phi*phi)
         case (3); d = (1.0_dp + 2.0_dp*z*z)/(phi**3)
         end select
      case (link_cloglog)
         z = log1p_v(-theta)
         select case (k)
         case (1)
            d = -1.0_dp/((1.0_dp - theta)*z)
         case (2)
            d = -(1.0_dp + z)/((1.0_dp - theta)*z)**2
         case (3)
            d = (1.0_dp/(1.0_dp - theta) - &
                2.0_dp*(1.0_dp + z)**2/((1.0_dp - theta)*z)) / &
                ((1.0_dp - theta)*z)**2
         end select
      case (link_clog)
         select case (k)
         case (1); d = 1.0_dp/(1.0_dp - theta)
         case (2); d = 1.0_dp/(1.0_dp - theta)**2
         case (3); d = 2.0_dp/(1.0_dp - theta)**3
         end select
      case (link_cauchit)
         z = pi*(theta - 0.5_dp)
         select case (k)
         case (1); d = pi/(cos(z)**2)
         case (2); d = 2.0_dp*pi*pi*tan(z)/(cos(z)**2)
         case (3)
            d = 2.0_dp*pi**3*(1.0_dp + 3.0_dp*tan(z)**2) / &
                (cos(z)**2)
         end select
      case (link_reciprocal)
         select case (k)
         case (1); d = -1.0_dp/theta**2
         case (2); d = 2.0_dp/theta**3
         case (3); d = -6.0_dp/theta**4
         end select
      case (link_negreciprocal)
         select case (k)
         case (1); d = 1.0_dp/theta**2
         case (2); d = -2.0_dp/theta**3
         case (3); d = 6.0_dp/theta**4
         end select
      case (link_sqrt)
         select case (k)
         case (1); d = 0.5_dp/sqrt(theta)
         case (2); d = -0.25_dp/theta**1.5_dp
         case (3); d = 0.375_dp/theta**2.5_dp
         end select
      case (link_fisherz)
         select case (k)
         case (1)
            d = 1.0_dp/(1.0_dp - theta*theta)
         case (2)
            d = 2.0_dp*theta/(1.0_dp - theta*theta)**2
         case (3)
            d = 2.0_dp*(1.0_dp + 3.0_dp*theta*theta) / &
                (1.0_dp - theta*theta)**3
         end select
      case (link_loglog)
         z = log(theta)
         select case (k)
         case (1); d = 1.0_dp/(theta*z)
         case (2); d = -(1.0_dp + z)/(theta*z)**2
         case (3)
            d = (2.0_dp*(1.0_dp + z)**2/(theta*z) - 1.0_dp/theta) / &
                (theta*z)**2
         end select
      case default
         if (.not. present(power)) then
            d = qnan()
            return
         end if
         p = power
         if (abs(p) < epsilon(1.0_dp)) then
            select case (k)
            case (1); d = 1.0_dp/theta
            case (2); d = -1.0_dp/theta**2
            case (3); d = 2.0_dp/theta**3
            end select
         else
            select case (k)
            case (1); d = p*theta**(p - 1.0_dp)
            case (2); d = p*(p - 1.0_dp)*theta**(p - 2.0_dp)
            case (3)
               d = p*(p - 1.0_dp)*(p - 2.0_dp)*theta**(p - 3.0_dp)
            end select
         end if
      end select
   end function link_derivative

   elemental real(dp) function inverse_link_derivative(eta, link_id) result(d)
      real(dp), intent(in) :: eta
      integer, intent(in) :: link_id
      real(dp) :: mu
      mu = link_inverse(eta, link_id)
      d = 1.0_dp/link_derivative(mu, link_id)
   end function inverse_link_derivative

   elemental real(dp) function loglink(theta, inverse) result(ans)
      real(dp), intent(in) :: theta
      logical, intent(in), optional :: inverse
      if (present(inverse) .and. inverse) then
         ans = exp(theta)
      else
         ans = log(theta)
      end if
   end function loglink

   elemental real(dp) function logitlink(theta, inverse) result(ans)
      real(dp), intent(in) :: theta
      logical, intent(in), optional :: inverse
      if (present(inverse) .and. inverse) then
         ans = logistic(theta)
      else
         ans = logit(theta)
      end if
   end function logitlink

   elemental real(dp) function probitlink(theta, inverse) result(ans)
      real(dp), intent(in) :: theta
      logical, intent(in), optional :: inverse
      if (present(inverse) .and. inverse) then
         ans = normal_cdf(theta)
      else
         ans = normal_quantile(theta)
      end if
   end function probitlink

   elemental real(dp) function clogloglink(theta, inverse) result(ans)
      real(dp), intent(in) :: theta
      logical, intent(in), optional :: inverse
      if (present(inverse) .and. inverse) then
         ans = link_inverse(theta, link_cloglog)
      else
         ans = link_value(theta, link_cloglog)
      end if
   end function clogloglink

   elemental real(dp) function cauchitlink(theta, inverse) result(ans)
      real(dp), intent(in) :: theta
      logical, intent(in), optional :: inverse
      if (present(inverse) .and. inverse) then
         ans = link_inverse(theta, link_cauchit)
      else
         ans = link_value(theta, link_cauchit)
      end if
   end function cauchitlink

   elemental real(dp) function fisherzlink(theta, inverse) result(ans)
      real(dp), intent(in) :: theta
      logical, intent(in), optional :: inverse
      if (present(inverse) .and. inverse) then
         ans = tanh(theta)
      else
         ans = link_value(theta, link_fisherz)
      end if
   end function fisherzlink

end module vgam_links
