! Modern Fortran translation of src/famstr.cc from geepack 1.3-13.
! Upstream license: GPL (>= 3). See LICENSE, NOTICE.md, and PROVENANCE.md.
! Exact upstream copyright-holder attribution is preserved in NOTICE.md and upstream/DESCRIPTION.
module geepack_links
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use r_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: LINK_IDENTITY = 1
   integer, parameter, public :: LINK_LOGIT = 2
   integer, parameter, public :: LINK_PROBIT = 3
   integer, parameter, public :: LINK_CLOGLOG = 4
   integer, parameter, public :: LINK_LOG = 5
   integer, parameter, public :: LINK_RECIPROCAL = 6
   integer, parameter, public :: LINK_FISHERZ = 7
   integer, parameter, public :: LINK_LWYBC2 = 8
   integer, parameter, public :: LINK_LWYLOG = 9

   integer, parameter, public :: VAR_GAUSSIAN = 1
   integer, parameter, public :: VAR_BINOMIAL = 2
   integer, parameter, public :: VAR_POISSON = 3
   integer, parameter, public :: VAR_GAMMA = 4

   public :: link_function, link_inverse, link_derivative
   public :: variance_function, variance_derivative, valid_mean

contains

   pure elemental real(dp) function link_function(mu, link_code) result(eta)
      real(dp), intent(in) :: mu !! Mean value at which the link is evaluated.
      integer, intent(in) :: link_code !! Integer link identifier LINK_*.
      real(dp), parameter :: tiny_prob = epsilon(1.0_dp)
      real(dp) :: p

      p = min(1.0_dp - tiny_prob, max(tiny_prob, mu))
      select case (link_code)
      case (LINK_IDENTITY)
         eta = mu
      case (LINK_LOGIT)
         eta = log(p / (1.0_dp - p))
      case (LINK_PROBIT)
         eta = normal_quantile(p)
      case (LINK_CLOGLOG)
         eta = log(-log(1.0_dp - p))
      case (LINK_LOG)
         eta = log(max(tiny_prob, mu))
      case (LINK_RECIPROCAL)
         eta = 1.0_dp / mu
      case (LINK_FISHERZ)
         eta = log(2.0_dp / (1.0_dp - mu) - 1.0_dp)
      case (LINK_LWYBC2)
         eta = log(sqrt(mu + 1.0_dp) - 1.0_dp)
      case (LINK_LWYLOG)
         eta = log(exp(mu) - 1.0_dp)
      case default
         eta = mu
      end select
   end function link_function

   pure elemental real(dp) function link_inverse(eta, link_code) result(mu)
      real(dp), intent(in) :: eta !! Linear predictor value.
      integer, intent(in) :: link_code !! Integer link identifier LINK_*.
      real(dp), parameter :: tiny_prob = epsilon(1.0_dp)
      real(dp) :: e
      real(dp) :: threshold

      threshold = -log(tiny_prob)
      select case (link_code)
      case (LINK_IDENTITY)
         mu = eta
      case (LINK_LOGIT)
         e = min(threshold, max(-threshold, eta))
         if (e >= 0.0_dp) then
            mu = 1.0_dp / (1.0_dp + exp(-e))
         else
            mu = exp(e) / (1.0_dp + exp(e))
         end if
      case (LINK_PROBIT)
         mu = 0.5_dp * erfc(-eta / sqrt(2.0_dp))
         mu = min(1.0_dp - tiny_prob, max(tiny_prob, mu))
      case (LINK_CLOGLOG)
         e = min(700.0_dp, eta)
         mu = 1.0_dp - exp(-exp(e))
         mu = min(1.0_dp - tiny_prob, max(tiny_prob, mu))
      case (LINK_LOG)
         mu = max(tiny_prob, exp(min(700.0_dp, eta)))
      case (LINK_RECIPROCAL)
         mu = 1.0_dp / eta
      case (LINK_FISHERZ)
         e = min(threshold, max(-threshold, eta))
         mu = tanh(0.5_dp * e)
      case (LINK_LWYBC2)
         e = max(tiny_prob, exp(min(700.0_dp, eta)))
         mu = (1.0_dp + e) ** 2 - 1.0_dp
      case (LINK_LWYLOG)
         if (eta > 36.0_dp) then
            mu = eta
         else
            mu = log(1.0_dp + exp(eta))
         end if
      case default
         mu = eta
      end select
   end function link_inverse

   pure elemental real(dp) function link_derivative(eta, link_code) result(dmu_deta)
      real(dp), intent(in) :: eta !! Linear predictor value.
      integer, intent(in) :: link_code !! Integer link identifier LINK_*.
      real(dp), parameter :: tiny_prob = epsilon(1.0_dp)
      real(dp), parameter :: inv_sqrt_2pi = 0.3989422804014326779399460599343819_dp
      real(dp) :: e
      real(dp) :: ex
      real(dp) :: threshold

      threshold = -log(tiny_prob)
      select case (link_code)
      case (LINK_IDENTITY)
         dmu_deta = 1.0_dp
      case (LINK_LOGIT)
         if (abs(eta) >= threshold) then
            dmu_deta = tiny_prob
         else
            e = exp(eta)
            dmu_deta = max(tiny_prob, e / (1.0_dp + e) ** 2)
         end if
      case (LINK_PROBIT)
         dmu_deta = max(tiny_prob, inv_sqrt_2pi * exp(-0.5_dp * eta * eta))
      case (LINK_CLOGLOG)
         e = min(700.0_dp, eta)
         dmu_deta = max(tiny_prob, exp(e) * exp(-exp(e)))
      case (LINK_LOG)
         dmu_deta = max(tiny_prob, exp(min(700.0_dp, eta)))
      case (LINK_RECIPROCAL)
         dmu_deta = -1.0_dp / (eta * eta)
      case (LINK_FISHERZ)
         if (abs(eta) >= threshold) then
            dmu_deta = tiny_prob
         else
            e = exp(eta)
            dmu_deta = max(tiny_prob, 2.0_dp * e / (1.0_dp + e) ** 2)
         end if
      case (LINK_LWYBC2)
         ex = exp(min(700.0_dp, eta))
         dmu_deta = max(tiny_prob, 2.0_dp * (1.0_dp + ex) * ex)
      case (LINK_LWYLOG)
         if (eta >= 0.0_dp) then
            dmu_deta = 1.0_dp / (1.0_dp + exp(-eta))
         else
            ex = exp(eta)
            dmu_deta = ex / (1.0_dp + ex)
         end if
      case default
         dmu_deta = 1.0_dp
      end select
   end function link_derivative

   pure elemental real(dp) function variance_function(mu, variance_code) result(v)
      real(dp), intent(in) :: mu !! Mean value used by the variance function.
      integer, intent(in) :: variance_code !! Integer variance identifier VAR_*.

      select case (variance_code)
      case (VAR_GAUSSIAN)
         v = 1.0_dp
      case (VAR_BINOMIAL)
         v = mu * (1.0_dp - mu)
      case (VAR_POISSON)
         v = mu
      case (VAR_GAMMA)
         v = mu * mu
      case default
         v = 1.0_dp
      end select
   end function variance_function

   pure elemental real(dp) function variance_derivative(mu, variance_code) result(dv_dmu)
      real(dp), intent(in) :: mu !! Mean value used by the variance function.
      integer, intent(in) :: variance_code !! Integer variance identifier VAR_*.

      select case (variance_code)
      case (VAR_GAUSSIAN)
         dv_dmu = 0.0_dp
      case (VAR_BINOMIAL)
         dv_dmu = 1.0_dp - 2.0_dp * mu
      case (VAR_POISSON)
         dv_dmu = 1.0_dp
      case (VAR_GAMMA)
         dv_dmu = 2.0_dp * mu
      case default
         dv_dmu = 0.0_dp
      end select
   end function variance_derivative

   pure elemental logical function valid_mean(mu, variance_code) result(ok)
      real(dp), intent(in) :: mu !! Candidate mean value.
      integer, intent(in) :: variance_code !! Integer variance identifier VAR_*.

      if (.not. ieee_is_finite(mu)) then
         ok = .false.
         return
      end if
      select case (variance_code)
      case (VAR_BINOMIAL)
         ok = mu > 0.0_dp .and. mu < 1.0_dp
      case (VAR_POISSON, VAR_GAMMA)
         ok = mu > 0.0_dp
      case default
         ok = .true.
      end select
   end function valid_mean

   pure elemental real(dp) function normal_quantile(p) result(x)
      real(dp), intent(in) :: p !! Probability strictly between zero and one.
      real(dp), parameter :: a1 = -3.969683028665376e1_dp
      real(dp), parameter :: a2 = 2.209460984245205e2_dp
      real(dp), parameter :: a3 = -2.759285104469687e2_dp
      real(dp), parameter :: a4 = 1.383577518672690e2_dp
      real(dp), parameter :: a5 = -3.066479806614716e1_dp
      real(dp), parameter :: a6 = 2.506628277459239_dp
      real(dp), parameter :: b1 = -5.447609879822406e1_dp
      real(dp), parameter :: b2 = 1.615858368580409e2_dp
      real(dp), parameter :: b3 = -1.556989798598866e2_dp
      real(dp), parameter :: b4 = 6.680131188771972e1_dp
      real(dp), parameter :: b5 = -1.328068155288572e1_dp
      real(dp), parameter :: c1 = -7.784894002430293e-3_dp
      real(dp), parameter :: c2 = -3.223964580411365e-1_dp
      real(dp), parameter :: c3 = -2.400758277161838_dp
      real(dp), parameter :: c4 = -2.549732539343734_dp
      real(dp), parameter :: c5 = 4.374664141464968_dp
      real(dp), parameter :: c6 = 2.938163982698783_dp
      real(dp), parameter :: d1 = 7.784695709041462e-3_dp
      real(dp), parameter :: d2 = 3.224671290700398e-1_dp
      real(dp), parameter :: d3 = 2.445134137142996_dp
      real(dp), parameter :: d4 = 3.754408661907416_dp
      real(dp), parameter :: plow = 0.02425_dp
      real(dp), parameter :: phigh = 1.0_dp - plow
      real(dp) :: q
      real(dp) :: r

      if (p <= 0.0_dp) then
         x = -huge(1.0_dp)
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
      else if (p < plow) then
         q = sqrt(-2.0_dp * log(p))
         x = (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / &
            ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0_dp)
      else if (p <= phigh) then
         q = p - 0.5_dp
         r = q * q
         x = (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q / &
            (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1.0_dp)
      else
         q = sqrt(-2.0_dp * log(1.0_dp - p))
         x = -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / &
            ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0_dp)
      end if
   end function normal_quantile

end module geepack_links
