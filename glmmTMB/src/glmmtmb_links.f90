! SPDX-License-Identifier: AGPL-3.0-only
! Derived from glmmTMB 1.1.14 computational sources; see NOTICE.md.
module glmmtmb_links
   use glmmtmb_kinds, only: dp
   use glmmtmb_codes
   use glmmtmb_math, only: expm1_safe, invlogit, log1mexp, logaddexp, logsubexp
   use glmmtmb_distributions, only: dcompois2
   use tmb_distributions, only: pnorm, qnorm
   use, intrinsic :: ieee_arithmetic, only: ieee_positive_inf, ieee_quiet_nan, ieee_value
   implicit none
   private
   public :: calc_log_nzprob, inverse_linkfun, linkfun, log_inverse_linkfun
   public :: log1m_inverse_linkfun, logit_inverse_linkfun, trunc_family
contains
   pure elemental logical function trunc_family(family) result(ans)
      integer, intent(in) :: family !! glmmTMB family code to classify as zero-truncated or ordinary.
      ans = family == truncated_poisson_family .or. family == truncated_genpois_family .or. &
         family == truncated_compois_family .or. family == truncated_nbinom1_family .or. &
         family == truncated_nbinom2_family
   end function trunc_family

   pure elemental real(dp) function inverse_linkfun(eta, link) result(ans)
      real(dp), intent(in) :: eta !! Linear predictor on the selected link scale.
      integer, intent(in) :: link !! glmmTMB link code controlling the inverse transformation.
      select case (link)
      case (log_link)
         ans = exp(eta)
      case (identity_link)
         ans = eta
      case (logit_link)
         ans = invlogit(eta)
      case (probit_link)
         ans = pnorm(eta, 0.0_dp, 1.0_dp)
      case (cloglog_link)
         ans = -expm1_safe(-exp(eta))
      case (inverse_link)
         ans = 1.0_dp / eta
      case (sqrt_link)
         ans = eta * eta
      case (lambertw_link)
         ans = exp(eta) * exp(exp(eta))
      case default
         ans = ieee_value(ans, ieee_quiet_nan)
      end select
   end function inverse_linkfun

   pure elemental real(dp) function linkfun(mu, link) result(ans)
      real(dp), intent(in) :: mu !! Mean-scale value to transform to the requested link scale.
      integer, intent(in) :: link !! glmmTMB link code controlling the forward transformation.
      select case (link)
      case (log_link)
         ans = log(mu)
      case (identity_link)
         ans = mu
      case (logit_link)
         ans = log(mu) - log(1.0_dp - mu)
      case (probit_link)
         ans = qnorm(mu, 0.0_dp, 1.0_dp)
      case (cloglog_link)
         ans = log(-log(1.0_dp - mu))
      case (inverse_link)
         ans = 1.0_dp / mu
      case (sqrt_link)
         ans = sqrt(mu)
      case default
         ans = ieee_value(ans, ieee_quiet_nan)
      end select
   end function linkfun

   pure elemental real(dp) function logit_pnorm_local(x) result(ans)
      real(dp), intent(in) :: x !! Standard-normal quantile whose CDF is returned on the logit scale.
      real(dp) :: lower, upper
      lower = pnorm(x, 0.0_dp, 1.0_dp)
      upper = pnorm(-x, 0.0_dp, 1.0_dp)
      if (lower == 0.0_dp) then
         ans = -ieee_value(ans, ieee_positive_inf)
      else if (upper == 0.0_dp) then
         ans = ieee_value(ans, ieee_positive_inf)
      else
         ans = log(lower) - log(upper)
      end if
   end function logit_pnorm_local

   pure elemental real(dp) function logit_invcloglog_local(x) result(ans)
      real(dp), intent(in) :: x !! Complementary-log-log predictor to transform to log odds.
      real(dp) :: ex
      ex = exp(x)
      if (ex > 40.0_dp) then
         ans = ex
      else
         ans = log(expm1_safe(ex))
      end if
   end function logit_invcloglog_local

   pure elemental real(dp) function logit_inverse_linkfun(eta, link) result(ans)
      real(dp), intent(in) :: eta !! Linear predictor to transform to log odds of the inverse-link mean.
      integer, intent(in) :: link !! glmmTMB link code used for the mean transformation.
      real(dp) :: p
      select case (link)
      case (logit_link)
         ans = eta
      case (probit_link)
         ans = logit_pnorm_local(eta)
      case (cloglog_link)
         ans = logit_invcloglog_local(eta)
      case default
         p = inverse_linkfun(eta, link)
         if (p <= 0.0_dp) then
            ans = -ieee_value(ans, ieee_positive_inf)
         else if (p >= 1.0_dp) then
            ans = ieee_value(ans, ieee_positive_inf)
         else
            ans = log(p) - log(1.0_dp - p)
         end if
      end select
   end function logit_inverse_linkfun

   pure elemental real(dp) function log_inverse_linkfun(eta, link) result(ans)
      real(dp), intent(in) :: eta !! Linear predictor whose inverse-link value is needed on the log scale.
      integer, intent(in) :: link !! glmmTMB link code used for the mean transformation.
      select case (link)
      case (log_link)
         ans = eta
      case (logit_link)
         ans = -logaddexp(0.0_dp, -eta)
      case default
         ans = log(inverse_linkfun(eta, link))
      end select
   end function log_inverse_linkfun

   pure elemental real(dp) function log1m_inverse_linkfun(eta, link) result(ans)
      real(dp), intent(in) :: eta !! Linear predictor whose one-minus inverse-link value is needed on the log scale.
      integer, intent(in) :: link !! glmmTMB link code used for the mean transformation.
      real(dp) :: logp
      select case (link)
      case (log_link)
         ans = logsubexp(0.0_dp, eta)
      case (logit_link)
         ans = -logaddexp(0.0_dp, eta)
      case default
         logp = log(inverse_linkfun(eta, link))
         ans = logsubexp(0.0_dp, logp)
      end select
   end function log1m_inverse_linkfun

   pure real(dp) function calc_log_nzprob(mu, phi, eta, etadisp, family, link) result(ans)
      real(dp), intent(in) :: mu !! Positive conditional mean on the data scale.
      real(dp), intent(in) :: phi !! Positive dispersion value exp(etadisp).
      real(dp), intent(in) :: eta !! Conditional linear predictor corresponding to mu.
      real(dp), intent(in) :: etadisp !! Log-dispersion linear predictor corresponding to phi.
      integer, intent(in) :: family !! glmmTMB family code, normally one of the zero-truncated families.
      integer, intent(in) :: link !! glmmTMB conditional-link code.
      real(dp) :: s1, s2, logp0
      select case (family)
      case (truncated_nbinom1_family)
         s2 = logaddexp(0.0_dp, etadisp)
         logp0 = -mu * s2 / phi
         ans = log1mexp(logp0)
      case (truncated_nbinom2_family)
         s1 = log_inverse_linkfun(eta, link)
         s2 = logaddexp(0.0_dp, s1 - etadisp)
         logp0 = -phi * s2
         ans = log1mexp(logp0)
      case (truncated_poisson_family)
         ans = log1mexp(-mu)
      case (truncated_genpois_family)
         s1 = mu / sqrt(phi)
         ans = log1mexp(-s1)
      case (truncated_compois_family)
         logp0 = dcompois2(0.0_dp, mu, 1.0_dp / phi, .true.)
         ans = log1mexp(logp0)
      case default
         ans = 0.0_dp
      end select
   end function calc_log_nzprob
end module glmmtmb_links
