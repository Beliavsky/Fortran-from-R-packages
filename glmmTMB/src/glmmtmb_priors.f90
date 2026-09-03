! SPDX-License-Identifier: AGPL-3.0-only
! Derived from glmmTMB 1.1.14 computational sources; see NOTICE.md.
module glmmtmb_priors
   use glmmtmb_kinds, only: dp
   use glmmtmb_codes, only: cauchy_prior, gamma_prior, lkj_prior, normal_prior, t_prior
   use glmmtmb_distributions, only: dcauchy, dlkj
   use tmb_distributions, only: dgamma, dnorm, dt_density
   use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
   implicit none
   private
   public :: correlation_prior_log_density, scalar_prior_log_density
contains
   pure real(dp) function scalar_prior_log_density(value, prior_code, params) result(ans)
      real(dp), intent(in) :: value !! Parameter value on the fitted scale to which the prior is applied.
      integer, intent(in) :: prior_code !! glmmTMB scalar prior code: normal, gamma, Student-t, or Cauchy.
      real(dp), intent(in) :: params(:) !! Prior hyperparameters in glmmTMB ordering for the selected distribution.
      real(dp) :: shape, scale, standardized
      ans = ieee_value(ans, ieee_quiet_nan)
      select case (prior_code)
      case (normal_prior)
         if (size(params) < 2) return
         ans = dnorm(value, params(1), params(2), .true.)
      case (gamma_prior)
         if (size(params) < 2 .or. params(2) <= 0.0_dp) return
         shape = params(2)
         scale = params(1) / params(2)
         ans = dgamma(exp(value), shape, scale, .true.)
      case (t_prior)
         if (size(params) < 3 .or. params(2) <= 0.0_dp .or. params(3) <= 0.0_dp) return
         standardized = (value - params(1)) / params(2)
         ans = dt_density(standardized, params(3), .true.) - log(params(2))
      case (cauchy_prior)
         if (size(params) < 2) return
         ans = dcauchy(value, params(1), params(2), .true.)
      end select
   end function scalar_prior_log_density

   pure real(dp) function correlation_prior_log_density(values, prior_code, params) result(ans)
      real(dp), intent(in) :: values(:) !! Vector of unconstrained correlation parameters receiving a multivariate prior.
      integer, intent(in) :: prior_code !! glmmTMB multivariate prior code, currently lkj_prior is supported.
      real(dp), intent(in) :: params(:) !! Prior hyperparameters, LKJ uses eta as params(1).
      ans = ieee_value(ans, ieee_quiet_nan)
      if (prior_code == lkj_prior .and. size(params) >= 1) then
         ans = dlkj(values, params(1), .true.)
      end if
   end function correlation_prior_log_density
end module glmmtmb_priors
