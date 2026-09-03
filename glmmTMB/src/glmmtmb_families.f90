! SPDX-License-Identifier: AGPL-3.0-only
! Derived from glmmTMB 1.1.14 computational sources; see NOTICE.md.
module glmmtmb_families
   use glmmtmb_kinds, only: dp
   use glmmtmb_codes
   use glmmtmb_math, only: invlogit, lambert_w0, logaddexp, logsubexp
   use glmmtmb_distributions
   use glmmtmb_links, only: calc_log_nzprob, inverse_linkfun, log1m_inverse_linkfun, &
      log_inverse_linkfun, logit_inverse_linkfun
   use tmb_distributions, only: dbeta, dbinom_robust, dgamma, dnorm, dt_density
   use, intrinsic :: ieee_arithmetic, only: ieee_positive_inf, ieee_quiet_nan, ieee_value
   implicit none
   private
   public :: family_variance, observation_loglik
contains
   pure real(dp) function observation_loglik(y, size_n, eta, etadisp, etazi, family, link, psi, zi_enabled) result(loglik)
      real(dp), intent(in) :: y !! Observed response value, counts are represented as real-valued integers.
      real(dp), intent(in) :: size_n !! Binomial or beta-binomial trial count, ignored by other families.
      real(dp), intent(in) :: eta !! Conditional linear predictor before applying the selected link.
      real(dp), intent(in) :: etadisp !! Dispersion linear predictor, with phi=exp(etadisp).
      real(dp), intent(in) :: etazi !! Zero-inflation logit predictor, ignored when zi_enabled is false.
      integer, intent(in) :: family !! glmmTMB response-family code.
      integer, intent(in) :: link !! glmmTMB conditional-link code.
      real(dp), intent(in) :: psi(:) !! Extra family parameters on their fitted scales, in upstream ordering.
      logical, intent(in) :: zi_enabled !! Apply the glmmTMB zero-inflation mixture when true.
      real(dp) :: mu, phi, s1, s2, s3, log_nzprob, log_pz, log_1mpz, btheta
      loglik = ieee_value(loglik, ieee_quiet_nan)
      mu = inverse_linkfun(eta, link)
      phi = exp(etadisp)
      log_nzprob = calc_log_nzprob(mu, phi, eta, etadisp, family, link)

      select case (family)
      case (gaussian_family)
         loglik = dnorm(y, mu, phi, .true.)
      case (skewnormal_family)
         if (size(psi) < 1) return
         loglik = dskewnorm(y, mu, phi, psi(1), .true.)
      case (poisson_family)
         loglik = dpois_glmmtmb(y, mu, .true.)
      case (binomial_family)
         s1 = logit_inverse_linkfun(eta, link)
         loglik = dbinom_robust(y, size_n, s1, .true.)
      case (gamma_family)
         if (y == 0.0_dp .and. zi_enabled) then
            loglik = -ieee_value(loglik, ieee_positive_inf)
         else
            loglik = dgamma(y, phi, mu / phi, .true.)
         end if
      case (beta_family)
         if (y == 0.0_dp .and. zi_enabled) then
            loglik = -ieee_value(loglik, ieee_positive_inf)
         else
            loglik = dbeta(y, mu * phi, (1.0_dp - mu) * phi, .true.)
         end if
      case (ordbeta_family)
         if (size(psi) < 2) return
         if (y == 0.0_dp) then
            loglik = log1m_inverse_linkfun(eta - psi(1), logit_link)
         else if (y == 1.0_dp) then
            loglik = log_inverse_linkfun(eta - psi(2), logit_link)
         else
            s1 = mu * phi
            s2 = (1.0_dp - mu) * phi
            s3 = logsubexp(log_inverse_linkfun(eta - psi(1), logit_link), &
               log_inverse_linkfun(eta - psi(2), logit_link))
            loglik = s3 + dbeta(y, s1, s2, .true.)
         end if
      case (betabinomial_family)
         s3 = logit_inverse_linkfun(eta, link)
         s1 = log_inverse_linkfun(s3, logit_link) + log(phi)
         s2 = log_inverse_linkfun(-s3, logit_link) + log(phi)
         loglik = dbetabinom_robust(y, s1, s2, size_n, .true.)
      case (nbinom1_family, truncated_nbinom1_family)
         s1 = log_inverse_linkfun(eta, link)
         s2 = s1 + etadisp
         loglik = dnbinom_robust(y, s1, s2, .true.)
         if (family == truncated_nbinom1_family) then
            loglik = loglik - log_nzprob
            if (y < 0.001_dp) loglik = -ieee_value(loglik, ieee_positive_inf)
         end if
      case (nbinom2_family, truncated_nbinom2_family)
         s1 = log_inverse_linkfun(eta, link)
         s2 = 2.0_dp * s1 - etadisp
         loglik = dnbinom_robust(y, s1, s2, .true.)
         if (family == truncated_nbinom2_family) then
            loglik = loglik - log_nzprob
            if (y < 0.001_dp) loglik = -ieee_value(loglik, ieee_positive_inf)
         end if
      case (nbinom12_family)
         if (size(psi) < 1) return
         s1 = log_inverse_linkfun(eta, link)
         s2 = s1 + logaddexp(etadisp, s1 - psi(1))
         loglik = dnbinom_robust(y, s1, s2, .true.)
      case (truncated_poisson_family)
         loglik = dpois_glmmtmb(y, mu, .true.) - log_nzprob
         if (y < 0.001_dp) loglik = -ieee_value(loglik, ieee_positive_inf)
      case (genpois_family, truncated_genpois_family)
         s1 = mu / sqrt(phi)
         s2 = 1.0_dp - 1.0_dp / sqrt(phi)
         loglik = dgenpois(y, s1, s2, .true.)
         if (family == truncated_genpois_family) then
            loglik = loglik - log_nzprob
            if (y < 0.001_dp) loglik = -ieee_value(loglik, ieee_positive_inf)
         end if
      case (compois_family, truncated_compois_family)
         loglik = dcompois2(y, mu, 1.0_dp / phi, .true.)
         if (family == truncated_compois_family) then
            loglik = loglik - log_nzprob
            if (y < 0.001_dp) loglik = -ieee_value(loglik, ieee_positive_inf)
         end if
      case (tweedie_family)
         if (size(psi) < 1) return
         s3 = invlogit(psi(1)) + 1.0_dp
         loglik = dtweedie_compound(y, mu, phi, s3, .true.)
      case (lognormal_family)
         if (y == 0.0_dp .and. zi_enabled) then
            loglik = -ieee_value(loglik, ieee_positive_inf)
         else if (y <= 0.0_dp) then
            loglik = -ieee_value(loglik, ieee_positive_inf)
         else
            s1 = logaddexp(2.0_dp * (log(phi) - log(mu)), 0.0_dp)
            s2 = log(mu) - 0.5_dp * s1
            s3 = sqrt(s1)
            loglik = dnorm(log(y), s2, s3, .true.) - log(y)
         end if
      case (t_family)
         if (size(psi) < 1) return
         s1 = (y - mu) / phi
         s2 = exp(psi(1))
         loglik = dt_density(s1, s2, .true.) - etadisp
      case (bell_family)
         btheta = lambert_w0(mu)
         loglik = dbell(y, btheta, .true.)
      case default
         return
      end select

      if (zi_enabled) then
         log_pz = -logaddexp(0.0_dp, -etazi)
         log_1mpz = -logaddexp(0.0_dp, etazi)
         if (y == 0.0_dp) then
            loglik = logaddexp(log_pz, log_1mpz + loglik)
         else
            loglik = log_1mpz + loglik
         end if
      end if
   end function observation_loglik

   pure real(dp) function family_variance(mu, phi, family, psi) result(ans)
      real(dp), intent(in) :: mu !! Conditional mean on the response scale.
      real(dp), intent(in) :: phi !! Positive dispersion parameter on the response-family scale.
      integer, intent(in) :: family !! glmmTMB response-family code.
      real(dp), intent(in) :: psi(:) !! Extra family parameters in the same scale/order used by observation_loglik.
      real(dp) :: theta, p0, mu_star, power, shape, delta
      select case (family)
      case (gaussian_family, lognormal_family)
         ans = phi * phi
      case (poisson_family)
         ans = mu
      case (nbinom1_family)
         ans = mu * (1.0_dp + phi)
      case (nbinom2_family)
         ans = mu * (1.0_dp + mu / phi)
      case (nbinom12_family)
         if (size(psi) < 1) then
            ans = ieee_value(ans, ieee_quiet_nan)
         else
            ans = mu * (1.0_dp + phi + mu / exp(psi(1)))
         end if
      case (compois_family)
         ans = compois_variance(mu, 1.0_dp / phi)
      case (genpois_family)
         ans = mu * phi
      case (truncated_poisson_family)
         if (mu <= 0.0_dp) then
            ans = 0.0_dp
         else
            p0 = exp(-mu)
            ans = (mu + mu * mu) / (1.0_dp - p0) - mu * mu / (1.0_dp - p0)**2
         end if
      case (truncated_nbinom2_family)
         theta = phi
         if (mu <= 0.0_dp .or. theta <= 0.0_dp) then
            ans = ieee_value(ans, ieee_quiet_nan)
         else
            p0 = exp(theta * (log(theta) - log(theta + mu)))
            mu_star = mu / (1.0_dp - p0)
            ans = mu_star + mu_star * mu * (1.0_dp + 1.0_dp / theta) - mu_star * mu_star
         end if
      case (beta_family)
         ans = mu * (1.0_dp - mu) / (1.0_dp + phi)
      case (betabinomial_family, ordbeta_family)
         ans = mu * (1.0_dp - mu)
      case (tweedie_family)
         if (size(psi) < 1) then
            ans = ieee_value(ans, ieee_quiet_nan)
         else
            power = invlogit(psi(1)) + 1.0_dp
            ans = phi * mu**power
         end if
      case (skewnormal_family)
         if (size(psi) < 1) then
            ans = ieee_value(ans, ieee_quiet_nan)
         else
            shape = psi(1)
            if (abs(shape) >= 1.0_dp) then
               ans = ieee_value(ans, ieee_quiet_nan)
            else
               delta = shape / sqrt(1.0_dp - shape * shape)
               ans = phi * phi * (1.0_dp - 2.0_dp * delta * delta / acos(-1.0_dp))
            end if
         end if
      case (t_family)
         ans = phi
      case (bell_family)
         ans = mu * (1.0_dp + lambert_w0(mu))
      case default
         ans = ieee_value(ans, ieee_quiet_nan)
      end select
   end function family_variance
end module glmmtmb_families
