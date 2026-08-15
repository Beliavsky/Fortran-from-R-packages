! Box-Cox distribution families translated from gamlss.dist 6.1-1.
! Original package: Copyright (C) gamlss.dist authors, GPL-2 | GPL-3.
! Translation: 2026. SPDX-License-Identifier: GPL-3.0-only
module gamlss_boxcox
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf
   use gamlss_kinds, only : dp
   use gamlss_base, only : pnorm_v, qnorm_v, pgamma_v, qgamma_v
   use gamlss_student_t, only : student_t_pdf, student_t_cdf, student_t_quantile
   implicit none
   private
   public :: dBCCG, pBCCG, qBCCG, rBCCG
   public :: dBCT, pBCT, qBCT, rBCT
   public :: dBCPE, pBCPE, qBCPE, rBCPE

contains

   elemental real(dp) function nanv() result(x)
      x = ieee_value(0.0_dp, ieee_quiet_nan)
   end function nanv

   elemental real(dp) function infv() result(x)
      x = ieee_value(0.0_dp, ieee_positive_inf)
   end function infv

   elemental logical function flag_value(flag, default_value) result(value)
      logical, intent(in), optional :: flag
      logical, intent(in) :: default_value
      value = default_value
      if (present(flag)) value = flag
   end function flag_value

   elemental real(dp) function finish_probability(probability, lower_tail, log_p) result(value)
      real(dp), intent(in) :: probability
      logical, intent(in), optional :: lower_tail, log_p
      logical :: lower_value, log_value
      value = max(0.0_dp, min(1.0_dp, probability))
      lower_value = flag_value(lower_tail, .true.)
      log_value = flag_value(log_p, .false.)
      if (.not. lower_value) value = 1.0_dp - value
      if (log_value) then
         if (value <= 0.0_dp) then
            value = -infv()
         else
            value = log(value)
         end if
      end if
   end function finish_probability

   elemental real(dp) function input_probability(probability, lower_tail, log_p) result(value)
      real(dp), intent(in) :: probability
      logical, intent(in), optional :: lower_tail, log_p
      value = probability
      if (flag_value(log_p, .false.)) value = exp(value)
      if (.not. flag_value(lower_tail, .true.)) value = 1.0_dp - value
   end function input_probability

   elemental real(dp) function boxcox_z(x, mu, sigma, nu) result(z)
      real(dp), intent(in) :: x, mu, sigma, nu
      if (x <= 0.0_dp .or. mu <= 0.0_dp .or. sigma <= 0.0_dp) then
         z = nanv()
      else if (abs(nu) <= 1.0e-12_dp) then
         z = log(x / mu) / sigma
      else
         z = ((x / mu)**nu - 1.0_dp) / (nu * sigma)
      end if
   end function boxcox_z

   elemental real(dp) function inv_boxcox_z(z, mu, sigma, nu) result(x)
      real(dp), intent(in) :: z, mu, sigma, nu
      real(dp) :: base
      if (mu <= 0.0_dp .or. sigma <= 0.0_dp) then
         x = nanv()
      else if (abs(nu) <= 1.0e-12_dp) then
         x = mu * exp(sigma * z)
      else
         base = 1.0_dp + nu * sigma * z
         if (base <= 0.0_dp) then
            x = nanv()
         else
            x = mu * base**(1.0_dp / nu)
         end if
      end if
   end function inv_boxcox_z

   elemental real(dp) function bccg_norm(sigma, nu) result(value)
      real(dp), intent(in) :: sigma, nu
      if (abs(nu) <= 1.0e-12_dp) then
         value = 1.0_dp
      else
         value = pnorm_v(1.0_dp / (sigma * abs(nu)), 0.0_dp, 1.0_dp)
      end if
   end function bccg_norm

   elemental real(dp) function dBCCG(x, mu, sigma, nu, log_density) result(value)
      real(dp), intent(in) :: x, mu, sigma, nu
      logical, intent(in), optional :: log_density
      real(dp) :: z, log_density_value, norm_value
      logical :: log_value
      log_value = flag_value(log_density, .false.)
      if (mu <= 0.0_dp .or. sigma <= 0.0_dp) then
         value = nanv()
         return
      end if
      if (x <= 0.0_dp) then
         value = merge(-infv(), 0.0_dp, log_value)
         return
      end if
      z = boxcox_z(x, mu, sigma, nu)
      norm_value = bccg_norm(sigma, nu)
      log_density_value = nu * log(x / mu) - log(sigma) - 0.5_dp * z * z &
         - log(x) - 0.5_dp * log(2.0_dp * acos(-1.0_dp)) - log(norm_value)
      value = merge(log_density_value, exp(log_density_value), log_value)
   end function dBCCG

   elemental real(dp) function pBCCG(q, mu, sigma, nu, lower_tail, log_p) result(value)
      real(dp), intent(in) :: q, mu, sigma, nu
      logical, intent(in), optional :: lower_tail, log_p
      real(dp) :: z, lower_cut, norm_value, probability
      if (mu <= 0.0_dp .or. sigma <= 0.0_dp) then
         value = nanv()
         return
      end if
      if (q <= 0.0_dp) then
         value = finish_probability(0.0_dp, lower_tail, log_p)
         return
      end if
      z = boxcox_z(q, mu, sigma, nu)
      norm_value = bccg_norm(sigma, nu)
      lower_cut = 0.0_dp
      if (nu > 1.0e-12_dp) then
         lower_cut = pnorm_v(-1.0_dp / (sigma * abs(nu)), 0.0_dp, 1.0_dp)
      end if
      probability = (pnorm_v(z, 0.0_dp, 1.0_dp) - lower_cut) / norm_value
      value = finish_probability(probability, lower_tail, log_p)
   end function pBCCG

   real(dp) function qBCCG(p, mu, sigma, nu, lower_tail, log_p) result(value)
      real(dp), intent(in) :: p, mu, sigma, nu
      logical, intent(in), optional :: lower_tail, log_p
      real(dp) :: probability, z, norm_value, adjusted_probability
      if (mu <= 0.0_dp .or. sigma <= 0.0_dp) then
         value = nanv()
         return
      end if
      probability = input_probability(p, lower_tail, log_p)
      if (probability < 0.0_dp .or. probability > 1.0_dp) then
         value = nanv()
         return
      end if
      if (probability == 0.0_dp) then
         value = 0.0_dp
         return
      else if (probability == 1.0_dp) then
         value = infv()
         return
      end if
      if (abs(nu) <= 1.0e-12_dp) then
         z = qnorm_v(probability, 0.0_dp, 1.0_dp)
      else
         norm_value = bccg_norm(sigma, nu)
         if (nu < 0.0_dp) then
            adjusted_probability = probability * norm_value
         else
            adjusted_probability = 1.0_dp - (1.0_dp - probability) * norm_value
         end if
         z = qnorm_v(adjusted_probability, 0.0_dp, 1.0_dp)
      end if
      value = inv_boxcox_z(z, mu, sigma, nu)
   end function qBCCG

   real(dp) function rBCCG(mu, sigma, nu) result(value)
      real(dp), intent(in) :: mu, sigma, nu
      real(dp) :: u
      call random_number(u)
      value = qBCCG(u, mu, sigma, nu)
   end function rBCCG

   elemental real(dp) function bct_norm(sigma, nu, tau) result(value)
      real(dp), intent(in) :: sigma, nu, tau
      if (abs(nu) <= 1.0e-12_dp) then
         value = 1.0_dp
      else
         value = student_t_cdf(1.0_dp / (sigma * abs(nu)), tau)
      end if
   end function bct_norm

   elemental real(dp) function dBCT(x, mu, sigma, nu, tau, log_density) result(value)
      real(dp), intent(in) :: x, mu, sigma, nu, tau
      logical, intent(in), optional :: log_density
      real(dp) :: z, log_density_value, norm_value
      logical :: log_value
      log_value = flag_value(log_density, .false.)
      if (mu <= 0.0_dp .or. sigma <= 0.0_dp .or. tau <= 0.0_dp) then
         value = nanv()
         return
      end if
      if (x <= 0.0_dp) then
         value = merge(-infv(), 0.0_dp, log_value)
         return
      end if
      if (tau > 1.0e6_dp) then
         value = dBCCG(x, mu, sigma, nu, log_density)
         return
      end if
      z = boxcox_z(x, mu, sigma, nu)
      norm_value = bct_norm(sigma, nu, tau)
      log_density_value = (nu - 1.0_dp) * log(x) - nu * log(mu) - log(sigma) &
         + log(student_t_pdf(z, tau)) - log(norm_value)
      value = merge(log_density_value, exp(log_density_value), log_value)
   end function dBCT

   elemental real(dp) function pBCT(q, mu, sigma, nu, tau, lower_tail, log_p) result(value)
      real(dp), intent(in) :: q, mu, sigma, nu, tau
      logical, intent(in), optional :: lower_tail, log_p
      real(dp) :: z, lower_cut, norm_value, probability
      if (mu <= 0.0_dp .or. sigma <= 0.0_dp .or. tau <= 0.0_dp) then
         value = nanv()
         return
      end if
      if (q <= 0.0_dp) then
         value = finish_probability(0.0_dp, lower_tail, log_p)
         return
      end if
      if (tau > 1.0e6_dp) then
         value = pBCCG(q, mu, sigma, nu, lower_tail, log_p)
         return
      end if
      z = boxcox_z(q, mu, sigma, nu)
      norm_value = bct_norm(sigma, nu, tau)
      lower_cut = 0.0_dp
      if (nu > 1.0e-12_dp) then
         lower_cut = student_t_cdf(-1.0_dp / (sigma * abs(nu)), tau)
      end if
      probability = (student_t_cdf(z, tau) - lower_cut) / norm_value
      value = finish_probability(probability, lower_tail, log_p)
   end function pBCT

   real(dp) function qBCT(p, mu, sigma, nu, tau, lower_tail, log_p) result(value)
      real(dp), intent(in) :: p, mu, sigma, nu, tau
      logical, intent(in), optional :: lower_tail, log_p
      real(dp) :: probability, z, norm_value, adjusted_probability
      if (mu <= 0.0_dp .or. sigma <= 0.0_dp .or. tau <= 0.0_dp) then
         value = nanv()
         return
      end if
      probability = input_probability(p, lower_tail, log_p)
      if (probability < 0.0_dp .or. probability > 1.0_dp) then
         value = nanv()
         return
      end if
      if (probability == 0.0_dp) then
         value = 0.0_dp
         return
      else if (probability == 1.0_dp) then
         value = infv()
         return
      end if
      if (tau > 1.0e6_dp) then
         value = qBCCG(probability, mu, sigma, nu)
         return
      end if
      if (abs(nu) <= 1.0e-12_dp) then
         z = student_t_quantile(probability, tau)
      else
         norm_value = bct_norm(sigma, nu, tau)
         if (nu < 0.0_dp) then
            adjusted_probability = probability * norm_value
         else
            adjusted_probability = 1.0_dp - (1.0_dp - probability) * norm_value
         end if
         z = student_t_quantile(adjusted_probability, tau)
      end if
      value = inv_boxcox_z(z, mu, sigma, nu)
   end function qBCT

   real(dp) function rBCT(mu, sigma, nu, tau) result(value)
      real(dp), intent(in) :: mu, sigma, nu, tau
      real(dp) :: u
      call random_number(u)
      value = qBCT(u, mu, sigma, nu, tau)
   end function rBCT

   elemental real(dp) function pe_scale(tau) result(cvalue)
      real(dp), intent(in) :: tau
      cvalue = exp(0.5_dp * (-(2.0_dp / tau) * log(2.0_dp) &
         + log_gamma(1.0_dp / tau) - log_gamma(3.0_dp / tau)))
   end function pe_scale

   elemental real(dp) function pe_logpdf_standard(x, tau) result(value)
      real(dp), intent(in) :: x, tau
      real(dp) :: cvalue
      if (tau <= 0.0_dp) then
         value = nanv()
         return
      end if
      cvalue = pe_scale(tau)
      value = log(tau) - log(cvalue) - 0.5_dp * abs(x / cvalue)**tau &
         - (1.0_dp + 1.0_dp / tau) * log(2.0_dp) - log_gamma(1.0_dp / tau)
   end function pe_logpdf_standard

   elemental real(dp) function pe_cdf_standard(x, tau) result(value)
      real(dp), intent(in) :: x, tau
      real(dp) :: cvalue, shape, svalue
      if (tau <= 0.0_dp) then
         value = nanv()
         return
      end if
      if (x == 0.0_dp) then
         value = 0.5_dp
         return
      end if
      cvalue = pe_scale(tau)
      shape = 1.0_dp / tau
      svalue = 0.5_dp * abs(x / cvalue)**tau
      value = 0.5_dp * (1.0_dp + sign(1.0_dp, x) * pgamma_v(svalue, shape, 1.0_dp))
   end function pe_cdf_standard

   real(dp) function pe_quantile_standard(probability, tau) result(value)
      real(dp), intent(in) :: probability, tau
      real(dp) :: cvalue, shape, gamma_probability, svalue
      if (tau <= 0.0_dp .or. probability < 0.0_dp .or. probability > 1.0_dp) then
         value = nanv()
         return
      end if
      if (probability == 0.0_dp) then
         value = -infv()
         return
      else if (probability == 1.0_dp) then
         value = infv()
         return
      else if (probability == 0.5_dp) then
         value = 0.0_dp
         return
      end if
      cvalue = pe_scale(tau)
      shape = 1.0_dp / tau
      gamma_probability = abs(2.0_dp * probability - 1.0_dp)
      svalue = qgamma_v(gamma_probability, shape, 1.0_dp)
      value = sign(1.0_dp, probability - 0.5_dp) * (2.0_dp * svalue)**(1.0_dp / tau) * cvalue
   end function pe_quantile_standard

   elemental real(dp) function bcpe_norm(sigma, nu, tau) result(value)
      real(dp), intent(in) :: sigma, nu, tau
      if (abs(nu) <= 1.0e-12_dp) then
         value = 1.0_dp
      else
         value = pe_cdf_standard(1.0_dp / (sigma * abs(nu)), tau)
      end if
   end function bcpe_norm

   elemental real(dp) function dBCPE(x, mu, sigma, nu, tau, log_density) result(value)
      real(dp), intent(in) :: x, mu, sigma, nu, tau
      logical, intent(in), optional :: log_density
      real(dp) :: z, log_density_value, norm_value
      logical :: log_value
      log_value = flag_value(log_density, .false.)
      if (mu <= 0.0_dp .or. sigma <= 0.0_dp .or. tau <= 0.0_dp) then
         value = nanv()
         return
      end if
      if (x <= 0.0_dp) then
         value = merge(-infv(), 0.0_dp, log_value)
         return
      end if
      z = boxcox_z(x, mu, sigma, nu)
      norm_value = bcpe_norm(sigma, nu, tau)
      log_density_value = (nu - 1.0_dp) * log(x) - nu * log(mu) - log(sigma) &
         + pe_logpdf_standard(z, tau) - log(norm_value)
      value = merge(log_density_value, exp(log_density_value), log_value)
   end function dBCPE

   elemental real(dp) function pBCPE(q, mu, sigma, nu, tau, lower_tail, log_p) result(value)
      real(dp), intent(in) :: q, mu, sigma, nu, tau
      logical, intent(in), optional :: lower_tail, log_p
      real(dp) :: z, lower_cut, norm_value, probability
      if (mu <= 0.0_dp .or. sigma <= 0.0_dp .or. tau <= 0.0_dp) then
         value = nanv()
         return
      end if
      if (q <= 0.0_dp) then
         value = finish_probability(0.0_dp, lower_tail, log_p)
         return
      end if
      z = boxcox_z(q, mu, sigma, nu)
      norm_value = bcpe_norm(sigma, nu, tau)
      lower_cut = 0.0_dp
      if (nu > 1.0e-12_dp) then
         lower_cut = pe_cdf_standard(-1.0_dp / (sigma * abs(nu)), tau)
      end if
      probability = (pe_cdf_standard(z, tau) - lower_cut) / norm_value
      value = finish_probability(probability, lower_tail, log_p)
   end function pBCPE

   real(dp) function qBCPE(p, mu, sigma, nu, tau, lower_tail, log_p) result(value)
      real(dp), intent(in) :: p, mu, sigma, nu, tau
      logical, intent(in), optional :: lower_tail, log_p
      real(dp) :: probability, z, norm_value, adjusted_probability
      if (mu <= 0.0_dp .or. sigma <= 0.0_dp .or. tau <= 0.0_dp) then
         value = nanv()
         return
      end if
      probability = input_probability(p, lower_tail, log_p)
      if (probability < 0.0_dp .or. probability > 1.0_dp) then
         value = nanv()
         return
      end if
      if (probability == 0.0_dp) then
         value = 0.0_dp
         return
      else if (probability == 1.0_dp) then
         value = infv()
         return
      end if
      if (abs(nu) <= 1.0e-12_dp) then
         z = pe_quantile_standard(probability, tau)
      else
         norm_value = bcpe_norm(sigma, nu, tau)
         if (nu < 0.0_dp) then
            adjusted_probability = probability * norm_value
         else
            adjusted_probability = 1.0_dp - (1.0_dp - probability) * norm_value
         end if
         z = pe_quantile_standard(adjusted_probability, tau)
      end if
      value = inv_boxcox_z(z, mu, sigma, nu)
   end function qBCPE

   real(dp) function rBCPE(mu, sigma, nu, tau) result(value)
      real(dp), intent(in) :: mu, sigma, nu, tau
      real(dp) :: u
      call random_number(u)
      value = qBCPE(u, mu, sigma, nu, tau)
   end function rBCPE

end module gamlss_boxcox
