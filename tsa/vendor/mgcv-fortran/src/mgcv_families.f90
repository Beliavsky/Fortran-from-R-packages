! SPDX-License-Identifier: GPL-2.0-or-later
module mgcv_families
   use mgcv_kinds, only : dp
   use mgcv_distributions, only : tweedie_deviance
   use mgcv_utils, only : log1p_stable
   implicit none
   private

   integer, parameter, public :: family_gaussian = 1
   integer, parameter, public :: family_binomial = 2
   integer, parameter, public :: family_poisson = 3
   integer, parameter, public :: family_gamma = 4
   integer, parameter, public :: family_inverse_gaussian = 5
   integer, parameter, public :: family_negative_binomial = 6
   integer, parameter, public :: family_tweedie = 7

   type, public :: family_t
      integer :: id = family_gaussian
      real(dp) :: theta = 1.0_dp
      real(dp) :: p = 1.5_dp
   end type family_t

   public :: family_name, link_inverse, link_function, mu_eta
   public :: variance_function, initialize_mu, deviance_sum, log_likelihood

contains

   function family_name(family) result(name)
      type(family_t), intent(in) :: family
      character(len=:), allocatable :: name
      select case (family%id)
      case (family_gaussian); name = 'gaussian(identity)'
      case (family_binomial); name = 'binomial(logit)'
      case (family_poisson); name = 'poisson(log)'
      case (family_gamma); name = 'Gamma(log)'
      case (family_inverse_gaussian); name = 'inverse.gaussian(log)'
      case (family_negative_binomial); name = 'negative.binomial(log)'
      case (family_tweedie); name = 'Tweedie(log)'
      case default; name = 'unknown'
      end select
   end function family_name

   elemental function link_inverse(eta, family) result(mu)
      real(dp), intent(in) :: eta
      type(family_t), intent(in) :: family
      real(dp) :: mu, e
      select case (family%id)
      case (family_gaussian)
         mu = eta
      case (family_binomial)
         if (eta >= 0.0_dp) then
            e = exp(-min(eta, 700.0_dp)); mu = 1.0_dp / (1.0_dp + e)
         else
            e = exp(max(eta, -700.0_dp)); mu = e / (1.0_dp + e)
         end if
         mu = min(1.0_dp - 1.0e-12_dp, max(1.0e-12_dp, mu))
      case default
         mu = exp(min(100.0_dp, max(-100.0_dp, eta)))
         mu = max(mu, 1.0e-12_dp)
      end select
   end function link_inverse

   elemental function link_function(mu, family) result(eta)
      real(dp), intent(in) :: mu
      type(family_t), intent(in) :: family
      real(dp) :: eta, m
      select case (family%id)
      case (family_gaussian)
         eta = mu
      case (family_binomial)
         m = min(1.0_dp - 1.0e-12_dp, max(1.0e-12_dp, mu))
         eta = log(m / (1.0_dp - m))
      case default
         eta = log(max(mu, 1.0e-12_dp))
      end select
   end function link_function

   elemental function mu_eta(eta, family) result(value)
      real(dp), intent(in) :: eta
      type(family_t), intent(in) :: family
      real(dp) :: value, mu
      mu = link_inverse(eta, family)
      select case (family%id)
      case (family_gaussian); value = 1.0_dp
      case (family_binomial); value = max(1.0e-12_dp, mu * (1.0_dp - mu))
      case default; value = max(1.0e-12_dp, mu)
      end select
   end function mu_eta

   elemental function variance_function(mu, family) result(value)
      real(dp), intent(in) :: mu
      type(family_t), intent(in) :: family
      real(dp) :: value, m
      m = max(mu, 1.0e-12_dp)
      select case (family%id)
      case (family_gaussian); value = 1.0_dp
      case (family_binomial); value = max(1.0e-12_dp, m * (1.0_dp - m))
      case (family_poisson); value = m
      case (family_gamma); value = m * m
      case (family_inverse_gaussian); value = m**3
      case (family_negative_binomial); value = m + m * m / max(family%theta, 1.0e-12_dp)
      case (family_tweedie); value = m**family%p
      case default; value = 1.0_dp
      end select
   end function variance_function

   subroutine initialize_mu(y, family, mu)
      real(dp), intent(in) :: y(:)
      type(family_t), intent(in) :: family
      real(dp), allocatable, intent(out) :: mu(:)
      real(dp) :: avg
      integer :: i
      allocate(mu(size(y)))
      avg = sum(y) / real(max(1, size(y)), dp)
      select case (family%id)
      case (family_gaussian)
         mu = y
      case (family_binomial)
         do i = 1, size(y)
            mu(i) = (0.5_dp + y(i)) / 2.0_dp
            mu(i) = min(0.99_dp, max(0.01_dp, mu(i)))
         end do
      case default
         avg = max(avg, 0.1_dp)
         do i = 1, size(y)
            mu(i) = max(0.1_dp, 0.5_dp * (max(y(i), 0.0_dp) + avg))
         end do
      end select
   end subroutine initialize_mu

   function deviance_sum(y, mu, weights, family) result(dev)
      real(dp), intent(in) :: y(:), mu(:), weights(:)
      type(family_t), intent(in) :: family
      real(dp) :: dev, yi, mi, term
      integer :: i
      dev = 0.0_dp
      do i = 1, size(y)
         yi = y(i); mi = max(mu(i), 1.0e-12_dp)
         select case (family%id)
         case (family_gaussian)
            term = (yi - mi)**2
         case (family_binomial)
            term = 0.0_dp
            if (yi > 0.0_dp) term = term + yi * log(yi / mi)
            if (yi < 1.0_dp) term = term + (1.0_dp - yi) * log((1.0_dp - yi) / (1.0_dp - mi))
            term = 2.0_dp * term
         case (family_poisson)
            if (yi <= tiny(1.0_dp)) then
               term = 2.0_dp * mi
            else
               term = 2.0_dp * (yi * log(yi / mi) - (yi - mi))
            end if
         case (family_gamma)
            if (yi <= 0.0_dp) then
               term = huge(1.0_dp)
            else
               term = 2.0_dp * ((yi - mi) / mi - log(yi / mi))
            end if
         case (family_inverse_gaussian)
            if (yi <= 0.0_dp) then
               term = huge(1.0_dp)
            else
               term = (yi - mi)**2 / (yi * mi * mi)
            end if
         case (family_negative_binomial)
            if (yi <= tiny(1.0_dp)) then
               term = 2.0_dp * family%theta * log1p_stable(mi / family%theta)
            else
               term = 2.0_dp * (yi * log(yi / mi) - &
                  (yi + family%theta) * log((yi + family%theta) / (mi + family%theta)))
            end if
         case (family_tweedie)
            term = tweedie_deviance(yi, mi, family%p)
         case default
            term = (yi - mi)**2
         end select
         dev = dev + weights(i) * term
      end do
   end function deviance_sum

   function log_likelihood(y, mu, weights, family, scale) result(value)
      real(dp), intent(in) :: y(:), mu(:), weights(:)
      type(family_t), intent(in) :: family
      real(dp), intent(in), optional :: scale
      real(dp) :: value, phi, yi, mi, th
      integer :: i
      phi = 1.0_dp; if (present(scale)) phi = max(scale, 1.0e-12_dp)
      value = 0.0_dp
      do i = 1, size(y)
         yi = y(i); mi = max(mu(i), 1.0e-12_dp)
         select case (family%id)
         case (family_gaussian)
            value = value - 0.5_dp * weights(i) * &
               (log(2.0_dp * acos(-1.0_dp) * phi) + (yi - mi)**2 / phi)
         case (family_binomial)
            value = value + weights(i) * (yi * log(mi) + (1.0_dp - yi) * log(1.0_dp - mi))
         case (family_poisson)
            value = value + weights(i) * (yi * log(mi) - mi - log_gamma(yi + 1.0_dp))
         case (family_gamma)
            th = 1.0_dp / phi
            if (yi > 0.0_dp) value = value + weights(i) * &
               (th * log(th / mi) + (th - 1.0_dp) * log(yi) - th * yi / mi - log_gamma(th))
         case (family_negative_binomial)
            th = max(family%theta, 1.0e-12_dp)
            value = value + weights(i) * (log_gamma(yi + th) - log_gamma(th) - log_gamma(yi + 1.0_dp) + &
               th * log(th / (th + mi)) + yi * log(mi / (th + mi)))
         case default
            value = value - 0.5_dp * weights(i) * deviance_sum([yi], [mi], [1.0_dp], family) / phi
         end select
      end do
   end function log_likelihood

end module mgcv_families
