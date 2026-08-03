! SPDX-License-Identifier: Artistic-2.0
module ldhmm_distribution
   use ldhmm_kinds, only : dp
   use ldhmm_math, only : regularized_gamma_p, gamma_random, quiet_nan
   use ldhmm_types, only : ecld_type
   implicit none
   private

   interface ecld_pdf
      module procedure ecld_pdf_scalar
      module procedure ecld_pdf_vector
   end interface ecld_pdf

   interface ecld_cdf
      module procedure ecld_cdf_scalar
      module procedure ecld_cdf_vector
   end interface ecld_cdf

   interface ecld_ccdf
      module procedure ecld_ccdf_scalar
      module procedure ecld_ccdf_vector
   end interface ecld_ccdf

   public :: ecld_create, ecld_pdf, ecld_cdf, ecld_ccdf
   public :: ecld_mean, ecld_variance, ecld_sd, ecld_skewness
   public :: ecld_kurtosis, ecld_random

contains

   function ecld_create(lambda, sigma, mu, status) result(distribution)
      real(dp), intent(in), optional :: lambda, sigma, mu
      integer, intent(out), optional :: status
      type(ecld_type) :: distribution

      if (present(lambda)) distribution%lambda = lambda
      if (present(sigma)) distribution%sigma = sigma
      if (present(mu)) distribution%mu = mu
      if (present(status)) status = 0
      if (distribution%lambda <= 0.0_dp .or. distribution%sigma <= 0.0_dp) then
         if (present(status)) status = 1
      end if
   end function ecld_create

   elemental real(dp) function ecld_pdf_scalar(distribution, x) result(density)
      type(ecld_type), intent(in) :: distribution
      real(dp), intent(in) :: x
      real(dp) :: beta, z, log_density

      if (distribution%lambda <= 0.0_dp .or. distribution%sigma <= 0.0_dp) then
         density = quiet_nan()
         return
      end if
      beta = 2.0_dp / distribution%lambda
      z = abs((x-distribution%mu)/distribution%sigma)
      log_density = log(beta) - log(2.0_dp*distribution%sigma) - &
         log_gamma(1.0_dp/beta) - z**beta
      density = exp(log_density)
   end function ecld_pdf_scalar

   pure function ecld_pdf_vector(distribution, x) result(density)
      type(ecld_type), intent(in) :: distribution
      real(dp), intent(in) :: x(:)
      real(dp), allocatable :: density(:)
      integer :: i

      allocate(density(size(x)))
      do i = 1, size(x)
         density(i) = ecld_pdf_scalar(distribution, x(i))
      end do
   end function ecld_pdf_vector

   real(dp) function ecld_cdf_scalar(distribution, x) result(probability)
      type(ecld_type), intent(in) :: distribution
      real(dp), intent(in) :: x
      real(dp) :: beta, z, pgamma

      if (distribution%lambda <= 0.0_dp .or. distribution%sigma <= 0.0_dp) then
         probability = quiet_nan()
         return
      end if
      beta = 2.0_dp / distribution%lambda
      z = abs((x-distribution%mu)/distribution%sigma)**beta
      pgamma = regularized_gamma_p(1.0_dp/beta, z)
      if (x < distribution%mu) then
         probability = 0.5_dp * (1.0_dp-pgamma)
      else
         probability = 0.5_dp * (1.0_dp+pgamma)
      end if
      probability = min(1.0_dp, max(0.0_dp, probability))
   end function ecld_cdf_scalar

   function ecld_cdf_vector(distribution, x) result(probability)
      type(ecld_type), intent(in) :: distribution
      real(dp), intent(in) :: x(:)
      real(dp), allocatable :: probability(:)
      integer :: i

      allocate(probability(size(x)))
      do i = 1, size(x)
         probability(i) = ecld_cdf_scalar(distribution, x(i))
      end do
   end function ecld_cdf_vector

   real(dp) function ecld_ccdf_scalar(distribution, x) result(probability)
      type(ecld_type), intent(in) :: distribution
      real(dp), intent(in) :: x
      probability = 1.0_dp - ecld_cdf_scalar(distribution, x)
   end function ecld_ccdf_scalar

   function ecld_ccdf_vector(distribution, x) result(probability)
      type(ecld_type), intent(in) :: distribution
      real(dp), intent(in) :: x(:)
      real(dp), allocatable :: probability(:)
      probability = 1.0_dp - ecld_cdf_vector(distribution, x)
   end function ecld_ccdf_vector

   elemental real(dp) function ecld_mean(distribution) result(value)
      type(ecld_type), intent(in) :: distribution
      value = distribution%mu
   end function ecld_mean

   elemental real(dp) function ecld_variance(distribution) result(value)
      type(ecld_type), intent(in) :: distribution
      real(dp) :: lambda
      lambda = distribution%lambda
      if (lambda <= 0.0_dp .or. distribution%sigma <= 0.0_dp) then
         value = quiet_nan()
      else
         value = distribution%sigma**2 * &
            exp(log_gamma(1.5_dp*lambda)-log_gamma(0.5_dp*lambda))
      end if
   end function ecld_variance

   elemental real(dp) function ecld_sd(distribution) result(value)
      type(ecld_type), intent(in) :: distribution
      value = sqrt(ecld_variance(distribution))
   end function ecld_sd

   elemental real(dp) function ecld_skewness(distribution) result(value)
      type(ecld_type), intent(in) :: distribution
      if (distribution%lambda <= 0.0_dp .or. distribution%sigma <= 0.0_dp) then
         value = quiet_nan()
      else
         value = 0.0_dp
      end if
   end function ecld_skewness

   elemental real(dp) function ecld_kurtosis(distribution) result(value)
      type(ecld_type), intent(in) :: distribution
      real(dp) :: lambda
      lambda = distribution%lambda
      if (lambda <= 0.0_dp) then
         value = quiet_nan()
      else
         value = exp(log_gamma(0.5_dp*lambda) + log_gamma(2.5_dp*lambda) - &
            2.0_dp*log_gamma(1.5_dp*lambda))
      end if
   end function ecld_kurtosis

   real(dp) function ecld_random(distribution) result(value)
      type(ecld_type), intent(in) :: distribution
      real(dp) :: beta, magnitude, u, sign_value

      if (distribution%lambda <= 0.0_dp .or. distribution%sigma <= 0.0_dp) then
         value = quiet_nan()
         return
      end if
      beta = 2.0_dp / distribution%lambda
      magnitude = gamma_random(1.0_dp/beta)**(1.0_dp/beta)
      call random_number(u)
      if (u < 0.5_dp) then
         sign_value = -1.0_dp
      else
         sign_value = 1.0_dp
      end if
      value = distribution%mu + distribution%sigma*sign_value*magnitude
   end function ecld_random

end module ldhmm_distribution
