! SPDX-License-Identifier: GPL-3.0-only
module ufrisk_math
   use kind_mod, only : dp
   use stats_mod, only : normal_quantile
   use special_functions_mod, only : regularized_beta, regularized_gamma_q
   implicit none
   private
   real(dp), parameter :: pi = acos(-1.0_dp)
   public :: normal_cdf, normal_pdf, normal_quantile
   public :: student_cdf, student_pdf, student_quantile
   public :: binomial_cdf, chi_square_cdf, sample_standard_deviation
   public :: log_probability, finite_vector
contains
   pure elemental real(dp) function normal_cdf(x) result(value)
      real(dp), intent(in) :: x
      value = 0.5_dp*erfc(-x/sqrt(2.0_dp))
   end function normal_cdf

   pure elemental real(dp) function normal_pdf(x) result(value)
      real(dp), intent(in) :: x
      value = exp(-0.5_dp*x*x)/sqrt(2.0_dp*pi)
   end function normal_pdf

   pure elemental real(dp) function student_pdf(x, degrees_freedom) result(value)
      real(dp), intent(in) :: x, degrees_freedom
      real(dp) :: nu
      nu = max(degrees_freedom, tiny(1.0_dp))
      value = exp(log_gamma(0.5_dp*(nu+1.0_dp)) - log_gamma(0.5_dp*nu) - &
         0.5_dp*log(nu*pi) - 0.5_dp*(nu+1.0_dp)*log(1.0_dp+x*x/nu))
   end function student_pdf

   pure elemental real(dp) function student_cdf(x, degrees_freedom) result(value)
      real(dp), intent(in) :: x, degrees_freedom
      real(dp) :: nu, ib
      nu = max(degrees_freedom, tiny(1.0_dp))
      ib = regularized_beta(nu/(nu+x*x), 0.5_dp*nu, 0.5_dp)
      if (x >= 0.0_dp) then
         value = 1.0_dp - 0.5_dp*ib
      else
         value = 0.5_dp*ib
      end if
      value = max(0.0_dp, min(1.0_dp, value))
   end function student_cdf

   pure real(dp) function student_quantile(probability, degrees_freedom) result(value)
      real(dp), intent(in) :: probability, degrees_freedom
      real(dp) :: lower, upper, middle, cdf_value
      integer :: iteration
      if (probability <= 0.0_dp) then
         value = -huge(1.0_dp)
         return
      else if (probability >= 1.0_dp) then
         value = huge(1.0_dp)
         return
      else if (abs(probability-0.5_dp) <= epsilon(1.0_dp)) then
         value = 0.0_dp
         return
      end if
      lower = -1.0_dp
      upper = 1.0_dp
      do while (student_cdf(lower, degrees_freedom) > probability)
         lower = 2.0_dp*lower
         if (lower < -1.0e12_dp) exit
      end do
      do while (student_cdf(upper, degrees_freedom) < probability)
         upper = 2.0_dp*upper
         if (upper > 1.0e12_dp) exit
      end do
      do iteration = 1, 160
         middle = 0.5_dp*(lower+upper)
         cdf_value = student_cdf(middle, degrees_freedom)
         if (cdf_value < probability) then
            lower = middle
         else
            upper = middle
         end if
         if (upper-lower <= 4.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(middle))) exit
      end do
      value = 0.5_dp*(lower+upper)
   end function student_quantile

   pure real(dp) function binomial_cdf(successes, trials, probability) result(value)
      integer, intent(in) :: successes, trials
      real(dp), intent(in) :: probability
      integer :: k
      real(dp) :: term
      if (trials < 0 .or. successes < 0) then
         value = 0.0_dp
      else if (successes >= trials) then
         value = 1.0_dp
      else if (probability <= 0.0_dp) then
         value = 1.0_dp
      else if (probability >= 1.0_dp) then
         value = 0.0_dp
      else
         term = exp(real(trials,dp)*log(1.0_dp-probability))
         value = term
         do k = 0, successes-1
            term = term*real(trials-k,dp)/real(k+1,dp)*probability/(1.0_dp-probability)
            value = value + term
         end do
         value = max(0.0_dp,min(1.0_dp,value))
      end if
   end function binomial_cdf

   pure real(dp) function chi_square_cdf(x, degrees_freedom) result(value)
      real(dp), intent(in) :: x
      integer, intent(in) :: degrees_freedom
      if (x <= 0.0_dp) then
         value = 0.0_dp
      else
         value = 1.0_dp-regularized_gamma_q(0.5_dp*real(degrees_freedom,dp),0.5_dp*x)
      end if
      value = max(0.0_dp,min(1.0_dp,value))
   end function chi_square_cdf

   pure real(dp) function sample_standard_deviation(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: mean_x
      if (size(x) < 2) then
         value = 0.0_dp
      else
         mean_x = sum(x)/real(size(x),dp)
         value = sqrt(sum((x-mean_x)**2)/real(size(x)-1,dp))
      end if
   end function sample_standard_deviation

   pure real(dp) function log_probability(probability, count) result(value)
      real(dp), intent(in) :: probability
      integer, intent(in) :: count
      if (count == 0) then
         value = 0.0_dp
      else if (probability <= 0.0_dp) then
         value = -huge(1.0_dp)
      else
         value = real(count,dp)*log(probability)
      end if
   end function log_probability

   pure logical function finite_vector(x) result(value)
      real(dp), intent(in) :: x(:)
      value = all(abs(x) <= huge(1.0_dp))
   end function finite_vector
end module ufrisk_math
