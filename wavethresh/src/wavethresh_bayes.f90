! SPDX-License-Identifier: GPL-2.0-or-later
! Empirical-Bayes wavelet shrinkage translated from wavethresh 4.7.3.
module wavethresh_bayes
   use ieee_arithmetic, only : ieee_is_finite
   use wavethresh_types, only : dp, wd_t
   use wavethresh_transform_1d, only : wd, wr_wd
   use wavethresh_threshold, only : threshold_value
   implicit none
   private

   public :: bayes_thr

contains

   function bayes_thr(data, alpha, beta, filter_number, family, boundary, j0) result(reconstruction)
      real(dp), intent(in) :: data(:) !! Dyadic observations to denoise with the BAYES.THR empirical-Bayes rule.
      real(dp), intent(in), optional :: alpha !! Prior scale-decay exponent; default is 0.5 and must be nonnegative.
      real(dp), intent(in), optional :: beta !! Prior inclusion-decay exponent; default is one and must be nonnegative.
      real(dp), intent(in), optional :: filter_number !! Wavethresh filter number; default is eight.
      character(len=*), intent(in), optional :: family !! Wavethresh filter family; default is DaubLeAsymm.
      character(len=*), intent(in), optional :: boundary !! Boundary mode; the translated numerical path supports periodic only.
      integer, intent(in), optional :: j0 !! First zero-based detail level used for initial universal soft shrinkage;
                                         !! default is five.
      real(dp), allocatable :: reconstruction(:)
      type(wd_t), allocatable :: object
      type(wd_t), allocatable :: universal
      type(wd_t), allocatable :: bayesian
      real(dp), allocatable :: nsignal(:)
      real(dp), allocatable :: sum2(:)
      real(dp), allocatable :: scale_decay(:)
      real(dp), allocatable :: tau2(:)
      real(dp), allocatable :: inclusion_base(:)
      real(dp), allocatable :: prior_probability(:)
      real(dp), allocatable :: shrink_ratio(:)
      real(dp), allocatable :: coefficient(:)
      real(dp), allocatable :: median(:)
      real(dp) :: a
      real(dp) :: b
      real(dp) :: fnum
      real(dp) :: sigma2
      real(dp) :: sigma
      real(dp) :: universal_threshold
      real(dp) :: c1
      real(dp) :: c2
      real(dp) :: likelihood
      real(dp) :: best_likelihood
      real(dp) :: candidate
      real(dp) :: tail_argument
      real(dp) :: tail_probability
      real(dp) :: weight
      real(dp) :: z
      real(dp) :: quantile
      real(dp) :: denominator
      character(len=24) :: fam
      character(len=16) :: bc
      integer :: first_level
      integer :: nlevels
      integer :: level
      integer :: i
      integer :: ngrid
      integer :: ndata

      if (size(data) < 2) then
         allocate(reconstruction(0))
         return
      end if
      a = 0.5_dp
      if (present(alpha)) a = alpha
      b = 1.0_dp
      if (present(beta)) b = beta
      if (a < 0.0_dp .or. b < 0.0_dp) then
         allocate(reconstruction(0))
         return
      end if
      fnum = 8.0_dp
      if (present(filter_number)) fnum = filter_number
      fam = "DaubLeAsymm"
      if (present(family)) fam = family
      bc = "periodic"
      if (present(boundary)) bc = boundary
      first_level = 5
      if (present(j0)) first_level = j0

      allocate(object, universal, bayesian)
      object = wd(data, filter_number=fnum, family=trim(fam), boundary=trim(bc))
      if (.not. object%ok .or. object%nlevels < 1) then
         allocate(reconstruction(0))
         return
      end if
      nlevels = object%nlevels
      ndata = size(data)
      first_level = max(0, min(first_level, nlevels - 1))
      sigma2 = sample_variance(object%detail(nlevels - 1)%values)
      sigma = sqrt(max(sigma2, 0.0_dp))
      if (sigma <= tiny(1.0_dp)) then
         reconstruction = data
         return
      end if
      universal_threshold = sqrt(2.0_dp * log(real(ndata / 2, dp))) * sigma
      universal = object
      do level = first_level, nlevels - 1
         universal%detail(level)%values = threshold_value(universal%detail(level)%values, universal_threshold, "soft")
      end do

      allocate(nsignal(0:nlevels - 1), source=0.0_dp)
      allocate(sum2(0:nlevels - 1), source=0.0_dp)
      allocate(scale_decay(0:nlevels - 1), source=0.0_dp)
      do level = 0, nlevels - 1
         coefficient = universal%detail(level)%values
         nsignal(level) = real(count(abs(coefficient) > 0.0_dp), dp)
         if (nsignal(level) > 0.0_dp) then
            sum2(level) = sum(coefficient**2, mask=abs(coefficient) > 0.0_dp)
         end if
         scale_decay(level) = 2.0_dp**(-a * real(level, dp))
      end do

      best_likelihood = -huge(1.0_dp)
      c1 = 1000.0_dp
      ngrid = (15000 - 1000) / 50 + 1
      do i = 0, ngrid - 1
         candidate = 1000.0_dp + 50.0_dp * real(i, dp)
         likelihood = 0.0_dp
         do level = 0, nlevels - 1
            denominator = sigma2 + candidate * scale_decay(level)
            tail_argument = -sigma * sqrt(2.0_dp * log(real(ndata, dp))) / sqrt(denominator)
            tail_probability = max(2.0_dp * standard_normal_cdf(tail_argument), tiny(1.0_dp))
            likelihood = likelihood - nsignal(level) * (log(denominator) + 2.0_dp * log(tail_probability)) - &
               sum2(level) / (2.0_dp * denominator)
         end do
         likelihood = 0.5_dp * likelihood
         if (likelihood > best_likelihood) then
            best_likelihood = likelihood
            c1 = candidate
         end if
      end do

      allocate(tau2(0:nlevels - 1), source=0.0_dp)
      allocate(inclusion_base(0:nlevels - 1), source=0.0_dp)
      allocate(prior_probability(0:nlevels - 1), source=0.0_dp)
      allocate(shrink_ratio(0:nlevels - 1), source=0.0_dp)
      tau2 = c1 * scale_decay
      do level = 0, nlevels - 1
         tail_argument = -sigma * sqrt(2.0_dp * log(real(ndata, dp))) / sqrt(sigma2 + tau2(level))
         inclusion_base(level) = max(2.0_dp * standard_normal_cdf(tail_argument), tiny(1.0_dp))
      end do
      if (abs(b - 1.0_dp) <= 32.0_dp * epsilon(1.0_dp)) then
         c2 = sum(nsignal / inclusion_base) / real(nlevels, dp)
      else
         c2 = (1.0_dp - 2.0_dp**(1.0_dp - b)) / &
            (1.0_dp - 2.0_dp**((1.0_dp - b) * real(nlevels, dp))) * sum(nsignal / inclusion_base)
      end if
      do level = 0, nlevels - 1
         prior_probability(level) = min(1.0_dp, max(0.0_dp, c2 * 2.0_dp**(-b * real(level, dp))))
         shrink_ratio(level) = tau2(level) / (sigma2 + tau2(level))
      end do

      bayesian = object
      do level = 0, nlevels - 1
         coefficient = object%detail(level)%values
         allocate(median(size(coefficient)), source=0.0_dp)
         do i = 1, size(coefficient)
            if (prior_probability(level) <= tiny(1.0_dp)) cycle
            weight = (1.0_dp - prior_probability(level)) / prior_probability(level)
            weight = weight / sqrt((sigma2 * shrink_ratio(level)) / tau2(level))
            weight = weight * exp(-shrink_ratio(level) * coefficient(i)**2 / (2.0_dp * sigma2))
            z = 0.5_dp * (1.0_dp + min(weight, 1.0_dp))
            if (z >= 1.0_dp - epsilon(1.0_dp)) cycle
            quantile = standard_normal_quantile(z)
            if (.not. ieee_is_finite(quantile)) cycle
            median(i) = sign(1.0_dp, coefficient(i)) * max(0.0_dp, &
               shrink_ratio(level) * abs(coefficient(i)) - sigma * sqrt(shrink_ratio(level)) * quantile)
         end do
         bayesian%detail(level)%values = median
         deallocate(median)
      end do
      reconstruction = wr_wd(bayesian)
   end function bayes_thr

   pure function sample_variance(x) result(value)
      real(dp), intent(in) :: x(:) !! Sample whose unbiased variance with denominator n-1 is requested.
      real(dp) :: value
      real(dp) :: mean_value

      if (size(x) <= 1) then
         value = 0.0_dp
         return
      end if
      mean_value = sum(x) / real(size(x), dp)
      value = sum((x - mean_value)**2) / real(size(x) - 1, dp)
   end function sample_variance

   pure elemental function standard_normal_cdf(x) result(probability)
      real(dp), intent(in) :: x !! Standard-normal variate at which the cumulative probability is evaluated.
      real(dp) :: probability

      probability = 0.5_dp * erfc(-x / sqrt(2.0_dp))
   end function standard_normal_cdf

   pure elemental function standard_normal_quantile(probability) result(x)
      real(dp), intent(in) :: probability !! Standard-normal probability strictly between zero and one.
      real(dp) :: x
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

      if (probability <= 0.0_dp) then
         x = -huge(1.0_dp)
      else if (probability >= 1.0_dp) then
         x = huge(1.0_dp)
      else if (probability < plow) then
         q = sqrt(-2.0_dp * log(probability))
         x = (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / &
            ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0_dp)
      else if (probability <= phigh) then
         q = probability - 0.5_dp
         r = q * q
         x = (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q / &
            (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1.0_dp)
      else
         q = sqrt(-2.0_dp * log(1.0_dp - probability))
         x = -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / &
            ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0_dp)
      end if
   end function standard_normal_quantile

end module wavethresh_bayes
