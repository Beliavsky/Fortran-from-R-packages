! SPDX-License-Identifier: GPL-2.0-or-later
module skellam_distribution
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
   use skellam_kinds, only : dp, i8, pi, log_two_pi
   use skellam_special, only : log_bessel_i_integer, log_sum_exp, log_add_exp, &
      normal_cdf, random_poisson
   implicit none
   private

   interface dskellam
      module procedure dskellam_scalar
      module procedure dskellam_vector
   end interface dskellam

   interface pskellam
      module procedure pskellam_scalar
      module procedure pskellam_vector
   end interface pskellam

   interface qskellam
      module procedure qskellam_scalar
      module procedure qskellam_vector
   end interface qskellam

   public :: dskellam, pskellam, qskellam, rskellam
   public :: dskellam_sp, pskellam_sp
   public :: skellam_log_pmf
   public :: skellam_mean, skellam_variance, skellam_skewness, skellam_excess_kurtosis

contains

   pure real(dp) function skellam_mean(lambda1, lambda2) result(value)
      real(dp), intent(in) :: lambda1, lambda2
      value = lambda1 - lambda2
   end function skellam_mean

   pure real(dp) function skellam_variance(lambda1, lambda2) result(value)
      real(dp), intent(in) :: lambda1, lambda2
      value = lambda1 + lambda2
   end function skellam_variance

   pure real(dp) function skellam_skewness(lambda1, lambda2) result(value)
      real(dp), intent(in) :: lambda1, lambda2
      real(dp) :: variance
      variance = lambda1 + lambda2
      if (variance > 0.0_dp) then
         value = (lambda1 - lambda2)/variance**1.5_dp
      else
         value = 0.0_dp
      end if
   end function skellam_skewness

   pure real(dp) function skellam_excess_kurtosis(lambda1, lambda2) result(value)
      real(dp), intent(in) :: lambda1, lambda2
      real(dp) :: variance
      variance = lambda1 + lambda2
      if (variance > 0.0_dp) then
         value = 1.0_dp/variance
      else
         value = 0.0_dp
      end if
   end function skellam_excess_kurtosis

   real(dp) function skellam_log_pmf(k, lambda1, lambda2, status) result(logp)
      integer(i8), intent(in) :: k
      real(dp), intent(in) :: lambda1, lambda2
      integer, intent(out), optional :: status
      real(dp) :: argument

      if (present(status)) status = 0
      if (.not. valid_rates(lambda1, lambda2)) then
         logp = quiet_nan()
         if (present(status)) status = 1
         return
      end if

      if (lambda1 <= tiny(1.0_dp) .and. lambda2 <= tiny(1.0_dp)) then
         if (k == 0_i8) then
            logp = 0.0_dp
         else
            logp = -huge(1.0_dp)
         end if
      else if (lambda1 <= tiny(1.0_dp)) then
         if (k <= 0_i8) then
            logp = -lambda2 + real(-k, dp)*log(lambda2) - log_gamma(real(-k + 1_i8, dp))
         else
            logp = -huge(1.0_dp)
         end if
      else if (lambda2 <= tiny(1.0_dp)) then
         if (k >= 0_i8) then
            logp = -lambda1 + real(k, dp)*log(lambda1) - log_gamma(real(k + 1_i8, dp))
         else
            logp = -huge(1.0_dp)
         end if
      else
         argument = 2.0_dp*sqrt(lambda1*lambda2)
         logp = -(lambda1 + lambda2) + 0.5_dp*real(k, dp)*log(lambda1/lambda2) &
              + log_bessel_i_integer(abs(k), argument)
         if (.not. ieee_is_finite(logp)) logp = dskellam_sp(k, lambda1, lambda2, log_p=.true.)
      end if
   end function skellam_log_pmf

   real(dp) function dskellam_scalar(k, lambda1, lambda2, log_p, status) result(value)
      integer(i8), intent(in) :: k
      real(dp), intent(in) :: lambda1
      real(dp), intent(in), optional :: lambda2
      logical, intent(in), optional :: log_p
      integer, intent(out), optional :: status
      real(dp) :: lambda2_use, log_value
      logical :: logarithm

      lambda2_use = lambda1
      if (present(lambda2)) lambda2_use = lambda2
      logarithm = .false.
      if (present(log_p)) logarithm = log_p
      log_value = skellam_log_pmf(k, lambda1, lambda2_use, status)
      if (logarithm) then
         value = log_value
      else if (log_value < log(tiny(1.0_dp))) then
         value = 0.0_dp
      else
         value = exp(log_value)
      end if
   end function dskellam_scalar

   function dskellam_vector(k, lambda1, lambda2, log_p, status) result(value)
      integer(i8), intent(in) :: k(:)
      real(dp), intent(in) :: lambda1
      real(dp), intent(in), optional :: lambda2
      logical, intent(in), optional :: log_p
      integer, intent(out), optional :: status
      real(dp) :: value(size(k))
      integer :: i, local_status
      real(dp) :: lambda2_use

      lambda2_use = lambda1
      if (present(lambda2)) lambda2_use = lambda2
      if (present(status)) status = 0
      do i = 1, size(k)
         value(i) = dskellam_scalar(k(i), lambda1, lambda2_use, log_p, local_status)
         if (present(status)) status = max(status, local_status)
      end do
   end function dskellam_vector

   real(dp) function dskellam_sp(k, lambda1, lambda2, log_p, status) result(value)
      integer(i8), intent(in) :: k
      real(dp), intent(in) :: lambda1, lambda2
      logical, intent(in), optional :: log_p
      integer, intent(out), optional :: status
      logical :: logarithm
      real(dp) :: x, s, cumulant, cumulant2, ratio, correction, log_density

      logarithm = .false.
      if (present(log_p)) logarithm = log_p
      if (present(status)) status = 0
      if (.not. valid_rates(lambda1, lambda2)) then
         value = quiet_nan()
         if (present(status)) status = 1
         return
      end if
      if (lambda1 <= tiny(1.0_dp) .or. lambda2 <= tiny(1.0_dp)) then
         value = dskellam_scalar(k, lambda1, lambda2, logarithm, status)
         return
      end if

      x = real(k, dp)
      s = log((x + sqrt(x*x + 4.0_dp*lambda1*lambda2))/(2.0_dp*lambda1))
      cumulant = lambda1*(exp(s) - 1.0_dp) + lambda2*(exp(-s) - 1.0_dp)
      cumulant2 = lambda1*exp(s) + lambda2*exp(-s)
      ratio = (lambda1*exp(s) - lambda2*exp(-s))/cumulant2
      correction = 1.0_dp + (1.0_dp - (5.0_dp/3.0_dp)*ratio*ratio)/(8.0_dp*cumulant2)
      log_density = cumulant - x*s - 0.5_dp*(log_two_pi + log(cumulant2)) &
         + log(max(0.5_dp*(1.0_dp + correction), tiny(1.0_dp)))
      if (logarithm) then
         value = log_density
      else
         value = exp(log_density)
      end if
   end function dskellam_sp

   real(dp) function pskellam_scalar(q, lambda1, lambda2, lower_tail, log_p, status) result(value)
      real(dp), intent(in) :: q
      real(dp), intent(in) :: lambda1
      real(dp), intent(in), optional :: lambda2
      logical, intent(in), optional :: lower_tail, log_p
      integer, intent(out), optional :: status
      real(dp) :: lambda2_use, log_total, log_part
      real(dp), allocatable :: log_mass(:)
      integer(i8) :: k, lo, hi, iq
      integer :: n, first_index, last_index
      logical :: lower, logarithm

      lambda2_use = lambda1
      if (present(lambda2)) lambda2_use = lambda2
      lower = .true.
      logarithm = .false.
      if (present(lower_tail)) lower = lower_tail
      if (present(log_p)) logarithm = log_p
      if (present(status)) status = 0

      if (.not. valid_rates(lambda1, lambda2_use) .or. .not. ieee_is_finite(q)) then
         value = quiet_nan()
         if (present(status)) status = 1
         return
      end if
      iq = int(floor(q), i8)
      call distribution_support(lambda1, lambda2_use, 1.0e-16_dp, lo, hi)
      if (iq < lo) then
         if (lower) then
            value = merge(-huge(1.0_dp), 0.0_dp, logarithm)
         else
            value = merge(0.0_dp, 1.0_dp, logarithm)
         end if
         return
      else if (iq >= hi) then
         if (lower) then
            value = merge(0.0_dp, 1.0_dp, logarithm)
         else
            value = merge(-huge(1.0_dp), 0.0_dp, logarithm)
         end if
         return
      end if

      n = int(hi - lo + 1_i8)
      allocate(log_mass(n))
      do k = lo, hi
         log_mass(int(k - lo + 1_i8)) = skellam_log_pmf(k, lambda1, lambda2_use)
      end do
      log_total = log_sum_exp(log_mass)
      if (lower) then
         first_index = 1
         last_index = int(iq - lo + 1_i8)
      else
         first_index = int(iq - lo + 2_i8)
         last_index = n
      end if
      log_part = log_sum_exp(log_mass(first_index:last_index)) - log_total
      if (logarithm) then
         value = min(0.0_dp, log_part)
      else
         value = min(1.0_dp, max(0.0_dp, exp(log_part)))
      end if
   end function pskellam_scalar

   function pskellam_vector(q, lambda1, lambda2, lower_tail, log_p, status) result(value)
      real(dp), intent(in) :: q(:)
      real(dp), intent(in) :: lambda1
      real(dp), intent(in), optional :: lambda2
      logical, intent(in), optional :: lower_tail, log_p
      integer, intent(out), optional :: status
      real(dp) :: value(size(q))
      real(dp) :: lambda2_use
      integer :: i, local_status

      lambda2_use = lambda1
      if (present(lambda2)) lambda2_use = lambda2
      if (present(status)) status = 0
      do i = 1, size(q)
         value(i) = pskellam_scalar(q(i), lambda1, lambda2_use, lower_tail, log_p, local_status)
         if (present(status)) status = max(status, local_status)
      end do
   end function pskellam_vector

   real(dp) function pskellam_sp(q, lambda1, lambda2, lower_tail, log_p, status) result(value)
      real(dp), intent(in) :: q, lambda1, lambda2
      logical, intent(in), optional :: lower_tail, log_p
      integer, intent(out), optional :: status
      logical :: lower, logarithm
      real(dp) :: p

      lower = .true.
      logarithm = .false.
      if (present(lower_tail)) lower = lower_tail
      if (present(log_p)) logarithm = log_p
      if (present(status)) status = 0
      if (.not. valid_rates(lambda1, lambda2) .or. .not. ieee_is_finite(q)) then
         value = quiet_nan()
         if (present(status)) status = 1
         return
      end if
      if (lambda1 <= tiny(1.0_dp) .or. lambda2 <= tiny(1.0_dp)) then
         value = pskellam_scalar(q, lambda1, lambda2, lower, logarithm, status)
         return
      end if

      if (lower) then
         p = saddlepoint_lower(q, lambda1, lambda2)
      else
         p = saddlepoint_lower(-q - 1.0_dp, lambda2, lambda1)
      end if
      p = min(1.0_dp, max(0.0_dp, p))
      if (logarithm) then
         if (p <= 0.0_dp) then
            value = -huge(1.0_dp)
         else
            value = log(p)
         end if
      else
         value = p
      end if
   end function pskellam_sp

   integer(i8) function qskellam_scalar(p, lambda1, lambda2, lower_tail, log_p, status) result(quantile)
      real(dp), intent(in) :: p, lambda1
      real(dp), intent(in), optional :: lambda2
      logical, intent(in), optional :: lower_tail, log_p
      integer, intent(out), optional :: status
      real(dp) :: lambda2_use, probability, target, tolerance, log_target
      real(dp), allocatable :: log_mass(:)
      real(dp) :: log_total, log_cumulative
      integer(i8) :: lo, hi, k
      integer :: i, n
      logical :: lower, logarithm

      lambda2_use = lambda1
      if (present(lambda2)) lambda2_use = lambda2
      lower = .true.
      logarithm = .false.
      if (present(lower_tail)) lower = lower_tail
      if (present(log_p)) logarithm = log_p
      if (present(status)) status = 0

      if (.not. valid_rates(lambda1, lambda2_use) .or. .not. ieee_is_finite(p)) then
         quantile = 0_i8
         if (present(status)) status = 1
         return
      end if
      if (logarithm) then
         if (p > 0.0_dp) then
            quantile = 0_i8
            if (present(status)) status = 2
            return
         end if
         probability = exp(p)
      else
         probability = p
      end if
      if (probability < 0.0_dp .or. probability > 1.0_dp) then
         quantile = 0_i8
         if (present(status)) status = 2
         return
      end if
      target = merge(probability, 1.0_dp - probability, lower)
      if (target <= 0.0_dp) then
         if (lambda2_use <= tiny(1.0_dp)) then
            quantile = 0_i8
         else
            quantile = -huge(1_i8)
         end if
         return
      else if (target >= 1.0_dp) then
         if (lambda1 <= tiny(1.0_dp)) then
            quantile = 0_i8
         else
            quantile = huge(1_i8)
         end if
         return
      end if

      tolerance = max(tiny(1.0_dp), min(1.0e-16_dp, 1.0e-4_dp*min(target, 1.0_dp - target)))
      call distribution_support(lambda1, lambda2_use, tolerance, lo, hi)
      n = int(hi - lo + 1_i8)
      allocate(log_mass(n))
      do k = lo, hi
         log_mass(int(k - lo + 1_i8)) = skellam_log_pmf(k, lambda1, lambda2_use)
      end do
      log_total = log_sum_exp(log_mass)
      log_target = log(target)
      log_cumulative = -huge(1.0_dp)
      do i = 1, n
         log_cumulative = log_add_exp(log_cumulative, log_mass(i))
         if (log_cumulative - log_total >= log_target - 32.0_dp*epsilon(1.0_dp)) then
            quantile = lo + int(i - 1, i8)
            return
         end if
      end do
      quantile = hi
   end function qskellam_scalar

   function qskellam_vector(p, lambda1, lambda2, lower_tail, log_p, status) result(quantile)
      real(dp), intent(in) :: p(:)
      real(dp), intent(in) :: lambda1
      real(dp), intent(in), optional :: lambda2
      logical, intent(in), optional :: lower_tail, log_p
      integer, intent(out), optional :: status
      integer(i8) :: quantile(size(p))
      real(dp) :: lambda2_use
      integer :: i, local_status

      lambda2_use = lambda1
      if (present(lambda2)) lambda2_use = lambda2
      if (present(status)) status = 0
      do i = 1, size(p)
         quantile(i) = qskellam_scalar(p(i), lambda1, lambda2_use, lower_tail, log_p, local_status)
         if (present(status)) status = max(status, local_status)
      end do
   end function qskellam_vector

   function rskellam(n, lambda1, lambda2, status) result(values)
      integer, intent(in) :: n
      real(dp), intent(in) :: lambda1
      real(dp), intent(in), optional :: lambda2
      integer, intent(out), optional :: status
      integer(i8), allocatable :: values(:)
      real(dp) :: lambda2_use
      integer :: i, status1, status2
      integer(i8) :: x1, x2

      lambda2_use = lambda1
      if (present(lambda2)) lambda2_use = lambda2
      if (present(status)) status = 0
      allocate(values(max(0, n)))
      if (n < 0 .or. .not. valid_rates(lambda1, lambda2_use)) then
         if (present(status)) status = 1
         return
      end if
      do i = 1, n
         x1 = random_poisson(lambda1, status1)
         x2 = random_poisson(lambda2_use, status2)
         values(i) = x1 - x2
         if (present(status)) status = max(status, max(status1, status2))
      end do
   end function rskellam

   subroutine distribution_support(lambda1, lambda2, tolerance, lo, hi)
      real(dp), intent(in) :: lambda1, lambda2, tolerance
      integer(i8), intent(out) :: lo, hi
      real(dp) :: mean, standard_deviation, multiplier, width, tol

      mean = lambda1 - lambda2
      standard_deviation = sqrt(lambda1 + lambda2)
      tol = max(tiny(1.0_dp), min(0.01_dp, tolerance))
      multiplier = sqrt(max(0.0_dp, -2.0_dp*log(tol))) + 6.0_dp
      width = max(20.0_dp, multiplier*standard_deviation + 10.0_dp)
      lo = int(floor(mean - width), i8)
      hi = int(ceiling(mean + width), i8)
      if (lambda2 <= tiny(1.0_dp)) lo = 0_i8
      if (lambda1 <= tiny(1.0_dp)) hi = 0_i8
   end subroutine distribution_support

   pure real(dp) function saddlepoint_lower(q, lambda1, lambda2) result(p)
      real(dp), intent(in) :: q, lambda1, lambda2
      real(dp) :: xm, s, cumulant, cumulant2, u, w, standardized, skewness, density

      xm = -floor(q) - 0.5_dp
      s = log((xm + sqrt(xm*xm + 4.0_dp*lambda2*lambda1))/(2.0_dp*lambda2))
      cumulant = lambda2*(exp(s) - 1.0_dp) + lambda1*(exp(-s) - 1.0_dp)
      cumulant2 = lambda2*exp(s) + lambda1*exp(-s)
      standardized = (xm + lambda1 - lambda2)/sqrt(lambda1 + lambda2)
      skewness = (lambda1 - lambda2)/(lambda1 + lambda2)**1.5_dp
      if (abs(standardized) < 1.0e-4_dp .or. abs(s) < 1.0e-8_dp) then
         density = exp(-0.5_dp*standardized*standardized)/sqrt(2.0_dp*pi)
         p = normal_cdf(-standardized) + density*skewness*(1.0_dp - standardized*standardized)/6.0_dp
      else
         u = 2.0_dp*sinh(0.5_dp*s)*sqrt(cumulant2)
         w = sign(sqrt(max(0.0_dp, 2.0_dp*(s*xm - cumulant))), s)
         density = exp(-0.5_dp*w*w)/sqrt(2.0_dp*pi)
         p = normal_cdf(-w) - density*(1.0_dp/w - 1.0_dp/u)
      end if
   end function saddlepoint_lower

   pure logical function valid_rates(lambda1, lambda2) result(valid)
      real(dp), intent(in) :: lambda1, lambda2
      valid = ieee_is_finite(lambda1) .and. ieee_is_finite(lambda2) &
         .and. lambda1 >= 0.0_dp .and. lambda2 >= 0.0_dp
   end function valid_rates

   pure real(dp) function quiet_nan() result(value)
      value = ieee_value(value, ieee_quiet_nan)
   end function quiet_nan

end module skellam_distribution
