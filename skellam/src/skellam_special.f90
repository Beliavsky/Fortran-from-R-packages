! SPDX-License-Identifier: GPL-2.0-or-later
module skellam_special
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
   use skellam_kinds, only : dp, i8, pi, sqrt_two, log_two_pi
   implicit none
   private

   public :: log_bessel_i_integer
   public :: normal_cdf, normal_survival, normal_log_cdf
   public :: log_add_exp, log_sum_exp
   public :: poisson_log_pmf, poisson_cdf
   public :: random_poisson, seed_random_number

contains

   pure real(dp) function log_add_exp(a, b) result(value)
      real(dp), intent(in) :: a, b
      real(dp) :: hi, lo

      hi = max(a, b)
      lo = min(a, b)
      if (hi < -0.5_dp*huge(1.0_dp)) then
         value = hi
      else if (hi - lo > 50.0_dp) then
         value = hi
      else
         value = hi + log(1.0_dp + exp(lo - hi))
      end if
   end function log_add_exp

   pure real(dp) function log_sum_exp(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: xmax

      if (size(x) == 0) then
         value = -huge(1.0_dp)
         return
      end if
      xmax = maxval(x)
      if (xmax < -0.5_dp*huge(1.0_dp)) then
         value = xmax
      else
         value = xmax + log(sum(exp(x - xmax)))
      end if
   end function log_sum_exp

   pure real(dp) function normal_cdf(x) result(p)
      real(dp), intent(in) :: x
      p = 0.5_dp*erfc(-x/sqrt_two)
   end function normal_cdf

   pure real(dp) function normal_survival(x) result(p)
      real(dp), intent(in) :: x
      p = 0.5_dp*erfc(x/sqrt_two)
   end function normal_survival

   pure real(dp) function normal_log_cdf(x) result(logp)
      real(dp), intent(in) :: x
      real(dp) :: invx2, series

      if (x > -8.0_dp) then
         logp = log(normal_cdf(x))
      else
         invx2 = 1.0_dp/(x*x)
         series = 1.0_dp - invx2 + 3.0_dp*invx2**2 - 15.0_dp*invx2**3 + 105.0_dp*invx2**4
         logp = -0.5_dp*x*x - log(-x) - 0.5_dp*log_two_pi + log(max(series, tiny(1.0_dp)))
      end if
   end function normal_log_cdf

   pure real(dp) function log_bessel_i_integer(n, x) result(log_i)
      integer(i8), intent(in) :: n
      real(dp), intent(in) :: x
      integer(i8) :: m
      integer(i8), parameter :: max_terms = 200000_i8
      real(dp) :: log_term, log_sum, ratio_log, mu, corr, inv8x
      logical :: descending

      if (n < 0_i8 .or. x < 0.0_dp) then
         log_i = ieee_nan()
         return
      end if
      if (x <= tiny(1.0_dp)) then
         if (n == 0_i8) then
            log_i = 0.0_dp
         else
            log_i = -huge(1.0_dp)
         end if
         return
      end if

      ! The fixed-order asymptotic expansion is fast and accurate when x is
      ! comfortably larger than the order. Otherwise the positive power
      ! series is accumulated entirely in log space.
      if (x > 5000.0_dp .and. real(n, dp)**2 < 0.10_dp*x) then
         mu = 4.0_dp*real(n, dp)**2
         inv8x = 1.0_dp/(8.0_dp*x)
         corr = 1.0_dp - (mu - 1.0_dp)*inv8x &
              + (mu - 1.0_dp)*(mu - 9.0_dp)*inv8x**2/2.0_dp &
              - (mu - 1.0_dp)*(mu - 9.0_dp)*(mu - 25.0_dp)*inv8x**3/6.0_dp &
              + (mu - 1.0_dp)*(mu - 9.0_dp)*(mu - 25.0_dp)*(mu - 49.0_dp)*inv8x**4/24.0_dp
         if (corr > 0.0_dp) then
            log_i = x - 0.5_dp*log(2.0_dp*pi*x) + log(corr)
            return
         end if
      end if

      log_term = real(n, dp)*log(0.5_dp*x) - log_gamma(real(n + 1_i8, dp))
      log_sum = log_term
      descending = .false.
      do m = 1_i8, max_terms
         ratio_log = 2.0_dp*log(0.5_dp*x) - log(real(m, dp)) - log(real(m + n, dp))
         log_term = log_term + ratio_log
         log_sum = log_add_exp(log_sum, log_term)
         if (ratio_log < 0.0_dp) descending = .true.
         if (descending .and. log_term - log_sum < -40.0_dp) exit
      end do
      log_i = log_sum
   end function log_bessel_i_integer

   pure real(dp) function poisson_log_pmf(k, lambda) result(logp)
      integer(i8), intent(in) :: k
      real(dp), intent(in) :: lambda

      if (lambda < 0.0_dp .or. .not. finite_number(lambda)) then
         logp = ieee_nan()
      else if (k < 0_i8) then
         logp = -huge(1.0_dp)
      else if (lambda <= tiny(1.0_dp)) then
         if (k == 0_i8) then
            logp = 0.0_dp
         else
            logp = -huge(1.0_dp)
         end if
      else
         logp = -lambda + real(k, dp)*log(lambda) - log_gamma(real(k + 1_i8, dp))
      end if
   end function poisson_log_pmf

   real(dp) function poisson_cdf(k, lambda, lower_tail, log_p) result(value)
      integer(i8), intent(in) :: k
      real(dp), intent(in) :: lambda
      logical, intent(in), optional :: lower_tail, log_p
      logical :: lower, logarithm
      integer(i8) :: j, hi
      real(dp), allocatable :: lp(:)
      real(dp) :: lsum, ltotal

      lower = .true.
      logarithm = .false.
      if (present(lower_tail)) lower = lower_tail
      if (present(log_p)) logarithm = log_p

      if (lambda < 0.0_dp .or. .not. finite_number(lambda)) then
         value = ieee_nan()
         return
      end if
      if (k < 0_i8) then
         if (lower) then
            value = merge(-huge(1.0_dp), 0.0_dp, logarithm)
         else
            value = merge(0.0_dp, 1.0_dp, logarithm)
         end if
         return
      end if
      if (lambda <= tiny(1.0_dp)) then
         if (lower) then
            value = merge(0.0_dp, 1.0_dp, logarithm)
         else
            value = merge(-huge(1.0_dp), 0.0_dp, logarithm)
         end if
         return
      end if

      hi = max(k, int(ceiling(lambda + 14.0_dp*sqrt(lambda) + 20.0_dp), i8))
      allocate(lp(0:hi))
      do j = 0_i8, hi
         lp(j) = poisson_log_pmf(j, lambda)
      end do
      ltotal = log_sum_exp(lp)
      if (lower) then
         lsum = log_sum_exp(lp(0:k)) - ltotal
      else if (k == hi) then
         lsum = -huge(1.0_dp)
      else
         lsum = log_sum_exp(lp(k + 1_i8:hi)) - ltotal
      end if
      if (logarithm) then
         value = lsum
      else
         value = exp(lsum)
      end if
   end function poisson_cdf

   subroutine seed_random_number(seed)
      integer, intent(in) :: seed
      integer :: n, i
      integer, allocatable :: put(:)
      integer(i8) :: state
      integer(i8), parameter :: modulus = 2147483647_i8

      call random_seed(size=n)
      allocate(put(n))
      state = modulo(abs(int(seed, i8)), modulus - 1_i8) + 1_i8
      do i = 1, n
         state = modulo(16807_i8*state, modulus)
         put(i) = int(state)
      end do
      call random_seed(put=put)
   end subroutine seed_random_number

   integer(i8) function random_poisson(lambda, status) result(k)
      real(dp), intent(in) :: lambda
      integer, intent(out), optional :: status
      real(dp) :: a, b, inv_alpha, v_r, u, v, us, p, product
      real(dp) :: sq, log_accept
      integer(i8) :: candidate

      if (present(status)) status = 0
      if (lambda < 0.0_dp .or. .not. finite_number(lambda)) then
         k = -1_i8
         if (present(status)) status = 1
         return
      end if
      if (lambda <= tiny(1.0_dp)) then
         k = 0_i8
         return
      end if

      if (lambda < 30.0_dp) then
         product = 1.0_dp
         p = exp(-lambda)
         k = -1_i8
         do
            k = k + 1_i8
            call random_number(u)
            product = product*u
            if (product <= p) exit
         end do
         return
      end if

      sq = sqrt(lambda)
      b = 0.931_dp + 2.53_dp*sq
      a = -0.059_dp + 0.02483_dp*b
      inv_alpha = 1.1239_dp + 1.1328_dp/(b - 3.4_dp)
      v_r = 0.9277_dp - 3.6224_dp/(b - 2.0_dp)
      do
         call random_number(u)
         call random_number(v)
         u = u - 0.5_dp
         us = 0.5_dp - abs(u)
         candidate = int(floor((2.0_dp*a/us + b)*u + lambda + 0.43_dp), i8)
         if (us >= 0.07_dp .and. v <= v_r .and. candidate >= 0_i8) then
            k = candidate
            return
         end if
         if (candidate < 0_i8 .or. (us < 0.013_dp .and. v > us)) cycle
         log_accept = log(v*inv_alpha/(a/(us*us) + b))
         if (log_accept <= poisson_log_pmf(candidate, lambda)) then
            k = candidate
            return
         end if
      end do
   end function random_poisson

   pure logical function finite_number(x) result(ok)
      real(dp), intent(in) :: x
      ok = ieee_is_finite(x)
   end function finite_number

   pure real(dp) function ieee_nan() result(x)
      x = ieee_value(x, ieee_quiet_nan)
   end function ieee_nan

end module skellam_special
