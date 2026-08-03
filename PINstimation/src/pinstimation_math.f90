! SPDX-License-Identifier: GPL-3.0-or-later
module pinstimation_math
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
   use pinstimation_kinds, only : dp, i8
   implicit none
   private
   public :: logistic, logit, softplus, inv_softplus, log_sum_exp, log_add_exp
   public :: poisson_log_pmf, random_poisson, seed_random_number, normal_cdf
   public :: sample_mean, sample_sd, quantile_sorted, sort_real, finite_number

contains

   pure real(dp) function logistic(x) result(y)
      real(dp), intent(in) :: x
      if (x >= 0.0_dp) then
         y = 1.0_dp/(1.0_dp + exp(-x))
      else
         y = exp(x)/(1.0_dp + exp(x))
      end if
   end function logistic

   pure real(dp) function logit(p) result(x)
      real(dp), intent(in) :: p
      real(dp) :: q
      q = min(1.0_dp - 1.0e-12_dp, max(1.0e-12_dp, p))
      x = log(q/(1.0_dp - q))
   end function logit

   pure real(dp) function softplus(x) result(y)
      real(dp), intent(in) :: x
      if (x > 30.0_dp) then
         y = x
      else if (x < -30.0_dp) then
         y = exp(x)
      else
         y = log(1.0_dp + exp(x))
      end if
   end function softplus

   pure real(dp) function inv_softplus(y) result(x)
      real(dp), intent(in) :: y
      real(dp) :: q
      q = max(y, 1.0e-12_dp)
      if (q > 30.0_dp) then
         x = q
      else
         x = log(exp(q) - 1.0_dp)
      end if
   end function inv_softplus

   pure real(dp) function log_add_exp(a, b) result(v)
      real(dp), intent(in) :: a, b
      real(dp) :: m
      if (a <= -huge(1.0_dp)/2.0_dp) then
         v = b
      else if (b <= -huge(1.0_dp)/2.0_dp) then
         v = a
      else
         m = max(a, b)
         v = m + log(exp(a - m) + exp(b - m))
      end if
   end function log_add_exp

   pure real(dp) function log_sum_exp(x) result(v)
      real(dp), intent(in) :: x(:)
      real(dp) :: m
      if (size(x) == 0) then
         v = -huge(1.0_dp)
         return
      end if
      m = maxval(x)
      if (.not. ieee_is_finite(m)) then
         v = m
      else
         v = m + log(sum(exp(x - m)))
      end if
   end function log_sum_exp

   pure real(dp) function poisson_log_pmf(k, lambda) result(logp)
      integer(i8), intent(in) :: k
      real(dp), intent(in) :: lambda
      if (k < 0_i8 .or. lambda < 0.0_dp .or. .not. ieee_is_finite(lambda)) then
         logp = -huge(1.0_dp)
      else if (lambda <= tiny(1.0_dp)) then
         logp = merge(0.0_dp, -huge(1.0_dp), k == 0_i8)
      else
         logp = -lambda + real(k, dp)*log(lambda) - log_gamma(real(k + 1_i8, dp))
      end if
   end function poisson_log_pmf

   pure real(dp) function normal_cdf(x) result(p)
      real(dp), intent(in) :: x
      p = 0.5_dp*erfc(-x/sqrt(2.0_dp))
   end function normal_cdf

   subroutine seed_random_number(seed)
      integer, intent(in) :: seed
      integer, allocatable :: put(:)
      integer(i8) :: state
      integer(i8), parameter :: modulus = 2147483647_i8
      integer :: n, i
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
      real(dp) :: a, b, inv_alpha, v_r, u, v, us, p, product, sq, log_accept
      integer(i8) :: candidate
      if (present(status)) status = 0
      if (lambda < 0.0_dp .or. .not. ieee_is_finite(lambda)) then
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
         if (us <= tiny(1.0_dp)) cycle
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

   pure real(dp) function sample_mean(x) result(value)
      real(dp), intent(in) :: x(:)
      if (size(x) > 0) then
         value = sum(x)/real(size(x), dp)
      else
         value = ieee_value(value, ieee_quiet_nan)
      end if
   end function sample_mean

   pure real(dp) function sample_sd(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: m
      if (size(x) > 1) then
         m = sum(x)/real(size(x), dp)
         value = sqrt(sum((x - m)**2)/real(size(x) - 1, dp))
      else
         value = 0.0_dp
      end if
   end function sample_sd

   subroutine sort_real(x)
      real(dp), intent(inout) :: x(:)
      integer :: i, j
      real(dp) :: key
      do i = 2, size(x)
         key = x(i)
         j = i - 1
         do while (j >= 1)
            if (x(j) <= key) exit
            x(j + 1) = x(j)
            j = j - 1
         end do
         x(j + 1) = key
      end do
   end subroutine sort_real

   real(dp) function quantile_sorted(x, probability) result(value)
      real(dp), intent(in) :: x(:), probability
      real(dp), allocatable :: work(:)
      real(dp) :: h, frac
      integer :: lo, hi
      if (size(x) == 0) then
         value = ieee_value(value, ieee_quiet_nan)
         return
      end if
      work = x
      call sort_real(work)
      h = 1.0_dp + (real(size(work), dp) - 1.0_dp)*min(1.0_dp, max(0.0_dp, probability))
      lo = int(floor(h))
      hi = int(ceiling(h))
      frac = h - real(lo, dp)
      value = (1.0_dp - frac)*work(lo) + frac*work(hi)
   end function quantile_sorted

end module pinstimation_math
