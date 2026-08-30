! SPDX-License-Identifier: GPL-2.0-or-later
! Computational translation derived from R package spdep 1.4-2.
! Translated/modified 2026-08-30 for modern free-form Fortran.
module spdep_math
   use spdep_kinds, only : dp
   use, intrinsic :: iso_fortran_env, only : int64
   use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
   implicit none
   private

   real(dp), parameter :: pi_dp = acos(-1.0_dp)
   real(dp), parameter :: sqrt2_dp = sqrt(2.0_dp)

   public :: mean_dp
   public :: variance_dp
   public :: normal_cdf
   public :: normal_two_sided_p
   public :: poisson_cdf
   public :: chi_square_sf
   public :: deterministic_seed
   public :: random_permutation
   public :: distance_metric
   public :: matrix_trace
   public :: symmetric_max_eigenvalue
   public :: safe_nan

contains

   pure real(dp) function safe_nan() result(x)
      x = ieee_value(0.0_dp, ieee_quiet_nan)
   end function safe_nan

   pure real(dp) function mean_dp(x) result(mu)
      real(dp), intent(in) :: x(:) !! Numeric observations whose arithmetic mean is required.

      if (size(x) == 0) then
         mu = safe_nan()
      else
         mu = sum(x) / real(size(x), dp)
      end if
   end function mean_dp

   pure real(dp) function variance_dp(x, sample) result(v)
      real(dp), intent(in) :: x(:) !! Numeric observations whose variance is required.
      logical, intent(in), optional :: sample !! If true, divide by n-1; otherwise divide by n.
      real(dp) :: mu
      integer :: denom

      if (size(x) == 0) then
         v = safe_nan()
         return
      end if
      mu = mean_dp(x)
      denom = size(x)
      if (present(sample)) then
         if (sample) denom = denom - 1
      end if
      if (denom <= 0) then
         v = safe_nan()
      else
         v = sum((x - mu) ** 2) / real(denom, dp)
      end if
   end function variance_dp

   pure elemental real(dp) function normal_cdf(x) result(p)
      real(dp), intent(in) :: x !! Standard-normal quantile at which the cumulative probability is evaluated.

      p = 0.5_dp * erfc(-x / sqrt2_dp)
   end function normal_cdf

   pure elemental real(dp) function normal_two_sided_p(z) result(p)
      real(dp), intent(in) :: z !! Standardized statistic for a two-sided normal approximation.

      p = erfc(abs(z) / sqrt2_dp)
   end function normal_two_sided_p

   pure real(dp) function poisson_cdf(k, lambda) result(p)
      integer, intent(in) :: k !! Largest nonnegative integer count included in the cumulative probability.
      real(dp), intent(in) :: lambda !! Poisson mean; must be nonnegative.
      integer :: j
      real(dp) :: term

      if (lambda < 0.0_dp) then
         p = safe_nan()
         return
      end if
      if (k < 0) then
         p = 0.0_dp
         return
      end if
      if (lambda == 0.0_dp) then
         p = 1.0_dp
         return
      end if
      term = exp(-lambda)
      p = term
      do j = 1, k
         term = term * lambda / real(j, dp)
         p = p + term
      end do
      p = min(1.0_dp, max(0.0_dp, p))
   end function poisson_cdf

   pure real(dp) function regularized_gamma_q(a, x) result(q)
      real(dp), intent(in) :: a !! Positive shape parameter of the regularized incomplete gamma function.
      real(dp), intent(in) :: x !! Nonnegative integration limit of the regularized incomplete gamma function.
      integer, parameter :: maxit = 200
      real(dp), parameter :: eps = 8.0_dp * epsilon(1.0_dp)
      real(dp), parameter :: fpmin = tiny(1.0_dp) / eps
      integer :: n
      real(dp) :: ap
      real(dp) :: del
      real(dp) :: sum_series
      real(dp) :: b
      real(dp) :: c
      real(dp) :: d
      real(dp) :: h
      real(dp) :: an
      real(dp) :: prefactor

      if (a <= 0.0_dp .or. x < 0.0_dp) then
         q = safe_nan()
         return
      end if
      if (x == 0.0_dp) then
         q = 1.0_dp
         return
      end if
      prefactor = exp(-x + a * log(x) - log_gamma(a))
      if (x < a + 1.0_dp) then
         ap = a
         del = 1.0_dp / a
         sum_series = del
         do n = 1, maxit
            ap = ap + 1.0_dp
            del = del * x / ap
            sum_series = sum_series + del
            if (abs(del) <= abs(sum_series) * eps) exit
         end do
         q = 1.0_dp - sum_series * prefactor
      else
         b = x + 1.0_dp - a
         c = 1.0_dp / fpmin
         d = 1.0_dp / b
         h = d
         do n = 1, maxit
            an = -real(n, dp) * (real(n, dp) - a)
            b = b + 2.0_dp
            d = an * d + b
            if (abs(d) < fpmin) d = fpmin
            c = b + an / c
            if (abs(c) < fpmin) c = fpmin
            d = 1.0_dp / d
            del = d * c
            h = h * del
            if (abs(del - 1.0_dp) <= eps) exit
         end do
         q = h * prefactor
      end if
      q = min(1.0_dp, max(0.0_dp, q))
   end function regularized_gamma_q

   pure real(dp) function chi_square_sf(x, df) result(p)
      real(dp), intent(in) :: x !! Chi-square statistic; values below zero are outside the distribution support.
      real(dp), intent(in) :: df !! Positive chi-square degrees of freedom.

      if (x < 0.0_dp .or. df <= 0.0_dp) then
         p = safe_nan()
      else
         p = regularized_gamma_q(0.5_dp * df, 0.5_dp * x)
      end if
   end function chi_square_sf

   subroutine deterministic_seed(seed)
      integer, intent(in) :: seed !! Integer seed mapped deterministically onto the processor random-number seed vector.
      integer :: n
      integer :: i
      integer, allocatable :: put(:)
      integer(int64) :: state

      call random_seed(size = n)
      allocate(put(n))
      state = int(seed, int64)
      if (state == 0_int64) state = 104729_int64
      do i = 1, n
         state = modulo(1664525_int64 * state + 1013904223_int64, 2147483647_int64)
         put(i) = int(max(1_int64, state))
      end do
      call random_seed(put = put)
   end subroutine deterministic_seed

   subroutine random_permutation(p)
      integer, intent(out) :: p(:) !! Output permutation containing each integer from 1 through size(p) exactly once.
      integer :: i
      integer :: j
      integer :: tmp
      real(dp) :: u

      do i = 1, size(p)
         p(i) = i
      end do
      do i = size(p), 2, -1
         call random_number(u)
         j = 1 + int(u * real(i, dp))
         j = min(i, j)
         tmp = p(i)
         p(i) = p(j)
         p(j) = tmp
      end do
   end subroutine random_permutation

   pure real(dp) function distance_metric(x, y, method, p) result(d)
      real(dp), intent(in) :: x(:) !! First feature vector; must have the same length as y.
      real(dp), intent(in) :: y(:) !! Second feature vector; must have the same length as x.
      character(len=*), intent(in), optional :: method !! Distance: euclidean, manhattan, maximum, canberra, binary, or minkowski.
      real(dp), intent(in), optional :: p !! Positive Minkowski power used only when method is minkowski; default is 2.
      character(len=16) :: m
      real(dp) :: pp
      real(dp) :: denom
      integer :: i
      integer :: mismatch
      integer :: active

      if (size(x) /= size(y)) then
         d = safe_nan()
         return
      end if
      m = "euclidean"
      if (present(method)) m = adjustl(method)
      select case (trim(m))
      case ("euclidean")
         d = sqrt(sum((x - y) ** 2))
      case ("manhattan")
         d = sum(abs(x - y))
      case ("maximum", "max")
         if (size(x) == 0) then
            d = 0.0_dp
         else
            d = maxval(abs(x - y))
         end if
      case ("canberra")
         d = 0.0_dp
         do i = 1, size(x)
            denom = abs(x(i)) + abs(y(i))
            if (denom > 0.0_dp) d = d + abs(x(i) - y(i)) / denom
         end do
      case ("binary")
         mismatch = 0
         active = 0
         do i = 1, size(x)
            if (x(i) /= 0.0_dp .or. y(i) /= 0.0_dp) then
               active = active + 1
               if ((x(i) /= 0.0_dp) .neqv. (y(i) /= 0.0_dp)) mismatch = mismatch + 1
            end if
         end do
         if (active == 0) then
            d = 0.0_dp
         else
            d = real(mismatch, dp) / real(active, dp)
         end if
      case ("minkowski")
         pp = 2.0_dp
         if (present(p)) pp = p
         if (pp <= 0.0_dp) then
            d = safe_nan()
         else
            d = sum(abs(x - y) ** pp) ** (1.0_dp / pp)
         end if
      case default
         d = safe_nan()
      end select
   end function distance_metric

   pure real(dp) function matrix_trace(a) result(tr)
      real(dp), intent(in) :: a(:, :) !! Square matrix whose diagonal sum is returned.
      integer :: i

      if (size(a, 1) /= size(a, 2)) then
         tr = safe_nan()
         return
      end if
      tr = 0.0_dp
      do i = 1, size(a, 1)
         tr = tr + a(i, i)
      end do
   end function matrix_trace

   pure real(dp) function symmetric_max_eigenvalue(a) result(lambda)
      real(dp), intent(in) :: a(:, :) !! Real symmetric matrix whose largest algebraic eigenvalue is approximated.
      integer, parameter :: maxit = 1000
      real(dp), parameter :: tol = 100.0_dp * epsilon(1.0_dp)
      real(dp), allocatable :: v(:)
      real(dp), allocatable :: w(:)
      real(dp) :: lambda_old
      real(dp) :: normw
      integer :: iter
      integer :: n

      n = size(a, 1)
      if (n == 0 .or. size(a, 2) /= n) then
         lambda = safe_nan()
         return
      end if
      allocate(v(n), w(n))
      v = 1.0_dp / sqrt(real(n, dp))
      lambda_old = huge(1.0_dp)
      lambda = 0.0_dp
      do iter = 1, maxit
         w = matmul(a, v)
         normw = sqrt(sum(w ** 2))
         if (normw <= tiny(1.0_dp)) then
            lambda = 0.0_dp
            return
         end if
         v = w / normw
         lambda = dot_product(v, matmul(a, v))
         if (abs(lambda - lambda_old) <= tol * max(1.0_dp, abs(lambda))) exit
         lambda_old = lambda
      end do
   end function symmetric_max_eigenvalue

end module spdep_math
