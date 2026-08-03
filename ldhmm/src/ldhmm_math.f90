! SPDX-License-Identifier: Artistic-2.0
module ldhmm_math
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_quiet_nan, ieee_value
   use ldhmm_kinds, only : dp
   implicit none
   private

   real(dp), parameter :: pi = acos(-1.0_dp)
   real(dp), parameter :: tiny_prob = 1.0e-300_dp

   public :: pi, quiet_nan, log_sum_exp, regularized_gamma_p
   public :: normal_quantile, seed_random, gamma_random, categorical_random
   public :: vector_mean, sample_sd, population_skewness, population_kurtosis
   public :: correlation, absolute_acf, simple_moving_average, remove_outliers
   public :: solve_linear_system, finite_count, finite_mean_abs

contains

   pure real(dp) function quiet_nan() result(x)
      x = ieee_value(0.0_dp, ieee_quiet_nan)
   end function quiet_nan

   real(dp) function log_sum_exp(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: xmax
      integer :: i

      if (size(x) == 0) then
         value = -huge(1.0_dp)
         return
      end if
      xmax = maxval(x)
      if (.not. ieee_is_finite(xmax)) then
         value = xmax
         return
      end if
      value = 0.0_dp
      do i = 1, size(x)
         value = value + exp(x(i) - xmax)
      end do
      value = xmax + log(max(value, tiny_prob))
   end function log_sum_exp

   real(dp) function regularized_gamma_p(a, x) result(p)
      real(dp), intent(in) :: a, x
      integer, parameter :: max_iter = 10000
      real(dp), parameter :: eps = 5.0e-15_dp
      real(dp), parameter :: fpmin = 1.0e-300_dp
      real(dp) :: ap, del, sum, b, c, d, h, an, q
      integer :: n

      if (a <= 0.0_dp .or. x < 0.0_dp) then
         p = quiet_nan()
         return
      end if
      if (x <= 0.0_dp) then
         p = 0.0_dp
         return
      end if

      if (x < a + 1.0_dp) then
         ap = a
         sum = 1.0_dp / a
         del = sum
         do n = 1, max_iter
            ap = ap + 1.0_dp
            del = del * x / ap
            sum = sum + del
            if (abs(del) <= abs(sum) * eps) exit
         end do
         p = sum * exp(-x + a * log(x) - log_gamma(a))
      else
         b = x + 1.0_dp - a
         c = 1.0_dp / fpmin
         d = 1.0_dp / max(abs(b), fpmin)
         if (b < 0.0_dp) d = -d
         h = d
         do n = 1, max_iter
            an = -real(n, dp) * (real(n, dp) - a)
            b = b + 2.0_dp
            d = an * d + b
            if (abs(d) < fpmin) d = sign(fpmin, d)
            c = b + an / c
            if (abs(c) < fpmin) c = sign(fpmin, c)
            d = 1.0_dp / d
            del = d * c
            h = h * del
            if (abs(del - 1.0_dp) <= eps) exit
         end do
         q = exp(-x + a * log(x) - log_gamma(a)) * h
         p = 1.0_dp - q
      end if
      p = min(1.0_dp, max(0.0_dp, p))
   end function regularized_gamma_p

   real(dp) function normal_quantile(p) result(x)
      real(dp), intent(in) :: p
      real(dp), parameter :: a1 = -3.969683028665376e+01_dp
      real(dp), parameter :: a2 =  2.209460984245205e+02_dp
      real(dp), parameter :: a3 = -2.759285104469687e+02_dp
      real(dp), parameter :: a4 =  1.383577518672690e+02_dp
      real(dp), parameter :: a5 = -3.066479806614716e+01_dp
      real(dp), parameter :: a6 =  2.506628277459239e+00_dp
      real(dp), parameter :: b1 = -5.447609879822406e+01_dp
      real(dp), parameter :: b2 =  1.615858368580409e+02_dp
      real(dp), parameter :: b3 = -1.556989798598866e+02_dp
      real(dp), parameter :: b4 =  6.680131188771972e+01_dp
      real(dp), parameter :: b5 = -1.328068155288572e+01_dp
      real(dp), parameter :: c1 = -7.784894002430293e-03_dp
      real(dp), parameter :: c2 = -3.223964580411365e-01_dp
      real(dp), parameter :: c3 = -2.400758277161838e+00_dp
      real(dp), parameter :: c4 = -2.549732539343734e+00_dp
      real(dp), parameter :: c5 =  4.374664141464968e+00_dp
      real(dp), parameter :: c6 =  2.938163982698783e+00_dp
      real(dp), parameter :: d1 =  7.784695709041462e-03_dp
      real(dp), parameter :: d2 =  3.224671290700398e-01_dp
      real(dp), parameter :: d3 =  2.445134137142996e+00_dp
      real(dp), parameter :: d4 =  3.754408661907416e+00_dp
      real(dp), parameter :: plow = 0.02425_dp
      real(dp), parameter :: phigh = 1.0_dp - plow
      real(dp) :: q, r

      if (p <= 0.0_dp) then
         x = -huge(1.0_dp)
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
      else if (p < plow) then
         q = sqrt(-2.0_dp * log(p))
         x = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / &
             ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
      else if (p <= phigh) then
         q = p - 0.5_dp
         r = q * q
         x = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6) * q / &
             (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0_dp)
      else
         q = sqrt(-2.0_dp * log(1.0_dp-p))
         x = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / &
              ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
      end if
   end function normal_quantile

   subroutine seed_random(seed)
      integer, intent(in), optional :: seed
      integer, allocatable :: values(:)
      integer :: n, i, base

      call random_seed(size=n)
      allocate(values(n))
      if (present(seed)) then
         base = seed
      else
         call system_clock(count=base)
      end if
      do i = 1, n
         values(i) = modulo(base + 104729*i + 8191*i*i, huge(1)-1)
         if (values(i) == 0) values(i) = i
      end do
      call random_seed(put=values)
   end subroutine seed_random

   real(dp) function standard_normal_random() result(z)
      real(dp) :: u1, u2
      call random_number(u1)
      call random_number(u2)
      u1 = max(u1, tiny_prob)
      z = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp*pi*u2)
   end function standard_normal_random

   recursive real(dp) function gamma_random(shape) result(g)
      real(dp), intent(in) :: shape
      real(dp) :: d, c, x, v, u

      if (shape <= 0.0_dp) then
         g = quiet_nan()
         return
      end if
      if (shape < 1.0_dp) then
         call random_number(u)
         g = gamma_random(shape + 1.0_dp) * max(u, tiny_prob)**(1.0_dp/shape)
         return
      end if

      d = shape - 1.0_dp/3.0_dp
      c = 1.0_dp / sqrt(9.0_dp*d)
      do
         do
            x = standard_normal_random()
            v = 1.0_dp + c*x
            if (v > 0.0_dp) exit
         end do
         v = v*v*v
         call random_number(u)
         if (u < 1.0_dp - 0.0331_dp*x**4) exit
         if (log(max(u, tiny_prob)) < 0.5_dp*x*x + d*(1.0_dp-v+log(v))) exit
      end do
      g = d*v
   end function gamma_random

   integer function categorical_random(probabilities) result(category)
      real(dp), intent(in) :: probabilities(:)
      real(dp) :: u, cumulative, total
      integer :: i

      total = sum(max(probabilities, 0.0_dp))
      if (total <= 0.0_dp) then
         category = 1
         return
      end if
      call random_number(u)
      u = u * total
      cumulative = 0.0_dp
      do i = 1, size(probabilities)
         cumulative = cumulative + max(probabilities(i), 0.0_dp)
         if (u <= cumulative) then
            category = i
            return
         end if
      end do
      category = size(probabilities)
   end function categorical_random

   integer function finite_count(x) result(n)
      real(dp), intent(in) :: x(:)
      integer :: i
      n = 0
      do i = 1, size(x)
         if (ieee_is_finite(x(i))) n = n + 1
      end do
   end function finite_count

   real(dp) function vector_mean(x) result(value)
      real(dp), intent(in) :: x(:)
      integer :: i, n
      value = 0.0_dp
      n = 0
      do i = 1, size(x)
         if (.not. ieee_is_finite(x(i))) cycle
         value = value + x(i)
         n = n + 1
      end do
      if (n == 0) then
         value = quiet_nan()
      else
         value = value / real(n, dp)
      end if
   end function vector_mean

   real(dp) function finite_mean_abs(x) result(value)
      real(dp), intent(in) :: x(:)
      integer :: i, n
      value = 0.0_dp
      n = 0
      do i = 1, size(x)
         if (.not. ieee_is_finite(x(i))) cycle
         value = value + abs(x(i))
         n = n + 1
      end do
      if (n == 0) then
         value = quiet_nan()
      else
         value = value / real(n, dp)
      end if
   end function finite_mean_abs

   real(dp) function sample_sd(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: mean_value, total
      integer :: i, n

      mean_value = vector_mean(x)
      if (.not. ieee_is_finite(mean_value)) then
         value = quiet_nan()
         return
      end if
      total = 0.0_dp
      n = 0
      do i = 1, size(x)
         if (.not. ieee_is_finite(x(i))) cycle
         total = total + (x(i)-mean_value)**2
         n = n + 1
      end do
      if (n < 2) then
         value = quiet_nan()
      else
         value = sqrt(total / real(n-1, dp))
      end if
   end function sample_sd

   real(dp) function population_skewness(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: mean_value, m2, m3, d
      integer :: i, n

      mean_value = vector_mean(x)
      m2 = 0.0_dp
      m3 = 0.0_dp
      n = 0
      do i = 1, size(x)
         if (.not. ieee_is_finite(x(i))) cycle
         d = x(i) - mean_value
         m2 = m2 + d*d
         m3 = m3 + d*d*d
         n = n + 1
      end do
      if (n == 0 .or. m2 <= 0.0_dp) then
         value = quiet_nan()
      else
         m2 = m2 / real(n, dp)
         m3 = m3 / real(n, dp)
         value = m3 / m2**1.5_dp
      end if
   end function population_skewness

   real(dp) function population_kurtosis(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: mean_value, m2, m4, d
      integer :: i, n

      mean_value = vector_mean(x)
      m2 = 0.0_dp
      m4 = 0.0_dp
      n = 0
      do i = 1, size(x)
         if (.not. ieee_is_finite(x(i))) cycle
         d = x(i) - mean_value
         m2 = m2 + d*d
         m4 = m4 + d**4
         n = n + 1
      end do
      if (n == 0 .or. m2 <= 0.0_dp) then
         value = quiet_nan()
      else
         m2 = m2 / real(n, dp)
         m4 = m4 / real(n, dp)
         value = m4 / (m2*m2)
      end if
   end function population_kurtosis

   real(dp) function correlation(x, y) result(value)
      real(dp), intent(in) :: x(:), y(:)
      real(dp) :: mx, my, sx, sy, sxy, dx, dy
      integer :: i, n

      if (size(x) /= size(y)) then
         value = quiet_nan()
         return
      end if
      mx = 0.0_dp
      my = 0.0_dp
      n = 0
      do i = 1, size(x)
         if (.not. ieee_is_finite(x(i)) .or. .not. ieee_is_finite(y(i))) cycle
         mx = mx + x(i)
         my = my + y(i)
         n = n + 1
      end do
      if (n < 2) then
         value = quiet_nan()
         return
      end if
      mx = mx / real(n, dp)
      my = my / real(n, dp)
      sx = 0.0_dp
      sy = 0.0_dp
      sxy = 0.0_dp
      do i = 1, size(x)
         if (.not. ieee_is_finite(x(i)) .or. .not. ieee_is_finite(y(i))) cycle
         dx = x(i) - mx
         dy = y(i) - my
         sx = sx + dx*dx
         sy = sy + dy*dy
         sxy = sxy + dx*dy
      end do
      if (sx <= 0.0_dp .or. sy <= 0.0_dp) then
         value = quiet_nan()
      else
         value = sxy / sqrt(sx*sy)
      end if
   end function correlation

   function absolute_acf(x, lag_max, drop) result(acf)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: lag_max
      integer, intent(in), optional :: drop
      real(dp), allocatable :: acf(:)
      real(dp), allocatable :: y(:), centered(:)
      real(dp) :: mean_value, denominator
      integer :: k, n, remove_n

      remove_n = 0
      if (present(drop)) remove_n = max(0, drop)
      y = remove_outliers(abs(x), remove_n)
      n = size(y)
      allocate(acf(max(0, lag_max)))
      if (lag_max <= 0) return
      if (n < 2) then
         acf = quiet_nan()
         return
      end if
      mean_value = vector_mean(y)
      allocate(centered(n))
      centered = y - mean_value
      denominator = sum(centered*centered)
      if (denominator <= 0.0_dp) then
         acf = quiet_nan()
         return
      end if
      do k = 1, lag_max
         if (k >= n) then
            acf(k) = quiet_nan()
         else
            acf(k) = sum(centered(1:n-k)*centered(1+k:n)) / denominator
         end if
      end do
   end function absolute_acf

   function simple_moving_average(x, order, na_backfill) result(ma)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: order
      logical, intent(in), optional :: na_backfill
      real(dp), allocatable :: ma(:), work(:)
      logical :: backfill
      integer :: i, j, n_used
      real(dp) :: total

      backfill = .true.
      if (present(na_backfill)) backfill = na_backfill
      allocate(ma(size(x)), work(size(x)))
      work = x
      if (size(x) == 0) return
      if (order <= 0) then
         ma = x
         return
      end if
      if (backfill) then
         do i = 2, size(work)
            if (.not. ieee_is_finite(work(i))) work(i) = work(i-1)
         end do
      end if
      do i = 1, size(work)
         total = 0.0_dp
         n_used = 0
         do j = max(1, i-order+1), i
            if (.not. ieee_is_finite(work(j))) cycle
            total = total + work(j)
            n_used = n_used + 1
         end do
         if (n_used == 0) then
            ma(i) = quiet_nan()
         else
            ma(i) = total / real(n_used, dp)
         end if
      end do
   end function simple_moving_average

   function remove_outliers(x, drop) result(y)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: drop
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: abs_values(:), sorted(:)
      logical, allocatable :: keep(:)
      real(dp) :: cutoff
      integer :: i, j, n_keep, d

      if (drop <= 0 .or. size(x) == 0) then
         y = x
         return
      end if
      allocate(abs_values(size(x)), sorted(size(x)), keep(size(x)))
      abs_values = abs(x)
      sorted = abs_values
      do i = 2, size(sorted)
         cutoff = sorted(i)
         j = i - 1
         do while (j >= 1)
            if (sorted(j) >= cutoff) exit
            sorted(j+1) = sorted(j)
            j = j - 1
         end do
         sorted(j+1) = cutoff
      end do
      d = min(drop, size(sorted))
      cutoff = sorted(d)
      keep = abs_values < cutoff
      n_keep = count(keep)
      allocate(y(n_keep))
      j = 0
      do i = 1, size(x)
         if (.not. keep(i)) cycle
         j = j + 1
         y(j) = x(i)
      end do
   end function remove_outliers

   subroutine solve_linear_system(a, b, x, status)
      real(dp), intent(in) :: a(:, :), b(:)
      real(dp), allocatable, intent(out) :: x(:)
      integer, intent(out) :: status
      real(dp), allocatable :: aug(:, :), row_temp(:)
      real(dp) :: factor, pivot_value
      integer :: i, k, pivot, n

      n = size(b)
      status = 0
      if (size(a,1) /= n .or. size(a,2) /= n) then
         allocate(x(0))
         status = 1
         return
      end if
      allocate(aug(n,n+1), row_temp(n+1), x(n))
      aug(:,1:n) = a
      aug(:,n+1) = b
      do k = 1, n
         pivot = k - 1 + maxloc(abs(aug(k:n,k)), dim=1)
         pivot_value = aug(pivot,k)
         if (abs(pivot_value) <= 100.0_dp*epsilon(1.0_dp)) then
            x = quiet_nan()
            status = 2
            return
         end if
         if (pivot /= k) then
            row_temp = aug(k,:)
            aug(k,:) = aug(pivot,:)
            aug(pivot,:) = row_temp
         end if
         aug(k,:) = aug(k,:) / aug(k,k)
         do i = 1, n
            if (i == k) cycle
            factor = aug(i,k)
            aug(i,:) = aug(i,:) - factor*aug(k,:)
         end do
      end do
      x = aug(:,n+1)
   end subroutine solve_linear_system

end module ldhmm_math
