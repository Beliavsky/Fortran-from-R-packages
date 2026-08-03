! SPDX-License-Identifier: GPL-2.0-or-later
module icsnp_special
   use icsnp_kinds, only : dp
   implicit none
   private
   public :: chi_square_survival, f_survival, normal_quantile
   public :: chi_square_quantile, rank_average, median_value

contains

   pure real(dp) function chi_square_survival(x, degrees_freedom) result(p)
      real(dp), intent(in) :: x, degrees_freedom
      if (x <= 0.0_dp) then
         p = 1.0_dp
      else if (degrees_freedom <= 0.0_dp) then
         p = 0.0_dp
      else
         p = regularized_gamma_q(0.5_dp * degrees_freedom, 0.5_dp * x)
      end if
   end function chi_square_survival

   pure real(dp) function f_survival(x, df1, df2) result(p)
      real(dp), intent(in) :: x, df1, df2
      real(dp) :: z
      if (x <= 0.0_dp) then
         p = 1.0_dp
      else if (df1 <= 0.0_dp .or. df2 <= 0.0_dp) then
         p = 0.0_dp
      else
         z = df2 / (df2 + df1 * x)
         p = regularized_beta(z, 0.5_dp * df2, 0.5_dp * df1)
      end if
      p = min(1.0_dp, max(0.0_dp, p))
   end function f_survival

   pure real(dp) function normal_quantile(probability) result(x)
      real(dp), intent(in) :: probability
      real(dp), parameter :: a1 = -3.969683028665376e1_dp
      real(dp), parameter :: a2 =  2.209460984245205e2_dp
      real(dp), parameter :: a3 = -2.759285104469687e2_dp
      real(dp), parameter :: a4 =  1.383577518672690e2_dp
      real(dp), parameter :: a5 = -3.066479806614716e1_dp
      real(dp), parameter :: a6 =  2.506628277459239_dp
      real(dp), parameter :: b1 = -5.447609879822406e1_dp
      real(dp), parameter :: b2 =  1.615858368580409e2_dp
      real(dp), parameter :: b3 = -1.556989798598866e2_dp
      real(dp), parameter :: b4 =  6.680131188771972e1_dp
      real(dp), parameter :: b5 = -1.328068155288572e1_dp
      real(dp), parameter :: c1 = -7.784894002430293e-3_dp
      real(dp), parameter :: c2 = -3.223964580411365e-1_dp
      real(dp), parameter :: c3 = -2.400758277161838_dp
      real(dp), parameter :: c4 = -2.549732539343734_dp
      real(dp), parameter :: c5 =  4.374664141464968_dp
      real(dp), parameter :: c6 =  2.938163982698783_dp
      real(dp), parameter :: d1 =  7.784695709041462e-3_dp
      real(dp), parameter :: d2 =  3.224671290700398e-1_dp
      real(dp), parameter :: d3 =  2.445134137142996_dp
      real(dp), parameter :: d4 =  3.754408661907416_dp
      real(dp), parameter :: p_low = 0.02425_dp
      real(dp), parameter :: p_high = 1.0_dp - p_low
      real(dp) :: p, q, r

      p = min(1.0_dp - epsilon(1.0_dp), max(epsilon(1.0_dp), probability))
      if (p < p_low) then
         q = sqrt(-2.0_dp * log(p))
         x = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / &
             ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
      else if (p <= p_high) then
         q = p - 0.5_dp
         r = q * q
         x = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6) * q / &
             (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0_dp)
      else
         q = sqrt(-2.0_dp * log(1.0_dp - p))
         x = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / &
              ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
      end if
   end function normal_quantile

   pure real(dp) function chi_square_quantile(probability, degrees_freedom) result(x)
      real(dp), intent(in) :: probability, degrees_freedom
      real(dp) :: lower, upper, middle, cdf, target, z
      integer :: i
      target = min(1.0_dp - 1.0e-14_dp, max(1.0e-14_dp, probability))
      if (degrees_freedom <= 0.0_dp) then
         x = 0.0_dp
         return
      end if
      z = normal_quantile(target)
      x = degrees_freedom * max(0.05_dp, 1.0_dp - 2.0_dp/(9.0_dp*degrees_freedom) + &
         z * sqrt(2.0_dp/(9.0_dp*degrees_freedom)))**3
      lower = 0.0_dp
      upper = max(1.0_dp, 2.0_dp * x)
      do while (1.0_dp - chi_square_survival(upper, degrees_freedom) < target)
         upper = 2.0_dp * upper
         if (upper > huge(1.0_dp) / 4.0_dp) exit
      end do
      do i = 1, 100
         middle = 0.5_dp * (lower + upper)
         cdf = 1.0_dp - chi_square_survival(middle, degrees_freedom)
         if (cdf < target) then
            lower = middle
         else
            upper = middle
         end if
      end do
      x = 0.5_dp * (lower + upper)
   end function chi_square_quantile

   subroutine rank_average(x, ranks)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: ranks(:)
      integer, allocatable :: order(:)
      integer :: n, i, j, key, left, right
      real(dp) :: average_rank, scale
      n = size(x)
      allocate(ranks(n), order(n))
      do i = 1, n
         order(i) = i
      end do
      do i = 2, n
         key = order(i)
         j = i - 1
         do while (j >= 1)
            if (x(order(j)) <= x(key)) exit
            order(j + 1) = order(j)
            j = j - 1
         end do
         order(j + 1) = key
      end do
      scale = max(1.0_dp, maxval(abs(x)))
      left = 1
      do while (left <= n)
         right = left
         do while (right < n)
            if (abs(x(order(right + 1)) - x(order(left))) > &
               16.0_dp * epsilon(1.0_dp) * scale) exit
            right = right + 1
         end do
         average_rank = 0.5_dp * real(left + right, dp)
         do i = left, right
            ranks(order(i)) = average_rank
         end do
         left = right + 1
      end do
   end subroutine rank_average

   real(dp) function median_value(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable :: work(:)
      integer :: n, i, j
      real(dp) :: key
      n = size(x)
      if (n == 0) then
         value = 0.0_dp
         return
      end if
      allocate(work(n))
      work = x
      do i = 2, n
         key = work(i)
         j = i - 1
         do while (j >= 1)
            if (work(j) <= key) exit
            work(j + 1) = work(j)
            j = j - 1
         end do
         work(j + 1) = key
      end do
      if (mod(n, 2) == 1) then
         value = work((n + 1) / 2)
      else
         value = 0.5_dp * (work(n / 2) + work(n / 2 + 1))
      end if
   end function median_value

   pure real(dp) function regularized_gamma_q(a, x) result(q)
      real(dp), intent(in) :: a, x
      real(dp) :: gln
      if (a <= 0.0_dp .or. x < 0.0_dp) then
         q = 0.0_dp
      else if (x <= tiny(1.0_dp)) then
         q = 1.0_dp
      else if (x < a + 1.0_dp) then
         gln = log_gamma(a)
         q = 1.0_dp - gamma_series(a, x, gln)
      else
         gln = log_gamma(a)
         q = gamma_continued_fraction(a, x, gln)
      end if
      q = min(1.0_dp, max(0.0_dp, q))
   end function regularized_gamma_q

   pure real(dp) function gamma_series(a, x, gln) result(p)
      real(dp), intent(in) :: a, x, gln
      integer, parameter :: itmax = 10000
      real(dp), parameter :: eps = 5.0e-15_dp
      real(dp) :: ap, del, sum_value
      integer :: n
      ap = a
      sum_value = 1.0_dp / a
      del = sum_value
      do n = 1, itmax
         ap = ap + 1.0_dp
         del = del * x / ap
         sum_value = sum_value + del
         if (abs(del) <= abs(sum_value) * eps) exit
      end do
      p = sum_value * exp(-x + a * log(x) - gln)
   end function gamma_series

   pure real(dp) function gamma_continued_fraction(a, x, gln) result(q)
      real(dp), intent(in) :: a, x, gln
      integer, parameter :: itmax = 10000
      real(dp), parameter :: eps = 5.0e-15_dp
      real(dp), parameter :: fpmin = tiny(1.0_dp) / eps
      real(dp) :: an, b, c, d, del, h
      integer :: i
      b = x + 1.0_dp - a
      c = 1.0_dp / fpmin
      d = 1.0_dp / max(abs(b), fpmin)
      if (b < 0.0_dp) d = -d
      h = d
      do i = 1, itmax
         an = -real(i, dp) * (real(i, dp) - a)
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
      q = exp(-x + a * log(x) - gln) * h
   end function gamma_continued_fraction

   pure real(dp) function regularized_beta(x, a, b) result(value)
      real(dp), intent(in) :: x, a, b
      real(dp) :: front
      if (x <= 0.0_dp) then
         value = 0.0_dp
         return
      else if (x >= 1.0_dp) then
         value = 1.0_dp
         return
      end if
      front = exp(log_gamma(a + b) - log_gamma(a) - log_gamma(b) + &
         a * log(x) + b * log(1.0_dp - x))
      if (x < (a + 1.0_dp) / (a + b + 2.0_dp)) then
         value = front * beta_fraction(x, a, b) / a
      else
         value = 1.0_dp - front * beta_fraction(1.0_dp - x, b, a) / b
      end if
      value = min(1.0_dp, max(0.0_dp, value))
   end function regularized_beta

   pure real(dp) function beta_fraction(x, a, b) result(h)
      real(dp), intent(in) :: x, a, b
      integer, parameter :: maxiter = 10000
      real(dp), parameter :: eps = 5.0e-15_dp
      real(dp), parameter :: fpmin = tiny(1.0_dp) / eps
      real(dp) :: qab, qap, qam, c, d, del, aa
      integer :: m, m2
      qab = a + b
      qap = a + 1.0_dp
      qam = a - 1.0_dp
      c = 1.0_dp
      d = 1.0_dp - qab * x / qap
      if (abs(d) < fpmin) d = sign(fpmin, d)
      d = 1.0_dp / d
      h = d
      do m = 1, maxiter
         m2 = 2 * m
         aa = real(m, dp) * (b - real(m, dp)) * x / &
              ((qam + real(m2, dp)) * (a + real(m2, dp)))
         d = 1.0_dp + aa * d
         if (abs(d) < fpmin) d = sign(fpmin, d)
         c = 1.0_dp + aa / c
         if (abs(c) < fpmin) c = sign(fpmin, c)
         d = 1.0_dp / d
         h = h * d * c
         aa = -(a + real(m, dp)) * (qab + real(m, dp)) * x / &
              ((a + real(m2, dp)) * (qap + real(m2, dp)))
         d = 1.0_dp + aa * d
         if (abs(d) < fpmin) d = sign(fpmin, d)
         c = 1.0_dp + aa / c
         if (abs(c) < fpmin) c = sign(fpmin, c)
         d = 1.0_dp / d
         del = d * c
         h = h * del
         if (abs(del - 1.0_dp) <= eps) exit
      end do
   end function beta_fraction

end module icsnp_special
