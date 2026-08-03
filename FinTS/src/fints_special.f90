! SPDX-License-Identifier: GPL-2.0-or-later
module fints_special
   use fints_kinds, only : dp
   implicit none
   private
   public :: chi_square_survival, rank_average

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

end module fints_special
