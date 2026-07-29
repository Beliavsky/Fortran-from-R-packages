! RiskPortfolios Fortran, derived from RiskPortfolios 2.1.7.
! Original code Copyright (C) 2013-2021 David Ardia.
! Original authors: David Ardia, Kris Boudt, Jean-Philippe Gagnon-Fleury.
! SPDX-License-Identifier: GPL-2.0-or-later
module riskportfolios_stats
   use riskportfolios_kinds, only : dp
   implicit none
   private
   public :: column_means, sample_covariance, population_covariance
   public :: standard_deviations, covariance_to_correlation
   public :: quantile_type7, median_value, frobenius_norm, symmetrize

contains

   pure function column_means(x) result(mu)
      real(dp), intent(in) :: x(:, :)
      real(dp) :: mu(size(x, 2))
      if (size(x, 1) > 0) then
         mu = sum(x, dim=1) / real(size(x, 1), dp)
      else
         mu = 0.0_dp
      end if
   end function column_means

   pure function sample_covariance(x) result(s)
      real(dp), intent(in) :: x(:, :)
      real(dp) :: s(size(x, 2), size(x, 2))
      real(dp) :: centered(size(x, 1), size(x, 2))
      real(dp) :: mu(size(x, 2))
      integer :: t

      t = size(x, 1)
      if (t <= 1) then
         s = 0.0_dp
         return
      end if
      mu = column_means(x)
      centered = x - spread(mu, 1, t)
      s = matmul(transpose(centered), centered) / real(t - 1, dp)
      s = symmetrize(s)
   end function sample_covariance

   pure function population_covariance(x) result(s)
      real(dp), intent(in) :: x(:, :)
      real(dp) :: s(size(x, 2), size(x, 2))
      real(dp) :: centered(size(x, 1), size(x, 2))
      real(dp) :: mu(size(x, 2))
      integer :: t

      t = size(x, 1)
      if (t <= 0) then
         s = 0.0_dp
         return
      end if
      mu = column_means(x)
      centered = x - spread(mu, 1, t)
      s = matmul(transpose(centered), centered) / real(t, dp)
      s = symmetrize(s)
   end function population_covariance

   pure function standard_deviations(x) result(sd)
      real(dp), intent(in) :: x(:, :)
      real(dp) :: sd(size(x, 2))
      real(dp) :: s(size(x, 2), size(x, 2))
      integer :: i
      s = sample_covariance(x)
      do i = 1, size(sd)
         sd(i) = sqrt(max(s(i, i), 0.0_dp))
      end do
   end function standard_deviations

   pure function covariance_to_correlation(sigma) result(rho)
      real(dp), intent(in) :: sigma(:, :)
      real(dp) :: rho(size(sigma, 1), size(sigma, 2))
      real(dp) :: sd(size(sigma, 1))
      integer :: i, j, n
      n = size(sigma, 1)
      sd = 0.0_dp
      do i = 1, n
         sd(i) = sqrt(max(sigma(i, i), 0.0_dp))
      end do
      do j = 1, n
         do i = 1, n
            if (sd(i) > 0.0_dp .and. sd(j) > 0.0_dp) then
               rho(i, j) = sigma(i, j) / (sd(i) * sd(j))
            else
               rho(i, j) = 0.0_dp
            end if
         end do
      end do
      do i = 1, n
         rho(i, i) = 1.0_dp
      end do
      rho = symmetrize(rho)
   end function covariance_to_correlation

   function quantile_type7(x, p) result(q)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: p
      real(dp) :: q
      real(dp), allocatable :: y(:)
      real(dp) :: h, g
      integer :: n, j
      n = size(x)
      if (n == 0) then
         q = 0.0_dp
         return
      end if
      allocate(y(n))
      y = x
      call sort_real(y)
      if (p <= 0.0_dp) then
         q = y(1)
      else if (p >= 1.0_dp) then
         q = y(n)
      else
         h = 1.0_dp + real(n - 1, dp) * p
         j = int(floor(h))
         g = h - real(j, dp)
         if (j >= n) then
            q = y(n)
         else
            q = (1.0_dp - g) * y(j) + g * y(j + 1)
         end if
      end if
   end function quantile_type7

   function median_value(x) result(m)
      real(dp), intent(in) :: x(:)
      real(dp) :: m
      m = quantile_type7(x, 0.5_dp)
   end function median_value

   pure function frobenius_norm(a) result(v)
      real(dp), intent(in) :: a(:, :)
      real(dp) :: v
      v = sqrt(sum(a * a))
   end function frobenius_norm

   pure function symmetrize(a) result(b)
      real(dp), intent(in) :: a(:, :)
      real(dp) :: b(size(a, 1), size(a, 2))
      b = 0.5_dp * (a + transpose(a))
   end function symmetrize

   subroutine sort_real(x)
      real(dp), intent(inout) :: x(:)
      if (size(x) > 1) call quicksort(x, 1, size(x))
   end subroutine sort_real

   recursive subroutine quicksort(x, left, right)
      real(dp), intent(inout) :: x(:)
      integer, intent(in) :: left, right
      integer :: i, j
      real(dp) :: pivot, tmp
      i = left
      j = right
      pivot = x((left + right) / 2)
      do
         do while (x(i) < pivot)
            i = i + 1
         end do
         do while (x(j) > pivot)
            j = j - 1
         end do
         if (i <= j) then
            tmp = x(i)
            x(i) = x(j)
            x(j) = tmp
            i = i + 1
            j = j - 1
         end if
         if (i > j) exit
      end do
      if (left < j) call quicksort(x, left, j)
      if (i < right) call quicksort(x, i, right)
   end subroutine quicksort

end module riskportfolios_stats
