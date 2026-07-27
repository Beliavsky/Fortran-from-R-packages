! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of garchx.
! Copyright (C) 2026 translation contributors.
! Original garchx package copyright (C) Genaro Sucarrat.
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 2 of the License, or
! (at your option) any later version.
module garchx_math
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use garchx_kinds, only : dp
   use garchx_linalg, only : cholesky_lower
   implicit none
   private
   public :: set_random_seed, random_normal, random_normal_vector
   public :: rmnorm, normal_cdf, normal_quantile, student_t_cdf
   public :: student_t_quantile, empirical_quantile, mean_value, variance_value
   public :: is_finite_vector, is_finite_matrix, sort_real
contains
   subroutine set_random_seed(seed)
      integer, intent(in) :: seed
      integer :: n, i
      integer, allocatable :: values(:)
      call random_seed(size=n)
      allocate(values(n))
      do i = 1, n
         values(i) = modulo(seed + 104729*i + 37*i*i, huge(1)-1)
         if (values(i) <= 0) values(i) = i + 17
      end do
      call random_seed(put=values)
   end subroutine set_random_seed

   real(dp) function random_normal() result(z)
      real(dp) :: u1, u2
      call random_number(u1)
      call random_number(u2)
      u1 = max(u1, tiny(1.0_dp))
      z = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*acos(-1.0_dp)*u2)
   end function random_normal

   subroutine random_normal_vector(z)
      real(dp), intent(out) :: z(:)
      integer :: i
      do i = 1, size(z)
         z(i) = random_normal()
      end do
   end subroutine random_normal_vector

   subroutine rmnorm(n, mean, vcov, draws, status)
      integer, intent(in) :: n
      real(dp), intent(in), optional :: mean(:)
      real(dp), intent(in) :: vcov(:, :)
      real(dp), allocatable, intent(out) :: draws(:, :)
      integer, intent(out) :: status
      integer :: d, i
      real(dp), allocatable :: l(:, :), z(:), mu(:)

      d = size(vcov, 1)
      if (size(vcov, 2) /= d .or. n < 0) then
         status = 1
         allocate(draws(0, 0))
         return
      end if
      allocate(mu(d))
      if (present(mean)) then
         if (size(mean) /= d) then
            status = 2
            allocate(draws(0, 0))
            return
         end if
         mu = mean
      else
         mu = 0.0_dp
      end if
      call cholesky_lower(vcov, l, status)
      if (status /= 0) then
         allocate(draws(0, 0))
         return
      end if
      allocate(draws(n, d), z(d))
      do i = 1, n
         call random_normal_vector(z)
         draws(i, :) = mu + matmul(l, z)
      end do
   end subroutine rmnorm

   pure real(dp) function normal_cdf(x) result(value)
      real(dp), intent(in) :: x
      value = 0.5_dp*erfc(-x/sqrt(2.0_dp))
   end function normal_cdf

   pure real(dp) function normal_quantile(p) result(x)
      real(dp), intent(in) :: p
      real(dp), parameter :: a(6) = [ &
         -3.969683028665376d1, 2.209460984245205d2, -2.759285104469687d2, &
          1.383577518672690d2, -3.066479806614716d1, 2.506628277459239d0]
      real(dp), parameter :: b(5) = [ &
         -5.447609879822406d1, 1.615858368580409d2, -1.556989798598866d2, &
          6.680131188771972d1, -1.328068155288572d1]
      real(dp), parameter :: c(6) = [ &
         -7.784894002430293d-3, -3.223964580411365d-1, -2.400758277161838d0, &
         -2.549732539343734d0, 4.374664141464968d0, 2.938163982698783d0]
      real(dp), parameter :: d(4) = [ &
          7.784695709041462d-3, 3.224671290700398d-1, 2.445134137142996d0, &
          3.754408661907416d0]
      real(dp), parameter :: plow = 0.02425_dp, phigh = 1.0_dp-plow
      real(dp) :: q, r

      if (p <= 0.0_dp) then
         x = -huge(1.0_dp)
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
      else if (p < plow) then
         q = sqrt(-2.0_dp*log(p))
         x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
             ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
      else if (p <= phigh) then
         q = p-0.5_dp
         r = q*q
         x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q/ &
             (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
      else
         q = sqrt(-2.0_dp*log(1.0_dp-p))
         x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
              ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
      end if
   end function normal_quantile

   pure real(dp) function log_beta(a, b) result(value)
      real(dp), intent(in) :: a, b
      value = log_gamma(a) + log_gamma(b) - log_gamma(a+b)
   end function log_beta

   pure real(dp) function beta_cont_frac(a, b, x) result(value)
      real(dp), intent(in) :: a, b, x
      integer, parameter :: max_iter = 300
      real(dp), parameter :: eps = 3.0e-14_dp, fpmin = 1.0e-300_dp
      integer :: m, m2
      real(dp) :: aa, c, d, del, h, qab, qam, qap

      qab = a+b
      qap = a+1.0_dp
      qam = a-1.0_dp
      c = 1.0_dp
      d = 1.0_dp-qab*x/qap
      if (abs(d) < fpmin) d = fpmin
      d = 1.0_dp/d
      h = d
      do m = 1, max_iter
         m2 = 2*m
         aa = real(m, dp)*(b-real(m, dp))*x/((qam+real(m2, dp))*(a+real(m2, dp)))
         d = 1.0_dp+aa*d
         if (abs(d) < fpmin) d = fpmin
         c = 1.0_dp+aa/c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp/d
         h = h*d*c
         aa = -(a+real(m, dp))*(qab+real(m, dp))*x/ &
              ((a+real(m2, dp))*(qap+real(m2, dp)))
         d = 1.0_dp+aa*d
         if (abs(d) < fpmin) d = fpmin
         c = 1.0_dp+aa/c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp/d
         del = d*c
         h = h*del
         if (abs(del-1.0_dp) <= eps) exit
      end do
      value = h
   end function beta_cont_frac

   pure real(dp) function regularized_beta(x, a, b) result(value)
      real(dp), intent(in) :: x, a, b
      real(dp) :: bt
      if (x <= 0.0_dp) then
         value = 0.0_dp
      else if (x >= 1.0_dp) then
         value = 1.0_dp
      else
         bt = exp(a*log(x)+b*log(1.0_dp-x)-log_beta(a, b))
         if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
            value = bt*beta_cont_frac(a, b, x)/a
         else
            value = 1.0_dp-bt*beta_cont_frac(b, a, 1.0_dp-x)/b
         end if
      end if
   end function regularized_beta

   pure real(dp) function student_t_cdf(x, df) result(value)
      real(dp), intent(in) :: x, df
      real(dp) :: z, ib
      if (df <= 0.0_dp) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      if (abs(x) <= tiny(1.0_dp)) then
         value = 0.5_dp
         return
      end if
      z = df/(df+x*x)
      ib = regularized_beta(z, 0.5_dp*df, 0.5_dp)
      if (x > 0.0_dp) then
         value = 1.0_dp-0.5_dp*ib
      else
         value = 0.5_dp*ib
      end if
   end function student_t_cdf

   real(dp) function student_t_quantile(p, df) result(x)
      real(dp), intent(in) :: p, df
      real(dp) :: lo, hi, mid
      integer :: i
      if (p <= 0.0_dp) then
         x = -huge(1.0_dp)
         return
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
         return
      end if
      lo = -1.0_dp
      hi = 1.0_dp
      do while (student_t_cdf(lo, df) > p)
         lo = 2.0_dp*lo
      end do
      do while (student_t_cdf(hi, df) < p)
         hi = 2.0_dp*hi
      end do
      do i = 1, 120
         mid = 0.5_dp*(lo+hi)
         if (student_t_cdf(mid, df) < p) then
            lo = mid
         else
            hi = mid
         end if
      end do
      x = 0.5_dp*(lo+hi)
   end function student_t_quantile

   subroutine sort_real(x)
      real(dp), intent(inout) :: x(:)
      integer :: i, j
      real(dp) :: key
      do i = 2, size(x)
         key = x(i)
         j = i-1
         do while (j >= 1)
            if (x(j) <= key) exit
            x(j+1) = x(j)
            j = j-1
         end do
         x(j+1) = key
      end do
   end subroutine sort_real

   real(dp) function empirical_quantile(x, p) result(value)
      real(dp), intent(in) :: x(:), p
      real(dp), allocatable :: work(:)
      real(dp) :: h, frac
      integer :: j, n
      n = size(x)
      if (n == 0) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      allocate(work(n))
      work = x
      call sort_real(work)
      if (p <= 0.0_dp) then
         value = work(1)
      else if (p >= 1.0_dp) then
         value = work(n)
      else
         h = 1.0_dp + real(n-1, dp)*p
         j = floor(h)
         frac = h-real(j, dp)
         if (j >= n) then
            value = work(n)
         else
            value = (1.0_dp-frac)*work(j)+frac*work(j+1)
         end if
      end if
   end function empirical_quantile

   pure real(dp) function mean_value(x) result(value)
      real(dp), intent(in) :: x(:)
      if (size(x) > 0) then
         value = sum(x)/real(size(x), dp)
      else
         value = ieee_value(0.0_dp, ieee_quiet_nan)
      end if
   end function mean_value

   pure real(dp) function variance_value(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: mu
      if (size(x) > 1) then
         mu = mean_value(x)
         value = sum((x-mu)**2)/real(size(x)-1, dp)
      else
         value = 0.0_dp
      end if
   end function variance_value

   pure logical function is_finite_vector(x) result(ok)
      use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
      real(dp), intent(in) :: x(:)
      ok = all(ieee_is_finite(x))
   end function is_finite_vector

   pure logical function is_finite_matrix(x) result(ok)
      use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
      real(dp), intent(in) :: x(:, :)
      ok = all(ieee_is_finite(x))
   end function is_finite_matrix
end module garchx_math
