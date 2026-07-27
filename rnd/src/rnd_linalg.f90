! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from RND 1.2, Copyright (C) 2017 Kam Hamidieh.
module rnd_linalg
   use rnd_kinds, only : dp
   implicit none
   private
   public :: solve_linear_system, quadratic_least_squares, simple_linear_regression

contains

   subroutine solve_linear_system(a, b, x, info)
      real(dp), intent(in) :: a(:, :), b(:)
      real(dp), intent(out) :: x(:)
      integer, intent(out) :: info
      real(dp), allocatable :: work(:, :), rhs(:), row(:)
      real(dp) :: factor, tmp
      integer :: i, k, n, pivot
      n = size(b)
      info = 0
      if (size(a,1) /= n .or. size(a,2) /= n .or. size(x) /= n) then
         info = -1
         x = 0.0_dp
         return
      end if
      allocate(work(n,n), rhs(n), row(n))
      work = a
      rhs = b
      do k = 1, n-1
         pivot = k-1+maxloc(abs(work(k:n,k)),dim=1)
         if (abs(work(pivot,k)) <= epsilon(1.0_dp)*max(1.0_dp,maxval(abs(work)))) then
            info = k
            x = 0.0_dp
            return
         end if
         if (pivot /= k) then
            row = work(k,:)
            work(k,:) = work(pivot,:)
            work(pivot,:) = row
            tmp = rhs(k)
            rhs(k) = rhs(pivot)
            rhs(pivot) = tmp
         end if
         do i = k+1, n
            factor = work(i,k)/work(k,k)
            work(i,k:n) = work(i,k:n)-factor*work(k,k:n)
            rhs(i) = rhs(i)-factor*rhs(k)
         end do
      end do
      if (abs(work(n,n)) <= epsilon(1.0_dp)*max(1.0_dp,maxval(abs(work)))) then
         info = n
         x = 0.0_dp
         return
      end if
      x(n) = rhs(n)/work(n,n)
      do i = n-1, 1, -1
         x(i) = (rhs(i)-sum(work(i,i+1:n)*x(i+1:n)))/work(i,i)
      end do
   end subroutine solve_linear_system

   subroutine quadratic_least_squares(k, volatility, coefficients, info)
      real(dp), intent(in) :: k(:), volatility(:)
      real(dp), intent(out) :: coefficients(3)
      integer, intent(out) :: info
      real(dp) :: normal(3,3), rhs(3)
      normal(1,1) = real(size(k),dp)
      normal(1,2) = sum(k)
      normal(1,3) = sum(k*k)
      normal(2,1) = normal(1,2)
      normal(2,2) = normal(1,3)
      normal(2,3) = sum(k**3)
      normal(3,1) = normal(1,3)
      normal(3,2) = normal(2,3)
      normal(3,3) = sum(k**4)
      rhs = [sum(volatility),sum(k*volatility),sum(k*k*volatility)]
      call solve_linear_system(normal,rhs,coefficients,info)
   end subroutine quadratic_least_squares

   subroutine simple_linear_regression(x, y, intercept, slope, info)
      real(dp), intent(in) :: x(:), y(:)
      real(dp), intent(out) :: intercept, slope
      integer, intent(out) :: info
      real(dp) :: xmean, ymean, denominator
      if (size(x) /= size(y) .or. size(x) < 2) then
         info = -1
         intercept = 0.0_dp
         slope = 0.0_dp
         return
      end if
      xmean = sum(x)/real(size(x),dp)
      ymean = sum(y)/real(size(y),dp)
      denominator = sum((x-xmean)**2)
      if (denominator <= epsilon(1.0_dp)) then
         info = 1
         intercept = 0.0_dp
         slope = 0.0_dp
         return
      end if
      slope = sum((x-xmean)*(y-ymean))/denominator
      intercept = ymean-slope*xmean
      info = 0
   end subroutine simple_linear_regression

end module rnd_linalg
