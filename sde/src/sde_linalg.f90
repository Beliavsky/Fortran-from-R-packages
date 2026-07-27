! This file is part of sde-fortran, a translation/adaptation of the R package
! sde 2.0.21 by Stefano Maria Iacus.
! Original package code Copyright (C) 2006 S. M. Iacus.
! SPDX-License-Identifier: GPL-2.0-or-later
module sde_linalg
   use sde_kinds, only : dp
   implicit none
   private

   public :: solve_linear_system
   public :: invert_matrix
   public :: quadratic_form
   public :: numerical_hessian

contains

   subroutine solve_linear_system(a, b, x, status)
      real(dp), intent(in) :: a(:, :)
      real(dp), intent(in) :: b(:)
      real(dp), intent(out) :: x(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: work(:, :), rhs(:), temp_row(:)
      real(dp) :: pivot_value, factor, temp_value
      integer :: n, i, j, pivot

      n = size(b)
      if (size(a, 1) /= n .or. size(a, 2) /= n .or. size(x) /= n) then
         error stop "solve_linear_system: incompatible dimensions"
      end if
      allocate(work(n, n), rhs(n), temp_row(n))
      work = a
      rhs = b
      if (present(status)) status = 0

      do i = 1, n
         pivot = i-1+maxloc(abs(work(i:n, i)), dim=1)
         pivot_value = work(pivot, i)
         if (abs(pivot_value) <= epsilon(1.0_dp)*max(1.0_dp, maxval(abs(work)))) then
            x = 0.0_dp
            if (present(status)) status = 1
            return
         end if
         if (pivot /= i) then
            temp_row = work(i, :)
            work(i, :) = work(pivot, :)
            work(pivot, :) = temp_row
            temp_value = rhs(i)
            rhs(i) = rhs(pivot)
            rhs(pivot) = temp_value
         end if
         do j = i+1, n
            factor = work(j, i)/work(i, i)
            work(j, i:n) = work(j, i:n)-factor*work(i, i:n)
            rhs(j) = rhs(j)-factor*rhs(i)
         end do
      end do

      x(n) = rhs(n)/work(n, n)
      do i = n-1, 1, -1
         x(i) = (rhs(i)-dot_product(work(i, i+1:n), x(i+1:n)))/work(i, i)
      end do
   end subroutine solve_linear_system

   subroutine invert_matrix(a, inverse, status, ridge)
      real(dp), intent(in) :: a(:, :)
      real(dp), intent(out) :: inverse(:, :)
      integer, intent(out), optional :: status
      real(dp), intent(in), optional :: ridge
      real(dp), allocatable :: adjusted(:, :), unit(:), column(:)
      real(dp) :: ridge_value
      integer :: n, i, local_status

      n = size(a, 1)
      if (size(a, 2) /= n .or. size(inverse, 1) /= n .or. size(inverse, 2) /= n) then
         error stop "invert_matrix: incompatible dimensions"
      end if
      ridge_value = 0.0_dp
      if (present(ridge)) ridge_value = ridge
      allocate(adjusted(n, n), unit(n), column(n))
      adjusted = a
      do i = 1, n
         adjusted(i, i) = adjusted(i, i)+ridge_value
      end do
      if (present(status)) status = 0
      do i = 1, n
         unit = 0.0_dp
         unit(i) = 1.0_dp
         call solve_linear_system(adjusted, unit, column, local_status)
         if (local_status /= 0) then
            inverse = 0.0_dp
            if (present(status)) status = local_status
            return
         end if
         inverse(:, i) = column
      end do
   end subroutine invert_matrix

   pure function quadratic_form(x, a) result(value)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: a(:, :)
      real(dp) :: value

      if (size(a, 1) /= size(x) .or. size(a, 2) /= size(x)) then
         error stop "quadratic_form: incompatible dimensions"
      end if
      value = dot_product(x, matmul(a, x))
   end function quadratic_form

   subroutine numerical_hessian(objective, x, hessian, relative_step)
      interface
         function objective(z) result(value)
            import dp
            real(dp), intent(in) :: z(:)
            real(dp) :: value
         end function objective
      end interface
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: hessian(:, :)
      real(dp), intent(in), optional :: relative_step
      real(dp), allocatable :: xp(:), xm(:), xpp(:), xpm(:), xmp(:), xmm(:)
      real(dp) :: step_scale, hi, hj, f0
      integer :: n, i, j

      n = size(x)
      if (size(hessian, 1) /= n .or. size(hessian, 2) /= n) then
         error stop "numerical_hessian: incompatible dimensions"
      end if
      step_scale = epsilon(1.0_dp)**0.25_dp
      if (present(relative_step)) step_scale = relative_step
      allocate(xp(n), xm(n), xpp(n), xpm(n), xmp(n), xmm(n))
      f0 = objective(x)
      do i = 1, n
         hi = step_scale*max(1.0_dp, abs(x(i)))
         xp = x
         xm = x
         xp(i) = xp(i)+hi
         xm(i) = xm(i)-hi
         hessian(i, i) = (objective(xp)-2.0_dp*f0+objective(xm))/(hi*hi)
         do j = i+1, n
            hj = step_scale*max(1.0_dp, abs(x(j)))
            xpp = x
            xpm = x
            xmp = x
            xmm = x
            xpp(i) = xpp(i)+hi
            xpp(j) = xpp(j)+hj
            xpm(i) = xpm(i)+hi
            xpm(j) = xpm(j)-hj
            xmp(i) = xmp(i)-hi
            xmp(j) = xmp(j)+hj
            xmm(i) = xmm(i)-hi
            xmm(j) = xmm(j)-hj
            hessian(i, j) = (objective(xpp)-objective(xpm)-objective(xmp)+objective(xmm))/(4.0_dp*hi*hj)
            hessian(j, i) = hessian(i, j)
         end do
      end do
   end subroutine numerical_hessian

end module sde_linalg
