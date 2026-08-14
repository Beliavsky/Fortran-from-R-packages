! Copyright (C) 1998 Douglas M. Bates and William N. Venables.
! Modern Fortran translation, 2026.
! SPDX-License-Identifier: GPL-2.0-or-later
module splines_basis
   use splines_kinds, only : dp
   use splines_core, only : spline_design, type7_quantile, sorted_copy
   use splines_linalg, only : nullspace_transform
   implicit none
   private
   public :: bs_basis, natural_spline_basis

contains

   subroutine bs_basis(x, basis, degree, knots, df, intercept, boundary_knots, &
                       interior_knots, status)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: basis(:, :)
      integer, intent(in), optional :: degree, df
      real(dp), intent(in), optional :: knots(:), boundary_knots(2)
      logical, intent(in), optional :: intercept
      real(dp), allocatable, intent(out), optional :: interior_knots(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: iknots(:), all_knots(:), raw(:, :), temp(:, :), pivot_x(:)
      real(dp) :: bounds(2), p, dx, factorial
      logical :: keep_intercept
      integer :: deg, ord, nik, i, j, ncoef, ierr
      integer, allocatable :: derivs(:)

      if (present(status)) status = 0
      deg = 3; if (present(degree)) deg = degree
      keep_intercept = .false.; if (present(intercept)) keep_intercept = intercept
      ord = deg + 1
      if (deg < 1 .or. size(x) == 0) then
         allocate(basis(0, 0)); if (present(status)) status = 1; return
      end if
      if (present(boundary_knots)) then
         bounds = boundary_knots
         if (bounds(2) < bounds(1)) bounds = bounds(2:1:-1)
      else
         bounds = [minval(x), maxval(x)]
      end if
      if (bounds(2) <= bounds(1)) then
         allocate(basis(0, 0)); if (present(status)) status = 2; return
      end if

      if (present(knots)) then
         allocate(iknots(size(knots))); iknots = sorted_copy(knots)
      else if (present(df)) then
         nik = df - ord + merge(0, 1, keep_intercept)
         nik = max(0, nik)
         allocate(iknots(nik))
         do i = 1, nik
            p = real(i, dp) / real(nik + 1, dp)
            iknots(i) = type7_quantile(pack(x, x >= bounds(1) .and. x <= bounds(2)), p)
         end do
      else
         allocate(iknots(0))
      end if
      if (present(interior_knots)) interior_knots = iknots

      allocate(all_knots(2 * ord + size(iknots)))
      all_knots(1:ord) = bounds(1)
      if (size(iknots) > 0) all_knots(ord + 1:ord + size(iknots)) = iknots
      all_knots(ord + size(iknots) + 1:) = bounds(2)
      ncoef = size(all_knots) - ord
      allocate(raw(size(x), ncoef)); raw = 0.0_dp

      do i = 1, size(x)
         if (x(i) >= bounds(1) .and. x(i) <= bounds(2)) then
            call spline_design(all_knots, [x(i)], ord, temp, status=ierr)
            if (ierr /= 0) then
               allocate(basis(0, 0)); if (present(status)) status = 3; return
            end if
            raw(i, :) = temp(1, :)
         else
            allocate(pivot_x(ord), derivs(ord))
            pivot_x = merge(bounds(1), bounds(2), x(i) < bounds(1))
            derivs = [(j - 1, j=1,ord)]
            call spline_design(all_knots, pivot_x, ord, temp, derivs, ierr)
            if (ierr /= 0) then
               allocate(basis(0, 0)); if (present(status)) status = 4; return
            end if
            dx = x(i) - pivot_x(1)
            factorial = 1.0_dp
            raw(i, :) = temp(1, :)
            do j = 2, ord
               factorial = factorial * real(j - 1, dp)
               raw(i, :) = raw(i, :) + dx**(j - 1) * temp(j, :) / factorial
            end do
            deallocate(pivot_x, derivs)
         end if
      end do

      if (keep_intercept) then
         call move_alloc(raw, basis)
      else
         allocate(basis(size(x), ncoef - 1))
         basis = raw(:, 2:ncoef)
      end if
   end subroutine bs_basis

   subroutine natural_spline_basis(x, basis, knots, df, intercept, boundary_knots, &
                                   interior_knots, status)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: basis(:, :)
      real(dp), intent(in), optional :: knots(:), boundary_knots(2)
      integer, intent(in), optional :: df
      logical, intent(in), optional :: intercept
      real(dp), allocatable, intent(out), optional :: interior_knots(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: iknots(:), all_knots(:), raw(:, :), temp(:, :)
      real(dp), allocatable :: constraints(:, :), reduced(:, :), pivot_x(:)
      real(dp) :: bounds(2), p, dx
      logical :: keep_intercept
      integer :: nik, i, ncoef, ierr
      integer :: derivs2(2)

      if (present(status)) status = 0
      keep_intercept = .false.; if (present(intercept)) keep_intercept = intercept
      if (size(x) == 0) then
         allocate(basis(0, 0)); if (present(status)) status = 1; return
      end if
      if (present(boundary_knots)) then
         bounds = boundary_knots
         if (bounds(2) < bounds(1)) bounds = bounds(2:1:-1)
      else
         bounds = [minval(x), maxval(x)]
      end if
      if (bounds(2) <= bounds(1)) then
         allocate(basis(0, 0)); if (present(status)) status = 2; return
      end if

      if (present(knots)) then
         allocate(iknots(size(knots))); iknots = sorted_copy(knots)
      else if (present(df)) then
         nik = max(0, df - 1 - merge(1, 0, keep_intercept))
         allocate(iknots(nik))
         do i = 1, nik
            p = real(i, dp) / real(nik + 1, dp)
            iknots(i) = type7_quantile(pack(x, x >= bounds(1) .and. x <= bounds(2)), p)
         end do
      else
         allocate(iknots(0))
      end if
      if (present(interior_knots)) interior_knots = iknots

      allocate(all_knots(8 + size(iknots)))
      all_knots(1:4) = bounds(1)
      if (size(iknots) > 0) all_knots(5:4 + size(iknots)) = iknots
      all_knots(5 + size(iknots):) = bounds(2)
      ncoef = size(all_knots) - 4
      allocate(raw(size(x), ncoef)); raw = 0.0_dp

      do i = 1, size(x)
         if (x(i) >= bounds(1) .and. x(i) <= bounds(2)) then
            call spline_design(all_knots, [x(i)], 4, temp, status=ierr)
            if (ierr /= 0) then
               allocate(basis(0, 0)); if (present(status)) status = 3; return
            end if
            raw(i, :) = temp(1, :)
         else
            allocate(pivot_x(2))
            pivot_x = merge(bounds(1), bounds(2), x(i) < bounds(1))
            call spline_design(all_knots, pivot_x, 4, temp, [0, 1], ierr)
            if (ierr /= 0) then
               allocate(basis(0, 0)); if (present(status)) status = 4; return
            end if
            dx = x(i) - pivot_x(1)
            raw(i, :) = temp(1, :) + dx * temp(2, :)
            deallocate(pivot_x)
         end if
      end do

      derivs2 = 2
      call spline_design(all_knots, bounds, 4, constraints, derivs2, ierr)
      if (ierr /= 0) then
         allocate(basis(0, 0)); if (present(status)) status = 5; return
      end if
      if (.not. keep_intercept) then
         reduced = raw(:, 2:ncoef)
         constraints = constraints(:, 2:ncoef)
      else
         reduced = raw
      end if
      call nullspace_transform(constraints, reduced, basis, ierr)
      if (ierr /= 0) then
         allocate(basis(0, 0)); if (present(status)) status = 6; return
      end if
   end subroutine natural_spline_basis

end module splines_basis
