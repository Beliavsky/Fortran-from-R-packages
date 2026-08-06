! Copyright (C) 1998 Douglas M. Bates and William N. Venables.
! Modern Fortran translation, 2026.
! SPDX-License-Identifier: GPL-2.0-or-later
module splines_core
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use splines_kinds, only : dp
   use splines_linalg, only : solve_linear
   implicit none
   private

   type, public :: b_spline_t
      real(dp), allocatable :: knots(:)
      real(dp), allocatable :: coefficients(:)
      integer :: order = 4
      logical :: natural = .false.
      logical :: periodic = .false.
      real(dp) :: period = 0.0_dp
   contains
      procedure :: evaluate_scalar => b_spline_evaluate_scalar
      procedure :: evaluate_vector => b_spline_evaluate_vector
      generic :: evaluate => evaluate_scalar, evaluate_vector
   end type b_spline_t

   type, public :: poly_spline_t
      real(dp), allocatable :: knots(:)
      real(dp), allocatable :: coefficients(:, :)
      logical :: natural = .false.
      logical :: periodic = .false.
      real(dp) :: period = 0.0_dp
   contains
      procedure :: evaluate_scalar => poly_spline_evaluate_scalar
      procedure :: evaluate_vector => poly_spline_evaluate_vector
      generic :: evaluate => evaluate_scalar, evaluate_vector
   end type poly_spline_t

   public :: spline_design, spline_basis_nonzero, linear_interp
   public :: fit_interpolating_spline, fit_periodic_spline
   public :: to_polynomial_spline, inverse_monotone_spline
   public :: sorted_copy, type7_quantile

contains

   pure function quiet_nan() result(x)
      real(dp) :: x
      x = ieee_value(0.0_dp, ieee_quiet_nan)
   end function quiet_nan

   pure subroutine insertion_sort(x)
      real(dp), intent(inout) :: x(:)
      real(dp) :: key
      integer :: i, j
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
   end subroutine insertion_sort

   function sorted_copy(x) result(y)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable :: y(:)
      allocate(y(size(x)))
      y = x
      call insertion_sort(y)
   end function sorted_copy

   function type7_quantile(x, probability, status) result(q)
      real(dp), intent(in) :: x(:), probability
      integer, intent(out), optional :: status
      real(dp) :: q, h, frac
      real(dp), allocatable :: sx(:)
      integer :: n, j

      if (present(status)) status = 0
      n = size(x)
      if (n == 0 .or. probability < 0.0_dp .or. probability > 1.0_dp) then
         q = quiet_nan()
         if (present(status)) status = 1
         return
      end if
      sx = sorted_copy(x)
      if (n == 1) then
         q = sx(1)
         return
      end if
      h = 1.0_dp + real(n - 1, dp) * probability
      j = int(floor(h))
      frac = h - real(j, dp)
      if (j >= n) then
         q = sx(n)
      else
         q = (1.0_dp - frac) * sx(j) + frac * sx(j + 1)
      end if
   end function type7_quantile

   pure subroutine difference_table(knots, ti, x, n, rdel, ldel)
      real(dp), intent(in) :: knots(:), x
      integer, intent(in) :: ti, n
      real(dp), intent(out) :: rdel(:), ldel(:)
      integer :: i
      do i = 1, n
         rdel(i) = knots(ti + i - 1) - x
         ldel(i) = x - knots(ti - i)
      end do
   end subroutine difference_table

   subroutine basis_functions(knots, ti, x, order, values)
      real(dp), intent(in) :: knots(:), x
      integer, intent(in) :: ti, order
      real(dp), intent(out) :: values(order)
      real(dp), allocatable :: rdel(:), ldel(:)
      real(dp) :: saved, term, denom
      integer :: j, r

      values = 0.0_dp
      values(1) = 1.0_dp
      if (order == 1) return
      allocate(rdel(order - 1), ldel(order - 1))
      call difference_table(knots, ti, x, order - 1, rdel, ldel)
      do j = 1, order - 1
         saved = 0.0_dp
         do r = 0, j - 1
            denom = rdel(r + 1) + ldel(j - r)
            if (abs(denom) <= tiny(1.0_dp)) then
               term = 0.0_dp
            else
               term = values(r + 1) / denom
            end if
            values(r + 1) = saved + rdel(r + 1) * term
            saved = ldel(j - r) * term
         end do
         values(j + 1) = saved
      end do
   end subroutine basis_functions

   function evaluate_local(knots, ti, x, coefficients, order, derivative) result(value)
      integer, intent(in) :: ti, order, derivative
      real(dp), intent(in) :: knots(:), x, coefficients(order)
      real(dp) :: value
      real(dp), allocatable :: a(:), rdel(:), ldel(:)
      real(dp) :: denom
      integer :: outer, inner, i, lidx

      if (derivative < 0 .or. derivative >= order) then
         value = quiet_nan()
         return
      end if
      allocate(a(order))
      a = coefficients
      outer = order - 1
      do i = 1, derivative
         do inner = 1, outer
            lidx = ti - outer + inner - 1
            denom = knots(lidx + outer) - knots(lidx)
            if (abs(denom) <= tiny(1.0_dp)) then
               a(inner) = 0.0_dp
            else
               a(inner) = real(outer, dp) * (a(inner + 1) - a(inner)) / denom
            end if
         end do
         outer = outer - 1
      end do
      if (outer == 0) then
         value = a(1)
         return
      end if
      allocate(rdel(outer), ldel(outer))
      call difference_table(knots, ti, x, outer, rdel, ldel)
      do i = outer, 1, -1
         do inner = 1, i
            denom = rdel(inner) + ldel(i + 1 - inner)
            if (abs(denom) <= tiny(1.0_dp)) then
               a(inner) = 0.0_dp
            else
               a(inner) = (a(inner + 1) * ldel(i + 1 - inner) + &
                           a(inner) * rdel(inner)) / denom
            end if
         end do
      end do
      value = a(1)
   end function evaluate_local

   subroutine spline_basis_nonzero(knots, x, order, derivative, values, offset, status)
      real(dp), intent(in) :: knots(:), x
      integer, intent(in) :: order, derivative
      real(dp), intent(out) :: values(:)
      integer, intent(out) :: offset
      integer, intent(out), optional :: status
      real(dp), allocatable :: unit(:)
      integer :: nk, ncoef, ti, i

      if (present(status)) status = 0
      values = 0.0_dp
      offset = 0
      nk = size(knots)
      ncoef = nk - order
      if (order < 1 .or. ncoef < order .or. size(values) /= order .or. &
          derivative < 0 .or. derivative >= order) then
         if (present(status)) status = 1
         return
      end if
      if (any(knots(2:nk) < knots(1:nk - 1))) then
         if (present(status)) status = 2
         return
      end if
      if (x < knots(order) .or. x > knots(ncoef + 1)) then
         if (present(status)) status = 3
         return
      end if

      ti = order + 1
      do while (ti < ncoef + 1 .and. knots(ti) <= x)
         ti = ti + 1
      end do
      offset = ti - (order + 1)
      if (derivative == 0) then
         call basis_functions(knots, ti, x, order, values)
      else
         allocate(unit(order))
         do i = 1, order
            unit = 0.0_dp
            unit(i) = 1.0_dp
            values(i) = evaluate_local(knots, ti, x, unit, order, derivative)
         end do
      end if
   end subroutine spline_basis_nonzero

   subroutine spline_design(knots, x, order, design, derivs, status)
      real(dp), intent(in) :: knots(:), x(:)
      integer, intent(in) :: order
      real(dp), allocatable, intent(out) :: design(:, :)
      integer, intent(in), optional :: derivs(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: local(:), sknots(:)
      integer :: ncoef, i, d, offset, ierr

      if (present(status)) status = 0
      ncoef = size(knots) - order
      if (order < 1 .or. ncoef < order) then
         allocate(design(0, 0))
         if (present(status)) status = 1
         return
      end if
      if (present(derivs)) then
         if (size(derivs) /= size(x)) then
            allocate(design(0, 0))
            if (present(status)) status = 2
            return
         end if
      end if
      sknots = sorted_copy(knots)
      allocate(design(size(x), ncoef), local(order))
      design = 0.0_dp
      do i = 1, size(x)
         d = 0
         if (present(derivs)) d = derivs(i)
         call spline_basis_nonzero(sknots, x(i), order, d, local, offset, ierr)
         if (ierr /= 0) then
            if (present(status)) status = 10 + ierr
            design = quiet_nan()
            return
         end if
         design(i, offset + 1:offset + order) = local
      end do
   end subroutine spline_design

   function wrap_periodic(x, lower, period) result(z)
      real(dp), intent(in) :: x, lower, period
      real(dp) :: z
      z = lower + modulo(x - lower, period)
      if (z < lower) z = z + period
   end function wrap_periodic

   recursive function b_spline_evaluate_scalar(self, x, derivative, status) result(y)
      class(b_spline_t), intent(in) :: self
      real(dp), intent(in) :: x
      integer, intent(in), optional :: derivative
      integer, intent(out), optional :: status
      real(dp) :: y, xx, lower, upper, slope
      real(dp), allocatable :: local(:)
      integer :: d, offset, ierr, ncoef

      if (present(status)) status = 0
      d = 0
      if (present(derivative)) d = derivative
      ncoef = size(self%coefficients)
      if (.not. allocated(self%knots) .or. ncoef == 0 .or. &
          size(self%knots) /= ncoef + self%order .or. d < 0 .or. d >= self%order) then
         y = quiet_nan()
         if (present(status)) status = 1
         return
      end if
      lower = self%knots(self%order)
      upper = self%knots(ncoef + 1)
      xx = x
      if (self%periodic) xx = wrap_periodic(x, lower, self%period)
      if (xx < lower .or. xx > upper) then
         if (.not. self%natural) then
            y = quiet_nan()
            if (present(status)) status = 2
            return
         end if
         if (d == 0) then
            if (xx < lower) then
               y = self%evaluate_scalar(lower, 0, ierr)
               slope = self%evaluate_scalar(lower, 1, ierr)
               y = y + slope * (xx - lower)
            else
               y = self%evaluate_scalar(upper, 0, ierr)
               slope = self%evaluate_scalar(upper, 1, ierr)
               y = y + slope * (xx - upper)
            end if
         else if (d == 1) then
            if (xx < lower) then
               y = self%evaluate_scalar(lower, 1, ierr)
            else
               y = self%evaluate_scalar(upper, 1, ierr)
            end if
         else
            y = 0.0_dp
         end if
         return
      end if
      allocate(local(self%order))
      call spline_basis_nonzero(self%knots, xx, self%order, d, local, offset, ierr)
      if (ierr /= 0) then
         y = quiet_nan()
         if (present(status)) status = 3
         return
      end if
      y = dot_product(local, self%coefficients(offset + 1:offset + self%order))
   end function b_spline_evaluate_scalar

   function b_spline_evaluate_vector(self, x, derivative, status) result(y)
      class(b_spline_t), intent(in) :: self
      real(dp), intent(in) :: x(:)
      integer, intent(in), optional :: derivative
      integer, intent(out), optional :: status
      real(dp), allocatable :: y(:)
      integer :: i, d, ierr
      d = 0
      if (present(derivative)) d = derivative
      if (present(status)) status = 0
      allocate(y(size(x)))
      do i = 1, size(x)
         y(i) = self%evaluate_scalar(x(i), d, ierr)
         if (ierr /= 0 .and. present(status)) status = ierr
      end do
   end function b_spline_evaluate_vector

   function poly_spline_evaluate_scalar(self, x, derivative, status) result(y)
      class(poly_spline_t), intent(in) :: self
      real(dp), intent(in) :: x
      integer, intent(in), optional :: derivative
      integer, intent(out), optional :: status
      real(dp) :: y, xx, dx
      real(dp), allocatable :: c(:)
      integer :: d, i, j, k, ord, nk

      if (present(status)) status = 0
      d = 0
      if (present(derivative)) d = derivative
      if (.not. allocated(self%knots) .or. .not. allocated(self%coefficients)) then
         y = quiet_nan(); if (present(status)) status = 1; return
      end if
      nk = size(self%knots)
      ord = size(self%coefficients, 2)
      if (size(self%coefficients, 1) /= nk .or. d < 0 .or. d >= ord) then
         y = quiet_nan(); if (present(status)) status = 1; return
      end if
      xx = x
      if (self%periodic) xx = wrap_periodic(x, self%knots(1), self%period)
      if (xx < self%knots(1) .or. xx > self%knots(nk)) then
         if (.not. self%natural) then
            y = quiet_nan(); if (present(status)) status = 2; return
         end if
         if (d == 0) then
            if (xx < self%knots(1)) then
               y = self%coefficients(1, 1) + self%coefficients(1, 2) * (xx - self%knots(1))
            else
               y = self%coefficients(nk, 1) + self%coefficients(nk, 2) * (xx - self%knots(nk))
            end if
         else if (d == 1) then
            if (xx < self%knots(1)) then
               y = self%coefficients(1, 2)
            else
               y = self%coefficients(nk, 2)
            end if
         else
            y = 0.0_dp
         end if
         return
      end if

      i = nk
      do j = 1, nk - 1
         if (xx < self%knots(j + 1)) then
            i = j
            exit
         end if
      end do
      allocate(c(ord))
      c = self%coefficients(i, :)
      do j = 1, d
         c(1:ord - j) = [(real(k, dp) * c(k + 1), k=1,ord-j)]
      end do
      ord = ord - d
      dx = xx - self%knots(i)
      y = c(ord)
      do j = ord - 1, 1, -1
         y = y * dx + c(j)
      end do
   end function poly_spline_evaluate_scalar

   function poly_spline_evaluate_vector(self, x, derivative, status) result(y)
      class(poly_spline_t), intent(in) :: self
      real(dp), intent(in) :: x(:)
      integer, intent(in), optional :: derivative
      integer, intent(out), optional :: status
      real(dp), allocatable :: y(:)
      integer :: i, d, ierr
      d = 0
      if (present(derivative)) d = derivative
      if (present(status)) status = 0
      allocate(y(size(x)))
      do i = 1, size(x)
         y(i) = self%evaluate_scalar(x(i), d, ierr)
         if (ierr /= 0 .and. present(status)) status = ierr
      end do
   end function poly_spline_evaluate_vector

   subroutine linear_interp(x, y, x0, y0, status)
      real(dp), intent(in) :: x(:), y(:), x0(:)
      real(dp), allocatable, intent(out) :: y0(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: sx(:), sy(:)
      integer, allocatable :: idx(:)
      integer :: n, i, j, k
      real(dp) :: tx, ty

      if (present(status)) status = 0
      n = size(x)
      allocate(y0(size(x0)))
      if (n < 1 .or. size(y) /= n) then
         y0 = quiet_nan(); if (present(status)) status = 1; return
      end if
      allocate(sx(n), sy(n), idx(n))
      sx = x; sy = y; idx = [(i, i=1,n)]
      do i = 2, n
         tx = sx(i); ty = sy(i); k = idx(i); j = i - 1
         do while (j >= 1)
            if (sx(j) <= tx) exit
            sx(j + 1) = sx(j); sy(j + 1) = sy(j); idx(j + 1) = idx(j)
            j = j - 1
         end do
         sx(j + 1) = tx; sy(j + 1) = ty; idx(j + 1) = k
      end do
      if (any(sx(2:n) <= sx(1:n - 1))) then
         y0 = quiet_nan(); if (present(status)) status = 2; return
      end if
      do i = 1, size(x0)
         if (x0(i) < sx(1) .or. x0(i) > sx(n)) then
            y0(i) = quiet_nan(); if (present(status)) status = 3
         else if (abs(x0(i) - sx(n)) <= epsilon(1.0_dp) * max(1.0_dp, abs(sx(n)))) then
            y0(i) = sy(n)
         else
            j = 1
            do while (j < n - 1 .and. sx(j + 1) <= x0(i))
               j = j + 1
            end do
            y0(i) = sy(j) + (sy(j + 1) - sy(j)) * (x0(i) - sx(j)) / (sx(j + 1) - sx(j))
         end if
      end do
   end subroutine linear_interp

   subroutine fit_interpolating_spline(x, y, spline, status)
      real(dp), intent(in) :: x(:), y(:)
      type(b_spline_t), intent(out) :: spline
      integer, intent(out), optional :: status
      real(dp), allocatable :: sx(:), sy(:), knots(:), xx(:), rhs(:), design(:, :), coeff(:)
      integer, allocatable :: derivs(:)
      real(dp) :: tx, ty
      integer :: n, i, j, ierr

      if (present(status)) status = 0
      n = size(x)
      if (n < 4 .or. size(y) /= n) then
         if (present(status)) status = 1
         return
      end if
      allocate(sx(n), sy(n))
      sx = x; sy = y
      do i = 2, n
         tx = sx(i); ty = sy(i); j = i - 1
         do while (j >= 1)
            if (sx(j) <= tx) exit
            sx(j + 1) = sx(j); sy(j + 1) = sy(j); j = j - 1
         end do
         sx(j + 1) = tx; sy(j + 1) = ty
      end do
      if (any(sx(2:n) <= sx(1:n - 1))) then
         if (present(status)) status = 2
         return
      end if
      allocate(knots(n + 6))
      knots(1:3) = sx(1:3) + sx(1) - sx(4)
      knots(4:n + 3) = sx
      knots(n + 4:n + 6) = sx(n - 2:n) + sx(n) - sx(n - 3)
      allocate(xx(n + 2), derivs(n + 2), rhs(n + 2))
      xx = [sx(1), sx, sx(n)]
      derivs = 0; derivs(1) = 2; derivs(n + 2) = 2
      rhs = [0.0_dp, sy, 0.0_dp]
      call spline_design(knots, xx, 4, design, derivs, ierr)
      if (ierr /= 0) then
         if (present(status)) status = 3
         return
      end if
      call solve_linear(design, rhs, coeff, ierr)
      if (ierr /= 0) then
         if (present(status)) status = 4
         return
      end if
      spline%knots = knots
      spline%coefficients = coeff
      spline%order = 4
      spline%natural = .true.
      spline%periodic = .false.
      spline%period = 0.0_dp
   end subroutine fit_interpolating_spline

   subroutine fit_periodic_spline(x, y, spline, period, order, knots, status)
      real(dp), intent(in) :: x(:), y(:)
      type(b_spline_t), intent(out) :: spline
      real(dp), intent(in), optional :: period
      integer, intent(in), optional :: order
      real(dp), intent(in), optional :: knots(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: sx(:), sy(:), aknots(:), design(:, :), system(:, :), coeff0(:), coeff(:)
      real(dp) :: p, tx, ty
      integer :: n, ord, i, j, ierr

      if (present(status)) status = 0
      n = size(x); ord = 4; if (present(order)) ord = order
      p = 2.0_dp * acos(-1.0_dp); if (present(period)) p = period
      if (n < ord .or. size(y) /= n .or. ord < 1) then
         if (present(status)) status = 1; return
      end if
      allocate(sx(n), sy(n)); sx = x; sy = y
      do i = 2, n
         tx = sx(i); ty = sy(i); j = i - 1
         do while (j >= 1)
            if (sx(j) <= tx) exit
            sx(j + 1) = sx(j); sy(j + 1) = sy(j); j = j - 1
         end do
         sx(j + 1) = tx; sy(j + 1) = ty
      end do
      if (any(sx(2:n) <= sx(1:n - 1))) then
         if (present(status)) status = 2; return
      end if
      if (present(knots)) then
         aknots = knots
         if (size(aknots) /= n + 2 * ord - 1) then
            if (present(status)) status = 3; return
         end if
         p = aknots(size(aknots) + 1 - ord) - aknots(1)
      else
         allocate(aknots(n + 2 * ord - 1))
         aknots(1:ord - 1) = sx(n - ord + 2:n) - p
         aknots(ord:ord + n - 1) = sx
         aknots(ord + n:) = sx(1:ord) + p
      end if
      if (sx(n) - sx(1) >= p) then
         if (present(status)) status = 4; return
      end if
      call spline_design(aknots, sx, ord, design, status=ierr)
      if (ierr /= 0) then
         if (present(status)) status = 5; return
      end if
      allocate(system(n, n))
      system = design(:, 1:n)
      system(:, 1:ord - 1) = system(:, 1:ord - 1) + design(:, n + 1:n + ord - 1)
      call solve_linear(system, sy, coeff0, ierr)
      if (ierr /= 0) then
         if (present(status)) status = 6; return
      end if
      allocate(coeff(n + ord - 1))
      coeff(1:n) = coeff0
      coeff(n + 1:) = coeff0(1:ord - 1)
      spline%knots = aknots
      spline%coefficients = coeff
      spline%order = ord
      spline%periodic = .true.
      spline%natural = .false.
      spline%period = p
   end subroutine fit_periodic_spline

   subroutine to_polynomial_spline(spline, poly, status)
      type(b_spline_t), intent(in) :: spline
      type(poly_spline_t), intent(out) :: poly
      integer, intent(out), optional :: status
      real(dp), allocatable :: breaks(:)
      real(dp) :: factorial
      integer :: ncoef, nbreak, i, d, ierr

      if (present(status)) status = 0
      ncoef = size(spline%coefficients)
      if (size(spline%knots) /= ncoef + spline%order) then
         if (present(status)) status = 1; return
      end if
      nbreak = size(spline%knots) + 2 - 2 * spline%order
      allocate(breaks(nbreak))
      breaks = spline%knots(spline%order:size(spline%knots) + 1 - spline%order)
      allocate(poly%coefficients(nbreak, spline%order))
      factorial = 1.0_dp
      do d = 0, spline%order - 1
         if (d > 0) factorial = factorial * real(d, dp)
         do i = 1, nbreak
            poly%coefficients(i, d + 1) = spline%evaluate_scalar(breaks(i), d, ierr) / factorial
         end do
      end do
      poly%knots = breaks
      poly%natural = spline%natural
      poly%periodic = spline%periodic
      poly%period = spline%period
      if (poly%natural .and. nbreak > 0 .and. spline%order >= 3) then
         poly%coefficients(1, 3) = 0.0_dp
         poly%coefficients(nbreak, 3) = 0.0_dp
      end if
   end subroutine to_polynomial_spline

   subroutine inverse_monotone_spline(spline, inverse, status)
      type(poly_spline_t), intent(in) :: spline
      type(poly_spline_t), intent(out) :: inverse
      integer, intent(out), optional :: status
      real(dp), allocatable :: adiff(:), kdiff(:), a(:, :), b(:), sol(:)
      integer :: n, i, ierr

      if (present(status)) status = 0
      n = size(spline%knots)
      if (n < 2 .or. size(spline%coefficients, 1) /= n .or. size(spline%coefficients, 2) < 4) then
         if (present(status)) status = 1; return
      end if
      kdiff = spline%knots(2:n) - spline%knots(1:n - 1)
      if (any(kdiff <= 0.0_dp)) then
         if (present(status)) status = 2; return
      end if
      inverse%knots = spline%coefficients(:, 1)
      adiff = inverse%knots(2:n) - inverse%knots(1:n - 1)
      if (.not. (all(adiff > 0.0_dp) .or. all(adiff < 0.0_dp))) then
         if (present(status)) status = 3; return
      end if
      allocate(inverse%coefficients(n, 4))
      inverse%coefficients = quiet_nan()
      inverse%coefficients(:, 1) = spline%knots
      inverse%coefficients(:, 2) = 1.0_dp / spline%coefficients(:, 2)
      allocate(a(2, 2), b(2))
      do i = 1, n - 1
         a = reshape([adiff(i)**2, 2.0_dp * adiff(i), adiff(i)**3, 3.0_dp * adiff(i)**2], [2, 2])
         b = [kdiff(i) - adiff(i) * inverse%coefficients(i, 2), &
              inverse%coefficients(i + 1, 2) - inverse%coefficients(i, 2)]
         call solve_linear(a, b, sol, ierr)
         if (ierr /= 0) then
            if (present(status)) status = 4; return
         end if
         inverse%coefficients(i, 3:4) = sol
      end do
      if (n > 2) then
         inverse%coefficients(1, 4) = 0.0_dp
         inverse%coefficients(n - 1, 4) = 0.0_dp
         a = reshape([adiff(1), 1.0_dp, adiff(1)**2, 2.0_dp * adiff(1)], [2, 2])
         b = [kdiff(1), 1.0_dp / spline%coefficients(2, 2)]
         call solve_linear(a, b, sol, ierr)
         inverse%coefficients(1, 2:3) = sol
         inverse%coefficients(n - 1, 3) = (kdiff(n - 1) - adiff(n - 1) * &
            inverse%coefficients(n - 1, 2)) / adiff(n - 1)**2
      end if
      if (inverse%coefficients(1, 3) > 0.0_dp) then
         inverse%coefficients(1, 3) = 0.0_dp
         inverse%coefficients(1, 2) = kdiff(1) / adiff(1)
      end if
      if (inverse%coefficients(n - 1, 3) < 0.0_dp) then
         inverse%coefficients(n - 1, 3) = 0.0_dp
         inverse%coefficients(n - 1, 2) = kdiff(n - 1) / adiff(n - 1)
      end if
      inverse%natural = .false.
      inverse%periodic = .false.
      inverse%period = 0.0_dp
   end subroutine inverse_monotone_spline

end module splines_core
