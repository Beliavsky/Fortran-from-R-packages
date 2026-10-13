! R package splines2 by Wenjie Wang and Jun Yan, Copyright (C) 2016-2025.
! Modern Fortran translation, 2026. SPDX-License-Identifier: GPL-3.0-or-later
module splines2_math
   use r_kinds, only : dp
   use splines_core, only : spline_basis_nonzero
   use splines_linalg, only : solve_linear
   implicit none
   private
   public :: clamped_knots, b_values, b_integrals, m_values, i_integrals
   public :: natural_transform, natural_values, periodic_values
contains
   pure function clamped_knots(knots, bounds, degree) result(t)
      real(dp), intent(in) :: knots(:) !! Sorted interior knots, strictly inside bounds.
      real(dp), intent(in) :: bounds(2) !! Increasing finite boundary knots.
      integer, intent(in) :: degree !! Nonnegative polynomial degree.
      real(dp), allocatable :: t(:)
      t = [spread(bounds(1), 1, degree + 1), knots, spread(bounds(2), 1, degree + 1)]
   end function clamped_knots

   function b_values(x, knots, bounds, degree, deriv) result(v)
      real(dp), intent(in) :: x !! Finite evaluation coordinate; endpoint polynomial is extended outside bounds.
      real(dp), intent(in) :: knots(:) !! Sorted interior knots.
      real(dp), intent(in) :: bounds(2) !! Increasing finite boundary knots.
      integer, intent(in) :: degree !! Nonnegative polynomial degree.
      integer, intent(in) :: deriv !! Nonnegative derivative order; orders above degree return zero.
      real(dp) :: v(size(knots) + degree + 1)
      real(dp) :: local(degree + 1), xx, factor
      real(dp), allocatable :: t(:)
      integer :: k, offset, stat
      v = 0.0_dp
      if (deriv > degree) return
      t = clamped_knots(knots, bounds, degree)
      xx = min(bounds(2), max(bounds(1), x))
      factor = 1.0_dp
      do k = deriv, degree
         if (k > deriv) factor = factor * (x - xx) / real(k - deriv, dp)
         if (factor == 0.0_dp) exit
         call spline_basis_nonzero(t, xx, degree + 1, k, local, offset, stat)
         if (stat /= 0) error stop 'Internal spline evaluator contract violation'
         v(offset + 1:offset + degree + 1) = v(offset + 1:offset + degree + 1) + factor * local
      end do
   end function b_values

   function b_integrals(x, knots, bounds, degree) result(v)
      real(dp), intent(in) :: x !! Finite evaluation point for integral from left boundary.
      real(dp), intent(in) :: knots(:) !! Sorted interior knots.
      real(dp), intent(in) :: bounds(2) !! Increasing finite boundary knots.
      integer, intent(in) :: degree !! Nonnegative degree of the original B-splines.
      real(dp) :: v(size(knots) + degree + 1)
      real(dp), allocatable :: elevated(:), t(:)
      real(dp) :: acc
      integer :: j
      elevated = b_values(x, knots, bounds, degree + 1, 0)
      t = clamped_knots(knots, bounds, degree)
      acc = 0.0_dp
      do j = size(v), 1, -1
         acc = acc + elevated(j + 1)
         v(j) = acc * (t(j + degree + 1) - t(j)) / real(degree + 1, dp)
      end do
   end function b_integrals

   function m_values(x, knots, bounds, degree, deriv) result(v)
      real(dp), intent(in) :: x !! Finite evaluation point.
      real(dp), intent(in) :: knots(:) !! Sorted interior knots.
      real(dp), intent(in) :: bounds(2) !! Increasing finite boundary knots.
      integer, intent(in) :: degree !! Nonnegative degree of the M-splines.
      integer, intent(in) :: deriv !! Derivative order; -1 denotes an integral from the left boundary.
      real(dp) :: v(size(knots) + degree + 1)
      real(dp), allocatable :: t(:)
      integer :: j
      t = clamped_knots(knots, bounds, degree)
      if (deriv == -1) then
         v = b_integrals(x, knots, bounds, degree)
      else
         v = b_values(x, knots, bounds, degree, deriv)
      end if
      do j = 1, size(v)
         v(j) = v(j) * real(degree + 1, dp) / (t(j + degree + 1) - t(j))
      end do
   end function m_values

   function i_integrals(x, knots, bounds, degree) result(v)
      real(dp), intent(in) :: x !! Finite evaluation point for integral of I-splines from left boundary.
      real(dp), intent(in) :: knots(:) !! Sorted interior knots.
      real(dp), intent(in) :: bounds(2) !! Increasing finite boundary knots.
      integer, intent(in) :: degree !! Degree of associated M-splines; I-splines have degree+1.
      real(dp) :: v(size(knots) + degree + 1)
      real(dp), allocatable :: elevated(:)
      real(dp) :: acc
      integer :: j
      elevated = b_integrals(x, knots, bounds, degree + 1)
      acc = 0.0_dp
      do j = size(v), 1, -1
         acc = acc + elevated(j + 1)
         v(j) = acc
      end do
   end function i_integrals

   subroutine natural_transform(knots, bounds, cardinal, h, status)
      real(dp), intent(in) :: knots(:) !! Distinct increasing interior knots for cubic natural splines.
      real(dp), intent(in) :: bounds(2) !! Increasing finite boundary knots.
      logical, intent(in) :: cardinal !! True selects nsk's unit-height basis at every knot.
      real(dp), allocatable, intent(out) :: h(:, :) !! B-spline-to-natural transform, shape (nknots+4,nknots+2).
      integer, intent(out) :: status !! Zero on success; nonzero indicates a singular cardinal transform.
      real(dp), allocatable :: kmat(:, :), row(:), temp(:, :), positions(:), bv(:)
      real(dp) :: w1, w2, w3
      integer :: k, n, i, j
      k = size(knots)
      n = k + 2
      allocate(h(n + 2, n))
      h = 0.0_dp
      status = 0
      if (k == 0) then
         h(1:3, 1) = [3.0_dp, 2.0_dp, 1.0_dp]
         h(2:4, 2) = [1.0_dp, 2.0_dp, 3.0_dp]
      else if (k == 1) then
         w1 = knots(1) - bounds(1)
         w2 = bounds(2) - knots(1)
         w3 = w1 + w2
         h(1, 1) = 1.0_dp + w1 / w3
         h(2, 1) = 1.0_dp
         h(2, 2) = 1.0_dp / (1.0_dp + w3 / w1)
         h(3, 2) = 1.0_dp
         h(4, 2) = 1.0_dp / (1.0_dp + w3 / w2)
         h(4, 3) = 1.0_dp
         h(5, 3) = 1.0_dp + w2 / w3
      else
         h(1:3, 1) = 1.0_dp
         h(n:n + 2, n) = 1.0_dp
         h(2, 2) = 1.0_dp
         h(3, 2) = 1.0_dp + (knots(2) - bounds(1)) / (knots(1) - bounds(1))
         h(n, n - 1) = 1.0_dp + (bounds(2) - knots(k - 1)) / (bounds(2) - knots(k))
         h(n + 1, n - 1) = 1.0_dp
         do j = 1, n - 4
            h(j + 3, j + 2) = 1.0_dp
         end do
      end if
      do j = 1, n
         h(:, j) = h(:, j) / sum(h(:, j))
      end do
      if (.not. cardinal) return
      positions = [bounds(1), knots, bounds(2)]
      allocate(kmat(n, n), temp(n + 2, n))
      do i = 1, n
         bv = b_values(positions(i), knots, bounds, 3, 0)
         kmat(i, :) = matmul(bv, h)
      end do
      do i = 1, n + 2
         call solve_linear(transpose(kmat), h(i, :), row, status)
         if (status /= 0) return
         temp(i, :) = row
      end do
      h = temp
   end subroutine natural_transform

   function natural_values(x, knots, bounds, deriv, h) result(v)
      real(dp), intent(in) :: x !! Finite evaluation point; linear continuation is used outside boundaries.
      real(dp), intent(in) :: knots(:) !! Distinct increasing interior knots.
      real(dp), intent(in) :: bounds(2) !! Increasing finite boundary knots.
      integer, intent(in) :: deriv !! Derivative order, or -1 for an integral (zero left of the left boundary).
      real(dp), intent(in) :: h(:, :) !! Natural-spline transform from natural_transform.
      real(dp) :: v(size(h, 2))
      real(dp), allocatable :: bv(:)
      real(dp) :: xx, dx
      xx = min(bounds(2), max(bounds(1), x))
      dx = x - xx
      if (deriv == -1) then
         bv = b_integrals(xx, knots, bounds, 3)
         if (x > bounds(2)) then
            bv = bv + dx * b_values(xx, knots, bounds, 3, 0) + &
               0.5_dp * dx**2 * b_values(xx, knots, bounds, 3, 1)
         end if
      else if (deriv == 0) then
         bv = b_values(xx, knots, bounds, 3, 0) + dx * b_values(xx, knots, bounds, 3, 1)
      else
         bv = b_values(xx, knots, bounds, 3, deriv)
         if (dx /= 0.0_dp .and. deriv > 1) bv = 0.0_dp
      end if
      v = matmul(bv, h)
   end function natural_values

   function periodic_values(x, knots, bounds, degree, deriv, mspline) result(v)
      real(dp), intent(in) :: x !! Finite coordinate, wrapped to the left-closed fundamental period.
      real(dp), intent(in) :: knots(:) !! Distinct increasing knots; size(knots)+1 must be at least degree.
      real(dp), intent(in) :: bounds(2) !! Increasing boundaries defining the fundamental period.
      integer, intent(in) :: degree !! Nonnegative polynomial degree.
      integer, intent(in) :: deriv !! Derivative order, or -1 for integral accumulated from left boundary.
      logical, intent(in) :: mspline !! True requests unit-period-integral M-splines instead of B-splines.
      real(dp) :: v(size(knots) + 1)
      real(dp), allocatable :: extended(:), raw(:), v0(:), v1(:)
      real(dp) :: period, shift, xx, eb(2)
      integer :: k, p, i, j, n
      k = size(knots)
      p = degree
      n = k + 1
      period = bounds(2) - bounds(1)
      xx = bounds(1) + modulo(x - bounds(1), period)
      shift = anint((x - xx) / period)
      allocate(extended(k + 2 * (p + 1)))
      extended(p + 1:p + k + 2) = [bounds(1), knots, bounds(2)]
      do i = 1, p
         extended(p - i + 1) = extended(p - i + 2) - &
            (extended(p + k - i + 3) - extended(p + k - i + 2))
         extended(p + k + i + 2) = extended(p + k + i + 1) + &
            (extended(p + i + 1) - extended(p + i))
      end do
      eb = [extended(1), extended(size(extended))]
      if (mspline) then
         raw = m_values(xx, extended(2:size(extended) - 1), eb, p, deriv)
      else if (deriv == -1) then
         raw = b_integrals(xx, extended(2:size(extended) - 1), eb, p)
      else
         raw = b_values(xx, extended(2:size(extended) - 1), eb, p, deriv)
      end if
      if (deriv == -1) then
         if (mspline) then
            v0 = m_values(bounds(1), extended(2:size(extended) - 1), eb, p, -1)
            v1 = m_values(bounds(2), extended(2:size(extended) - 1), eb, p, -1)
         else
            v0 = b_integrals(bounds(1), extended(2:size(extended) - 1), eb, p)
            v1 = b_integrals(bounds(2), extended(2:size(extended) - 1), eb, p)
         end if
         raw = raw - v0 + shift * (v1 - v0)
         if (shift < 0.0_dp) raw = 0.0_dp
      end if
      do j = 1, n - p
         v(j) = raw(2 * p + j)
      end do
      do j = 1, p
         v(n - p + j) = raw(p + j) + raw(p + n + j)
      end do
   end function periodic_values
end module splines2_math
