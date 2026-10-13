! R package splines2 by Wenjie Wang and Jun Yan, Copyright (C) 2016-2025.
! Modern Fortran translation, 2026. SPDX-License-Identifier: GPL-3.0-or-later
module splines2
   use r_kinds, only : dp
   use r_quantiles, only : r_quantile_type7
   use r_ordering, only : r_sort_values_in_place
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite, ieee_value, ieee_quiet_nan
   use splines2_math
   implicit none
   private
   public :: dp, spline_basis, b_spline, m_spline, i_spline, c_spline
   public :: natural_spline, nsk, bernstein_poly, dbs, ibs, predict_spline, derivative_spline
   public :: bsp, msp, isp, csp, nsp, bpoly
   integer, parameter, public :: SPLINES2_OK = 0, SPLINES2_INVALID_ARGUMENT = 1, SPLINES2_SINGULAR = 2
   type :: spline_basis
      real(dp), allocatable :: values(:, :) !! Rows correspond to x; columns are the requested spline basis.
      real(dp), allocatable :: x(:), knots(:), transform(:, :), scales(:)
      real(dp) :: boundary_knots(2) = [0.0_dp, 1.0_dp]
      integer :: family = 1, degree = 3, derivs = 0
      logical :: intercept = .false., periodic = .false., integral = .false., scale = .true.
      logical :: equidistant_fallback = .false. !! Quantile knots were replaced by equidistant knots.
      integer :: status = SPLINES2_OK
      character(:), allocatable :: message
   end type spline_basis
   interface bsp
      module procedure b_spline
   end interface bsp
   interface msp
      module procedure m_spline
   end interface msp
   interface isp
      module procedure i_spline
   end interface isp
   interface csp
      module procedure c_spline
   end interface csp
   interface nsp
      module procedure natural_spline
   end interface nsp
   interface bpoly
      module procedure bernstein_poly
   end interface bpoly
contains
   function b_spline(x, df, knots, degree, intercept, &
      boundary_knots, periodic, derivs, integral) result(out)
      real(dp), intent(in) :: x(:) !! Evaluation coordinates; NaN rows propagate, but all-NaN or infinite x is invalid.
      integer, intent(in), optional :: df !! Requested column count for automatic knots; zero or absence uses explicit knots.
      real(dp), intent(in), optional :: knots(:) !! Interior knots, sorted internally; nonempty knots override df.
      integer, intent(in), optional :: degree !! Nonnegative B/M polynomial degree, default 3; I/C use associated M degree.
      logical, intent(in), optional :: intercept !! Retain the first basis column; defaults false except I- and C-splines.
      real(dp), intent(in), optional :: boundary_knots(2) !! Finite distinct boundaries; default is the range of nonmissing x.
      logical, intent(in), optional :: periodic !! Wrap B/M basis periodically across boundaries; default false.
      integer, intent(in), optional :: derivs !! Nonnegative derivative order, default zero (dbs defaults to one).
      logical, intent(in), optional :: integral !! Integrate once from the left boundary, then apply derivs; default false.
      type(spline_basis) :: out
      out = construct_basis(x, 1, &
         df=df, knots=knots, degree=degree, &
         intercept=intercept, boundary_knots=boundary_knots, periodic=periodic, &
         derivs=derivs, integral=integral)
   end function b_spline

   function m_spline(x, df, knots, degree, intercept, &
      boundary_knots, periodic, derivs, integral) result(out)
      real(dp), intent(in) :: x(:) !! Evaluation coordinates; NaN rows propagate, but all-NaN or infinite x is invalid.
      integer, intent(in), optional :: df !! Requested column count for automatic knots; zero or absence uses explicit knots.
      real(dp), intent(in), optional :: knots(:) !! Interior knots, sorted internally; nonempty knots override df.
      integer, intent(in), optional :: degree !! Nonnegative B/M polynomial degree, default 3; I/C use associated M degree.
      logical, intent(in), optional :: intercept !! Retain the first basis column; defaults false except I- and C-splines.
      real(dp), intent(in), optional :: boundary_knots(2) !! Finite distinct boundaries; default is the range of nonmissing x.
      logical, intent(in), optional :: periodic !! Wrap B/M basis periodically across boundaries; default false.
      integer, intent(in), optional :: derivs !! Nonnegative derivative order, default zero (dbs defaults to one).
      logical, intent(in), optional :: integral !! Integrate once from the left boundary, then apply derivs; default false.
      type(spline_basis) :: out
      out = construct_basis(x, 2, &
         df=df, knots=knots, degree=degree, &
         intercept=intercept, boundary_knots=boundary_knots, periodic=periodic, &
         derivs=derivs, integral=integral)
   end function m_spline

   function i_spline(x, df, knots, degree, intercept, &
      boundary_knots, derivs) result(out)
      real(dp), intent(in) :: x(:) !! Evaluation coordinates; NaN rows propagate, but all-NaN or infinite x is invalid.
      integer, intent(in), optional :: df !! Requested column count for automatic knots; zero or absence uses explicit knots.
      real(dp), intent(in), optional :: knots(:) !! Interior knots, sorted internally; nonempty knots override df.
      integer, intent(in), optional :: degree !! Nonnegative B/M polynomial degree, default 3; I/C use associated M degree.
      logical, intent(in), optional :: intercept !! Retain the first basis column; defaults false except I- and C-splines.
      real(dp), intent(in), optional :: boundary_knots(2) !! Finite distinct boundaries; default is the range of nonmissing x.
      integer, intent(in), optional :: derivs !! Nonnegative derivative order, default zero (dbs defaults to one).
      type(spline_basis) :: out
      out = construct_basis(x, 3, &
         df=df, knots=knots, degree=degree, &
         intercept=intercept, boundary_knots=boundary_knots, derivs=derivs)
   end function i_spline

   function c_spline(x, df, knots, degree, intercept, &
      boundary_knots, derivs, scale) result(out)
      real(dp), intent(in) :: x(:) !! Evaluation coordinates; NaN rows propagate, but all-NaN or infinite x is invalid.
      integer, intent(in), optional :: df !! Requested column count for automatic knots; zero or absence uses explicit knots.
      real(dp), intent(in), optional :: knots(:) !! Interior knots, sorted internally; nonempty knots override df.
      integer, intent(in), optional :: degree !! Nonnegative B/M polynomial degree, default 3; I/C use associated M degree.
      logical, intent(in), optional :: intercept !! Retain the first basis column; defaults false except I- and C-splines.
      real(dp), intent(in), optional :: boundary_knots(2) !! Finite distinct boundaries; default is the range of nonmissing x.
      integer, intent(in), optional :: derivs !! Nonnegative derivative order, default zero (dbs defaults to one).
      logical, intent(in), optional :: scale !! Normalize each C-spline to one at the right boundary; default true.
      type(spline_basis) :: out
      out = construct_basis(x, 4, &
         df=df, knots=knots, degree=degree, &
         intercept=intercept, boundary_knots=boundary_knots, derivs=derivs, &
         scale=scale)
   end function c_spline

   function natural_spline(x, df, knots, intercept, boundary_knots, &
      trim, derivs, integral) result(out)
      real(dp), intent(in) :: x(:) !! Evaluation coordinates; NaN rows propagate, but all-NaN or infinite x is invalid.
      integer, intent(in), optional :: df !! Requested column count for automatic knots; zero or absence uses explicit knots.
      real(dp), intent(in), optional :: knots(:) !! Interior knots, sorted internally; nonempty knots override df.
      logical, intent(in), optional :: intercept !! Retain the first basis column; defaults false except I- and C-splines.
      real(dp), intent(in), optional :: boundary_knots(2) !! Finite distinct boundaries; default is the range of nonmissing x.
      real(dp), intent(in), optional :: trim !! Fraction in [0,0.5) trimmed from each end for default natural-spline boundaries.
      integer, intent(in), optional :: derivs !! Nonnegative derivative order, default zero (dbs defaults to one).
      logical, intent(in), optional :: integral !! Integrate once from the left boundary, then apply derivs; default false.
      type(spline_basis) :: out
      out = construct_basis(x, 5, &
         df=df, knots=knots, intercept=intercept, &
         boundary_knots=boundary_knots, trim=trim, derivs=derivs, &
         integral=integral)
   end function natural_spline

   function nsk(x, df, knots, intercept, boundary_knots, &
      trim, derivs, integral) result(out)
      real(dp), intent(in) :: x(:) !! Evaluation coordinates; NaN rows propagate, but all-NaN or infinite x is invalid.
      integer, intent(in), optional :: df !! Requested column count for automatic knots; zero or absence uses explicit knots.
      real(dp), intent(in), optional :: knots(:) !! Interior knots, sorted internally; nonempty knots override df.
      logical, intent(in), optional :: intercept !! Retain the first basis column; defaults false except I- and C-splines.
      real(dp), intent(in), optional :: boundary_knots(2) !! Finite distinct boundaries; default is the range of nonmissing x.
      real(dp), intent(in), optional :: trim !! Fraction in [0,0.5) trimmed from each end for default natural-spline boundaries.
      integer, intent(in), optional :: derivs !! Nonnegative derivative order, default zero (dbs defaults to one).
      logical, intent(in), optional :: integral !! Integrate once from the left boundary, then apply derivs; default false.
      type(spline_basis) :: out
      out = construct_basis(x, 6, &
         df=df, knots=knots, intercept=intercept, &
         boundary_knots=boundary_knots, trim=trim, derivs=derivs, &
         integral=integral)
   end function nsk

   pure function bernstein_poly(x, degree, intercept, boundary_knots, derivs, &
      integral) result(out)
      real(dp), intent(in) :: x(:) !! Evaluation coordinates; NaN rows propagate, but all-NaN or infinite x is invalid.
      integer, intent(in), optional :: degree !! Nonnegative B/M polynomial degree, default 3; I/C use associated M degree.
      logical, intent(in), optional :: intercept !! Retain the first basis column; defaults false except I- and C-splines.
      real(dp), intent(in), optional :: boundary_knots(2) !! Finite distinct boundaries; default is the range of nonmissing x.
      integer, intent(in), optional :: derivs !! Nonnegative derivative order, default zero (dbs defaults to one).
      logical, intent(in), optional :: integral !! Integrate once from the left boundary, then apply derivs; default false.
      type(spline_basis) :: out
      real(dp), allocatable :: finite_x(:), raw(:, :)
      integer :: i, d, first
      out%family = 7
      out%x = x
      allocate(out%values(0, 0), out%knots(0), out%transform(0, 0), out%scales(0))
      out%message = ''
      if (present(degree)) out%degree = degree
      if (present(derivs)) out%derivs = derivs
      if (present(integral)) out%integral = integral
      if (present(intercept)) out%intercept = intercept
      if (out%degree < 0 .or. out%derivs < 0) then
         call invalid(out, 'Degree and derivative order must be nonnegative.')
         return
      end if
      finite_x = pack(x, .not. ieee_is_nan(x))
      if (size(finite_x) == 0) then
         call invalid(out, 'At least one nonmissing coordinate is required.')
         return
      end if
      if (.not. all(ieee_is_finite(finite_x))) then
         call invalid(out, 'Infinite coordinates are not supported.')
         return
      end if
      out%boundary_knots = [minval(finite_x), maxval(finite_x)]
      if (present(boundary_knots)) out%boundary_knots = boundary_knots
      if (.not. all(ieee_is_finite(out%boundary_knots)) .or. &
         out%boundary_knots(1) >= out%boundary_knots(2)) then
         call invalid(out, 'Bernstein boundary knots must be finite, distinct, and increasing.')
         return
      end if
      if (any(finite_x < out%boundary_knots(1) .or. finite_x > out%boundary_knots(2))) then
         call invalid(out, 'Bernstein coordinates must lie inside the boundary knots.')
         return
      end if
      first = merge(1, 2, out%intercept)
      if (out%degree + 1 < first) then
         call invalid(out, 'No basis column remains without the intercept.')
         return
      end if
      allocate(raw(size(x), out%degree + 1))
      d = out%derivs - merge(1, 0, out%integral)
      do i = 1, size(x)
         if (ieee_is_nan(x(i))) then
            raw(i, :) = ieee_value(0.0_dp, ieee_quiet_nan)
         else
            raw(i, :) = bernstein_values(x(i), out%boundary_knots, out%degree, d)
         end if
      end do
      out%values = raw(:, first:)
   end function bernstein_poly

   function dbs(x, derivs, df, knots, degree, &
      intercept, boundary_knots, periodic) result(out)
      real(dp), intent(in) :: x(:) !! Evaluation coordinates; NaN rows propagate, but all-NaN or infinite x is invalid.
      integer, intent(in), optional :: derivs !! Nonnegative derivative order, default zero (dbs defaults to one).
      integer, intent(in), optional :: df !! Requested column count for automatic knots; zero or absence uses explicit knots.
      real(dp), intent(in), optional :: knots(:) !! Interior knots, sorted internally; nonempty knots override df.
      integer, intent(in), optional :: degree !! Nonnegative B/M polynomial degree, default 3; I/C use associated M degree.
      logical, intent(in), optional :: intercept !! Retain the first basis column; defaults false except I- and C-splines.
      real(dp), intent(in), optional :: boundary_knots(2) !! Finite distinct boundaries; default is the range of nonmissing x.
      logical, intent(in), optional :: periodic !! Wrap B/M basis periodically across boundaries; default false.
      type(spline_basis) :: out
      integer :: d
      d = 1
      if (present(derivs)) d = derivs
      out = construct_basis(x, 1, &
         derivs=d, df=df, knots=knots, &
         degree=degree, intercept=intercept, boundary_knots=boundary_knots, &
         periodic=periodic)
   end function dbs

   function ibs(x, df, knots, degree, intercept, &
      boundary_knots, periodic) result(out)
      real(dp), intent(in) :: x(:) !! Evaluation coordinates; NaN rows propagate, but all-NaN or infinite x is invalid.
      integer, intent(in), optional :: df !! Requested column count for automatic knots; zero or absence uses explicit knots.
      real(dp), intent(in), optional :: knots(:) !! Interior knots, sorted internally; nonempty knots override df.
      integer, intent(in), optional :: degree !! Nonnegative B/M polynomial degree, default 3; I/C use associated M degree.
      logical, intent(in), optional :: intercept !! Retain the first basis column; defaults false except I- and C-splines.
      real(dp), intent(in), optional :: boundary_knots(2) !! Finite distinct boundaries; default is the range of nonmissing x.
      logical, intent(in), optional :: periodic !! Wrap B/M basis periodically across boundaries; default false.
      type(spline_basis) :: out
      out = construct_basis(x, 1, &
         df=df, knots=knots, degree=degree, &
         intercept=intercept, boundary_knots=boundary_knots, periodic=periodic, &
         integral=.true.)
   end function ibs

   function construct_basis(x, family, df, knots, degree, intercept, boundary_knots, &
      periodic, derivs, integral, scale, trim) result(out)
      real(dp), intent(in) :: x(:) !! Coordinates, allowing isolated NaNs but no infinities.
      integer, intent(in) :: family !! Basis code: 1 B, 2 M, 3 I, 4 C, 5 natural, 6 nsk, 7 Bernstein.
      integer, intent(in), optional :: df !! Requested number of columns for automatic interior knots.
      real(dp), intent(in), optional :: knots(:) !! Interior knots, sorted before validation.
      integer, intent(in), optional :: degree !! Nonnegative degree; ignored for natural cubic splines.
      logical, intent(in), optional :: intercept !! Retain the first column, default true only for I/C.
      real(dp), intent(in), optional :: boundary_knots(2) !! Finite distinct boundaries; overrides trimming.
      logical, intent(in), optional :: periodic !! Use periodic B/M basis, default false.
      integer, intent(in), optional :: derivs !! Nonnegative derivative order, default zero.
      logical, intent(in), optional :: integral !! Integrate once before differentiating, default false.
      logical, intent(in), optional :: scale !! Normalize C-splines at the right boundary, default true.
      real(dp), intent(in), optional :: trim !! Fraction used to quantile-trim default natural boundaries.
      type(spline_basis) :: out
      real(dp), allocatable :: finite_x(:), inside(:), row(:), raw(:, :)
      real(dp) :: fraction, prob, a, b
      integer :: n, d, k, i, nr, first, stat, multiplicity
      logical :: fallback
      out%family = family
      out%x = x
      allocate(out%values(0, 0), out%knots(0), out%transform(0, 0), out%scales(0))
      out%message = ''
      out%intercept = family == 3 .or. family == 4
      if (present(intercept)) out%intercept = intercept
      if (present(degree)) out%degree = degree
      if (family == 5 .or. family == 6) out%degree = 3
      if (present(periodic)) out%periodic = periodic
      if (present(derivs)) out%derivs = derivs
      if (present(integral)) out%integral = integral
      if (present(scale)) out%scale = scale
      if (out%degree < 0 .or. out%derivs < 0) then
         call invalid(out, 'Degree and derivative order must be nonnegative.')
         return
      end if
      if (size(x) == 0 .or. all(ieee_is_nan(x))) then
         call invalid(out, 'At least one nonmissing coordinate is required.')
         return
      end if
      finite_x = pack(x, .not. ieee_is_nan(x))
      if (.not. all(ieee_is_finite(finite_x))) then
         call invalid(out, 'Infinite coordinates are not supported.')
         return
      end if
      fraction = 0.0_dp
      if (present(trim) .and. .not. present(boundary_knots)) fraction = trim
      if (.not. ieee_is_finite(fraction) .or. fraction < 0.0_dp .or. fraction >= 0.5_dp) then
         call invalid(out, 'Trim must be finite and in [0, 0.5).')
         return
      end if
      if (present(boundary_knots)) then
         out%boundary_knots = boundary_knots
         if (family /= 7) call r_sort_values_in_place(out%boundary_knots)
      else
         out%boundary_knots = [r_quantile_type7(finite_x, fraction), &
            r_quantile_type7(finite_x, 1.0_dp - fraction)]
      end if
      a = out%boundary_knots(1)
      b = out%boundary_knots(2)
      if (.not. all(ieee_is_finite(out%boundary_knots)) .or. a >= b) then
         call invalid(out, 'Boundary knots must be finite and distinct, increasing for Bernstein.')
         return
      end if
      if (family == 7 .and. any(finite_x < a .or. finite_x > b)) then
         call invalid(out, 'Bernstein coordinates must lie inside the boundary knots.')
         return
      end if
      if (present(knots)) then
         out%knots = knots
         if (.not. all(ieee_is_finite(out%knots))) then
            call invalid(out, 'Interior knots must be finite.')
            return
         end if
         call r_sort_values_in_place(out%knots)
      end if
      if (present(df)) then
         if (family == 2 .and. out%periodic .and. df < out%degree) then
            call invalid(out, 'Periodic M-splines require explicit df to be at least degree.')
            return
         end if
         if (df < 0 .or. ((family == 5 .or. family == 6) .and. df < 2)) then
            call invalid(out, 'Invalid requested degrees of freedom.')
            return
         end if
         if (size(out%knots) == 0 .and. df > 0) then
            k = df + merge(0, 1, out%intercept) - out%degree - 1
            if (family == 5 .or. family == 6) k = df + merge(0, 1, out%intercept) - 2
            if (out%periodic) k = df + merge(0, 1, out%intercept) - 1
            if (k < 0) then
               call invalid(out, 'Requested df is too small for the degree and intercept.')
               return
            end if
            if (out%periodic) then
               inside = a + modulo(finite_x - a, b - a)
            else
               inside = pack(finite_x, finite_x >= a .and. finite_x <= b)
            end if
            if (k > 0 .and. size(inside) == 0) then
               call invalid(out, 'Automatic knots require observations inside the boundaries.')
               return
            end if
            deallocate(out%knots)
            allocate(out%knots(k))
            do i = 1, k
               prob = real(i, dp) / real(k + 1, dp)
               out%knots(i) = r_quantile_type7(inside, prob)
            end do
            fallback = any(out%knots <= a .or. out%knots >= b)
            if (k > 1) fallback = fallback .or. any(out%knots(2:) <= out%knots(:k - 1))
            if (fallback) then
               do i = 1, k
                  out%knots(i) = a + (b - a) * real(i, dp) / real(k + 1, dp)
               end do
               out%equidistant_fallback = .true.
            end if
         end if
      end if
      k = size(out%knots)
      if (any(out%knots <= a .or. out%knots >= b)) then
         call invalid(out, 'Interior knots must lie strictly inside the boundary knots.')
         return
      end if
      multiplicity = 1
      do i = 2, k
         if (out%knots(i) == out%knots(i - 1)) then
            multiplicity = multiplicity + 1
         else
            multiplicity = 1
         end if
         if (multiplicity > out%degree + 1 .or. &
            ((family == 5 .or. family == 6 .or. out%periodic) .and. multiplicity > 1)) then
            call invalid(out, 'Unsupported interior knot multiplicity.')
            return
         end if
      end do
      n = k + out%degree + 1
      if (family == 5 .or. family == 6) n = k + 2
      if (out%periodic) then
         n = k + 1
         if (n < out%degree) then
            call invalid(out, 'Periodic splines need at least degree-1 interior knots.')
            return
         end if
      end if
      first = merge(1, 2, out%intercept)
      if (n < first) then
         call invalid(out, 'No basis column remains without the intercept.')
         return
      end if
      if (family == 5 .or. family == 6) then
         call natural_transform(out%knots, out%boundary_knots, family == 6, out%transform, stat)
         if (stat /= 0) then
            call invalid(out, 'The nsk transformation is numerically singular.')
            out%status = SPLINES2_SINGULAR
            return
         end if
      end if
      if (family == 4) then
         out%scales = i_integrals(b, out%knots, out%boundary_knots, out%degree)
         if (.not. out%scale) out%scales = 1.0_dp
      end if
      nr = size(x)
      allocate(raw(nr, n))
      d = out%derivs - merge(1, 0, out%integral)
      do i = 1, nr
         if (ieee_is_nan(x(i))) then
            raw(i, :) = ieee_value(0.0_dp, ieee_quiet_nan)
            cycle
         end if
         if (out%periodic) then
            row = periodic_values(x(i), out%knots, out%boundary_knots, out%degree, d, family == 2)
         else
            select case (family)
            case (1)
               if (d == -1) then
                  row = b_integrals(x(i), out%knots, out%boundary_knots, out%degree)
               else
                  row = b_values(x(i), out%knots, out%boundary_knots, out%degree, d)
               end if
            case (2)
               row = m_values(x(i), out%knots, out%boundary_knots, out%degree, d)
            case (3)
               row = m_values(x(i), out%knots, out%boundary_knots, out%degree, d - 1)
            case (4)
               if (d == 0) then
                  row = i_integrals(x(i), out%knots, out%boundary_knots, out%degree)
               else
                  row = m_values(x(i), out%knots, out%boundary_knots, out%degree, d - 2)
               end if
               row = row / out%scales
            case (5, 6)
               row = natural_values(x(i), out%knots, out%boundary_knots, d, out%transform)
            case (7)
               row = bernstein_values(x(i), out%boundary_knots, out%degree, d)
            end select
         end if
         raw(i, :) = row
      end do
      out%values = raw(:, first:)
   end function construct_basis

   pure subroutine invalid(out, message)
      type(spline_basis), intent(inout) :: out !! Result whose status and explanatory message are set.
      character(*), intent(in) :: message !! Description of the invalid argument or numerical failure.
      out%status = SPLINES2_INVALID_ARGUMENT
      out%message = message
   end subroutine invalid

   pure function bernstein_values(x, bounds, degree, deriv) result(v)
      real(dp), intent(in) :: x !! Finite coordinate in the closed boundary interval.
      real(dp), intent(in) :: bounds(2) !! Increasing finite boundary knots.
      integer, intent(in) :: degree !! Nonnegative Bernstein polynomial degree.
      integer, intent(in) :: deriv !! Nonnegative derivative order, or -1 for integral from left boundary.
      real(dp) :: v(degree + 1)
      real(dp) :: work(degree + 2), t, saved, term, width, acc
      integer :: p, i, j
      v = 0.0_dp
      if (deriv > degree) return
      width = bounds(2) - bounds(1)
      t = (x - bounds(1)) / width
      p = degree - deriv
      work = 0.0_dp
      work(1) = 1.0_dp
      do i = 1, p
         saved = 0.0_dp
         do j = 1, i
            term = work(j)
            work(j) = saved + (1.0_dp - t) * term
            saved = t * term
         end do
         work(i + 1) = saved
      end do
      if (deriv == -1) then
         acc = 0.0_dp
         do j = degree + 1, 1, -1
            acc = acc + work(j + 1)
            v(j) = width * acc / real(degree + 1, dp)
         end do
      else
         do i = 1, deriv
            saved = 0.0_dp
            do j = 1, p + i
               term = real(p + i, dp) * work(j) / width
               work(j) = saved - term
               saved = term
            end do
            work(p + i + 1) = saved
         end do
         v = work(:degree + 1)
      end if
   end function bernstein_values

   function predict_spline(object, newx) result(out)
      type(spline_basis), intent(in) :: object !! Successful basis object containing reusable knots and settings.
      real(dp), intent(in) :: newx(:) !! New coordinates; NaNs produce missing rows and infinities are invalid.
      type(spline_basis) :: out
      if (object%status /= SPLINES2_OK) then
         out = object
         return
      end if
      out = construct_basis(newx, object%family, knots=object%knots, degree=object%degree, &
         intercept=object%intercept, boundary_knots=object%boundary_knots, periodic=object%periodic, &
         derivs=object%derivs, integral=object%integral, scale=object%scale)
   end function predict_spline

   function derivative_spline(object, derivs) result(out)
      type(spline_basis), intent(in) :: object !! Successful basis object to differentiate at its stored x values.
      integer, intent(in), optional :: derivs !! Additional derivative order, default one, required nonnegative.
      type(spline_basis) :: out
      integer :: d
      d = 1
      if (present(derivs)) d = derivs
      if (object%status /= SPLINES2_OK) then
         out = object
         return
      end if
      if (d < 0) then
         out = object
         call invalid(out, 'Additional derivative order must be nonnegative.')
         return
      end if
      out = construct_basis(object%x, object%family, knots=object%knots, degree=object%degree, &
         intercept=object%intercept, boundary_knots=object%boundary_knots, periodic=object%periodic, &
         derivs=object%derivs + d, integral=object%integral, scale=object%scale)
   end function derivative_spline
end module splines2
