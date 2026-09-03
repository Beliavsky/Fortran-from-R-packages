! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from the computational code of R package fda 6.3.0.
module fda_basis
   use r_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: basis_constant = 1
   integer, parameter, public :: basis_bspline = 2
   integer, parameter, public :: basis_fourier = 3
   integer, parameter, public :: basis_monomial = 4
   integer, parameter, public :: basis_exponential = 5
   integer, parameter, public :: basis_power = 6
   integer, parameter, public :: basis_polygonal = 7

   type, public :: basis_type
      integer :: kind = basis_constant
      integer :: nbasis = 1
      integer :: norder = 1
      real(dp) :: rangeval(2) = [0.0_dp, 1.0_dp]
      real(dp) :: period = 1.0_dp
      real(dp) :: argtrans(2) = [0.0_dp, 1.0_dp]
      real(dp), allocatable :: params(:)
   end type basis_type

   public :: make_constant_basis
   public :: make_bspline_basis
   public :: make_fourier_basis
   public :: make_monomial_basis
   public :: make_exponential_basis
   public :: make_power_basis
   public :: make_polygonal_basis
   public :: eval_basis
   public :: basis_penalty
   public :: basis_gram
   public :: basis_breaks

contains

   pure subroutine make_constant_basis(rangeval, basis, info)
      real(dp), intent(in) :: rangeval(2) !! End points of the basis domain; the second value must exceed the first.
      type(basis_type), intent(out) :: basis !! Constructed single-function constant basis on `rangeval`.
      integer, intent(out) :: info !! Zero on success; nonzero when the supplied range is invalid.

      basis = basis_type()
      info = 0
      if (rangeval(2) <= rangeval(1)) then
         info = 1
         return
      end if
      basis%kind = basis_constant
      basis%nbasis = 1
      basis%norder = 1
      basis%rangeval = rangeval
      basis%period = rangeval(2) - rangeval(1)
   end subroutine make_constant_basis

   pure subroutine make_bspline_basis(breaks, norder, basis, info)
      real(dp), intent(in) :: breaks(:) !! Strictly increasing break points, including both domain end points.
      integer, intent(in) :: norder !! B-spline order, one greater than polynomial degree; valid values are 1 through 20.
      type(basis_type), intent(out) :: basis !! Constructed B-spline basis with endpoint multiplicity equal to `norder`.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid order, break count, or nonincreasing breaks.
      integer :: i

      basis = basis_type()
      info = 0
      if (size(breaks) < 2) then
         info = 1
         return
      end if
      if (norder < 1 .or. norder > 20) then
         info = 2
         return
      end if
      do i = 2, size(breaks)
         if (breaks(i) <= breaks(i - 1)) then
            info = 3
            return
         end if
      end do
      basis%kind = basis_bspline
      basis%norder = norder
      basis%nbasis = size(breaks) + norder - 2
      basis%rangeval = [breaks(1), breaks(size(breaks))]
      basis%period = basis%rangeval(2) - basis%rangeval(1)
      allocate(basis%params(size(breaks)))
      basis%params = breaks
   end subroutine make_bspline_basis

   pure subroutine make_fourier_basis(rangeval, nbasis, period, basis, info)
      real(dp), intent(in) :: rangeval(2) !! End points of the basis domain; the second value must exceed the first.
      integer, intent(in) :: nbasis !! Requested number of Fourier functions; an even positive value is increased by one.
      real(dp), intent(in) :: period !! Positive period used to normalize and evaluate the Fourier system.
      type(basis_type), intent(out) :: basis !! Constructed normalized Fourier basis.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid range, basis count, or period.
      integer :: nuse

      basis = basis_type()
      info = 0
      if (rangeval(2) <= rangeval(1)) then
         info = 1
         return
      end if
      if (nbasis <= 0) then
         info = 2
         return
      end if
      if (period <= 0.0_dp) then
         info = 3
         return
      end if
      nuse = nbasis
      if (mod(nuse, 2) == 0) nuse = nuse + 1
      basis%kind = basis_fourier
      basis%nbasis = nuse
      basis%norder = 1
      basis%rangeval = rangeval
      basis%period = period
   end subroutine make_fourier_basis

   pure subroutine make_monomial_basis(rangeval, exponents, argtrans, basis, info)
      real(dp), intent(in) :: rangeval(2) !! End points of the untransformed argument domain.
      real(dp), intent(in) :: exponents(:) !! Distinct nonnegative integer powers defining the monomial functions.
      real(dp), intent(in) :: argtrans(2) !! Shift and nonzero scale used as `(x-shift)/scale` before exponentiation.
      type(basis_type), intent(out) :: basis !! Constructed transformed monomial basis.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid range, scale, exponent, or duplicate exponent.
      integer :: i, j

      basis = basis_type()
      info = 0
      if (rangeval(2) <= rangeval(1) .or. size(exponents) < 1) then
         info = 1
         return
      end if
      if (abs(argtrans(2)) <= tiny(1.0_dp)) then
         info = 2
         return
      end if
      do i = 1, size(exponents)
         if (exponents(i) < 0.0_dp .or. abs(exponents(i) - nint(exponents(i))) > 16.0_dp * epsilon(1.0_dp)) then
            info = 3
            return
         end if
         do j = 1, i - 1
            if (abs(exponents(i) - exponents(j)) <= 16.0_dp * epsilon(1.0_dp)) then
               info = 4
               return
            end if
         end do
      end do
      basis%kind = basis_monomial
      basis%nbasis = size(exponents)
      basis%norder = 1
      basis%rangeval = rangeval
      basis%period = rangeval(2) - rangeval(1)
      basis%argtrans = argtrans
      allocate(basis%params(size(exponents)))
      basis%params = exponents
   end subroutine make_monomial_basis

   pure subroutine make_exponential_basis(rangeval, rates, basis, info)
      real(dp), intent(in) :: rangeval(2) !! End points of the argument domain.
      real(dp), intent(in) :: rates(:) !! Exponential rate constants multiplying the argument in `exp(rate*x)`.
      type(basis_type), intent(out) :: basis !! Constructed exponential basis.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid range or an empty rate vector.

      basis = basis_type()
      info = 0
      if (rangeval(2) <= rangeval(1) .or. size(rates) < 1) then
         info = 1
         return
      end if
      basis%kind = basis_exponential
      basis%nbasis = size(rates)
      basis%norder = 1
      basis%rangeval = rangeval
      basis%period = rangeval(2) - rangeval(1)
      allocate(basis%params(size(rates)))
      basis%params = rates
   end subroutine make_exponential_basis

   pure subroutine make_power_basis(rangeval, exponents, basis, info)
      real(dp), intent(in) :: rangeval(2) !! End points of the argument domain on which powers are evaluated.
      real(dp), intent(in) :: exponents(:) !! Real exponents defining the functions `x**exponent`.
      type(basis_type), intent(out) :: basis !! Constructed power basis.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid range or an empty exponent vector.

      basis = basis_type()
      info = 0
      if (rangeval(2) <= rangeval(1) .or. size(exponents) < 1) then
         info = 1
         return
      end if
      basis%kind = basis_power
      basis%nbasis = size(exponents)
      basis%norder = 1
      basis%rangeval = rangeval
      basis%period = rangeval(2) - rangeval(1)
      allocate(basis%params(size(exponents)))
      basis%params = exponents
   end subroutine make_power_basis

   pure subroutine make_polygonal_basis(breaks, basis, info)
      real(dp), intent(in) :: breaks(:) !! Strictly increasing knot locations at which polygonal basis functions peak.
      type(basis_type), intent(out) :: basis !! Constructed piecewise-linear interpolating basis.
      integer, intent(out) :: info !! Zero on success; nonzero for too few or nonincreasing knots.
      integer :: i

      basis = basis_type()
      info = 0
      if (size(breaks) < 2) then
         info = 1
         return
      end if
      do i = 2, size(breaks)
         if (breaks(i) <= breaks(i - 1)) then
            info = 2
            return
         end if
      end do
      basis%kind = basis_polygonal
      basis%nbasis = size(breaks)
      basis%norder = 2
      basis%rangeval = [breaks(1), breaks(size(breaks))]
      basis%period = basis%rangeval(2) - basis%rangeval(1)
      allocate(basis%params(size(breaks)))
      basis%params = breaks
   end subroutine make_polygonal_basis

   pure subroutine eval_basis(x, basis, nderiv, values, info)
      real(dp), intent(in) :: x(:) !! Argument values lying inside the basis domain, including its end points.
      type(basis_type), intent(in) :: basis !! Basis specification to evaluate.
      integer, intent(in) :: nderiv !! Nonnegative derivative order to evaluate.
      real(dp), allocatable, intent(out) :: values(:, :) !! Allocated matrix with `size(x)` rows and `basis%nbasis` columns.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid derivative order, domain values, or basis data.

      info = 0
      if (nderiv < 0) then
         allocate(values(0, 0))
         info = 1
         return
      end if
      if (any(x < basis%rangeval(1) - 64.0_dp * epsilon(1.0_dp)) .or. &
          any(x > basis%rangeval(2) + 64.0_dp * epsilon(1.0_dp))) then
         allocate(values(0, 0))
         info = 2
         return
      end if

      select case (basis%kind)
      case (basis_constant)
         call constant_matrix(x, nderiv, values)
      case (basis_bspline)
         call bspline_matrix(x, basis%params, basis%norder, nderiv, values, info)
      case (basis_fourier)
         call fourier_matrix(x, basis%nbasis, basis%period, nderiv, values)
      case (basis_monomial)
         call monomial_matrix(x, basis%params, basis%argtrans, nderiv, values)
      case (basis_exponential)
         call exponential_matrix(x, basis%params, nderiv, values)
      case (basis_power)
         call power_matrix(x, basis%params, nderiv, values, info)
      case (basis_polygonal)
         call polygonal_matrix(x, basis%params, nderiv, values, info)
      case default
         allocate(values(0, 0))
         info = 3
      end select
   end subroutine eval_basis

   pure subroutine constant_matrix(x, nderiv, values)
      real(dp), intent(in) :: x(:) !! Argument values; only their count matters for the constant basis.
      integer, intent(in) :: nderiv !! Derivative order; derivatives above zero vanish.
      real(dp), allocatable, intent(out) :: values(:, :) !! Allocated constant-basis evaluation matrix.

      allocate(values(size(x), 1))
      if (nderiv == 0) then
         values(:, 1) = 1.0_dp
      else
         values = 0.0_dp
      end if
   end subroutine constant_matrix

   pure subroutine fourier_matrix(x, nbasis, period, nderiv, values)
      real(dp), intent(in) :: x(:) !! Argument values at which the normalized Fourier functions are evaluated.
      integer, intent(in) :: nbasis !! Positive odd number of Fourier basis functions.
      real(dp), intent(in) :: period !! Positive Fourier period used for frequency and normalization.
      integer, intent(in) :: nderiv !! Nonnegative derivative order.
      real(dp), allocatable, intent(out) :: values(:, :) !! Allocated Fourier evaluation matrix.
      real(dp) :: angle, fac, omega, phase, scale
      integer :: i, j, k

      allocate(values(size(x), nbasis))
      values = 0.0_dp
      omega = 2.0_dp * acos(-1.0_dp) / period
      if (nderiv == 0) then
         values(:, 1) = 1.0_dp / sqrt(period)
      end if
      scale = sqrt(2.0_dp / period)
      do j = 2, nbasis - 1, 2
         k = j / 2
         fac = (real(k, dp) * omega)**nderiv
         phase = real(nderiv, dp) * acos(-1.0_dp) / 2.0_dp
         do i = 1, size(x)
            angle = real(k, dp) * omega * x(i) + phase
            values(i, j) = scale * fac * sin(angle)
            values(i, j + 1) = scale * fac * cos(angle)
         end do
      end do
   end subroutine fourier_matrix

   pure subroutine monomial_matrix(x, exponents, argtrans, nderiv, values)
      real(dp), intent(in) :: x(:) !! Untransformed argument values.
      real(dp), intent(in) :: exponents(:) !! Nonnegative integer exponents defining the basis functions.
      real(dp), intent(in) :: argtrans(2) !! Shift and scale applied to the argument before taking powers.
      integer, intent(in) :: nderiv !! Nonnegative derivative order with respect to the original argument.
      real(dp), allocatable, intent(out) :: values(:, :) !! Allocated monomial derivative matrix.
      real(dp) :: fac, z
      integer :: degree, i, j, k

      allocate(values(size(x), size(exponents)))
      values = 0.0_dp
      do j = 1, size(exponents)
         degree = nint(exponents(j))
         if (nderiv > degree) cycle
         fac = 1.0_dp
         do k = 0, nderiv - 1
            fac = fac * real(degree - k, dp) / argtrans(2)
         end do
         do i = 1, size(x)
            z = (x(i) - argtrans(1)) / argtrans(2)
            values(i, j) = fac * z**(degree - nderiv)
         end do
      end do
   end subroutine monomial_matrix

   pure subroutine exponential_matrix(x, rates, nderiv, values)
      real(dp), intent(in) :: x(:) !! Argument values at which exponentials are evaluated.
      real(dp), intent(in) :: rates(:) !! Exponential rate constants.
      integer, intent(in) :: nderiv !! Nonnegative derivative order.
      real(dp), allocatable, intent(out) :: values(:, :) !! Allocated exponential derivative matrix.
      integer :: i, j

      allocate(values(size(x), size(rates)))
      do j = 1, size(rates)
         do i = 1, size(x)
            values(i, j) = rates(j)**nderiv * exp(rates(j) * x(i))
         end do
      end do
   end subroutine exponential_matrix

   pure subroutine power_matrix(x, exponents, nderiv, values, info)
      real(dp), intent(in) :: x(:) !! Argument values used as bases of the real powers.
      real(dp), intent(in) :: exponents(:) !! Real exponents defining the power basis.
      integer, intent(in) :: nderiv !! Nonnegative derivative order.
      real(dp), allocatable, intent(out) :: values(:, :) !! Allocated power-basis derivative matrix.
      integer, intent(out) :: info !! Zero on success; nonzero when a required negative power is evaluated at a nonpositive value.
      real(dp) :: fac, power
      integer :: i, j, k

      allocate(values(size(x), size(exponents)))
      values = 0.0_dp
      info = 0
      do j = 1, size(exponents)
         power = exponents(j) - real(nderiv, dp)
         fac = 1.0_dp
         do k = 0, nderiv - 1
            fac = fac * (exponents(j) - real(k, dp))
         end do
         if (abs(fac) <= tiny(1.0_dp)) cycle
         do i = 1, size(x)
            if (x(i) <= 0.0_dp .and. power < 0.0_dp) then
               info = 1
               return
            end if
            values(i, j) = fac * x(i)**power
         end do
      end do
   end subroutine power_matrix

   pure subroutine polygonal_matrix(x, breaks, nderiv, values, info)
      real(dp), intent(in) :: x(:) !! Argument values inside the polygonal knot range.
      real(dp), intent(in) :: breaks(:) !! Strictly increasing polygonal knots.
      integer, intent(in) :: nderiv !! Derivative order; only orders zero and one are nonzero away from knots.
      real(dp), allocatable, intent(out) :: values(:, :) !! Allocated polygonal basis evaluation matrix.
      integer, intent(out) :: info !! Zero on success; nonzero for a negative derivative order.
      real(dp) :: h
      integer :: i, j

      allocate(values(size(x), size(breaks)))
      values = 0.0_dp
      info = 0
      if (nderiv < 0) then
         info = 1
         return
      end if
      if (nderiv >= 2) return
      do i = 1, size(x)
         if (x(i) >= breaks(size(breaks))) then
            j = size(breaks) - 1
         else
            j = 1
            do while (j < size(breaks) - 1 .and. x(i) >= breaks(j + 1))
               j = j + 1
            end do
         end if
         h = breaks(j + 1) - breaks(j)
         if (nderiv == 0) then
            values(i, j) = (breaks(j + 1) - x(i)) / h
            values(i, j + 1) = (x(i) - breaks(j)) / h
         else
            values(i, j) = -1.0_dp / h
            values(i, j + 1) = 1.0_dp / h
         end if
      end do
   end subroutine polygonal_matrix

   pure subroutine bspline_matrix(x, breaks, norder, nderiv, values, info)
      real(dp), intent(in) :: x(:) !! Argument values lying inside the break-point range.
      real(dp), intent(in) :: breaks(:) !! Strictly increasing B-spline break points including end points.
      integer, intent(in) :: norder !! B-spline order, one greater than polynomial degree.
      integer, intent(in) :: nderiv !! Nonnegative derivative order.
      real(dp), allocatable, intent(out) :: values(:, :) !! Allocated B-spline derivative matrix.
      integer, intent(out) :: info !! Zero on success; nonzero when the spline definition is invalid.
      real(dp), allocatable :: knots(:), work(:, :, :)
      real(dp) :: den1, den2, xx
      integer :: d, i, ix, k, nbase, nknots

      nbase = size(breaks) + norder - 2
      allocate(values(size(x), nbase))
      values = 0.0_dp
      info = 0
      if (nderiv >= norder) return
      if (norder < 1 .or. size(breaks) < 2) then
         info = 1
         return
      end if

      nknots = size(breaks) + 2 * (norder - 1)
      allocate(knots(nknots))
      knots(1:norder - 1) = breaks(1)
      knots(norder:norder + size(breaks) - 1) = breaks
      knots(norder + size(breaks):nknots) = breaks(size(breaks))
      allocate(work(nknots - 1, norder, nderiv + 1))

      do ix = 1, size(x)
         work = 0.0_dp
         xx = x(ix)
         if (xx >= breaks(size(breaks))) xx = nearest(breaks(size(breaks)), -1.0_dp)
         do i = 1, nknots - 1
            if (xx >= knots(i) .and. xx < knots(i + 1)) work(i, 1, 1) = 1.0_dp
         end do
         do k = 2, norder
            do i = 1, nknots - k
               den1 = knots(i + k - 1) - knots(i)
               den2 = knots(i + k) - knots(i + 1)
               if (den1 > 0.0_dp) work(i, k, 1) = work(i, k, 1) + &
                  (xx - knots(i)) * work(i, k - 1, 1) / den1
               if (den2 > 0.0_dp) work(i, k, 1) = work(i, k, 1) + &
                  (knots(i + k) - xx) * work(i + 1, k - 1, 1) / den2
               do d = 1, min(nderiv, k - 1)
                  if (den1 > 0.0_dp) work(i, k, d + 1) = work(i, k, d + 1) + &
                     real(k - 1, dp) * work(i, k - 1, d) / den1
                  if (den2 > 0.0_dp) work(i, k, d + 1) = work(i, k, d + 1) - &
                     real(k - 1, dp) * work(i + 1, k - 1, d) / den2
               end do
            end do
         end do
         values(ix, :) = work(1:nbase, norder, nderiv + 1)
      end do
   end subroutine bspline_matrix

   pure subroutine basis_breaks(basis, breaks)
      type(basis_type), intent(in) :: basis !! Basis whose natural integration break points are requested.
      real(dp), allocatable, intent(out) :: breaks(:) !! Sorted integration break points, including both domain ends.

      if (basis%kind == basis_bspline .or. basis%kind == basis_polygonal) then
         allocate(breaks(size(basis%params)))
         breaks = basis%params
      else
         allocate(breaks(2))
         breaks = basis%rangeval
      end if
   end subroutine basis_breaks

   subroutine basis_penalty(basis, nderiv, penalty, info, nquad)
      type(basis_type), intent(in) :: basis !! Basis whose derivative inner-product penalty matrix is required.
      integer, intent(in) :: nderiv !! Nonnegative derivative order used in the roughness penalty.
      real(dp), allocatable, intent(out) :: penalty(:, :) !! Allocated symmetric penalty matrix with shape `(nbasis, nbasis)`.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid derivative order or basis evaluation failure.
      integer, intent(in), optional :: nquad !! Odd Simpson point count per natural interval; defaults to  nine points.
      real(dp), allocatable :: breaks(:), pts(:), vals(:, :), wts(:)
      integer :: nq

      nq = 9
      if (present(nquad)) nq = max(5, nquad)
      if (mod(nq, 2) == 0) nq = nq + 1
      call basis_breaks(basis, breaks)
      call simpson_points(breaks, nq, pts, wts)
      call eval_basis(pts, basis, nderiv, vals, info)
      if (info /= 0) then
         allocate(penalty(0, 0))
         return
      end if
      allocate(penalty(basis%nbasis, basis%nbasis))
      penalty = matmul(transpose(vals), vals * spread(wts, 2, basis%nbasis))
      penalty = 0.5_dp * (penalty + transpose(penalty))
   end subroutine basis_penalty

   subroutine basis_gram(basis, gram, info, nquad)
      type(basis_type), intent(in) :: basis !! Basis whose ordinary L2 Gram matrix is required.
      real(dp), allocatable, intent(out) :: gram(:, :) !! Allocated matrix of pairwise basis inner products.
      integer, intent(out) :: info !! Zero on success; otherwise a basis-evaluation error code.
      integer, intent(in), optional :: nquad !! Odd Simpson point count per natural interval; defaults to nine.

      call basis_penalty(basis, 0, gram, info, nquad)
   end subroutine basis_gram

   pure subroutine simpson_points(breaks, nquad, points, weights)
      real(dp), intent(in) :: breaks(:) !! Strictly increasing interval boundaries for composite Simpson integration.
      integer, intent(in) :: nquad !! Odd number of Simpson points placed on every interval.
      real(dp), allocatable, intent(out) :: points(:) !! Concatenated quadrature locations, including interval boundaries.
      real(dp), allocatable, intent(out) :: weights(:) !! Matching Simpson integration weights.
      real(dp) :: h
      integer :: i, j, k, nint

      nint = size(breaks) - 1
      allocate(points(nint * nquad), weights(nint * nquad))
      k = 0
      do i = 1, nint
         h = (breaks(i + 1) - breaks(i)) / real(nquad - 1, dp)
         do j = 1, nquad
            k = k + 1
            points(k) = breaks(i) + real(j - 1, dp) * h
            if (j == 1 .or. j == nquad) then
               weights(k) = h / 3.0_dp
            else if (mod(j, 2) == 0) then
               weights(k) = 4.0_dp * h / 3.0_dp
            else
               weights(k) = 2.0_dp * h / 3.0_dp
            end if
         end do
      end do
   end subroutine simpson_points

end module fda_basis
