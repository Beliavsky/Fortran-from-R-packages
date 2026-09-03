! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from the computational code of R package fda 6.3.0.
module fda_numeric
   use r_kinds, only : dp
   use r_linalg, only : cholesky_factor, inverse_matrix, solve_spd, thin_svd
   implicit none
   private

   public :: trapz_mat
   public :: quadset
   public :: polint_matrix
   public :: symsolve
   public :: geigen
   public :: zero_find

contains

   pure subroutine trapz_mat(xmat, ymat, delta, weights, result, info)
      real(dp), intent(in) :: xmat(:, :) !! First sampled-function matrix with observations in rows.
      real(dp), intent(in) :: ymat(:, :) !! Second sampled-function matrix with the same number of rows as `xmat`.
      real(dp), intent(in) :: delta !! Positive equal spacing between successive argument values.
      real(dp), intent(in), optional :: weights(:) !! Optional nonnegative observation weights; endpoints are halved internally.
      real(dp), allocatable, intent(out) :: result(:, :) !! Allocated trapezoidal estimate of `transpose(xmat) W ymat`.
      integer, intent(out) :: info !! Zero on success; nonzero for incompatible shapes, spacing, or weight length.
      real(dp), allocatable :: wt(:), xweighted(:, :)

      info = 0
      if (size(ymat, 1) /= size(xmat, 1)) then
         allocate(result(0, 0))
         info = 1
         return
      end if
      if (delta <= 0.0_dp) then
         allocate(result(0, 0))
         info = 2
         return
      end if
      allocate(wt(size(xmat, 1)))
      wt = 1.0_dp
      if (present(weights)) then
         if (size(weights) /= size(wt)) then
            allocate(result(0, 0))
            info = 3
            return
         end if
         wt = weights
      end if
      if (size(wt) > 0) then
         wt(1) = 0.5_dp * wt(1)
         wt(size(wt)) = 0.5_dp * wt(size(wt))
      end if
      wt = delta * wt
      allocate(xweighted(size(xmat, 1), size(xmat, 2)))
      xweighted = xmat * spread(wt, 2, size(xmat, 2))
      allocate(result(size(xmat, 2), size(ymat, 2)))
      result = matmul(transpose(xweighted), ymat)
   end subroutine trapz_mat

   pure subroutine quadset(breaks, nquad, points, weights, info)
      real(dp), intent(in) :: breaks(:) !! Strictly increasing interval boundaries used by composite Simpson quadrature.
      integer, intent(in) :: nquad !! Requested point count per interval; values below five or even values are increased.
      real(dp), allocatable, intent(out) :: points(:) !! Concatenated Simpson quadrature locations for all intervals.
      real(dp), allocatable, intent(out) :: weights(:) !! Matching Simpson quadrature weights.
      integer, intent(out) :: info !! Zero on success; nonzero when fewer than two or nonincreasing breaks are supplied.
      real(dp) :: h
      integer :: i, j, k, nq

      info = 0
      if (size(breaks) < 2) then
         allocate(points(0), weights(0))
         info = 1
         return
      end if
      do i = 2, size(breaks)
         if (breaks(i) <= breaks(i - 1)) then
            allocate(points(0), weights(0))
            info = 2
            return
         end if
      end do
      nq = max(5, nquad)
      if (mod(nq, 2) == 0) nq = nq + 1
      allocate(points((size(breaks) - 1) * nq), weights((size(breaks) - 1) * nq))
      k = 0
      do i = 1, size(breaks) - 1
         h = (breaks(i + 1) - breaks(i)) / real(nq - 1, dp)
         do j = 1, nq
            k = k + 1
            points(k) = breaks(i) + real(j - 1, dp) * h
            if (j == 1 .or. j == nq) then
               weights(k) = h / 3.0_dp
            else if (mod(j, 2) == 0) then
               weights(k) = 4.0_dp * h / 3.0_dp
            else
               weights(k) = 2.0_dp * h / 3.0_dp
            end if
         end do
      end do
   end subroutine quadset

   pure subroutine polint_matrix(xa, ya, x, y, dy, info)
      real(dp), intent(in) :: xa(:) !! Distinct abscissas of the polynomial interpolation sequence.
      real(dp), intent(in) :: ya(:, :) !! Sequence values with first dimension matching `xa` and arbitrary columns.
      real(dp), intent(in) :: x !! Abscissa at which the interpolating polynomial is evaluated.
      real(dp), allocatable, intent(out) :: y(:) !! Interpolated values, one for each column of `ya`.
      real(dp), allocatable, intent(out) :: dy(:) !! Last Neville correction, usable as a local error indicator.
      integer, intent(out) :: info !! Zero on success; nonzero for shape mismatch or duplicate abscissas.
      real(dp), allocatable :: cs(:, :), ds(:, :), difx(:), w(:)
      real(dp) :: ho, hp
      integer :: i, m, n, ns

      n = size(xa)
      info = 0
      if (size(ya, 1) /= n .or. n < 1) then
         allocate(y(0), dy(0))
         info = 1
         return
      end if
      do i = 1, n
         if (count(abs(xa - xa(i)) <= 16.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(xa(i)))) > 1) then
            allocate(y(0), dy(0))
            info = 2
            return
         end if
      end do
      allocate(cs(n, size(ya, 2)), ds(n, size(ya, 2)), difx(n))
      cs = ya
      ds = ya
      difx = xa - x
      ns = minloc(abs(difx), dim=1)
      allocate(y(size(ya, 2)), dy(size(ya, 2)), w(size(ya, 2)))
      y = ya(ns, :)
      dy = 0.0_dp
      ns = ns - 1
      do m = 1, n - 1
         do i = 1, n - m
            ho = difx(i)
            hp = difx(i + m)
            w = (cs(i + 1, :) - ds(i, :)) / (ho - hp)
            ds(i, :) = hp * w
            cs(i, :) = ho * w
         end do
         if (2 * ns < n - m) then
            dy = cs(ns + 1, :)
         else
            dy = ds(ns, :)
            ns = ns - 1
         end if
         y = y + dy
      end do
   end subroutine polint_matrix

   pure subroutine symsolve(asym, bmat, xmat, info)
      real(dp), intent(in) :: asym(:, :) !! Symmetric positive-definite coefficient matrix.
      real(dp), intent(in) :: bmat(:, :) !! Right-hand-side matrix with the same row count as `asym`.
      real(dp), allocatable, intent(out) :: xmat(:, :) !! Allocated solution matrix satisfying `asym*xmat=bmat`.
      integer, intent(out) :: info !! Zero on success; nonzero for shape, symmetry, or factorization failure.
      real(dp) :: scale

      if (size(asym, 1) /= size(asym, 2) .or. size(bmat, 1) /= size(asym, 1)) then
         allocate(xmat(0, 0))
         info = 1
         return
      end if
      scale = max(1.0_dp, maxval(abs(asym)))
      if (maxval(abs(asym - transpose(asym))) > 1.0e-10_dp * scale) then
         allocate(xmat(0, 0))
         info = 2
         return
      end if
      allocate(xmat(size(asym, 1), size(bmat, 2)))
      call solve_spd(0.5_dp * (asym + transpose(asym)), bmat, xmat, info)
   end subroutine symsolve

   subroutine geigen(amat, bmat, cmat, values, lmat, mmat, info)
      real(dp), intent(in) :: amat(:, :) !! Rectangular cross-matrix in the generalized singular-value problem.
      real(dp), intent(in) :: bmat(:, :) !! Symmetric positive-definite row metric with order `size(amat,1)`.
      real(dp), intent(in) :: cmat(:, :) !! Symmetric positive-definite column metric with order `size(amat,2)`.
      real(dp), allocatable, intent(out) :: values(:) !! Generalized singular values in descending order.
      real(dp), allocatable, intent(out) :: lmat(:, :) !! Left generalized vectors normalized in the `bmat` metric.
      real(dp), allocatable, intent(out) :: mmat(:, :) !! Right generalized vectors normalized in the `cmat` metric.
      integer, intent(out) :: info !! Zero on success; otherwise a shape, Cholesky, inverse, or SVD error code.
      real(dp), allocatable :: bfac(:, :), binv(:, :), cfac(:, :), cinv(:, :)
      real(dp), allocatable :: dmat(:, :), u(:, :), vt(:, :)

      info = 0
      if (size(bmat, 1) /= size(bmat, 2) .or. size(cmat, 1) /= size(cmat, 2) .or. &
          size(amat, 1) /= size(bmat, 1) .or. size(amat, 2) /= size(cmat, 1)) then
         allocate(values(0), lmat(0, 0), mmat(0, 0))
         info = 1
         return
      end if
      call cholesky_factor(bmat, bfac, info, upper=.true.)
      if (info /= 0) then
         allocate(values(0), lmat(0, 0), mmat(0, 0))
         return
      end if
      call cholesky_factor(cmat, cfac, info, upper=.true.)
      if (info /= 0) then
         allocate(values(0), lmat(0, 0), mmat(0, 0))
         return
      end if
      call inverse_matrix(bfac, binv, info)
      if (info /= 0) then
         allocate(values(0), lmat(0, 0), mmat(0, 0))
         return
      end if
      call inverse_matrix(cfac, cinv, info)
      if (info /= 0) then
         allocate(values(0), lmat(0, 0), mmat(0, 0))
         return
      end if
      allocate(dmat(size(amat, 1), size(amat, 2)))
      dmat = matmul(transpose(binv), matmul(amat, cinv))
      call thin_svd(dmat, u, values, vt, info)
      if (info /= 0) then
         allocate(lmat(0, 0), mmat(0, 0))
         return
      end if
      allocate(lmat(size(amat, 1), size(values)), mmat(size(amat, 2), size(values)))
      lmat = matmul(binv, u)
      mmat = matmul(cinv, transpose(vt))
   end subroutine geigen

   pure logical function zero_find(fmat) result(has_zero)
      real(dp), intent(in) :: fmat(:) !! Numeric values whose range is tested for containing zero.

      if (size(fmat) == 0) then
         has_zero = .false.
      else
         has_zero = minval(fmat) <= 0.0_dp .and. maxval(fmat) >= 0.0_dp
      end if
   end function zero_find


end module fda_numeric
