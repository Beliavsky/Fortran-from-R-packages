! Modern Fortran translation of correlation structures in src/famstr.cc.
! Upstream license: GPL (>= 3). See LICENSE, NOTICE.md, and PROVENANCE.md.
! Exact upstream copyright-holder attribution is preserved in NOTICE.md and upstream/DESCRIPTION.
module geepack_correlations
   use r_kinds, only : dp
   use geepack_status, only : GEE_OK, GEE_ERR_ARGUMENT, GEE_ERR_CORRELATION
   implicit none
   private

   integer, parameter, public :: COR_INDEPENDENCE = 1
   integer, parameter, public :: COR_EXCHANGEABLE = 2
   integer, parameter, public :: COR_AR1 = 3
   integer, parameter, public :: COR_UNSTRUCTURED = 4
   integer, parameter, public :: COR_USERDEFINED = 5
   integer, parameter, public :: COR_FIXED = 6

   public :: working_correlation, correlation_rho_derivative
   public :: upper_triangle, pair_products

contains

   pure subroutine working_correlation(rho, wave, corstr, rmat, status)
      real(dp), intent(in) :: rho(:) !! Correlation coefficients or pairwise fixed correlations.
      integer, intent(in) :: wave(:) !! Positive wave/visit indices for observations in one cluster.
      integer, intent(in) :: corstr !! Correlation structure identifier COR_*.
      real(dp), intent(out) :: rmat(:, :) !! Working correlation matrix, shape size(wave) by size(wave).
      integer, intent(out) :: status !! GEE_OK on success or an error code.
      integer :: i
      integer :: j
      integer :: k
      integer :: n
      integer :: maxwave
      integer :: idx
      real(dp) :: lag

      n = size(wave)
      rmat = 0.0_dp
      status = GEE_OK
      if (size(rmat, 1) /= n .or. size(rmat, 2) /= n) then
         status = GEE_ERR_ARGUMENT
         return
      end if
      do i = 1, n
         rmat(i, i) = 1.0_dp
      end do
      if (n <= 1) return

      select case (corstr)
      case (COR_INDEPENDENCE)
         return
      case (COR_EXCHANGEABLE)
         if (size(rho) < 1) then
            status = GEE_ERR_ARGUMENT
            return
         end if
         do i = 1, n - 1
            do j = i + 1, n
               rmat(i, j) = rho(1)
               rmat(j, i) = rho(1)
            end do
         end do
      case (COR_AR1)
         if (size(rho) < 1) then
            status = GEE_ERR_ARGUMENT
            return
         end if
         do i = 1, n - 1
            do j = i + 1, n
               lag = real(abs(wave(j) - wave(i)), dp)
               rmat(i, j) = rho(1) ** lag
               rmat(j, i) = rmat(i, j)
            end do
         end do
      case (COR_UNSTRUCTURED, COR_FIXED)
         maxwave = infer_maxwave(size(rho))
         if (maxwave < maxval(wave) .or. size(rho) < maxwave * (maxwave - 1) / 2) then
            status = GEE_ERR_ARGUMENT
            return
         end if
         do i = 1, n - 1
            do j = i + 1, n
               idx = pair_index(min(wave(i), wave(j)), max(wave(i), wave(j)), maxwave)
               rmat(i, j) = rho(idx)
               rmat(j, i) = rho(idx)
            end do
         end do
      case (COR_USERDEFINED)
         if (size(rho) /= n * (n - 1) / 2) then
            status = GEE_ERR_ARGUMENT
            return
         end if
         k = 0
         do i = 1, n - 1
            do j = i + 1, n
               k = k + 1
               rmat(i, j) = rho(k)
               rmat(j, i) = rho(k)
            end do
         end do
      case default
         status = GEE_ERR_ARGUMENT
         return
      end select
      if (any(abs(rmat) > 1.0_dp + 100.0_dp * epsilon(1.0_dp))) status = GEE_ERR_CORRELATION
   end subroutine working_correlation

   pure subroutine correlation_rho_derivative(rho, wave, corstr, deriv, status)
      real(dp), intent(in) :: rho(:) !! Correlation coefficients at which derivatives are evaluated.
      integer, intent(in) :: wave(:) !! Positive wave/visit indices for one cluster.
      integer, intent(in) :: corstr !! Correlation structure identifier COR_*.
      real(dp), intent(out) :: deriv(:, :) !! Pair-by-rho derivative matrix in upper-triangle row order.
      integer, intent(out) :: status !! GEE_OK on success or an error code.
      integer :: i
      integer :: j
      integer :: k
      integer :: n
      integer :: npair
      integer :: maxwave
      integer :: idx
      integer :: ilag
      real(dp) :: lag

      n = size(wave)
      npair = n * (n - 1) / 2
      deriv = 0.0_dp
      status = GEE_OK
      if (size(deriv, 1) /= npair) then
         status = GEE_ERR_ARGUMENT
         return
      end if
      select case (corstr)
      case (COR_INDEPENDENCE, COR_FIXED)
         return
      case (COR_EXCHANGEABLE)
         if (size(deriv, 2) < 1) then
            status = GEE_ERR_ARGUMENT
            return
         end if
         if (npair > 0) deriv(:, 1) = 1.0_dp
      case (COR_AR1)
         if (size(rho) < 1 .or. size(deriv, 2) < 1) then
            status = GEE_ERR_ARGUMENT
            return
         end if
         k = 0
         do i = 1, n - 1
            do j = i + 1, n
               k = k + 1
               ilag = abs(wave(j) - wave(i))
               lag = real(ilag, dp)
               if (ilag == 1) then
                  deriv(k, 1) = 1.0_dp
               else if (abs(rho(1)) <= tiny(1.0_dp) .and. ilag > 1) then
                  deriv(k, 1) = 0.0_dp
               else
                  deriv(k, 1) = lag * rho(1) ** (ilag - 1)
               end if
            end do
         end do
      case (COR_UNSTRUCTURED)
         maxwave = infer_maxwave(size(rho))
         if (maxwave < maxval(wave) .or. size(deriv, 2) < maxwave * (maxwave - 1) / 2) then
            status = GEE_ERR_ARGUMENT
            return
         end if
         k = 0
         do i = 1, n - 1
            do j = i + 1, n
               k = k + 1
               idx = pair_index(min(wave(i), wave(j)), max(wave(i), wave(j)), maxwave)
               deriv(k, idx) = 1.0_dp
            end do
         end do
      case (COR_USERDEFINED)
         if (size(deriv, 2) /= npair) then
            status = GEE_ERR_ARGUMENT
            return
         end if
         do k = 1, npair
            deriv(k, k) = 1.0_dp
         end do
      case default
         status = GEE_ERR_ARGUMENT
      end select
   end subroutine correlation_rho_derivative

   pure subroutine upper_triangle(matrix, values)
      real(dp), intent(in) :: matrix(:, :) !! Square matrix whose strict upper triangle is extracted.
      real(dp), intent(out) :: values(:) !! Output values in row-major pair order i<j.
      integer :: i
      integer :: j
      integer :: k
      integer :: n

      n = size(matrix, 1)
      k = 0
      do i = 1, n - 1
         do j = i + 1, n
            k = k + 1
            if (k <= size(values)) values(k) = matrix(i, j)
         end do
      end do
   end subroutine upper_triangle

   pure subroutine pair_products(values, products)
      real(dp), intent(in) :: values(:) !! Values whose i<j products are required.
      real(dp), intent(out) :: products(:) !! Pair products in row-major upper-triangle order.
      integer :: i
      integer :: j
      integer :: k

      k = 0
      do i = 1, size(values) - 1
         do j = i + 1, size(values)
            k = k + 1
            if (k <= size(products)) products(k) = values(i) * values(j)
         end do
      end do
   end subroutine pair_products

   pure integer function infer_maxwave(npair) result(n)
      integer, intent(in) :: npair !! Number of unique strict-upper-triangle correlation parameters.

      n = 1
      do while (n * (n - 1) / 2 < npair)
         n = n + 1
      end do
      if (n * (n - 1) / 2 /= npair) n = 0
   end function infer_maxwave

   pure integer function pair_index(i, j, n) result(idx)
      integer, intent(in) :: i !! Smaller one-based wave index.
      integer, intent(in) :: j !! Larger one-based wave index.
      integer, intent(in) :: n !! Maximum wave index represented by the correlation vector.
      integer :: a

      idx = 0
      do a = 1, i - 1
         idx = idx + n - a
      end do
      idx = idx + j - i
   end function pair_index

end module geepack_correlations
