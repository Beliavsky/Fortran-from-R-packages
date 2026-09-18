! SPDX-License-Identifier: GPL-3.0-only
! Scalar and diagonal BEKK kernels translated from dccmidas R/functions.R.
module dccmidas_bekk
   use r_kinds, only : dp
   use dccmidas_types, only : DCCMIDAS_SUCCESS, DCCMIDAS_INVALID_INPUT
   use dccmidas_matrix, only : outer_product, sample_covariance_rows, logdet_quadratic
   implicit none
   private

   public :: sBEKK_loglik, sBEKK_mat_est, dBEKK_loglik, dBEKK_mat_est

contains

   pure subroutine unpack_sbekk(param, k, c, a, b, status)
      real(dp), intent(in) :: param(:) !! Scalar-BEKK parameter vector in the upstream lower-triangle ordering.
      integer, intent(in) :: k !! Number of assets; upstream supports two through five.
      real(dp), intent(out) :: c(:, :) !! Lower-triangular intercept factor with shape `(k, k)`.
      real(dp), intent(out) :: a !! Scalar ARCH coefficient.
      real(dp), intent(out) :: b !! Scalar GARCH coefficient.
      integer, intent(out) :: status !! Zero on success or invalid-input code for incompatible dimensions.
      integer :: i, j, m, p

      c = 0.0_dp
      a = 0.0_dp
      b = 0.0_dp
      m = k * (k + 1) / 2
      if (k < 2 .or. k > 5 .or. size(c, 1) /= k .or. size(c, 2) /= k .or. size(param) /= m + 2) then
         status = DCCMIDAS_INVALID_INPUT
         return
      end if
      p = 0
      do i = 1, k
         do j = 1, i
            p = p + 1
            c(i, j) = param(p)
         end do
      end do
      a = param(m + 1)
      b = param(m + 2)
      status = DCCMIDAS_SUCCESS
   end subroutine unpack_sbekk

   pure subroutine unpack_dbekk(param, k, c, a, b, status)
      real(dp), intent(in) :: param(:) !! Diagonal-BEKK parameter vector in the upstream lower-triangle/diagonal ordering.
      integer, intent(in) :: k !! Number of assets; upstream supports two through five.
      real(dp), intent(out) :: c(:, :) !! Lower-triangular intercept factor with shape `(k, k)`.
      real(dp), intent(out) :: a(:, :) !! Diagonal ARCH matrix with shape `(k, k)`.
      real(dp), intent(out) :: b(:, :) !! Diagonal GARCH matrix with shape `(k, k)`.
      integer, intent(out) :: status !! Zero on success or invalid-input code for incompatible dimensions.
      integer :: i, j, m, p

      c = 0.0_dp
      a = 0.0_dp
      b = 0.0_dp
      m = k * (k + 1) / 2
      if (k < 2 .or. k > 5 .or. size(param) /= m + 2 * k) then
         status = DCCMIDAS_INVALID_INPUT
         return
      end if
      if (size(c, 1) /= k .or. size(c, 2) /= k .or. size(a, 1) /= k .or. size(a, 2) /= k) then
         status = DCCMIDAS_INVALID_INPUT
         return
      end if
      if (size(b, 1) /= k .or. size(b, 2) /= k) then
         status = DCCMIDAS_INVALID_INPUT
         return
      end if
      p = 0
      do i = 1, k
         do j = 1, i
            p = p + 1
            c(i, j) = param(p)
         end do
      end do
      do i = 1, k
         a(i, i) = param(m + i)
         b(i, i) = param(m + k + i)
      end do
      status = DCCMIDAS_SUCCESS
   end subroutine unpack_dbekk

   pure subroutine sBEKK_mat_est(param, ret, h_t, status)
      real(dp), intent(in) :: param(:) !! Scalar-BEKK parameter vector.
      real(dp), intent(in) :: ret(:, :) !! Daily returns with shape `(time, assets)`.
      real(dp), allocatable, intent(out) :: h_t(:, :, :) !! Conditional covariance matrices with shape `(assets, assets, time)`.
      integer, intent(out) :: status !! Zero on success or a dccmidas error code.
      real(dp) :: c(size(ret, 2), size(ret, 2)), w(size(ret, 2), size(ret, 2))
      real(dp) :: a, b, cross(size(ret, 2), size(ret, 2))
      integer :: k, tt, stat

      k = size(ret, 2)
      allocate(h_t(k, k, size(ret, 1)))
      h_t = 0.0_dp
      call unpack_sbekk(param, k, c, a, b, stat)
      if (stat /= DCCMIDAS_SUCCESS .or. size(ret, 1) < 2) then
         status = DCCMIDAS_INVALID_INPUT
         return
      end if
      call sample_covariance_rows(ret, h_t(:, :, 1), stat)
      if (stat /= DCCMIDAS_SUCCESS) then
         status = stat
         return
      end if
      w = matmul(c, transpose(c))
      do tt = 2, size(ret, 1)
         cross = outer_product(ret(tt - 1, :))
         h_t(:, :, tt) = w + a ** 2 * cross + b ** 2 * h_t(:, :, tt - 1)
      end do
      status = DCCMIDAS_SUCCESS
   end subroutine sBEKK_mat_est

   pure subroutine sBEKK_loglik(param, ret, ll, status)
      real(dp), intent(in) :: param(:) !! Scalar-BEKK parameter vector.
      real(dp), intent(in) :: ret(:, :) !! Daily returns with shape `(time, assets)`.
      real(dp), allocatable, intent(out) :: ll(:) !! Per-observation upstream likelihood terms without constants or half-factor.
      integer, intent(out) :: status !! Zero when all recursions are evaluable; a linear-algebra error otherwise.
      real(dp), allocatable :: h_t(:, :, :)
      real(dp) :: logdet, quad
      integer :: tt, stat

      allocate(ll(size(ret, 1)))
      ll = 0.0_dp
      call sBEKK_mat_est(param, ret, h_t, stat)
      if (stat /= DCCMIDAS_SUCCESS) then
         status = stat
         return
      end if
      do tt = 2, size(ret, 1)
         call logdet_quadratic(h_t(:, :, tt), ret(tt, :), logdet, quad, stat)
         if (stat /= DCCMIDAS_SUCCESS) then
            ll(tt) = -huge(1.0_dp)
            status = stat
            return
         end if
         ll(tt) = -(logdet + quad)
      end do
      status = DCCMIDAS_SUCCESS
   end subroutine sBEKK_loglik

   pure subroutine dBEKK_mat_est(param, ret, h_t, status)
      real(dp), intent(in) :: param(:) !! Diagonal-BEKK parameter vector.
      real(dp), intent(in) :: ret(:, :) !! Daily returns with shape `(time, assets)`.
      real(dp), allocatable, intent(out) :: h_t(:, :, :) !! Conditional covariance matrices with shape `(assets, assets, time)`.
      integer, intent(out) :: status !! Zero on success or a dccmidas error code.
      real(dp) :: c(size(ret, 2), size(ret, 2)), a(size(ret, 2), size(ret, 2))
      real(dp) :: b(size(ret, 2), size(ret, 2)), w(size(ret, 2), size(ret, 2))
      real(dp) :: cross(size(ret, 2), size(ret, 2))
      integer :: k, tt, stat

      k = size(ret, 2)
      allocate(h_t(k, k, size(ret, 1)))
      h_t = 0.0_dp
      call unpack_dbekk(param, k, c, a, b, stat)
      if (stat /= DCCMIDAS_SUCCESS .or. size(ret, 1) < 2) then
         status = DCCMIDAS_INVALID_INPUT
         return
      end if
      call sample_covariance_rows(ret, h_t(:, :, 1), stat)
      if (stat /= DCCMIDAS_SUCCESS) then
         status = stat
         return
      end if
      w = matmul(c, transpose(c))
      do tt = 2, size(ret, 1)
         cross = outer_product(ret(tt - 1, :))
         h_t(:, :, tt) = w + matmul(a, matmul(cross, transpose(a))) + &
            matmul(b, matmul(h_t(:, :, tt - 1), transpose(b)))
      end do
      status = DCCMIDAS_SUCCESS
   end subroutine dBEKK_mat_est

   pure subroutine dBEKK_loglik(param, ret, ll, status)
      real(dp), intent(in) :: param(:) !! Diagonal-BEKK parameter vector.
      real(dp), intent(in) :: ret(:, :) !! Daily returns with shape `(time, assets)`.
      real(dp), allocatable, intent(out) :: ll(:) !! Per-observation upstream likelihood terms without constants or half-factor.
      integer, intent(out) :: status !! Zero when all recursions are evaluable; a linear-algebra error otherwise.
      real(dp), allocatable :: h_t(:, :, :)
      real(dp) :: logdet, quad
      integer :: tt, stat

      allocate(ll(size(ret, 1)))
      ll = 0.0_dp
      call dBEKK_mat_est(param, ret, h_t, stat)
      if (stat /= DCCMIDAS_SUCCESS) then
         status = stat
         return
      end if
      do tt = 2, size(ret, 1)
         call logdet_quadratic(h_t(:, :, tt), ret(tt, :), logdet, quad, stat)
         if (stat /= DCCMIDAS_SUCCESS) then
            ll(tt) = -huge(1.0_dp)
            status = stat
            return
         end if
         ll(tt) = -(logdet + quad)
      end do
      status = DCCMIDAS_SUCCESS
   end subroutine dBEKK_loglik

end module dccmidas_bekk
