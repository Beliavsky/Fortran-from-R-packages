! Multivariate normal and Wishart-family draws used by jomo.
! Upstream jomo 2.7-6 by Matteo Quartagno and James Carpenter; License: GPL-2.
! Upstream jomo includes LGPL Wishart/multivariate-normal support derived from
! John Burkardt and earlier Brown/Lovato algorithms. This file independently
! implements the same mathematical distributions for the translated kernels.
! Modern Fortran translation, 2026. Distributed under GPL-2.0-only.
module jomo_distributions
   use jomo_kinds, only : dp
   use jomo_rng, only : rng_state, rng_normal, rng_chisq
   use jomo_linalg, only : chol_lower, inverse_spd, symmetrize
   implicit none
   private

   public :: mvnormal_sample
   public :: matrix_normal_sample
   public :: wishart_sample
   public :: invwishart_sample

contains

   subroutine mvnormal_sample(rng, mean, covariance, draw, info)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for the multivariate-normal draw.
      real(dp), intent(in) :: mean(:) !! Mean vector of length p.
      real(dp), intent(in) :: covariance(:, :) !! Symmetric positive-definite covariance matrix, p by p.
      real(dp), intent(out) :: draw(:) !! Sampled vector of length p.
      integer, intent(out) :: info !! Zero on success; positive when covariance is not positive definite.
      real(dp), allocatable :: l(:, :)
      real(dp), allocatable :: z(:)
      integer :: i

      if (size(covariance, 1) /= size(mean) .or. size(covariance, 2) /= size(mean)) &
         error stop "mvnormal_sample: covariance shape mismatch"
      if (size(draw) /= size(mean)) error stop "mvnormal_sample: draw shape mismatch"
      allocate(l(size(mean), size(mean)), z(size(mean)))
      call chol_lower(covariance, l, info)
      if (info /= 0) then
         draw = mean
         return
      end if
      do i = 1, size(mean)
         z(i) = rng_normal(rng, 0.0_dp, 1.0_dp)
      end do
      draw = mean + matmul(l, z)
   end subroutine mvnormal_sample

   subroutine matrix_normal_sample(rng, mean, row_covariance, col_covariance, draw, info)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for the matrix-normal draw.
      real(dp), intent(in) :: mean(:, :) !! Mean matrix with shape q by p.
      real(dp), intent(in) :: row_covariance(:, :) !! Positive-definite q by q covariance across matrix rows.
      real(dp), intent(in) :: col_covariance(:, :) !! Positive-definite p by p covariance across matrix columns.
      real(dp), intent(out) :: draw(:, :) !! Sampled q by p matrix.
      integer, intent(out) :: info !! Zero on success; positive when either covariance matrix is not positive definite.
      real(dp), allocatable :: lr(:, :)
      real(dp), allocatable :: lc(:, :)
      real(dp), allocatable :: z(:, :)
      integer :: i
      integer :: j
      integer :: info2

      if (size(mean, 1) /= size(row_covariance, 1) .or. size(row_covariance, 1) /= size(row_covariance, 2)) &
         error stop "matrix_normal_sample: row covariance shape mismatch"
      if (size(mean, 2) /= size(col_covariance, 1) .or. size(col_covariance, 1) /= size(col_covariance, 2)) &
         error stop "matrix_normal_sample: column covariance shape mismatch"
      if (any(shape(draw) /= shape(mean))) error stop "matrix_normal_sample: draw shape mismatch"
      allocate(lr(size(row_covariance, 1), size(row_covariance, 2)))
      allocate(lc(size(col_covariance, 1), size(col_covariance, 2)))
      allocate(z(size(mean, 1), size(mean, 2)))
      call chol_lower(row_covariance, lr, info)
      if (info /= 0) then
         draw = mean
         return
      end if
      call chol_lower(col_covariance, lc, info2)
      if (info2 /= 0) then
         info = info2
         draw = mean
         return
      end if
      do j = 1, size(z, 2)
         do i = 1, size(z, 1)
            z(i, j) = rng_normal(rng, 0.0_dp, 1.0_dp)
         end do
      end do
      draw = mean + matmul(matmul(lr, z), transpose(lc))
      info = 0
   end subroutine matrix_normal_sample

   subroutine wishart_sample(rng, df, scale, draw, info)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for the Wishart draw.
      real(dp), intent(in) :: df !! Wishart degrees of freedom, required to exceed p-1.
      real(dp), intent(in) :: scale(:, :) !! Symmetric positive-definite Wishart scale matrix, p by p.
      real(dp), intent(out) :: draw(:, :) !! Sampled p by p Wishart matrix.
      integer, intent(out) :: info !! Zero on success; positive when the scale matrix is not positive definite.
      real(dp), allocatable :: l(:, :)
      real(dp), allocatable :: a(:, :)
      real(dp), allocatable :: b(:, :)
      integer :: i
      integer :: j
      integer :: p

      p = size(scale, 1)
      if (size(scale, 2) /= p) error stop "wishart_sample: scale must be square"
      if (any(shape(draw) /= shape(scale))) error stop "wishart_sample: draw shape mismatch"
      if (df <= real(p - 1, dp)) error stop "wishart_sample: df must exceed p-1"
      allocate(l(p, p), a(p, p), b(p, p))
      call chol_lower(scale, l, info)
      if (info /= 0) then
         draw = 0.0_dp
         return
      end if
      a = 0.0_dp
      do i = 1, p
         a(i, i) = sqrt(rng_chisq(rng, df - real(i - 1, dp)))
         do j = 1, i - 1
            a(i, j) = rng_normal(rng, 0.0_dp, 1.0_dp)
         end do
      end do
      b = matmul(l, a)
      draw = matmul(b, transpose(b))
      call symmetrize(draw)
   end subroutine wishart_sample

   subroutine invwishart_sample(rng, df, scale, draw, info)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for the inverse-Wishart draw.
      real(dp), intent(in) :: df !! Inverse-Wishart degrees of freedom, required to exceed p-1.
      real(dp), intent(in) :: scale(:, :) !! Symmetric positive-definite inverse-Wishart scale matrix, p by p.
      real(dp), intent(out) :: draw(:, :) !! Sampled p by p inverse-Wishart covariance matrix.
      integer, intent(out) :: info !! Zero on success; positive when an SPD factorization fails.
      real(dp), allocatable :: scale_inv(:, :)
      real(dp), allocatable :: precision(:, :)
      integer :: info2

      allocate(scale_inv(size(scale, 1), size(scale, 2)))
      allocate(precision(size(scale, 1), size(scale, 2)))
      call inverse_spd(scale, scale_inv, info)
      if (info /= 0) then
         draw = scale
         return
      end if
      call wishart_sample(rng, df, scale_inv, precision, info)
      if (info /= 0) then
         draw = scale
         return
      end if
      call inverse_spd(precision, draw, info2)
      if (info2 /= 0) then
         info = info2
         draw = scale
      else
         info = 0
      end if
   end subroutine invwishart_sample

end module jomo_distributions
