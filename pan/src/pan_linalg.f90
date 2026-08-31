! SPDX-License-Identifier: GPL-3.0-only
! Derived from computational code in R package pan 2.0.
! Upstream authorship/maintenance: Joseph L. Schafer and Jing Hua Zhao.
! Dense linear-algebra and Gaussian/Wishart kernels for the pan translation.
module pan_linalg
   use pan_kinds, only : dp
   use pan_rng, only : rng_state, rng_normal, rng_gamma
   implicit none
   private

   public :: chol_lower
   public :: spd_inverse
   public :: spd_solve_vec
   public :: spd_solve_mat
   public :: spd_logdet
   public :: mvn_draw
   public :: matrix_normal_draw
   public :: invwishart_draw
   public :: symmetrize
   public :: is_spd

contains

   pure subroutine chol_lower(a, l, info)
      real(dp), intent(in) :: a(:, :) !! Symmetric positive-definite matrix to factor.
      real(dp), intent(out) :: l(:, :) !! Lower-triangular Cholesky factor satisfying A = L L^T.
      integer, intent(out) :: info !! Zero on success; positive pivot index if A is not positive definite.

      integer :: i
      integer :: j
      integer :: k
      integer :: n
      real(dp) :: s

      n = size(a, 1)
      l = 0.0_dp
      info = 0

      if (size(a, 2) /= n .or. size(l, 1) /= n .or. size(l, 2) /= n) then
         info = -1
         return
      end if

      do i = 1, n
         do j = 1, i
            s = a(i, j)
            do k = 1, j - 1
               s = s - l(i, k) * l(j, k)
            end do

            if (i == j) then
               if (s <= 0.0_dp) then
                  info = i
                  return
               end if
               l(i, j) = sqrt(s)
            else
               l(i, j) = s / l(j, j)
            end if
         end do
      end do
   end subroutine chol_lower

   pure subroutine forward_solve(l, b, x, info)
      real(dp), intent(in) :: l(:, :) !! Nonsingular lower-triangular coefficient matrix.
      real(dp), intent(in) :: b(:) !! Right-hand-side vector.
      real(dp), intent(out) :: x(:) !! Solution of L x = b.
      integer, intent(out) :: info !! Zero on success; positive index for a zero diagonal.

      integer :: i
      integer :: j
      integer :: n
      real(dp) :: s

      n = size(b)
      x = 0.0_dp
      info = 0

      do i = 1, n
         if (abs(l(i, i)) <= tiny(1.0_dp)) then
            info = i
            return
         end if
         s = b(i)
         do j = 1, i - 1
            s = s - l(i, j) * x(j)
         end do
         x(i) = s / l(i, i)
      end do
   end subroutine forward_solve

   pure subroutine backward_solve_lt(l, b, x, info)
      real(dp), intent(in) :: l(:, :) !! Nonsingular lower-triangular matrix whose transpose is solved against.
      real(dp), intent(in) :: b(:) !! Right-hand-side vector.
      real(dp), intent(out) :: x(:) !! Solution of L^T x = b.
      integer, intent(out) :: info !! Zero on success; positive index for a zero diagonal.

      integer :: i
      integer :: j
      integer :: n
      real(dp) :: s

      n = size(b)
      x = 0.0_dp
      info = 0

      do i = n, 1, -1
         if (abs(l(i, i)) <= tiny(1.0_dp)) then
            info = i
            return
         end if
         s = b(i)
         do j = i + 1, n
            s = s - l(j, i) * x(j)
         end do
         x(i) = s / l(i, i)
      end do
   end subroutine backward_solve_lt

   pure subroutine spd_solve_vec(a, b, x, info)
      real(dp), intent(in) :: a(:, :) !! Symmetric positive-definite coefficient matrix.
      real(dp), intent(in) :: b(:) !! Right-hand-side vector.
      real(dp), intent(out) :: x(:) !! Solution of A x = b.
      integer, intent(out) :: info !! Zero on success; nonzero if dimensions or factorization fail.

      real(dp) :: l(size(a, 1), size(a, 1))
      real(dp) :: y(size(b))
      integer :: stat

      call chol_lower(a, l, info)
      if (info /= 0) return

      call forward_solve(l, b, y, stat)
      if (stat /= 0) then
         info = stat
         return
      end if

      call backward_solve_lt(l, y, x, stat)
      if (stat /= 0) info = stat
   end subroutine spd_solve_vec

   pure subroutine spd_solve_mat(a, b, x, info)
      real(dp), intent(in) :: a(:, :) !! Symmetric positive-definite coefficient matrix.
      real(dp), intent(in) :: b(:, :) !! Matrix of right-hand sides, one per column.
      real(dp), intent(out) :: x(:, :) !! Solution matrix satisfying A X = B.
      integer, intent(out) :: info !! Zero on success; nonzero if dimensions or factorization fail.

      integer :: j
      integer :: stat

      info = 0
      do j = 1, size(b, 2)
         call spd_solve_vec(a, b(:, j), x(:, j), stat)
         if (stat /= 0) then
            info = stat
            return
         end if
      end do
   end subroutine spd_solve_mat

   pure subroutine spd_inverse(a, ainv, info)
      real(dp), intent(in) :: a(:, :) !! Symmetric positive-definite matrix to invert.
      real(dp), intent(out) :: ainv(:, :) !! Symmetric inverse of A.
      integer, intent(out) :: info !! Zero on success; nonzero if dimensions or factorization fail.

      integer :: i
      integer :: n
      real(dp) :: eye(size(a, 1), size(a, 1))

      n = size(a, 1)
      eye = 0.0_dp
      do i = 1, n
         eye(i, i) = 1.0_dp
      end do

      call spd_solve_mat(a, eye, ainv, info)
      if (info == 0) call symmetrize(ainv)
   end subroutine spd_inverse

   pure subroutine spd_logdet(a, logdet, info)
      real(dp), intent(in) :: a(:, :) !! Symmetric positive-definite matrix.
      real(dp), intent(out) :: logdet !! Natural logarithm of det(A).
      integer, intent(out) :: info !! Zero on success; nonzero if factorization fails.

      integer :: i
      real(dp) :: l(size(a, 1), size(a, 1))

      call chol_lower(a, l, info)
      if (info /= 0) then
         logdet = 0.0_dp
         return
      end if

      logdet = 0.0_dp
      do i = 1, size(a, 1)
         logdet = logdet + 2.0_dp * log(l(i, i))
      end do
   end subroutine spd_logdet

   pure subroutine symmetrize(a)
      real(dp), intent(inout) :: a(:, :) !! Square matrix replaced by its averaged symmetric part.

      integer :: i
      integer :: j
      real(dp) :: s

      do j = 1, size(a, 2)
         do i = j + 1, size(a, 1)
            s = 0.5_dp * (a(i, j) + a(j, i))
            a(i, j) = s
            a(j, i) = s
         end do
      end do
   end subroutine symmetrize

   pure logical function is_spd(a) result(ok)
      real(dp), intent(in) :: a(:, :) !! Symmetric matrix tested for positive definiteness.

      real(dp) :: l(size(a, 1), size(a, 1))
      integer :: info

      call chol_lower(a, l, info)
      ok = info == 0
   end function is_spd

   subroutine mvn_draw(rng, mean, cov, x, info)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for the Gaussian draw.
      real(dp), intent(in) :: mean(:) !! Mean vector of the multivariate normal distribution.
      real(dp), intent(in) :: cov(:, :) !! Symmetric positive-definite covariance matrix.
      real(dp), intent(out) :: x(:) !! Random draw from N(mean, cov).
      integer, intent(out) :: info !! Zero on success; nonzero if the covariance is not positive definite.

      integer :: i
      real(dp) :: l(size(cov, 1), size(cov, 1))
      real(dp) :: z(size(mean))

      call chol_lower(cov, l, info)
      if (info /= 0) then
         x = mean
         return
      end if

      do i = 1, size(mean)
         z(i) = rng_normal(rng)
      end do
      x = mean + matmul(l, z)
   end subroutine mvn_draw

   subroutine matrix_normal_draw(rng, mean, row_cov, col_cov, x, info)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for the matrix-normal draw.
      real(dp), intent(in) :: mean(:, :) !! Mean matrix with dimensions matching the requested draw.
      real(dp), intent(in) :: row_cov(:, :) !! Positive-definite covariance across matrix rows.
      real(dp), intent(in) :: col_cov(:, :) !! Positive-definite covariance across matrix columns.
      real(dp), intent(out) :: x(:, :) !! Matrix-normal random draw.
      integer, intent(out) :: info !! Zero on success; nonzero if either covariance is not positive definite.

      integer :: i
      integer :: j
      integer :: stat
      real(dp) :: lr(size(row_cov, 1), size(row_cov, 1))
      real(dp) :: lc(size(col_cov, 1), size(col_cov, 1))
      real(dp) :: z(size(mean, 1), size(mean, 2))

      call chol_lower(row_cov, lr, info)
      if (info /= 0) return

      call chol_lower(col_cov, lc, stat)
      if (stat /= 0) then
         info = stat
         return
      end if

      do j = 1, size(z, 2)
         do i = 1, size(z, 1)
            z(i, j) = rng_normal(rng)
         end do
      end do

      x = mean + matmul(matmul(lr, z), transpose(lc))
      info = 0
   end subroutine matrix_normal_draw

   subroutine invwishart_draw(rng, scale, df, x, info)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for the inverse-Wishart draw.
      real(dp), intent(in) :: scale(:, :) !! Positive-definite inverse-Wishart scale matrix.
      real(dp), intent(in) :: df !! Degrees of freedom, required to exceed dimension minus one.
      real(dp), intent(out) :: x(:, :) !! Random inverse-Wishart covariance draw.
      integer, intent(out) :: info !! Zero on success; negative for invalid df, positive for factorization failure.

      integer :: i
      integer :: j
      integer :: n
      integer :: stat
      real(dp) :: a(size(scale, 1), size(scale, 1))
      real(dp) :: lscale(size(scale, 1), size(scale, 1))
      real(dp) :: winv(size(scale, 1), size(scale, 1))
      real(dp) :: w(size(scale, 1), size(scale, 1))

      n = size(scale, 1)
      if (df <= real(n - 1, dp)) then
         info = -1
         x = 0.0_dp
         return
      end if

      call chol_lower(scale, lscale, info)
      if (info /= 0) return

      a = 0.0_dp
      do i = 1, n
         a(i, i) = sqrt(2.0_dp * rng_gamma(rng, 0.5_dp * (df - real(i, dp) + 1.0_dp)))
         do j = 1, i - 1
            a(i, j) = rng_normal(rng)
         end do
      end do

      w = matmul(a, transpose(a))
      call spd_inverse(w, winv, stat)
      if (stat /= 0) then
         info = stat
         return
      end if

      x = matmul(matmul(lscale, winv), transpose(lscale))
      call symmetrize(x)
      info = 0
   end subroutine invwishart_draw

end module pan_linalg
