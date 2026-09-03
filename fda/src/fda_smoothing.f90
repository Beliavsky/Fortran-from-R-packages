! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from the computational code of R package fda 6.3.0.
module fda_smoothing
   use r_kinds, only : dp
   use r_linalg, only : solve_spd, symmetric_eigen
   use fda_basis, only : basis_type, basis_penalty, eval_basis
   use fda_fd, only : fd_type, make_fd
   implicit none
   private

   type, public :: smooth_result_type
      type(fd_type) :: fd
      real(dp) :: df = 0.0_dp
      real(dp) :: sse = 0.0_dp
      real(dp), allocatable :: gcv(:)
      real(dp), allocatable :: yhat(:, :)
      real(dp), allocatable :: y2c_map(:, :)
   end type smooth_result_type

   public :: smooth_basis
   public :: lambda_to_df
   public :: df_to_lambda
   public :: lambda_to_gcv
   public :: project_basis

contains

   subroutine smooth_basis(argvals, y, basis, lambda, nderiv, result, info, weights)
      real(dp), intent(in) :: argvals(:) !! Observation argument values, one per row of `y`.
      real(dp), intent(in) :: y(:, :) !! Observed curves with argument values in rows and replications in columns.
      type(basis_type), intent(in) :: basis !! Basis in which the smooth is represented.
      real(dp), intent(in) :: lambda !! Nonnegative roughness-penalty multiplier; zero gives unpenalized least squares.
      integer, intent(in) :: nderiv !! Nonnegative derivative order used to construct the roughness penalty.
      type(smooth_result_type), intent(out) :: result !! Smooth, df, GCV, fitted values, SSE, and observation map.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid shapes, weights, penalty, or linear algebra failure.
      real(dp), intent(in), optional :: weights(:) !! Optional nonnegative observation weights; defaults to one for every row.
      real(dp), allocatable :: amat(:, :), basismat(:, :), bmat(:, :), coefs(:, :), penalty(:, :)
      real(dp), allocatable :: rhs(:, :), wt(:), weighted_basis(:, :)
      integer :: i, n

      info = 0
      n = size(argvals)
      if (size(y, 1) /= n .or. n < 1 .or. size(y, 2) < 1 .or. lambda < 0.0_dp) then
         info = 1
         return
      end if
      allocate(wt(n))
      wt = 1.0_dp
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) then
            info = 2
            return
         end if
         wt = weights
      end if
      call eval_basis(argvals, basis, 0, basismat, info)
      if (info /= 0) return
      if (n < basis%nbasis .and. lambda <= 0.0_dp) then
         info = 3
         return
      end if
      call basis_penalty(basis, nderiv, penalty, info)
      if (info /= 0) return

      allocate(weighted_basis(n, basis%nbasis))
      weighted_basis = basismat * spread(wt, 2, basis%nbasis)
      allocate(bmat(basis%nbasis, basis%nbasis))
      bmat = matmul(transpose(basismat), weighted_basis)
      allocate(amat(basis%nbasis, basis%nbasis))
      amat = 0.5_dp * (bmat + transpose(bmat)) + lambda * penalty
      allocate(rhs(basis%nbasis, size(y, 2)))
      rhs = matmul(transpose(weighted_basis), y)
      call symmetric_solve_with_pinv(amat, rhs, coefs, info)
      if (info /= 0) return
      call make_fd(coefs, basis, result%fd, info)
      if (info /= 0) return

      allocate(result%yhat(n, size(y, 2)))
      result%yhat = matmul(basismat, coefs)
      result%sse = sum((y - result%yhat)**2)
      call smoothing_map(amat, basismat, wt, result%y2c_map, info)
      if (info /= 0) return
      result%df = matrix_trace(matmul(result%y2c_map, basismat))
      allocate(result%gcv(size(y, 2)))
      if (result%df < real(n, dp)) then
         do i = 1, size(y, 2)
            result%gcv(i) = (sum((y(:, i) - result%yhat(:, i))**2) / real(n, dp)) / &
               ((real(n, dp) - result%df) / real(n, dp))**2
         end do
      else
         result%gcv = huge(1.0_dp)
      end if
   end subroutine smooth_basis

   subroutine smoothing_map(amat, basismat, weights, y2c_map, info)
      real(dp), intent(in) :: amat(:, :) !! Penalized normal-equation matrix `B'WB + lambda*R`.
      real(dp), intent(in) :: basismat(:, :) !! Unweighted basis evaluation matrix `B`.
      real(dp), intent(in) :: weights(:) !! Observation weights forming diagonal matrix `W`.
      real(dp), allocatable, intent(out) :: y2c_map(:, :) !! Linear map from observations to basis coefficients.
      integer, intent(out) :: info !! Zero on success; otherwise a symmetric solve failure code.
      real(dp), allocatable :: rhs(:, :)

      allocate(rhs(size(basismat, 2), size(basismat, 1)))
      rhs = transpose(basismat) * spread(weights, 1, size(basismat, 2))
      call symmetric_solve_with_pinv(amat, rhs, y2c_map, info)
   end subroutine smoothing_map

   subroutine lambda_to_df(argvals, basis, nderiv, lambda, df, info, weights)
      real(dp), intent(in) :: argvals(:) !! Observation argument values used to construct the smoother matrix.
      type(basis_type), intent(in) :: basis !! Basis defining both the model matrix and roughness penalty.
      integer, intent(in) :: nderiv !! Nonnegative derivative order used for the penalty matrix.
      real(dp), intent(in) :: lambda !! Nonnegative smoothing parameter whose effective degrees of freedom are requested.
      real(dp), intent(out) :: df !! Trace of the smoothing hat matrix for the supplied `lambda`.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid lambda, weights, or numerical failure.
      real(dp), intent(in), optional :: weights(:) !! Optional nonnegative observation weights; defaults to one.
      real(dp), allocatable :: amat(:, :), basismat(:, :), bmat(:, :), penalty(:, :), solved(:, :), wt(:)

      info = 0
      df = 0.0_dp
      if (lambda < 0.0_dp .or. size(argvals) < 1) then
         info = 1
         return
      end if
      allocate(wt(size(argvals)))
      wt = 1.0_dp
      if (present(weights)) then
         if (size(weights) /= size(argvals) .or. any(weights < 0.0_dp)) then
            info = 2
            return
         end if
         wt = weights
      end if
      call eval_basis(argvals, basis, 0, basismat, info)
      if (info /= 0) return
      allocate(bmat(basis%nbasis, basis%nbasis))
      bmat = matmul(transpose(basismat), basismat * spread(wt, 2, basis%nbasis))
      if (lambda <= 0.0_dp) then
         df = real(min(basis%nbasis, numerical_rank_symmetric(bmat)), dp)
         return
      end if
      call basis_penalty(basis, nderiv, penalty, info)
      if (info /= 0) return
      allocate(amat(basis%nbasis, basis%nbasis))
      amat = 0.5_dp * (bmat + transpose(bmat)) + lambda * penalty
      call symmetric_solve_with_pinv(amat, bmat, solved, info)
      if (info /= 0) return
      df = matrix_trace(solved)
   end subroutine lambda_to_df

   subroutine df_to_lambda(argvals, basis, nderiv, target_df, lambda, info, weights)
      real(dp), intent(in) :: argvals(:) !! Observation argument values defining the smoother matrix.
      type(basis_type), intent(in) :: basis !! Basis defining the model matrix and roughness penalty.
      integer, intent(in) :: nderiv !! Nonnegative derivative order used for the penalty matrix.
      real(dp), intent(in) :: target_df !! Desired effective degrees of freedom between the penalty nullity and unpenalized rank.
      real(dp), intent(out) :: lambda !! Nonnegative smoothing parameter whose effective df approximately matches `target_df`.
      integer, intent(out) :: info !! Zero on success; nonzero when the target is unattainable or a df evaluation fails.
      real(dp), intent(in), optional :: weights(:) !! Optional nonnegative observation weights; defaults to one.
      real(dp) :: df_hi, df_lo, df_mid, hi, lo, mid
      integer :: iter

      info = 0
      lambda = 0.0_dp
      call lambda_to_df(argvals, basis, nderiv, 0.0_dp, df_hi, info, weights)
      if (info /= 0) return
      if (target_df >= df_hi) return
      lo = -16.0_dp
      hi = 16.0_dp
      call lambda_to_df(argvals, basis, nderiv, 10.0_dp**hi, df_lo, info, weights)
      if (info /= 0) return
      if (target_df < df_lo - 1.0e-8_dp) then
         info = 3
         lambda = 10.0_dp**hi
         return
      end if
      do iter = 1, 100
         mid = 0.5_dp * (lo + hi)
         call lambda_to_df(argvals, basis, nderiv, 10.0_dp**mid, df_mid, info, weights)
         if (info /= 0) return
         if (abs(df_mid - target_df) <= 1.0e-8_dp * max(1.0_dp, target_df)) exit
         if (df_mid > target_df) then
            lo = mid
         else
            hi = mid
         end if
      end do
      lambda = 10.0_dp**mid
   end subroutine df_to_lambda

   subroutine lambda_to_gcv(log10lambda, argvals, y, basis, nderiv, gcv, info, weights)
      real(dp), intent(in) :: log10lambda !! Base-10 logarithm of the nonnegative smoothing parameter.
      real(dp), intent(in) :: argvals(:) !! Observation argument values, one per row of `y`.
      real(dp), intent(in) :: y(:, :) !! Observed curves with rows corresponding to `argvals`.
      type(basis_type), intent(in) :: basis !! Basis used for penalized smoothing.
      integer, intent(in) :: nderiv !! Nonnegative derivative order defining the roughness penalty.
      real(dp), allocatable, intent(out) :: gcv(:) !! Per-replication generalized cross-validation values.
      integer, intent(out) :: info !! Zero on success; otherwise the error code returned by `smooth_basis`.
      real(dp), intent(in), optional :: weights(:) !! Optional nonnegative observation weights; defaults to one.
      type(smooth_result_type) :: fit

      call smooth_basis(argvals, y, basis, 10.0_dp**log10lambda, nderiv, fit, info, weights)
      if (info /= 0) then
         allocate(gcv(0))
         return
      end if
      allocate(gcv(size(fit%gcv)))
      gcv = fit%gcv
   end subroutine lambda_to_gcv

   subroutine project_basis(y, argvals, basis, coefs, info, penalize, nderiv)
      real(dp), intent(in) :: y(:, :) !! Sampled curves to project, with argument values in rows.
      real(dp), intent(in) :: argvals(:) !! Argument values matching the rows of `y`.
      type(basis_type), intent(in) :: basis !! Basis onto which the sampled curves are projected.
      real(dp), allocatable, intent(out) :: coefs(:, :) !! Allocated coefficient matrix with basis functions in rows.
      integer, intent(out) :: info !! Zero on success; nonzero for shape, evaluation, penalty, or solve failure.
      logical, intent(in), optional :: penalize !! When true, add the upstream small stabilizing roughness penalty.
      integer, intent(in), optional :: nderiv !! Derivative order for the optional stabilizing penalty; defaults to two.
      real(dp), allocatable :: amat(:, :), basismat(:, :), penalty(:, :), rhs(:, :)
      real(dp) :: lambda
      integer :: deriv_order, i
      logical :: use_penalty

      info = 0
      if (size(y, 1) /= size(argvals)) then
         allocate(coefs(0, 0))
         info = 1
         return
      end if
      use_penalty = .false.
      if (present(penalize)) use_penalty = penalize
      deriv_order = 2
      if (present(nderiv)) deriv_order = nderiv
      call eval_basis(argvals, basis, 0, basismat, info)
      if (info /= 0) then
         allocate(coefs(0, 0))
         return
      end if
      allocate(amat(basis%nbasis, basis%nbasis))
      amat = matmul(transpose(basismat), basismat)
      if (use_penalty) then
         call basis_penalty(basis, deriv_order, penalty, info)
         if (info /= 0) then
            allocate(coefs(0, 0))
            return
         end if
         if (sum([(amat(i, i), i=1, basis%nbasis)]) > 0.0_dp .and. &
             sum([(penalty(i, i), i=1, basis%nbasis)]) > 0.0_dp) then
            lambda = 1.0e-4_dp * sum([(amat(i, i), i=1, basis%nbasis)]) / &
               sum([(penalty(i, i), i=1, basis%nbasis)])
            amat = amat + lambda * penalty
         end if
      end if
      allocate(rhs(basis%nbasis, size(y, 2)))
      rhs = matmul(transpose(basismat), y)
      call symmetric_solve_with_pinv(0.5_dp * (amat + transpose(amat)), rhs, coefs, info)
   end subroutine project_basis

   subroutine symmetric_solve_with_pinv(a, b, x, info)
      real(dp), intent(in) :: a(:, :) !! Symmetric positive-semidefinite matrix to solve or pseudoinvert.
      real(dp), intent(in) :: b(:, :) !! Right-hand-side matrix with the same row count as `a`.
      real(dp), allocatable, intent(out) :: x(:, :) !! Allocated minimum-norm solution matrix.
      integer, intent(out) :: info !! Zero on success; nonzero for incompatible shapes or eigensolver failure.
      real(dp), allocatable :: temp(:, :), values(:), vectors(:, :)
      real(dp) :: threshold
      integer :: i, local_info

      if (size(a, 1) /= size(a, 2) .or. size(b, 1) /= size(a, 1)) then
         allocate(x(0, 0))
         info = 1
         return
      end if
      allocate(x(size(a, 1), size(b, 2)))
      call solve_spd(a, b, x, local_info)
      if (local_info == 0) then
         info = 0
         return
      end if
      call symmetric_eigen(a, values, vectors, info, descending=.true.)
      if (info /= 0) return
      threshold = real(max(1, size(a, 1)), dp) * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(values)))
      allocate(temp(size(a, 1), size(b, 2)))
      temp = matmul(transpose(vectors), b)
      do i = 1, size(values)
         if (values(i) > threshold) then
            temp(i, :) = temp(i, :) / values(i)
         else
            temp(i, :) = 0.0_dp
         end if
      end do
      x = matmul(vectors, temp)
      info = 0
   end subroutine symmetric_solve_with_pinv

   function numerical_rank_symmetric(a) result(rank)
      real(dp), intent(in) :: a(:, :) !! Symmetric matrix whose numerical rank is estimated from its eigenvalues.
      integer :: rank
      real(dp), allocatable :: values(:), vectors(:, :)
      real(dp) :: threshold
      integer :: info

      call symmetric_eigen(a, values, vectors, info, descending=.true.)
      if (info /= 0 .or. size(values) == 0) then
         rank = 0
         return
      end if
      threshold = real(max(1, size(a, 1)), dp) * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(values)))
      rank = count(values > threshold)
   end function numerical_rank_symmetric

   pure real(dp) function matrix_trace(a) result(value)
      real(dp), intent(in) :: a(:, :) !! Matrix whose main diagonal is summed up to its smaller dimension.
      integer :: i

      value = 0.0_dp
      do i = 1, min(size(a, 1), size(a, 2))
         value = value + a(i, i)
      end do
   end function matrix_trace

end module fda_smoothing
