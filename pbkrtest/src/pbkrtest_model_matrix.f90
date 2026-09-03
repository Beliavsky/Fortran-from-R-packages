! SPDX-License-Identifier: GPL-2.0-or-later
module pbkrtest_model_matrix
   use r_kinds, only : dp
   use r_linalg, only : full_svd, numerical_rank
   use pbkrtest_types, only : pbkr_invalid_argument, pbkr_invalid_shape, pbkr_linalg_failure, pbkr_success
   implicit none
   private
   public :: compare_column_space
   public :: force_full_rank
   public :: make_model_matrix
   public :: make_restriction_matrix
   public :: orthogonal_complement

contains

   subroutine compare_column_space(x1, x2, relationship, status)
      real(dp), intent(in) :: x1(:, :) !! First matrix whose column space is compared.
      real(dp), intent(in) :: x2(:, :) !! Second matrix; it must have the same number of rows as `x1`.
      integer, intent(out) :: relationship !! `1` if `C(x2)` is within `C(x1)`, `0` for the reverse, otherwise `-1`.
      integer, intent(out) :: status !! `pbkr_success` on success or a package error code.
      real(dp), allocatable :: both(:, :)
      integer :: info
      integer :: r1
      integer :: r2
      integer :: rboth

      status = pbkr_success
      relationship = -1
      if (size(x1, 1) /= size(x2, 1)) then
         status = pbkr_invalid_shape
         return
      end if
      allocate(both(size(x1, 1), size(x1, 2) + size(x2, 2)))
      both(:, 1:size(x1, 2)) = x1
      both(:, size(x1, 2) + 1:) = x2
      call numerical_rank(x1, r1, info)
      if (info /= 0) then
         status = pbkr_linalg_failure
         return
      end if
      call numerical_rank(x2, r2, info)
      if (info /= 0) then
         status = pbkr_linalg_failure
         return
      end if
      call numerical_rank(both, rboth, info)
      if (info /= 0) then
         status = pbkr_linalg_failure
         return
      end if
      if (rboth == max(r1, r2)) then
         if (r2 < r1) then
            relationship = 1
         else if (r2 > r1) then
            relationship = 0
         end if
      end if
   end subroutine compare_column_space

   subroutine orthogonal_complement(w, basis, status, tolerance)
      real(dp), intent(in) :: w(:, :) !! Matrix whose column space is complemented in the ambient row dimension.
      real(dp), allocatable, intent(out) :: basis(:, :) !! Orthonormal basis for the orthogonal complement of `C(w)`.
      integer, intent(out) :: status !! `pbkr_success` on success or a package error code.
      real(dp), intent(in), optional :: tolerance !! Absolute singular-value threshold used to determine rank.
      real(dp), allocatable :: s(:)
      real(dp), allocatable :: u(:, :)
      real(dp), allocatable :: vt(:, :)
      real(dp) :: threshold
      integer :: i
      integer :: info
      integer :: m
      integer :: n
      integer :: rank_value

      status = pbkr_success
      m = size(w, 1)
      n = size(w, 2)
      if (m <= 0) then
         status = pbkr_invalid_shape
         allocate(basis(0, 0))
         return
      end if
      if (n == 0) then
         allocate(basis(m, m))
         basis = 0.0_dp
         do i = 1, m
            basis(i, i) = 1.0_dp
         end do
         return
      end if
      call full_svd(w, u, s, vt, info)
      if (info /= 0) then
         status = pbkr_linalg_failure
         allocate(basis(0, 0))
         return
      end if
      if (size(s) == 0) then
         rank_value = 0
      else
         threshold = real(max(m, n), dp) * epsilon(1.0_dp) * max(1.0_dp, s(1))
         if (present(tolerance)) threshold = tolerance
         rank_value = count(s > threshold)
      end if
      allocate(basis(m, m - rank_value))
      if (m - rank_value > 0) basis = u(:, rank_value + 1:m)
   end subroutine orthogonal_complement

   subroutine force_full_rank(l, l_full, status, tolerance)
      real(dp), intent(in) :: l(:, :) !! Restriction matrix whose row space must be represented without redundant rows.
      real(dp), allocatable, intent(out) :: l_full(:, :) !! Full-row-rank matrix spanning the same row space as `l`.
      integer, intent(out) :: status !! `pbkr_success` on success or a package error code.
      real(dp), intent(in), optional :: tolerance !! Absolute singular-value threshold used to determine rank.
      real(dp), allocatable :: s(:)
      real(dp), allocatable :: u(:, :)
      real(dp), allocatable :: vt(:, :)
      real(dp) :: threshold
      integer :: info
      integer :: rank_value

      status = pbkr_success
      if (size(l, 1) <= 0 .or. size(l, 2) <= 0) then
         status = pbkr_invalid_shape
         allocate(l_full(0, 0))
         return
      end if
      call full_svd(l, u, s, vt, info)
      if (info /= 0) then
         status = pbkr_linalg_failure
         allocate(l_full(0, 0))
         return
      end if
      threshold = real(max(size(l, 1), size(l, 2)), dp) * epsilon(1.0_dp) * max(1.0_dp, s(1))
      if (present(tolerance)) threshold = tolerance
      rank_value = count(s > threshold)
      if (rank_value == size(l, 1)) then
         allocate(l_full(size(l, 1), size(l, 2)))
         l_full = l
      else
         allocate(l_full(rank_value, size(l, 2)))
         if (rank_value > 0) l_full = vt(1:rank_value, :)
      end if
   end subroutine force_full_rank

   subroutine make_model_matrix(x, l, x2, status)
      real(dp), intent(in) :: x(:, :) !! Large-model fixed-effect matrix, shape `(n,p)`.
      real(dp), intent(in) :: l(:, :) !! Restriction matrix with exactly `p` columns.
      real(dp), allocatable, intent(out) :: x2(:, :) !! Model matrix spanning `{X b : L b = 0}`.
      integer, intent(out) :: status !! `pbkr_success` on success or a package error code.
      real(dp), allocatable :: basis(:, :)

      status = pbkr_success
      if (size(x, 2) /= size(l, 2)) then
         status = pbkr_invalid_shape
         allocate(x2(0, 0))
         return
      end if
      call orthogonal_complement(transpose(l), basis, status)
      if (status /= pbkr_success) then
         allocate(x2(0, 0))
         return
      end if
      allocate(x2(size(x, 1), size(basis, 2)))
      x2 = matmul(x, basis)
   end subroutine make_model_matrix

   subroutine make_restriction_matrix(x, x2, l, status)
      real(dp), intent(in) :: x(:, :) !! Large-model matrix whose column space must contain `C(x2)`.
      real(dp), intent(in) :: x2(:, :) !! Small-model matrix with the same number of rows as `x`.
      real(dp), allocatable, intent(out) :: l(:, :) !! Full-row-rank restriction matrix defining the same nested mean subspace.
      integer, intent(out) :: status !! `pbkr_success` on success or a package error code.
      real(dp), allocatable :: combined(:, :)
      real(dp), allocatable :: complement(:, :)
      real(dp), allocatable :: l_raw(:, :)
      real(dp), allocatable :: residual(:, :)
      real(dp), allocatable :: s(:)
      real(dp), allocatable :: s2(:)
      real(dp), allocatable :: u(:, :)
      real(dp), allocatable :: u2(:, :)
      real(dp), allocatable :: vt(:, :)
      real(dp), allocatable :: vt2(:, :)
      real(dp) :: threshold
      integer :: info
      integer :: rank_combined
      integer :: rank_residual
      integer :: rank_x
      integer :: rank_x2

      status = pbkr_success
      if (size(x, 1) /= size(x2, 1) .or. size(x, 1) <= 0 .or. size(x, 2) <= 0) then
         status = pbkr_invalid_shape
         allocate(l(0, 0))
         return
      end if
      allocate(combined(size(x, 1), size(x, 2) + size(x2, 2)))
      combined(:, 1:size(x2, 2)) = x2
      combined(:, size(x2, 2) + 1:) = x
      call numerical_rank(x, rank_x, info)
      if (info /= 0) then
         status = pbkr_linalg_failure
         allocate(l(0, 0))
         return
      end if
      call numerical_rank(x2, rank_x2, info)
      if (info /= 0) then
         status = pbkr_linalg_failure
         allocate(l(0, 0))
         return
      end if
      call numerical_rank(combined, rank_combined, info)
      if (info /= 0) then
         status = pbkr_linalg_failure
         allocate(l(0, 0))
         return
      end if
      if (rank_combined > rank_x) then
         status = pbkr_invalid_argument
         allocate(l(0, 0))
         return
      end if
      if (rank_x2 >= rank_x) then
         allocate(l(0, size(x, 2)))
         return
      end if

      call full_svd(x2, u2, s2, vt2, info)
      if (info /= 0) then
         status = pbkr_linalg_failure
         allocate(l(0, 0))
         return
      end if
      allocate(residual(size(x, 1), size(x, 2)))
      residual = x
      if (rank_x2 > 0) then
         residual = residual - matmul(u2(:, 1:rank_x2), &
            matmul(transpose(u2(:, 1:rank_x2)), x))
      end if
      call full_svd(residual, u, s, vt, info)
      if (info /= 0) then
         status = pbkr_linalg_failure
         allocate(l(0, 0))
         return
      end if
      threshold = real(max(size(residual, 1), size(residual, 2)), dp) * &
         epsilon(1.0_dp) * max(1.0_dp, s(1))
      rank_residual = count(s > threshold)
      if (rank_residual /= rank_x - rank_x2) then
         status = pbkr_linalg_failure
         allocate(l(0, 0))
         return
      end if
      allocate(complement(size(x, 1), rank_residual))
      complement = u(:, 1:rank_residual)
      allocate(l_raw(rank_residual, size(x, 2)))
      l_raw = matmul(transpose(complement), x)

      call full_svd(l_raw, u, s, vt, info)
      if (info /= 0) then
         status = pbkr_linalg_failure
         allocate(l(0, 0))
         return
      end if
      allocate(l(rank_residual, size(x, 2)))
      l = vt(1:rank_residual, :)
      where (abs(l) < 32.0_dp * epsilon(1.0_dp)) l = 0.0_dp
   end subroutine make_restriction_matrix

end module pbkrtest_model_matrix
