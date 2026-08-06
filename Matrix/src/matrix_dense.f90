! SPDX-License-Identifier: GPL-3.0-only
module matrix_dense
   use matrix_kinds, only : dp
   use matrix_status, only : matrix_success, matrix_err_shape, matrix_err_invalid
   implicit none
   private
   public :: eye, diag_matrix, diag_values, band_matrix, block_diag
   public :: kronecker_product, khatri_rao, symmpart, skewpart
   public :: row_sums, col_sums, row_means, col_means, nnzero_dense
   public :: is_symmetric, is_triangular, is_diagonal
   public :: pack_triangular, unpack_triangular, dense_permute
   public :: trace_matrix, frobenius_norm, one_norm, infinity_norm, max_abs_norm
   public :: scale_rows, scale_cols, crossprod, tcrossprod

contains

   function eye(n) result(a)
      integer, intent(in) :: n
      real(dp), allocatable :: a(:,:)
      integer :: i
      allocate(a(max(0, n), max(0, n)), source=0.0_dp)
      do i = 1, n
         a(i, i) = 1.0_dp
      end do
   end function eye

   function diag_matrix(d, nrow, ncol) result(a)
      real(dp), intent(in) :: d(:)
      integer, intent(in), optional :: nrow, ncol
      real(dp), allocatable :: a(:,:)
      integer :: nr, nc, i, k
      nr = size(d)
      nc = size(d)
      if (present(nrow)) nr = nrow
      if (present(ncol)) nc = ncol
      allocate(a(max(0, nr), max(0, nc)), source=0.0_dp)
      k = min(size(d), min(nr, nc))
      do i = 1, k
         a(i, i) = d(i)
      end do
   end function diag_matrix

   function diag_values(a, k) result(d)
      real(dp), intent(in) :: a(:,:)
      integer, intent(in), optional :: k
      real(dp), allocatable :: d(:)
      integer :: kk, n, i, row0, col0
      kk = 0
      if (present(k)) kk = k
      if (kk >= 0) then
         row0 = 1
         col0 = 1 + kk
         n = min(size(a, 1), size(a, 2) - kk)
      else
         row0 = 1 - kk
         col0 = 1
         n = min(size(a, 1) + kk, size(a, 2))
      end if
      n = max(0, n)
      allocate(d(n))
      do i = 1, n
         d(i) = a(row0 + i - 1, col0 + i - 1)
      end do
   end function diag_values

   function band_matrix(a, lower, upper) result(b)
      real(dp), intent(in) :: a(:,:)
      integer, intent(in) :: lower, upper
      real(dp), allocatable :: b(:,:)
      integer :: i, j
      allocate(b(size(a, 1), size(a, 2)), source=0.0_dp)
      do j = 1, size(a, 2)
         do i = max(1, j - upper), min(size(a, 1), j - lower)
            b(i, j) = a(i, j)
         end do
      end do
   end function band_matrix

   function block_diag(blocks, info) result(a)
      real(dp), intent(in) :: blocks(:,:,:)
      integer, intent(out), optional :: info
      real(dp), allocatable :: a(:,:)
      integer :: n, m, nb, k, r0, c0
      n = size(blocks, 1)
      m = size(blocks, 2)
      nb = size(blocks, 3)
      allocate(a(n * nb, m * nb), source=0.0_dp)
      do k = 1, nb
         r0 = (k - 1) * n
         c0 = (k - 1) * m
         a(r0 + 1:r0 + n, c0 + 1:c0 + m) = blocks(:, :, k)
      end do
      if (present(info)) info = matrix_success
   end function block_diag

   function kronecker_product(a, b) result(c)
      real(dp), intent(in) :: a(:,:), b(:,:)
      real(dp), allocatable :: c(:,:)
      integer :: i, j, nr, nc
      nr = size(b, 1)
      nc = size(b, 2)
      allocate(c(size(a, 1) * nr, size(a, 2) * nc))
      do j = 1, size(a, 2)
         do i = 1, size(a, 1)
            c((i - 1) * nr + 1:i * nr, (j - 1) * nc + 1:j * nc) = a(i, j) * b
         end do
      end do
   end function kronecker_product

   function khatri_rao(a, b, info) result(c)
      real(dp), intent(in) :: a(:,:), b(:,:)
      integer, intent(out), optional :: info
      real(dp), allocatable :: c(:,:)
      integer :: j, istat
      istat = matrix_success
      if (size(a, 2) /= size(b, 2)) then
         allocate(c(0, 0))
         istat = matrix_err_shape
      else
         allocate(c(size(a, 1) * size(b, 1), size(a, 2)))
         do j = 1, size(a, 2)
            c(:, j) = reshape(spread(b(:, j), 2, size(a, 1)) * &
               spread(a(:, j), 1, size(b, 1)), [size(c, 1)])
         end do
      end if
      if (present(info)) info = istat
   end function khatri_rao

   function symmpart(a, info) result(b)
      real(dp), intent(in) :: a(:,:)
      integer, intent(out), optional :: info
      real(dp), allocatable :: b(:,:)
      integer :: istat
      if (size(a, 1) /= size(a, 2)) then
         allocate(b(0, 0))
         istat = matrix_err_shape
      else
         allocate(b(size(a, 1), size(a, 2)))
         b = 0.5_dp * (a + transpose(a))
         istat = matrix_success
      end if
      if (present(info)) info = istat
   end function symmpart

   function skewpart(a, info) result(b)
      real(dp), intent(in) :: a(:,:)
      integer, intent(out), optional :: info
      real(dp), allocatable :: b(:,:)
      integer :: istat
      if (size(a, 1) /= size(a, 2)) then
         allocate(b(0, 0))
         istat = matrix_err_shape
      else
         allocate(b(size(a, 1), size(a, 2)))
         b = 0.5_dp * (a - transpose(a))
         istat = matrix_success
      end if
      if (present(info)) info = istat
   end function skewpart

   function row_sums(a) result(x)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable :: x(:)
      allocate(x(size(a, 1)))
      x = sum(a, dim=2)
   end function row_sums

   function col_sums(a) result(x)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable :: x(:)
      allocate(x(size(a, 2)))
      x = sum(a, dim=1)
   end function col_sums

   function row_means(a) result(x)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable :: x(:)
      allocate(x(size(a, 1)))
      if (size(a, 2) == 0) then
         x = 0.0_dp
      else
         x = sum(a, dim=2) / real(size(a, 2), dp)
      end if
   end function row_means

   function col_means(a) result(x)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable :: x(:)
      allocate(x(size(a, 2)))
      if (size(a, 1) == 0) then
         x = 0.0_dp
      else
         x = sum(a, dim=1) / real(size(a, 1), dp)
      end if
   end function col_means

   integer function nnzero_dense(a, tol) result(n)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(in), optional :: tol
      real(dp) :: eps
      eps = 0.0_dp
      if (present(tol)) eps = max(0.0_dp, tol)
      n = count(abs(a) > eps)
   end function nnzero_dense

   logical function is_symmetric(a, tol) result(ok)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(in), optional :: tol
      real(dp) :: eps
      if (size(a, 1) /= size(a, 2)) then
         ok = .false.
         return
      end if
      eps = 100.0_dp * epsilon(1.0_dp)
      if (present(tol)) eps = max(0.0_dp, tol)
      if (size(a) == 0) then
         ok = .true.
      else
         ok = maxval(abs(a - transpose(a))) <= eps
      end if
   end function is_symmetric

   logical function is_triangular(a, upper, tol) result(ok)
      real(dp), intent(in) :: a(:,:)
      logical, intent(in), optional :: upper
      real(dp), intent(in), optional :: tol
      logical :: use_upper
      real(dp) :: eps
      integer :: i, j
      use_upper = .true.
      if (present(upper)) use_upper = upper
      eps = 100.0_dp * epsilon(1.0_dp)
      if (present(tol)) eps = max(0.0_dp, tol)
      ok = .true.
      if (use_upper) then
         do j = 1, size(a, 2)
            do i = j + 1, size(a, 1)
               if (abs(a(i, j)) > eps) then
                  ok = .false.
                  return
               end if
            end do
         end do
      else
         do j = 1, size(a, 2)
            do i = 1, min(j - 1, size(a, 1))
               if (abs(a(i, j)) > eps) then
                  ok = .false.
                  return
               end if
            end do
         end do
      end if
   end function is_triangular

   logical function is_diagonal(a, tol) result(ok)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(in), optional :: tol
      real(dp) :: eps
      integer :: i, j
      eps = 100.0_dp * epsilon(1.0_dp)
      if (present(tol)) eps = max(0.0_dp, tol)
      ok = .true.
      do j = 1, size(a, 2)
         do i = 1, size(a, 1)
            if (i /= j .and. abs(a(i, j)) > eps) then
               ok = .false.
               return
            end if
         end do
      end do
   end function is_diagonal

   function pack_triangular(a, upper, unit_diag, info) result(x)
      real(dp), intent(in) :: a(:,:)
      logical, intent(in), optional :: upper, unit_diag
      integer, intent(out), optional :: info
      real(dp), allocatable :: x(:)
      logical :: use_upper, omit_diag
      integer :: n, i, j, k, nout, istat
      use_upper = .true.
      omit_diag = .false.
      if (present(upper)) use_upper = upper
      if (present(unit_diag)) omit_diag = unit_diag
      if (size(a, 1) /= size(a, 2)) then
         allocate(x(0))
         istat = matrix_err_shape
      else
         n = size(a, 1)
         nout = n * (n + 1) / 2
         if (omit_diag) nout = n * (n - 1) / 2
         allocate(x(nout))
         k = 0
         if (use_upper) then
            do j = 1, n
               do i = 1, merge(j - 1, j, omit_diag)
                  k = k + 1
                  x(k) = a(i, j)
               end do
            end do
         else
            do j = 1, n
               do i = merge(j + 1, j, omit_diag), n
                  k = k + 1
                  x(k) = a(i, j)
               end do
            end do
         end if
         istat = matrix_success
      end if
      if (present(info)) info = istat
   end function pack_triangular

   function unpack_triangular(x, n, upper, unit_diag, info) result(a)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: n
      logical, intent(in), optional :: upper, unit_diag
      integer, intent(out), optional :: info
      real(dp), allocatable :: a(:,:)
      logical :: use_upper, set_unit
      integer :: expected, i, j, k, istat
      use_upper = .true.
      set_unit = .false.
      if (present(upper)) use_upper = upper
      if (present(unit_diag)) set_unit = unit_diag
      expected = n * (n + 1) / 2
      if (set_unit) expected = n * (n - 1) / 2
      if (n < 0 .or. size(x) /= expected) then
         allocate(a(0, 0))
         istat = matrix_err_invalid
      else
         allocate(a(n, n), source=0.0_dp)
         if (set_unit) then
            do i = 1, n
               a(i, i) = 1.0_dp
            end do
         end if
         k = 0
         if (use_upper) then
            do j = 1, n
               do i = 1, merge(j - 1, j, set_unit)
                  k = k + 1
                  a(i, j) = x(k)
               end do
            end do
         else
            do j = 1, n
               do i = merge(j + 1, j, set_unit), n
                  k = k + 1
                  a(i, j) = x(k)
               end do
            end do
         end if
         istat = matrix_success
      end if
      if (present(info)) info = istat
   end function unpack_triangular

   function dense_permute(a, row_perm, col_perm, info) result(b)
      real(dp), intent(in) :: a(:,:)
      integer, intent(in), optional :: row_perm(:), col_perm(:)
      integer, intent(out), optional :: info
      real(dp), allocatable :: b(:,:)
      integer, allocatable :: rp(:), cp(:)
      integer :: i, istat
      allocate(rp(size(a, 1)), cp(size(a, 2)))
      rp = [(i, i = 1, size(a, 1))]
      cp = [(i, i = 1, size(a, 2))]
      if (present(row_perm)) then
         if (.not. valid_permutation(row_perm, size(a, 1))) then
            allocate(b(0, 0))
            if (present(info)) info = matrix_err_invalid
            return
         end if
         rp = row_perm
      end if
      if (present(col_perm)) then
         if (.not. valid_permutation(col_perm, size(a, 2))) then
            allocate(b(0, 0))
            if (present(info)) info = matrix_err_invalid
            return
         end if
         cp = col_perm
      end if
      allocate(b(size(a, 1), size(a, 2)))
      b = a(rp, cp)
      istat = matrix_success
      if (present(info)) info = istat
   end function dense_permute

   pure logical function valid_permutation(p, n) result(ok)
      integer, intent(in) :: p(:), n
      logical, allocatable :: seen(:)
      integer :: i
      if (size(p) /= n .or. any(p < 1) .or. any(p > n)) then
         ok = .false.
         return
      end if
      allocate(seen(n), source=.false.)
      do i = 1, n
         if (seen(p(i))) then
            ok = .false.
            return
         end if
         seen(p(i)) = .true.
      end do
      ok = .true.
   end function valid_permutation

   real(dp) function trace_matrix(a) result(x)
      real(dp), intent(in) :: a(:,:)
      integer :: i
      x = 0.0_dp
      do i = 1, min(size(a, 1), size(a, 2))
         x = x + a(i, i)
      end do
   end function trace_matrix

   real(dp) function frobenius_norm(a) result(x)
      real(dp), intent(in) :: a(:,:)
      x = sqrt(sum(a * a))
   end function frobenius_norm

   real(dp) function one_norm(a) result(x)
      real(dp), intent(in) :: a(:,:)
      if (size(a, 2) == 0) then
         x = 0.0_dp
      else
         x = maxval(sum(abs(a), dim=1))
      end if
   end function one_norm

   real(dp) function infinity_norm(a) result(x)
      real(dp), intent(in) :: a(:,:)
      if (size(a, 1) == 0) then
         x = 0.0_dp
      else
         x = maxval(sum(abs(a), dim=2))
      end if
   end function infinity_norm

   real(dp) function max_abs_norm(a) result(x)
      real(dp), intent(in) :: a(:,:)
      if (size(a) == 0) then
         x = 0.0_dp
      else
         x = maxval(abs(a))
      end if
   end function max_abs_norm

   function scale_rows(a, s, info) result(b)
      real(dp), intent(in) :: a(:,:), s(:)
      integer, intent(out), optional :: info
      real(dp), allocatable :: b(:,:)
      integer :: i, istat
      if (size(s) /= size(a, 1)) then
         allocate(b(0, 0))
         istat = matrix_err_shape
      else
         allocate(b(size(a, 1), size(a, 2)))
         do i = 1, size(a, 1)
            b(i, :) = s(i) * a(i, :)
         end do
         istat = matrix_success
      end if
      if (present(info)) info = istat
   end function scale_rows

   function scale_cols(a, s, info) result(b)
      real(dp), intent(in) :: a(:,:), s(:)
      integer, intent(out), optional :: info
      real(dp), allocatable :: b(:,:)
      integer :: j, istat
      if (size(s) /= size(a, 2)) then
         allocate(b(0, 0))
         istat = matrix_err_shape
      else
         allocate(b(size(a, 1), size(a, 2)))
         do j = 1, size(a, 2)
            b(:, j) = s(j) * a(:, j)
         end do
         istat = matrix_success
      end if
      if (present(info)) info = istat
   end function scale_cols

   function crossprod(a, b, info) result(c)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(in), optional :: b(:,:)
      integer, intent(out), optional :: info
      real(dp), allocatable :: c(:,:)
      integer :: istat
      if (present(b)) then
         if (size(a, 1) /= size(b, 1)) then
            allocate(c(0, 0))
            istat = matrix_err_shape
         else
            c = matmul(transpose(a), b)
            istat = matrix_success
         end if
      else
         c = matmul(transpose(a), a)
         istat = matrix_success
      end if
      if (present(info)) info = istat
   end function crossprod

   function tcrossprod(a, b, info) result(c)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(in), optional :: b(:,:)
      integer, intent(out), optional :: info
      real(dp), allocatable :: c(:,:)
      integer :: istat
      if (present(b)) then
         if (size(a, 2) /= size(b, 2)) then
            allocate(c(0, 0))
            istat = matrix_err_shape
         else
            c = matmul(a, transpose(b))
            istat = matrix_success
         end if
      else
         c = matmul(a, transpose(a))
         istat = matrix_success
      end if
      if (present(info)) info = istat
   end function tcrossprod

end module matrix_dense
