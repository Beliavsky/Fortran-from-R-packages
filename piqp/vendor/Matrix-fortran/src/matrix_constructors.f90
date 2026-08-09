! SPDX-License-Identifier: GPL-3.0-only
module matrix_constructors
   use matrix_kinds, only : dp
   use matrix_status, only : matrix_success, matrix_err_invalid
   use matrix_sparse, only : csr_matrix, csr_from_triplet
   implicit none
   private
   public :: toeplitz_matrix, hilbert_matrix, permutation_matrix
   public :: sparse_identity, sparse_diagonal, random_sparse_matrix
   public :: companion_matrix

contains

   function toeplitz_matrix(first_col, first_row, info) result(a)
      real(dp), intent(in) :: first_col(:)
      real(dp), intent(in), optional :: first_row(:)
      integer, intent(out), optional :: info
      real(dp), allocatable :: a(:,:), row(:)
      integer :: i, j, nrow, ncol, istat
      nrow = size(first_col)
      if (present(first_row)) then
         ncol = size(first_row)
         row = first_row
         if (nrow > 0 .and. ncol > 0) then
            if (abs(first_col(1) - row(1)) > 100.0_dp * epsilon(1.0_dp)) then
               allocate(a(0, 0))
               istat = matrix_err_invalid
               if (present(info)) info = istat
               return
            end if
         end if
      else
         ncol = nrow
         allocate(row(ncol))
         if (ncol > 0) row(1) = first_col(1)
         do j = 2, ncol
            row(j) = first_col(j)
         end do
      end if
      allocate(a(nrow, ncol))
      do j = 1, ncol
         do i = 1, nrow
            if (i >= j) then
               a(i, j) = first_col(i - j + 1)
            else
               a(i, j) = row(j - i + 1)
            end if
         end do
      end do
      istat = matrix_success
      if (present(info)) info = istat
   end function toeplitz_matrix

   function hilbert_matrix(n, m) result(a)
      integer, intent(in) :: n
      integer, intent(in), optional :: m
      real(dp), allocatable :: a(:,:)
      integer :: nr, nc, i, j
      nr = max(0, n)
      nc = nr
      if (present(m)) nc = max(0, m)
      allocate(a(nr, nc))
      do j = 1, nc
         do i = 1, nr
            a(i, j) = 1.0_dp / real(i + j - 1, dp)
         end do
      end do
   end function hilbert_matrix

   function permutation_matrix(p, info) result(a)
      integer, intent(in) :: p(:)
      integer, intent(out), optional :: info
      real(dp), allocatable :: a(:,:)
      logical, allocatable :: seen(:)
      integer :: i, istat, n
      n = size(p)
      if (any(p < 1) .or. any(p > n)) then
         allocate(a(0, 0))
         istat = matrix_err_invalid
      else
         allocate(seen(n), source=.false.)
         istat = matrix_success
         do i = 1, n
            if (seen(p(i))) then
               istat = matrix_err_invalid
               exit
            end if
            seen(p(i)) = .true.
         end do
         if (istat == matrix_success) then
            allocate(a(n, n), source=0.0_dp)
            do i = 1, n
               a(i, p(i)) = 1.0_dp
            end do
         else
            allocate(a(0, 0))
         end if
      end if
      if (present(info)) info = istat
   end function permutation_matrix

   subroutine sparse_identity(n, a)
      integer, intent(in) :: n
      type(csr_matrix), intent(out) :: a
      integer, allocatable :: idx(:)
      real(dp), allocatable :: values(:)
      integer :: i, info
      allocate(idx(max(0, n)))
      allocate(values(max(0, n)), source=1.0_dp)
      idx = [(i, i = 1, n)]
      call csr_from_triplet(max(0, n), max(0, n), idx, idx, values, a, info)
   end subroutine sparse_identity

   subroutine sparse_diagonal(d, a, nrow, ncol)
      real(dp), intent(in) :: d(:)
      type(csr_matrix), intent(out) :: a
      integer, intent(in), optional :: nrow, ncol
      integer, allocatable :: idx(:)
      integer :: nr, nc, n, i, info
      nr = size(d)
      nc = size(d)
      if (present(nrow)) nr = nrow
      if (present(ncol)) nc = ncol
      n = min(size(d), min(nr, nc))
      allocate(idx(max(0, n)))
      idx = [(i, i = 1, n)]
      call csr_from_triplet(nr, nc, idx, idx, d(:n), a, info)
   end subroutine sparse_diagonal

   subroutine random_sparse_matrix(nrow, ncol, density, a, info, symmetric)
      integer, intent(in) :: nrow, ncol
      real(dp), intent(in) :: density
      type(csr_matrix), intent(out) :: a
      integer, intent(out) :: info
      logical, intent(in), optional :: symmetric
      integer, allocatable :: rows(:), cols(:)
      real(dp), allocatable :: vals(:)
      real(dp) :: u, v
      integer :: i, j, n, cap
      logical :: sym
      sym = .false.
      if (present(symmetric)) sym = symmetric
      if (nrow < 0 .or. ncol < 0 .or. density < 0.0_dp .or. density > 1.0_dp .or. &
         (sym .and. nrow /= ncol)) then
         call csr_from_triplet(0, 0, [integer ::], [integer ::], [real(dp) ::], a, info)
         info = matrix_err_invalid
         return
      end if
      cap = max(1, ceiling(density * real(nrow * ncol, dp) * merge(2.0_dp, 1.2_dp, sym)))
      allocate(rows(cap), cols(cap), vals(cap))
      n = 0
      do j = 1, ncol
         do i = 1, nrow
            if (sym .and. i < j) cycle
            call random_number(u)
            if (u <= density) then
               call random_number(v)
               if (n == cap) call grow(rows, cols, vals, cap)
               n = n + 1
               rows(n) = i
               cols(n) = j
               vals(n) = 2.0_dp * v - 1.0_dp
               if (sym .and. i /= j) then
                  if (n == cap) call grow(rows, cols, vals, cap)
                  n = n + 1
                  rows(n) = j
                  cols(n) = i
                  vals(n) = vals(n - 1)
               end if
            end if
         end do
      end do
      call csr_from_triplet(nrow, ncol, rows(:n), cols(:n), vals(:n), a, info)
   end subroutine random_sparse_matrix

   subroutine grow(rows, cols, vals, cap)
      integer, allocatable, intent(inout) :: rows(:), cols(:)
      real(dp), allocatable, intent(inout) :: vals(:)
      integer, intent(inout) :: cap
      integer, allocatable :: r(:), c(:)
      real(dp), allocatable :: v(:)
      integer :: old
      old = cap
      cap = 2 * cap
      allocate(r(cap), c(cap), v(cap))
      r(:old) = rows
      c(:old) = cols
      v(:old) = vals
      call move_alloc(r, rows)
      call move_alloc(c, cols)
      call move_alloc(v, vals)
   end subroutine grow

   function companion_matrix(coefficients) result(a)
      real(dp), intent(in) :: coefficients(:)
      real(dp), allocatable :: a(:,:)
      integer :: n, i
      n = size(coefficients)
      allocate(a(n, n), source=0.0_dp)
      if (n == 0) return
      a(1, :) = -coefficients
      do i = 2, n
         a(i, i - 1) = 1.0_dp
      end do
   end function companion_matrix

end module matrix_constructors
