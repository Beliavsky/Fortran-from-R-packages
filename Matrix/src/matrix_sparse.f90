! SPDX-License-Identifier: GPL-3.0-only
module matrix_sparse
   use matrix_kinds, only : dp
   use matrix_status, only : matrix_success, matrix_err_shape, matrix_err_invalid
   implicit none
   private

   type, public :: csr_matrix
      integer :: nrow = 0
      integer :: ncol = 0
      integer, allocatable :: row_ptr(:)
      integer, allocatable :: col_ind(:)
      real(dp), allocatable :: values(:)
   contains
      procedure :: nnz => csr_nnz
      procedure :: valid => csr_valid
   end type csr_matrix

   type, public :: csc_matrix
      integer :: nrow = 0
      integer :: ncol = 0
      integer, allocatable :: col_ptr(:)
      integer, allocatable :: row_ind(:)
      real(dp), allocatable :: values(:)
   contains
      procedure :: nnz => csc_nnz
      procedure :: valid => csc_valid
   end type csc_matrix

   public :: csr_from_dense, csr_from_triplet, csr_to_dense
   public :: csc_from_csr, csr_from_csc
   public :: csr_transpose, csr_matvec, csr_matmat
   public :: csr_add, csr_multiply, csr_scale, csr_drop0
   public :: csr_band, csr_diagonal, csr_kronecker, csr_permute
   public :: csr_is_symmetric, csr_equal

contains

   integer function csr_nnz(self) result(n)
      class(csr_matrix), intent(in) :: self
      if (allocated(self%values)) then
         n = size(self%values)
      else
         n = 0
      end if
   end function csr_nnz

   integer function csc_nnz(self) result(n)
      class(csc_matrix), intent(in) :: self
      if (allocated(self%values)) then
         n = size(self%values)
      else
         n = 0
      end if
   end function csc_nnz

   logical function csr_valid(self) result(ok)
      class(csr_matrix), intent(in) :: self
      integer :: i, k0, k1
      ok = self%nrow >= 0 .and. self%ncol >= 0
      if (.not. allocated(self%row_ptr) .or. .not. allocated(self%col_ind) .or. &
         .not. allocated(self%values)) then
         ok = .false.
         return
      end if
      if (size(self%row_ptr) /= self%nrow + 1 .or. size(self%col_ind) /= size(self%values)) then
         ok = .false.
         return
      end if
      if (self%row_ptr(1) /= 1 .or. self%row_ptr(self%nrow + 1) /= size(self%values) + 1) then
         ok = .false.
         return
      end if
      if (any(self%row_ptr(2:) < self%row_ptr(:self%nrow))) then
         ok = .false.
         return
      end if
      if (any(self%col_ind < 1) .or. any(self%col_ind > self%ncol)) then
         ok = .false.
         return
      end if
      do i = 1, self%nrow
         k0 = self%row_ptr(i)
         k1 = self%row_ptr(i + 1) - 1
         if (k1 > k0) then
            if (any(self%col_ind(k0 + 1:k1) <= self%col_ind(k0:k1 - 1))) then
               ok = .false.
               return
            end if
         end if
      end do
   end function csr_valid

   logical function csc_valid(self) result(ok)
      class(csc_matrix), intent(in) :: self
      integer :: j, k0, k1
      ok = self%nrow >= 0 .and. self%ncol >= 0
      if (.not. allocated(self%col_ptr) .or. .not. allocated(self%row_ind) .or. &
         .not. allocated(self%values)) then
         ok = .false.
         return
      end if
      if (size(self%col_ptr) /= self%ncol + 1 .or. size(self%row_ind) /= size(self%values)) then
         ok = .false.
         return
      end if
      if (self%col_ptr(1) /= 1 .or. self%col_ptr(self%ncol + 1) /= size(self%values) + 1) then
         ok = .false.
         return
      end if
      if (any(self%col_ptr(2:) < self%col_ptr(:self%ncol))) then
         ok = .false.
         return
      end if
      if (any(self%row_ind < 1) .or. any(self%row_ind > self%nrow)) then
         ok = .false.
         return
      end if
      do j = 1, self%ncol
         k0 = self%col_ptr(j)
         k1 = self%col_ptr(j + 1) - 1
         if (k1 > k0) then
            if (any(self%row_ind(k0 + 1:k1) <= self%row_ind(k0:k1 - 1))) then
               ok = .false.
               return
            end if
         end if
      end do
   end function csc_valid

   subroutine csr_from_dense(x, a, tol)
      real(dp), intent(in) :: x(:,:)
      type(csr_matrix), intent(out) :: a
      real(dp), intent(in), optional :: tol
      real(dp) :: eps
      integer :: i, j, k, n
      eps = 0.0_dp
      if (present(tol)) eps = max(0.0_dp, tol)
      a%nrow = size(x, 1)
      a%ncol = size(x, 2)
      n = count(abs(x) > eps)
      allocate(a%row_ptr(a%nrow + 1), a%col_ind(n), a%values(n))
      k = 0
      a%row_ptr(1) = 1
      do i = 1, a%nrow
         do j = 1, a%ncol
            if (abs(x(i, j)) > eps) then
               k = k + 1
               a%col_ind(k) = j
               a%values(k) = x(i, j)
            end if
         end do
         a%row_ptr(i + 1) = k + 1
      end do
   end subroutine csr_from_dense

   subroutine csr_from_triplet(nrow, ncol, rows, cols, vals, a, info, tol)
      integer, intent(in) :: nrow, ncol, rows(:), cols(:)
      real(dp), intent(in) :: vals(:)
      type(csr_matrix), intent(out) :: a
      integer, intent(out) :: info
      real(dp), intent(in), optional :: tol
      integer, allocatable :: r(:), c(:), row_counts(:)
      real(dp), allocatable :: v(:), cv(:)
      real(dp) :: eps, sumv
      integer :: n, i, k, out_n
      if (nrow < 0 .or. ncol < 0 .or. size(rows) /= size(cols) .or. &
         size(rows) /= size(vals)) then
         call empty_csr(a)
         info = matrix_err_shape
         return
      end if
      if (any(rows < 1) .or. any(rows > nrow) .or. any(cols < 1) .or. any(cols > ncol)) then
         call empty_csr(a)
         info = matrix_err_invalid
         return
      end if
      eps = 0.0_dp
      if (present(tol)) eps = max(0.0_dp, tol)
      n = size(vals)
      allocate(r(n), c(n), v(n))
      r = rows
      c = cols
      v = vals
      call sort_triplets(r, c, v)
      allocate(cv(n), source=0.0_dp)
      allocate(row_counts(nrow), source=0)
      allocate(a%col_ind(n))
      out_n = 0
      i = 1
      do while (i <= n)
         sumv = v(i)
         k = i + 1
         do while (k <= n)
            if (r(k) /= r(i) .or. c(k) /= c(i)) exit
            sumv = sumv + v(k)
            k = k + 1
         end do
         if (abs(sumv) > eps) then
            out_n = out_n + 1
            a%col_ind(out_n) = c(i)
            cv(out_n) = sumv
            row_counts(r(i)) = row_counts(r(i)) + 1
         end if
         i = k
      end do
      a%nrow = nrow
      a%ncol = ncol
      allocate(a%row_ptr(nrow + 1))
      a%row_ptr(1) = 1
      do i = 1, nrow
         a%row_ptr(i + 1) = a%row_ptr(i) + row_counts(i)
      end do
      allocate(a%values(out_n))
      if (out_n > 0) a%values = cv(:out_n)
      if (out_n < n) a%col_ind = a%col_ind(:out_n)
      info = matrix_success
   end subroutine csr_from_triplet

   subroutine sort_triplets(r, c, v)
      integer, intent(inout) :: r(:), c(:)
      real(dp), intent(inout) :: v(:)
      integer :: gap, i, j, tr, tc
      real(dp) :: tv
      gap = size(r) / 2
      do while (gap > 0)
         do i = gap + 1, size(r)
            tr = r(i)
            tc = c(i)
            tv = v(i)
            j = i
            do while (j > gap)
               if (r(j - gap) < tr) exit
               if (r(j - gap) == tr .and. c(j - gap) <= tc) exit
               r(j) = r(j - gap)
               c(j) = c(j - gap)
               v(j) = v(j - gap)
               j = j - gap
            end do
            r(j) = tr
            c(j) = tc
            v(j) = tv
         end do
         gap = gap / 2
      end do
   end subroutine sort_triplets

   subroutine empty_csr(a)
      type(csr_matrix), intent(out) :: a
      a%nrow = 0
      a%ncol = 0
      allocate(a%row_ptr(1), a%col_ind(0), a%values(0))
      a%row_ptr = 1
   end subroutine empty_csr

   function csr_to_dense(a) result(x)
      type(csr_matrix), intent(in) :: a
      real(dp), allocatable :: x(:,:)
      integer :: i, k
      allocate(x(a%nrow, a%ncol), source=0.0_dp)
      do i = 1, a%nrow
         do k = a%row_ptr(i), a%row_ptr(i + 1) - 1
            x(i, a%col_ind(k)) = a%values(k)
         end do
      end do
   end function csr_to_dense

   subroutine csc_from_csr(a, b)
      type(csr_matrix), intent(in) :: a
      type(csc_matrix), intent(out) :: b
      integer, allocatable :: count_col(:), next(:)
      integer :: i, j, k, pos
      b%nrow = a%nrow
      b%ncol = a%ncol
      allocate(count_col(b%ncol), source=0)
      do k = 1, a%nnz()
         count_col(a%col_ind(k)) = count_col(a%col_ind(k)) + 1
      end do
      allocate(b%col_ptr(b%ncol + 1))
      b%col_ptr(1) = 1
      do j = 1, b%ncol
         b%col_ptr(j + 1) = b%col_ptr(j) + count_col(j)
      end do
      allocate(b%row_ind(a%nnz()), b%values(a%nnz()), next(b%ncol))
      next = b%col_ptr(:b%ncol)
      do i = 1, a%nrow
         do k = a%row_ptr(i), a%row_ptr(i + 1) - 1
            j = a%col_ind(k)
            pos = next(j)
            b%row_ind(pos) = i
            b%values(pos) = a%values(k)
            next(j) = next(j) + 1
         end do
      end do
   end subroutine csc_from_csr

   subroutine csr_from_csc(a, b)
      type(csc_matrix), intent(in) :: a
      type(csr_matrix), intent(out) :: b
      integer, allocatable :: count_row(:), next(:)
      integer :: i, j, k, pos
      b%nrow = a%nrow
      b%ncol = a%ncol
      allocate(count_row(b%nrow), source=0)
      do k = 1, a%nnz()
         count_row(a%row_ind(k)) = count_row(a%row_ind(k)) + 1
      end do
      allocate(b%row_ptr(b%nrow + 1))
      b%row_ptr(1) = 1
      do i = 1, b%nrow
         b%row_ptr(i + 1) = b%row_ptr(i) + count_row(i)
      end do
      allocate(b%col_ind(a%nnz()), b%values(a%nnz()), next(b%nrow))
      next = b%row_ptr(:b%nrow)
      do j = 1, a%ncol
         do k = a%col_ptr(j), a%col_ptr(j + 1) - 1
            i = a%row_ind(k)
            pos = next(i)
            b%col_ind(pos) = j
            b%values(pos) = a%values(k)
            next(i) = next(i) + 1
         end do
      end do
   end subroutine csr_from_csc

   subroutine csr_transpose(a, at)
      type(csr_matrix), intent(in) :: a
      type(csr_matrix), intent(out) :: at
      type(csc_matrix) :: csc
      call csc_from_csr(a, csc)
      at%nrow = a%ncol
      at%ncol = a%nrow
      at%row_ptr = csc%col_ptr
      at%col_ind = csc%row_ind
      at%values = csc%values
   end subroutine csr_transpose

   recursive subroutine csr_matvec(a, x, y, info, transpose_a)
      type(csr_matrix), intent(in) :: a
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: y(:)
      integer, intent(out) :: info
      logical, intent(in), optional :: transpose_a
      logical :: trans
      type(csr_matrix) :: at
      integer :: i, k
      trans = .false.
      if (present(transpose_a)) trans = transpose_a
      if (trans) then
         call csr_transpose(a, at)
         call csr_matvec(at, x, y, info)
         return
      end if
      if (size(x) /= a%ncol) then
         allocate(y(0))
         info = matrix_err_shape
         return
      end if
      allocate(y(a%nrow), source=0.0_dp)
      do i = 1, a%nrow
         do k = a%row_ptr(i), a%row_ptr(i + 1) - 1
            y(i) = y(i) + a%values(k) * x(a%col_ind(k))
         end do
      end do
      info = matrix_success
   end subroutine csr_matvec

   subroutine csr_matmat(a, x, y, info)
      type(csr_matrix), intent(in) :: a
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable, intent(out) :: y(:,:)
      integer, intent(out) :: info
      integer :: i, k
      if (size(x, 1) /= a%ncol) then
         allocate(y(0, 0))
         info = matrix_err_shape
         return
      end if
      allocate(y(a%nrow, size(x, 2)), source=0.0_dp)
      do i = 1, a%nrow
         do k = a%row_ptr(i), a%row_ptr(i + 1) - 1
            y(i, :) = y(i, :) + a%values(k) * x(a%col_ind(k), :)
         end do
      end do
      info = matrix_success
   end subroutine csr_matmat

   subroutine csr_add(a, b, c, info, alpha, beta, tol)
      type(csr_matrix), intent(in) :: a, b
      type(csr_matrix), intent(out) :: c
      integer, intent(out) :: info
      real(dp), intent(in), optional :: alpha, beta, tol
      integer, allocatable :: rows(:), cols(:)
      real(dp), allocatable :: vals(:)
      real(dp) :: sa, sb, eps
      integer :: i, ka, kb, enda, endb, n
      if (a%nrow /= b%nrow .or. a%ncol /= b%ncol) then
         call empty_csr(c)
         info = matrix_err_shape
         return
      end if
      sa = 1.0_dp
      sb = 1.0_dp
      eps = 0.0_dp
      if (present(alpha)) sa = alpha
      if (present(beta)) sb = beta
      if (present(tol)) eps = max(0.0_dp, tol)
      allocate(rows(a%nnz() + b%nnz()), cols(a%nnz() + b%nnz()), vals(a%nnz() + b%nnz()))
      n = 0
      do i = 1, a%nrow
         ka = a%row_ptr(i)
         kb = b%row_ptr(i)
         enda = a%row_ptr(i + 1) - 1
         endb = b%row_ptr(i + 1) - 1
         do while (ka <= enda .or. kb <= endb)
            if (ka > enda) then
               n = n + 1
               rows(n) = i
               cols(n) = b%col_ind(kb)
               vals(n) = sb * b%values(kb)
               kb = kb + 1
            else if (kb > endb) then
               n = n + 1
               rows(n) = i
               cols(n) = a%col_ind(ka)
               vals(n) = sa * a%values(ka)
               ka = ka + 1
            else if (a%col_ind(ka) < b%col_ind(kb)) then
               n = n + 1
               rows(n) = i
               cols(n) = a%col_ind(ka)
               vals(n) = sa * a%values(ka)
               ka = ka + 1
            else if (b%col_ind(kb) < a%col_ind(ka)) then
               n = n + 1
               rows(n) = i
               cols(n) = b%col_ind(kb)
               vals(n) = sb * b%values(kb)
               kb = kb + 1
            else
               n = n + 1
               rows(n) = i
               cols(n) = a%col_ind(ka)
               vals(n) = sa * a%values(ka) + sb * b%values(kb)
               ka = ka + 1
               kb = kb + 1
            end if
         end do
      end do
      call csr_from_triplet(a%nrow, a%ncol, rows(:n), cols(:n), vals(:n), c, info, tol=eps)
   end subroutine csr_add

   subroutine csr_multiply(a, b, c, info, tol)
      type(csr_matrix), intent(in) :: a, b
      type(csr_matrix), intent(out) :: c
      integer, intent(out) :: info
      real(dp), intent(in), optional :: tol
      integer, allocatable :: rows(:), cols(:), marker(:), used(:)
      real(dp), allocatable :: vals(:), accum(:)
      real(dp) :: eps
      integer :: i, ka, kb, j, col, n, nused, cap, q
      if (a%ncol /= b%nrow) then
         call empty_csr(c)
         info = matrix_err_shape
         return
      end if
      eps = 0.0_dp
      if (present(tol)) eps = max(0.0_dp, tol)
      cap = max(16, min(max(1, a%nrow * b%ncol), max(16, a%nnz() + b%nnz())))
      allocate(rows(cap), cols(cap), vals(cap))
      allocate(marker(b%ncol), source=0)
      allocate(accum(b%ncol), source=0.0_dp)
      allocate(used(b%ncol))
      n = 0
      do i = 1, a%nrow
         nused = 0
         do ka = a%row_ptr(i), a%row_ptr(i + 1) - 1
            j = a%col_ind(ka)
            do kb = b%row_ptr(j), b%row_ptr(j + 1) - 1
               col = b%col_ind(kb)
               if (marker(col) /= i) then
                  marker(col) = i
                  nused = nused + 1
                  used(nused) = col
                  accum(col) = a%values(ka) * b%values(kb)
               else
                  accum(col) = accum(col) + a%values(ka) * b%values(kb)
               end if
            end do
         end do
         call sort_integer(used(:nused))
         do q = 1, nused
            col = used(q)
            if (abs(accum(col)) > eps) then
               if (n == cap) call grow_triplet_arrays(rows, cols, vals, cap)
               n = n + 1
               rows(n) = i
               cols(n) = col
               vals(n) = accum(col)
            end if
            accum(col) = 0.0_dp
         end do
      end do
      call csr_from_triplet(a%nrow, b%ncol, rows(:n), cols(:n), vals(:n), c, info, tol=eps)
   end subroutine csr_multiply

   subroutine grow_triplet_arrays(rows, cols, vals, cap)
      integer, allocatable, intent(inout) :: rows(:), cols(:)
      real(dp), allocatable, intent(inout) :: vals(:)
      integer, intent(inout) :: cap
      integer, allocatable :: ri(:), ci(:)
      real(dp), allocatable :: vi(:)
      integer :: old
      old = cap
      cap = 2 * cap
      allocate(ri(cap), ci(cap), vi(cap))
      ri(:old) = rows
      ci(:old) = cols
      vi(:old) = vals
      call move_alloc(ri, rows)
      call move_alloc(ci, cols)
      call move_alloc(vi, vals)
   end subroutine grow_triplet_arrays

   subroutine sort_integer(x)
      integer, intent(inout) :: x(:)
      integer :: i, j, key
      do i = 2, size(x)
         key = x(i)
         j = i - 1
         do while (j >= 1)
            if (x(j) <= key) exit
            x(j + 1) = x(j)
            j = j - 1
         end do
         x(j + 1) = key
      end do
   end subroutine sort_integer

   subroutine csr_scale(a, scalar, b)
      type(csr_matrix), intent(in) :: a
      real(dp), intent(in) :: scalar
      type(csr_matrix), intent(out) :: b
      b%nrow = a%nrow
      b%ncol = a%ncol
      b%row_ptr = a%row_ptr
      b%col_ind = a%col_ind
      b%values = scalar * a%values
   end subroutine csr_scale

   subroutine csr_drop0(a, b, tol)
      type(csr_matrix), intent(in) :: a
      type(csr_matrix), intent(out) :: b
      real(dp), intent(in), optional :: tol
      integer, allocatable :: rows(:), cols(:)
      real(dp), allocatable :: vals(:)
      real(dp) :: eps
      integer :: i, k, n, info
      eps = 0.0_dp
      if (present(tol)) eps = max(0.0_dp, tol)
      allocate(rows(a%nnz()), cols(a%nnz()), vals(a%nnz()))
      n = 0
      do i = 1, a%nrow
         do k = a%row_ptr(i), a%row_ptr(i + 1) - 1
            if (abs(a%values(k)) > eps) then
               n = n + 1
               rows(n) = i
               cols(n) = a%col_ind(k)
               vals(n) = a%values(k)
            end if
         end do
      end do
      call csr_from_triplet(a%nrow, a%ncol, rows(:n), cols(:n), vals(:n), b, info)
   end subroutine csr_drop0

   subroutine csr_band(a, lower, upper, b)
      type(csr_matrix), intent(in) :: a
      integer, intent(in) :: lower, upper
      type(csr_matrix), intent(out) :: b
      integer, allocatable :: rows(:), cols(:)
      real(dp), allocatable :: vals(:)
      integer :: i, j, k, n, info
      allocate(rows(a%nnz()), cols(a%nnz()), vals(a%nnz()))
      n = 0
      do i = 1, a%nrow
         do k = a%row_ptr(i), a%row_ptr(i + 1) - 1
            j = a%col_ind(k)
            if (i - j >= lower .and. i - j <= upper) then
               n = n + 1
               rows(n) = i
               cols(n) = j
               vals(n) = a%values(k)
            end if
         end do
      end do
      call csr_from_triplet(a%nrow, a%ncol, rows(:n), cols(:n), vals(:n), b, info)
   end subroutine csr_band

   function csr_diagonal(a, offset) result(d)
      type(csr_matrix), intent(in) :: a
      integer, intent(in), optional :: offset
      real(dp), allocatable :: d(:)
      integer :: kdiag, n, row0, col0, q, i, k
      kdiag = 0
      if (present(offset)) kdiag = offset
      if (kdiag >= 0) then
         row0 = 1
         col0 = 1 + kdiag
         n = min(a%nrow, a%ncol - kdiag)
      else
         row0 = 1 - kdiag
         col0 = 1
         n = min(a%nrow + kdiag, a%ncol)
      end if
      n = max(0, n)
      allocate(d(n), source=0.0_dp)
      do q = 1, n
         i = row0 + q - 1
         do k = a%row_ptr(i), a%row_ptr(i + 1) - 1
            if (a%col_ind(k) == col0 + q - 1) then
               d(q) = a%values(k)
               exit
            end if
         end do
      end do
   end function csr_diagonal

   subroutine csr_kronecker(a, b, c, info, tol)
      type(csr_matrix), intent(in) :: a, b
      type(csr_matrix), intent(out) :: c
      integer, intent(out) :: info
      real(dp), intent(in), optional :: tol
      integer, allocatable :: rows(:), cols(:)
      real(dp), allocatable :: vals(:)
      integer :: ia, ib, ka, kb, n
      n = a%nnz() * b%nnz()
      allocate(rows(n), cols(n), vals(n))
      n = 0
      do ia = 1, a%nrow
         do ka = a%row_ptr(ia), a%row_ptr(ia + 1) - 1
            do ib = 1, b%nrow
               do kb = b%row_ptr(ib), b%row_ptr(ib + 1) - 1
                  n = n + 1
                  rows(n) = (ia - 1) * b%nrow + ib
                  cols(n) = (a%col_ind(ka) - 1) * b%ncol + b%col_ind(kb)
                  vals(n) = a%values(ka) * b%values(kb)
               end do
            end do
         end do
      end do
      call csr_from_triplet(a%nrow * b%nrow, a%ncol * b%ncol, rows, cols, vals, c, info, tol)
   end subroutine csr_kronecker

   subroutine csr_permute(a, b, info, row_perm, col_perm)
      type(csr_matrix), intent(in) :: a
      type(csr_matrix), intent(out) :: b
      integer, intent(out) :: info
      integer, intent(in), optional :: row_perm(:), col_perm(:)
      integer, allocatable :: inv_r(:), inv_c(:), rows(:), cols(:)
      real(dp), allocatable :: vals(:)
      integer :: i, k, n, new_i, new_j
      allocate(inv_r(a%nrow), inv_c(a%ncol))
      inv_r = [(i, i = 1, a%nrow)]
      inv_c = [(i, i = 1, a%ncol)]
      if (present(row_perm)) then
         if (.not. valid_permutation(row_perm, a%nrow)) then
            call empty_csr(b)
            info = matrix_err_invalid
            return
         end if
         do i = 1, a%nrow
            inv_r(row_perm(i)) = i
         end do
      end if
      if (present(col_perm)) then
         if (.not. valid_permutation(col_perm, a%ncol)) then
            call empty_csr(b)
            info = matrix_err_invalid
            return
         end if
         do i = 1, a%ncol
            inv_c(col_perm(i)) = i
         end do
      end if
      n = a%nnz()
      allocate(rows(n), cols(n), vals(n))
      n = 0
      do i = 1, a%nrow
         new_i = inv_r(i)
         do k = a%row_ptr(i), a%row_ptr(i + 1) - 1
            new_j = inv_c(a%col_ind(k))
            n = n + 1
            rows(n) = new_i
            cols(n) = new_j
            vals(n) = a%values(k)
         end do
      end do
      call csr_from_triplet(a%nrow, a%ncol, rows, cols, vals, b, info)
   end subroutine csr_permute

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

   logical function csr_equal(a, b, tol) result(ok)
      type(csr_matrix), intent(in) :: a, b
      real(dp), intent(in), optional :: tol
      real(dp) :: eps
      eps = 0.0_dp
      if (present(tol)) eps = max(0.0_dp, tol)
      ok = a%nrow == b%nrow .and. a%ncol == b%ncol
      if (.not. ok) return
      if (a%nnz() /= b%nnz()) then
         ok = .false.
         return
      end if
      ok = all(a%row_ptr == b%row_ptr) .and. all(a%col_ind == b%col_ind)
      if (ok) then
         if (a%nnz() > 0) ok = maxval(abs(a%values - b%values)) <= eps
      end if
   end function csr_equal

   logical function csr_is_symmetric(a, tol) result(ok)
      type(csr_matrix), intent(in) :: a
      real(dp), intent(in), optional :: tol
      type(csr_matrix) :: at
      if (a%nrow /= a%ncol) then
         ok = .false.
         return
      end if
      call csr_transpose(a, at)
      ok = csr_equal(a, at, tol)
   end function csr_is_symmetric

end module matrix_sparse
