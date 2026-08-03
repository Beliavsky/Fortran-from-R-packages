! SPDX-License-Identifier: Apache-2.0
module osqp_sparse
   use osqp_kinds, only : dp, osqp_int
   implicit none
   private

   type, public :: osqp_sparse_matrix
      integer(osqp_int) :: nrow = 0
      integer(osqp_int) :: ncol = 0
      integer(osqp_int), allocatable :: col_ptr(:)
      integer(osqp_int), allocatable :: row_index(:)
      real(dp), allocatable :: value(:)
   contains
      procedure :: valid => sparse_valid
      procedure :: nnz => sparse_nnz
      procedure :: to_dense => sparse_to_dense
   end type osqp_sparse_matrix

   public :: osqp_csc_from_dense, osqp_csc_from_triplet, osqp_empty_matrix

contains

   pure integer(osqp_int) function sparse_nnz(self) result(n)
      class(osqp_sparse_matrix), intent(in) :: self
      if (allocated(self%value)) then
         n = int(size(self%value), osqp_int)
      else
         n = 0
      end if
   end function sparse_nnz

   pure logical function sparse_valid(self) result(ok)
      class(osqp_sparse_matrix), intent(in) :: self
      integer :: j, k, first, last
      ok = self%nrow >= 0 .and. self%ncol >= 0
      if (.not. ok) return
      ok = allocated(self%col_ptr) .and. allocated(self%row_index) .and. allocated(self%value)
      if (.not. ok) return
      ok = size(self%col_ptr) == self%ncol + 1
      ok = ok .and. size(self%row_index) == size(self%value)
      if (.not. ok) return
      ok = self%col_ptr(1) == 1 .and. self%col_ptr(size(self%col_ptr)) == size(self%value) + 1
      if (.not. ok) return
      if (any(self%col_ptr(2:) < self%col_ptr(:size(self%col_ptr)-1))) then
         ok = .false.
         return
      end if
      if (size(self%row_index) > 0) then
         if (any(self%row_index < 1) .or. any(self%row_index > self%nrow)) then
            ok = .false.
            return
         end if
      end if
      do j = 1, self%ncol
         first = self%col_ptr(j)
         last = self%col_ptr(j+1) - 1
         do k = first + 1, last
            if (self%row_index(k) <= self%row_index(k-1)) then
               ok = .false.
               return
            end if
         end do
      end do
   end function sparse_valid

   function sparse_to_dense(self) result(a)
      class(osqp_sparse_matrix), intent(in) :: self
      real(dp), allocatable :: a(:,:)
      integer :: j, k
      allocate(a(self%nrow, self%ncol), source=0.0_dp)
      if (.not. self%valid()) return
      do j = 1, self%ncol
         do k = self%col_ptr(j), self%col_ptr(j+1) - 1
            a(self%row_index(k), j) = self%value(k)
         end do
      end do
   end function sparse_to_dense

   function osqp_empty_matrix(nrow, ncol) result(mat)
      integer, intent(in) :: nrow, ncol
      type(osqp_sparse_matrix) :: mat
      mat%nrow = nrow
      mat%ncol = ncol
      allocate(mat%col_ptr(ncol+1), source=1_osqp_int)
      allocate(mat%row_index(0), mat%value(0))
   end function osqp_empty_matrix

   function osqp_csc_from_dense(a, upper_only, tolerance) result(mat)
      real(dp), intent(in) :: a(:,:)
      logical, intent(in), optional :: upper_only
      real(dp), intent(in), optional :: tolerance
      type(osqp_sparse_matrix) :: mat
      logical :: upper
      real(dp) :: tol
      integer :: i, j, k, nz

      upper = .false.
      if (present(upper_only)) upper = upper_only
      tol = 0.0_dp
      if (present(tolerance)) tol = max(0.0_dp, tolerance)
      mat%nrow = size(a, 1)
      mat%ncol = size(a, 2)
      nz = 0
      do j = 1, size(a,2)
         do i = 1, size(a,1)
            if (upper .and. i > j) cycle
            if (abs(a(i,j)) > tol) nz = nz + 1
         end do
      end do
      allocate(mat%col_ptr(mat%ncol+1), mat%row_index(nz), mat%value(nz))
      k = 1
      mat%col_ptr(1) = 1
      do j = 1, mat%ncol
         do i = 1, mat%nrow
            if (upper .and. i > j) cycle
            if (abs(a(i,j)) <= tol) cycle
            mat%row_index(k) = i
            mat%value(k) = a(i,j)
            k = k + 1
         end do
         mat%col_ptr(j+1) = k
      end do
   end function osqp_csc_from_dense

   function osqp_csc_from_triplet(nrow, ncol, row, col, x, upper_only, tolerance, status) result(mat)
      integer, intent(in) :: nrow, ncol
      integer(osqp_int), intent(in) :: row(:), col(:)
      real(dp), intent(in) :: x(:)
      logical, intent(in), optional :: upper_only
      real(dp), intent(in), optional :: tolerance
      integer(osqp_int), intent(out), optional :: status
      type(osqp_sparse_matrix) :: mat
      integer(osqp_int), allocatable :: r(:), c(:)
      real(dp), allocatable :: v(:)
      logical :: upper
      real(dp) :: tol, sumv
      integer :: k, n, kept, outn, j, start

      if (present(status)) status = 1
      mat = osqp_empty_matrix(max(0,nrow), max(0,ncol))
      if (nrow < 0 .or. ncol < 0) return
      if (size(row) /= size(col) .or. size(row) /= size(x)) return
      if (size(row) > 0) then
         if (any(row < 1) .or. any(row > nrow) .or. any(col < 1) .or. any(col > ncol)) return
      end if
      upper = .false.
      if (present(upper_only)) upper = upper_only
      tol = 0.0_dp
      if (present(tolerance)) tol = max(0.0_dp, tolerance)

      kept = 0
      do k = 1, size(x)
         if (upper .and. row(k) > col(k)) cycle
         kept = kept + 1
      end do
      allocate(r(kept), c(kept), v(kept))
      n = 0
      do k = 1, size(x)
         if (upper .and. row(k) > col(k)) cycle
         n = n + 1
         r(n) = row(k)
         c(n) = col(k)
         v(n) = x(k)
      end do
      if (n > 1) call sort_triplets(r, c, v, 1, n)

      outn = 0
      k = 1
      do while (k <= n)
         sumv = v(k)
         start = k
         do while (k < n)
            if (r(k+1) /= r(start) .or. c(k+1) /= c(start)) exit
            k = k + 1
            sumv = sumv + v(k)
         end do
         if (abs(sumv) > tol) outn = outn + 1
         k = k + 1
      end do

      mat%nrow = nrow
      mat%ncol = ncol
      if (allocated(mat%col_ptr)) deallocate(mat%col_ptr, mat%row_index, mat%value)
      allocate(mat%col_ptr(ncol+1), mat%row_index(outn), mat%value(outn))
      mat%col_ptr = 1
      k = 1
      outn = 0
      do j = 1, ncol
         mat%col_ptr(j) = outn + 1
         do while (k <= n)
            if (c(k) >= j) exit
            k = k + 1
         end do
         do while (k <= n)
            if (c(k) /= j) exit
            sumv = v(k)
            start = k
            do while (k < n)
               if (r(k+1) /= r(start) .or. c(k+1) /= c(start)) exit
               k = k + 1
               sumv = sumv + v(k)
            end do
            if (abs(sumv) > tol) then
               outn = outn + 1
               mat%row_index(outn) = r(start)
               mat%value(outn) = sumv
            end if
            k = k + 1
         end do
         mat%col_ptr(j+1) = outn + 1
      end do
      if (present(status)) status = 0
   end function osqp_csc_from_triplet

   recursive subroutine sort_triplets(r, c, v, left, right)
      integer(osqp_int), intent(inout) :: r(:), c(:)
      real(dp), intent(inout) :: v(:)
      integer, intent(in) :: left, right
      integer :: i, j
      integer(osqp_int) :: pr, pc, tr, tc
      real(dp) :: tv
      if (left >= right) return
      i = left
      j = right
      pc = c((left+right)/2)
      pr = r((left+right)/2)
      do
         do while (triplet_less(c(i), r(i), pc, pr))
            i = i + 1
         end do
         do while (triplet_less(pc, pr, c(j), r(j)))
            j = j - 1
         end do
         if (i <= j) then
            tc = c(i); c(i) = c(j); c(j) = tc
            tr = r(i); r(i) = r(j); r(j) = tr
            tv = v(i); v(i) = v(j); v(j) = tv
            i = i + 1
            j = j - 1
         end if
         if (i > j) exit
      end do
      if (left < j) call sort_triplets(r, c, v, left, j)
      if (i < right) call sort_triplets(r, c, v, i, right)
   end subroutine sort_triplets

   pure logical function triplet_less(c1, r1, c2, r2) result(less)
      integer(osqp_int), intent(in) :: c1, r1, c2, r2
      less = c1 < c2 .or. (c1 == c2 .and. r1 < r2)
   end function triplet_less

end module osqp_sparse
