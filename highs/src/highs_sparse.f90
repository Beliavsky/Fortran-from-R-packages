! SPDX-License-Identifier: GPL-2.0-or-later
module highs_sparse
   use highs_kinds, only : dp, highs_int
   use highs_constants, only : highs_matrix_colwise, highs_matrix_rowwise, &
      highs_hessian_triangular
   implicit none
   private

   type, public :: highs_sparse_matrix
      integer(highs_int) :: nrow = 0
      integer(highs_int) :: ncol = 0
      integer(highs_int) :: format = highs_matrix_colwise
      integer(highs_int), allocatable :: start(:)
      integer(highs_int), allocatable :: index(:)
      real(dp), allocatable :: value(:)
   contains
      procedure :: nnz => sparse_nnz
      procedure :: valid => sparse_valid
      procedure :: to_dense => sparse_to_dense
   end type highs_sparse_matrix

   public :: highs_csc_from_dense, highs_csr_from_dense
   public :: highs_csc_from_triplets, highs_empty_matrix, highs_hessian_from_dense

contains

   pure integer function sparse_nnz(self) result(n)
      class(highs_sparse_matrix), intent(in) :: self
      if (allocated(self%value)) then
         n = size(self%value)
      else
         n = 0
      end if
   end function sparse_nnz

   pure logical function sparse_valid(self) result(ok)
      class(highs_sparse_matrix), intent(in) :: self
      integer :: major
      ok = self%nrow >= 0 .and. self%ncol >= 0
      if (.not. ok) return
      major = merge(self%ncol, self%nrow, self%format == highs_matrix_colwise)
      if (.not. allocated(self%start)) then
         ok = major == 0
         return
      end if
      ok = size(self%start) == major + 1
      if (.not. ok) return
      if (.not. allocated(self%index) .or. .not. allocated(self%value)) then
         ok = self%start(major + 1) == 0
         return
      end if
      ok = size(self%index) == size(self%value)
      if (.not. ok) return
      ok = self%start(1) == 0 .and. self%start(major + 1) == size(self%value)
      if (.not. ok) return
      ok = all(self%start(2:) >= self%start(:major))
      if (.not. ok) return
      if (size(self%index) == 0) return
      if (self%format == highs_matrix_colwise) then
         ok = minval(self%index) >= 0 .and. maxval(self%index) < self%nrow
      else
         ok = minval(self%index) >= 0 .and. maxval(self%index) < self%ncol
      end if
   end function sparse_valid

   function sparse_to_dense(self) result(a)
      class(highs_sparse_matrix), intent(in) :: self
      real(dp), allocatable :: a(:,:)
      integer :: i, j, k, first, last
      allocate(a(self%nrow, self%ncol), source=0.0_dp)
      if (.not. allocated(self%start)) return
      if (self%format == highs_matrix_colwise) then
         do j = 1, self%ncol
            first = self%start(j) + 1
            last = self%start(j + 1)
            do k = first, last
               i = self%index(k) + 1
               a(i,j) = a(i,j) + self%value(k)
            end do
         end do
      else
         do i = 1, self%nrow
            first = self%start(i) + 1
            last = self%start(i + 1)
            do k = first, last
               j = self%index(k) + 1
               a(i,j) = a(i,j) + self%value(k)
            end do
         end do
      end if
   end function sparse_to_dense

   function highs_empty_matrix(nrow, ncol, format) result(s)
      integer, intent(in) :: nrow, ncol
      integer, intent(in), optional :: format
      type(highs_sparse_matrix) :: s
      integer :: major
      s%nrow = nrow
      s%ncol = ncol
      if (present(format)) s%format = format
      major = merge(ncol, nrow, s%format == highs_matrix_colwise)
      allocate(s%start(major + 1), source=0_highs_int)
      allocate(s%index(0), s%value(0))
   end function highs_empty_matrix

   function highs_csc_from_dense(a, tolerance) result(s)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(in), optional :: tolerance
      type(highs_sparse_matrix) :: s
      real(dp) :: tol
      integer :: i, j, k, nz
      tol = 0.0_dp
      if (present(tolerance)) tol = max(0.0_dp, tolerance)
      s%nrow = size(a,1)
      s%ncol = size(a,2)
      s%format = highs_matrix_colwise
      nz = count(abs(a) > tol)
      allocate(s%start(s%ncol + 1), s%index(nz), s%value(nz))
      k = 0
      s%start(1) = 0
      do j = 1, s%ncol
         do i = 1, s%nrow
            if (abs(a(i,j)) > tol) then
               k = k + 1
               s%index(k) = i - 1
               s%value(k) = a(i,j)
            end if
         end do
         s%start(j + 1) = k
      end do
   end function highs_csc_from_dense

   function highs_csr_from_dense(a, tolerance) result(s)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(in), optional :: tolerance
      type(highs_sparse_matrix) :: s
      real(dp) :: tol
      integer :: i, j, k, nz
      tol = 0.0_dp
      if (present(tolerance)) tol = max(0.0_dp, tolerance)
      s%nrow = size(a,1)
      s%ncol = size(a,2)
      s%format = highs_matrix_rowwise
      nz = count(abs(a) > tol)
      allocate(s%start(s%nrow + 1), s%index(nz), s%value(nz))
      k = 0
      s%start(1) = 0
      do i = 1, s%nrow
         do j = 1, s%ncol
            if (abs(a(i,j)) > tol) then
               k = k + 1
               s%index(k) = j - 1
               s%value(k) = a(i,j)
            end if
         end do
         s%start(i + 1) = k
      end do
   end function highs_csr_from_dense

   function highs_csc_from_triplets(nrow, ncol, row, col, val, tolerance) result(s)
      integer, intent(in) :: nrow, ncol
      integer, intent(in) :: row(:), col(:)
      real(dp), intent(in) :: val(:)
      real(dp), intent(in), optional :: tolerance
      type(highs_sparse_matrix) :: s
      real(dp), allocatable :: a(:,:)
      real(dp) :: tol
      integer :: k
      if (size(row) /= size(col) .or. size(row) /= size(val)) error stop "triplet sizes differ"
      allocate(a(nrow,ncol), source=0.0_dp)
      do k = 1, size(val)
         if (row(k) < 1 .or. row(k) > nrow .or. col(k) < 1 .or. col(k) > ncol) &
            error stop "triplet index out of range"
         a(row(k), col(k)) = a(row(k), col(k)) + val(k)
      end do
      tol = 0.0_dp
      if (present(tolerance)) tol = tolerance
      s = highs_csc_from_dense(a, tol)
   end function highs_csc_from_triplets

   function highs_hessian_from_dense(q, tolerance) result(s)
      real(dp), intent(in) :: q(:,:)
      real(dp), intent(in), optional :: tolerance
      type(highs_sparse_matrix) :: s
      real(dp) :: tol
      integer :: i, j, k, nz, n
      if (size(q,1) /= size(q,2)) error stop "Hessian must be square"
      tol = 0.0_dp
      if (present(tolerance)) tol = max(0.0_dp, tolerance)
      n = size(q,1)
      s%nrow = n
      s%ncol = n
      s%format = highs_hessian_triangular
      nz = 0
      do j = 1, n
         do i = 1, j
            if (abs(q(i,j)) > tol) nz = nz + 1
         end do
      end do
      allocate(s%start(n + 1), s%index(nz), s%value(nz))
      k = 0
      s%start(1) = 0
      do j = 1, n
         do i = 1, j
            if (abs(q(i,j)) > tol) then
               k = k + 1
               s%index(k) = i - 1
               s%value(k) = q(i,j)
            end if
         end do
         s%start(j + 1) = k
      end do
   end function highs_hessian_from_dense

end module highs_sparse
