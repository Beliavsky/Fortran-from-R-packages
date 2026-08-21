! SPDX-License-Identifier: GPL-3.0-only
module matrixextra_conversions
   use matrix_kinds, only : dp
   use matrix_status, only : matrix_success, matrix_err_shape, matrix_err_invalid
   use matrix_sparse, only : csr_matrix, csc_matrix, csr_from_dense, csr_to_dense, &
      csc_from_csr, csr_from_csc, csr_from_triplet
   use matrixextra_types, only : coo_matrix, sparse_vector
   implicit none
   private
   public :: coo_from_dense, coo_to_dense, coo_from_csr, csr_from_coo
   public :: coo_from_csc, csc_from_coo
   public :: sparse_vector_from_dense, sparse_vector_to_dense
   public :: csr_from_sparse_vector, coo_from_sparse_vector
   public :: csr_transpose_shallow, csc_transpose_shallow
   public :: csr_to_csc, csc_to_csr

contains

   subroutine coo_from_dense(x, a, tol)
      real(dp), intent(in) :: x(:,:)
      type(coo_matrix), intent(out) :: a
      real(dp), intent(in), optional :: tol
      real(dp) :: eps
      integer :: i, j, k, nz
      eps = 0.0_dp
      if (present(tol)) eps = max(0.0_dp, tol)
      a%nrow = size(x,1)
      a%ncol = size(x,2)
      nz = count(abs(x) > eps)
      allocate(a%row_ind(nz), a%col_ind(nz), a%values(nz))
      k = 0
      do i = 1, a%nrow
         do j = 1, a%ncol
            if (abs(x(i,j)) > eps) then
               k = k + 1
               a%row_ind(k) = i
               a%col_ind(k) = j
               a%values(k) = x(i,j)
            end if
         end do
      end do
   end subroutine coo_from_dense

   function coo_to_dense(a) result(x)
      type(coo_matrix), intent(in) :: a
      real(dp), allocatable :: x(:,:)
      integer :: k
      allocate(x(a%nrow,a%ncol), source=0.0_dp)
      do k = 1, a%nnz()
         x(a%row_ind(k),a%col_ind(k)) = x(a%row_ind(k),a%col_ind(k)) + a%values(k)
      end do
   end function coo_to_dense

   subroutine coo_from_csr(a, b)
      type(csr_matrix), intent(in) :: a
      type(coo_matrix), intent(out) :: b
      integer :: i, k, p
      b%nrow = a%nrow
      b%ncol = a%ncol
      allocate(b%row_ind(a%nnz()), b%col_ind(a%nnz()), b%values(a%nnz()))
      p = 0
      do i = 1, a%nrow
         do k = a%row_ptr(i), a%row_ptr(i+1)-1
            p = p + 1
            b%row_ind(p) = i
            b%col_ind(p) = a%col_ind(k)
            b%values(p) = a%values(k)
         end do
      end do
   end subroutine coo_from_csr

   subroutine csr_from_coo(a, b, info, tol)
      type(coo_matrix), intent(in) :: a
      type(csr_matrix), intent(out) :: b
      integer, intent(out), optional :: info
      real(dp), intent(in), optional :: tol
      integer :: istat
      if (.not. a%valid()) then
         call csr_from_triplet(0,0,[integer ::],[integer ::],[real(dp) ::],b,istat)
         istat = matrix_err_invalid
      else
         if (present(tol)) then
            call csr_from_triplet(a%nrow,a%ncol,a%row_ind,a%col_ind,a%values,b,istat,tol)
         else
            call csr_from_triplet(a%nrow,a%ncol,a%row_ind,a%col_ind,a%values,b,istat)
         end if
      end if
      if (present(info)) info = istat
   end subroutine csr_from_coo

   subroutine coo_from_csc(a, b)
      type(csc_matrix), intent(in) :: a
      type(coo_matrix), intent(out) :: b
      type(csr_matrix) :: r
      call csr_from_csc(a,r)
      call coo_from_csr(r,b)
   end subroutine coo_from_csc

   subroutine csc_from_coo(a,b,info,tol)
      type(coo_matrix), intent(in) :: a
      type(csc_matrix), intent(out) :: b
      integer, intent(out), optional :: info
      real(dp), intent(in), optional :: tol
      type(csr_matrix) :: r
      integer :: istat
      if (present(tol)) then
         call csr_from_coo(a,r,istat,tol)
      else
         call csr_from_coo(a,r,istat)
      end if
      if (istat == matrix_success) call csc_from_csr(r,b)
      if (present(info)) info = istat
   end subroutine csc_from_coo

   subroutine sparse_vector_from_dense(x, v, tol)
      real(dp), intent(in) :: x(:)
      type(sparse_vector), intent(out) :: v
      real(dp), intent(in), optional :: tol
      real(dp) :: eps
      integer :: i, k, nz
      eps = 0.0_dp
      if (present(tol)) eps = max(0.0_dp,tol)
      v%n = size(x)
      nz = count(abs(x) > eps)
      allocate(v%index(nz),v%values(nz))
      k = 0
      do i = 1, size(x)
         if (abs(x(i)) > eps) then
            k = k + 1
            v%index(k) = i
            v%values(k) = x(i)
         end if
      end do
   end subroutine sparse_vector_from_dense

   function sparse_vector_to_dense(v) result(x)
      type(sparse_vector), intent(in) :: v
      real(dp), allocatable :: x(:)
      integer :: k
      allocate(x(v%n), source=0.0_dp)
      do k = 1, v%nnz()
         x(v%index(k)) = x(v%index(k)) + v%values(k)
      end do
   end function sparse_vector_to_dense

   subroutine csr_from_sparse_vector(v, a, as_column, info)
      type(sparse_vector), intent(in) :: v
      type(csr_matrix), intent(out) :: a
      logical, intent(in), optional :: as_column
      integer, intent(out), optional :: info
      logical :: col
      integer, allocatable :: rows(:), cols(:)
      integer :: istat
      col = .true.
      if (present(as_column)) col = as_column
      allocate(rows(v%nnz()),cols(v%nnz()))
      if (col) then
         rows = v%index
         cols = 1
         call csr_from_triplet(v%n,1,rows,cols,v%values,a,istat)
      else
         rows = 1
         cols = v%index
         call csr_from_triplet(1,v%n,rows,cols,v%values,a,istat)
      end if
      if (present(info)) info = istat
   end subroutine csr_from_sparse_vector

   subroutine coo_from_sparse_vector(v,a,as_column)
      type(sparse_vector), intent(in) :: v
      type(coo_matrix), intent(out) :: a
      logical, intent(in), optional :: as_column
      logical :: col
      col = .true.
      if (present(as_column)) col = as_column
      if (col) then
         a%nrow = v%n; a%ncol = 1
         allocate(a%row_ind(v%nnz()),a%col_ind(v%nnz()),a%values(v%nnz()))
         a%row_ind = v%index; a%col_ind = 1; a%values = v%values
      else
         a%nrow = 1; a%ncol = v%n
         allocate(a%row_ind(v%nnz()),a%col_ind(v%nnz()),a%values(v%nnz()))
         a%row_ind = 1; a%col_ind = v%index; a%values = v%values
      end if
   end subroutine coo_from_sparse_vector

   subroutine csr_transpose_shallow(a,b)
      type(csr_matrix), intent(in) :: a
      type(csc_matrix), intent(out) :: b
      b%nrow = a%ncol
      b%ncol = a%nrow
      b%col_ptr = a%row_ptr
      b%row_ind = a%col_ind
      b%values = a%values
   end subroutine csr_transpose_shallow

   subroutine csc_transpose_shallow(a,b)
      type(csc_matrix), intent(in) :: a
      type(csr_matrix), intent(out) :: b
      b%nrow = a%ncol
      b%ncol = a%nrow
      b%row_ptr = a%col_ptr
      b%col_ind = a%row_ind
      b%values = a%values
   end subroutine csc_transpose_shallow

   subroutine csr_to_csc(a,b)
      type(csr_matrix), intent(in) :: a
      type(csc_matrix), intent(out) :: b
      call csc_from_csr(a,b)
   end subroutine csr_to_csc

   subroutine csc_to_csr(a,b)
      type(csc_matrix), intent(in) :: a
      type(csr_matrix), intent(out) :: b
      call csr_from_csc(a,b)
   end subroutine csc_to_csr

end module matrixextra_conversions
