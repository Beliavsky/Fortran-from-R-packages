! SPDX-License-Identifier: GPL-3.0-only
module matrixextra_bind
   use matrix_kinds, only : dp
   use matrix_sparse, only : csr_matrix, csr_from_triplet
   implicit none
   private
   public :: csr_rbind, csr_cbind

contains

   subroutine csr_rbind(a,b,c,info)
      type(csr_matrix), intent(in) :: a,b
      type(csr_matrix), intent(out) :: c
      integer, intent(out), optional :: info
      integer, allocatable :: rows(:),cols(:)
      real(dp), allocatable :: vals(:)
      integer :: i,k,p,istat
      if (a%ncol/=b%ncol) then
         call csr_from_triplet(0,0,[integer ::],[integer ::],[real(dp) ::],c,istat)
         if (present(info)) info=1
         return
      end if
      allocate(rows(a%nnz()+b%nnz()),cols(a%nnz()+b%nnz()),vals(a%nnz()+b%nnz()))
      p=0
      do i=1,a%nrow
         do k=a%row_ptr(i),a%row_ptr(i+1)-1
            p=p+1; rows(p)=i; cols(p)=a%col_ind(k); vals(p)=a%values(k)
         end do
      end do
      do i=1,b%nrow
         do k=b%row_ptr(i),b%row_ptr(i+1)-1
            p=p+1; rows(p)=a%nrow+i; cols(p)=b%col_ind(k); vals(p)=b%values(k)
         end do
      end do
      call csr_from_triplet(a%nrow+b%nrow,a%ncol,rows,cols,vals,c,istat)
      if (present(info)) info=istat
   end subroutine csr_rbind

   subroutine csr_cbind(a,b,c,info)
      type(csr_matrix), intent(in) :: a,b
      type(csr_matrix), intent(out) :: c
      integer, intent(out), optional :: info
      integer, allocatable :: rows(:),cols(:)
      real(dp), allocatable :: vals(:)
      integer :: i,k,p,istat
      if (a%nrow/=b%nrow) then
         call csr_from_triplet(0,0,[integer ::],[integer ::],[real(dp) ::],c,istat)
         if (present(info)) info=1
         return
      end if
      allocate(rows(a%nnz()+b%nnz()),cols(a%nnz()+b%nnz()),vals(a%nnz()+b%nnz()))
      p=0
      do i=1,a%nrow
         do k=a%row_ptr(i),a%row_ptr(i+1)-1
            p=p+1; rows(p)=i; cols(p)=a%col_ind(k); vals(p)=a%values(k)
         end do
         do k=b%row_ptr(i),b%row_ptr(i+1)-1
            p=p+1; rows(p)=i; cols(p)=a%ncol+b%col_ind(k); vals(p)=b%values(k)
         end do
      end do
      call csr_from_triplet(a%nrow,a%ncol+b%ncol,rows,cols,vals,c,istat)
      if (present(info)) info=istat
   end subroutine csr_cbind

end module matrixextra_bind
