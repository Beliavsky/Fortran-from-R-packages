! SPDX-License-Identifier: GPL-3.0-only
module matrixextra_pattern
   use matrix_kinds, only : dp
   use matrix_sparse, only : csr_matrix, csr_from_triplet
   implicit none
   private
   public :: csr_logical_and, csr_logical_or, csr_patternize

contains

   subroutine csr_patternize(a)
      type(csr_matrix), intent(inout) :: a
      if (allocated(a%values)) a%values=1.0_dp
   end subroutine csr_patternize

   subroutine csr_logical_and(a,b,c,info)
      type(csr_matrix), intent(in) :: a,b
      type(csr_matrix), intent(out) :: c
      integer, intent(out), optional :: info
      integer, allocatable :: rows(:),cols(:)
      real(dp), allocatable :: vals(:)
      integer :: i,ka,kb,nz,istat
      if (a%nrow/=b%nrow .or. a%ncol/=b%ncol) then
         call csr_from_triplet(0,0,[integer ::],[integer ::],[real(dp) ::],c,istat)
         if (present(info)) info=1
         return
      end if
      allocate(rows(min(a%nnz(),b%nnz())),cols(min(a%nnz(),b%nnz())), &
         vals(min(a%nnz(),b%nnz())))
      nz=0
      do i=1,a%nrow
         ka=a%row_ptr(i); kb=b%row_ptr(i)
         do while (ka<a%row_ptr(i+1) .and. kb<b%row_ptr(i+1))
            if (a%col_ind(ka)==b%col_ind(kb)) then
               nz=nz+1; rows(nz)=i; cols(nz)=a%col_ind(ka); vals(nz)=1.0_dp
               ka=ka+1; kb=kb+1
            else if (a%col_ind(ka)<b%col_ind(kb)) then
               ka=ka+1
            else
               kb=kb+1
            end if
         end do
      end do
      call csr_from_triplet(a%nrow,a%ncol,rows(:nz),cols(:nz),vals(:nz),c,istat)
      if (present(info)) info=istat
   end subroutine csr_logical_and

   subroutine csr_logical_or(a,b,c,info)
      type(csr_matrix), intent(in) :: a,b
      type(csr_matrix), intent(out) :: c
      integer, intent(out), optional :: info
      integer, allocatable :: rows(:),cols(:)
      real(dp), allocatable :: vals(:)
      integer :: i,ka,kb,nz,istat,col
      if (a%nrow/=b%nrow .or. a%ncol/=b%ncol) then
         call csr_from_triplet(0,0,[integer ::],[integer ::],[real(dp) ::],c,istat)
         if (present(info)) info=1
         return
      end if
      allocate(rows(a%nnz()+b%nnz()),cols(a%nnz()+b%nnz()),vals(a%nnz()+b%nnz()))
      nz=0
      do i=1,a%nrow
         ka=a%row_ptr(i); kb=b%row_ptr(i)
         do while (ka<a%row_ptr(i+1) .or. kb<b%row_ptr(i+1))
            if (ka>=a%row_ptr(i+1)) then
               col=b%col_ind(kb); kb=kb+1
            else if (kb>=b%row_ptr(i+1)) then
               col=a%col_ind(ka); ka=ka+1
            else if (a%col_ind(ka)==b%col_ind(kb)) then
               col=a%col_ind(ka); ka=ka+1; kb=kb+1
            else if (a%col_ind(ka)<b%col_ind(kb)) then
               col=a%col_ind(ka); ka=ka+1
            else
               col=b%col_ind(kb); kb=kb+1
            end if
            nz=nz+1; rows(nz)=i; cols(nz)=col; vals(nz)=1.0_dp
         end do
      end do
      call csr_from_triplet(a%nrow,a%ncol,rows(:nz),cols(:nz),vals(:nz),c,istat)
      if (present(info)) info=istat
   end subroutine csr_logical_or

end module matrixextra_pattern
