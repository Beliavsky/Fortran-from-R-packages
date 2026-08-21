! SPDX-License-Identifier: GPL-3.0-only
module matrixextra_slice
   use matrix_kinds, only : dp
   use matrix_sparse, only : csr_matrix, csr_from_triplet
   use matrixextra_types, only : coo_matrix
   use matrixextra_conversions, only : coo_from_csr, csr_from_coo
   implicit none
   private
   public :: csr_get, csr_slice, coo_slice
   public :: csr_set_value, csr_set_row_constant, csr_set_col_constant
   public :: csr_set_block_constant, csr_assign_dense_block

contains

   real(dp) function csr_get(a,row,col) result(v)
      type(csr_matrix), intent(in) :: a
      integer, intent(in) :: row,col
      integer :: k
      v=0.0_dp
      if (row<1 .or. row>a%nrow .or. col<1 .or. col>a%ncol) return
      do k=a%row_ptr(row),a%row_ptr(row+1)-1
         if (a%col_ind(k)==col) then
            v=a%values(k); return
         else if (a%col_ind(k)>col) then
            return
         end if
      end do
   end function csr_get

   subroutine csr_slice(a,rows,cols,b)
      type(csr_matrix), intent(in) :: a
      integer, intent(in) :: rows(:), cols(:)
      type(csr_matrix), intent(out) :: b
      integer :: rr,cc,k,nz,p,info
      integer, allocatable :: ro(:),co(:)
      real(dp), allocatable :: va(:)
      if (size(rows)==0 .or. size(cols)==0) then
         call csr_from_triplet(size(rows),size(cols),[integer ::],[integer ::], &
            [real(dp) ::],b,info)
         return
      end if
      if (any(rows<1) .or. any(rows>a%nrow) .or. any(cols<1) .or. any(cols>a%ncol)) then
         call csr_from_triplet(0,0,[integer ::],[integer ::],[real(dp) ::],b,info)
         return
      end if
      nz=0
      do rr=1,size(rows)
         do k=a%row_ptr(rows(rr)),a%row_ptr(rows(rr)+1)-1
            nz=nz+count(cols==a%col_ind(k))
         end do
      end do
      allocate(ro(nz),co(nz),va(nz))
      p=0
      do rr=1,size(rows)
         do k=a%row_ptr(rows(rr)),a%row_ptr(rows(rr)+1)-1
            do cc=1,size(cols)
               if (cols(cc)==a%col_ind(k)) then
                  p=p+1; ro(p)=rr; co(p)=cc; va(p)=a%values(k)
               end if
            end do
         end do
      end do
      call csr_from_triplet(size(rows),size(cols),ro,co,va,b,info)
   end subroutine csr_slice

   subroutine coo_slice(a,rows,cols,b)
      type(coo_matrix), intent(in) :: a
      integer, intent(in) :: rows(:),cols(:)
      type(coo_matrix), intent(out) :: b
      type(csr_matrix) :: r,s
      integer :: info
      call csr_from_coo(a,r,info)
      call csr_slice(r,rows,cols,s)
      call coo_from_csr(s,b)
   end subroutine coo_slice

   subroutine csr_set_value(a,row,col,value)
      type(csr_matrix), intent(inout) :: a
      integer, intent(in) :: row,col
      real(dp), intent(in) :: value
      integer, allocatable :: rr(:),cc(:)
      real(dp), allocatable :: vv(:)
      integer :: i,k,p,info,nr,nc
      nr=a%nrow; nc=a%ncol
      if (row<1 .or. row>a%nrow .or. col<1 .or. col>a%ncol) return
      allocate(rr(a%nnz()+1),cc(a%nnz()+1),vv(a%nnz()+1))
      p=0
      do i=1,a%nrow
         do k=a%row_ptr(i),a%row_ptr(i+1)-1
            if (i==row .and. a%col_ind(k)==col) cycle
            p=p+1; rr(p)=i; cc(p)=a%col_ind(k); vv(p)=a%values(k)
         end do
      end do
      p=p+1; rr(p)=row; cc(p)=col; vv(p)=value
      call csr_from_triplet(nr,nc,rr(:p),cc(:p),vv(:p),a,info)
   end subroutine csr_set_value

   subroutine csr_set_row_constant(a,row,value)
      type(csr_matrix), intent(inout) :: a
      integer, intent(in) :: row
      real(dp), intent(in) :: value
      integer, allocatable :: rr(:),cc(:)
      real(dp), allocatable :: vv(:)
      integer :: i,j,k,p,info,cap,nr,nc
      nr=a%nrow; nc=a%ncol
      if (row<1 .or. row>nr) return
      cap=a%nnz()+nc
      allocate(rr(cap),cc(cap),vv(cap)); p=0
      do i=1,nr
         if (i==row) cycle
         do k=a%row_ptr(i),a%row_ptr(i+1)-1
            p=p+1; rr(p)=i; cc(p)=a%col_ind(k); vv(p)=a%values(k)
         end do
      end do
      do j=1,nc
         p=p+1; rr(p)=row; cc(p)=j; vv(p)=value
      end do
      call csr_from_triplet(nr,nc,rr(:p),cc(:p),vv(:p),a,info)
   end subroutine csr_set_row_constant

   subroutine csr_set_col_constant(a,col,value)
      type(csr_matrix), intent(inout) :: a
      integer, intent(in) :: col
      real(dp), intent(in) :: value
      integer, allocatable :: rr(:),cc(:)
      real(dp), allocatable :: vv(:)
      integer :: i,k,p,info,cap,nr,nc
      nr=a%nrow; nc=a%ncol
      if (col<1 .or. col>nc) return
      cap=a%nnz()+nr
      allocate(rr(cap),cc(cap),vv(cap)); p=0
      do i=1,nr
         do k=a%row_ptr(i),a%row_ptr(i+1)-1
            if (a%col_ind(k)==col) cycle
            p=p+1; rr(p)=i; cc(p)=a%col_ind(k); vv(p)=a%values(k)
         end do
         p=p+1; rr(p)=i; cc(p)=col; vv(p)=value
      end do
      call csr_from_triplet(nr,nc,rr(:p),cc(:p),vv(:p),a,info)
   end subroutine csr_set_col_constant

   subroutine csr_set_block_constant(a,rows,cols,value)
      type(csr_matrix), intent(inout) :: a
      integer, intent(in) :: rows(:),cols(:)
      real(dp), intent(in) :: value
      real(dp), allocatable :: block(:,:)
      allocate(block(size(rows),size(cols)),source=value)
      call csr_assign_dense_block(a,rows,cols,block)
   end subroutine csr_set_block_constant

   subroutine csr_assign_dense_block(a,rows,cols,value)
      type(csr_matrix), intent(inout) :: a
      integer, intent(in) :: rows(:),cols(:)
      real(dp), intent(in) :: value(:,:)
      integer, allocatable :: rr(:),cc(:)
      real(dp), allocatable :: vv(:)
      integer :: i,j,k,p,info,cap,nr,nc
      logical :: replace
      nr=a%nrow; nc=a%ncol
      if (size(value,1)/=size(rows) .or. size(value,2)/=size(cols)) return
      if (any(rows<1) .or. any(rows>nr) .or. any(cols<1) .or. any(cols>nc)) return
      cap=a%nnz()+size(rows)*size(cols)
      allocate(rr(cap),cc(cap),vv(cap)); p=0
      do i=1,nr
         do k=a%row_ptr(i),a%row_ptr(i+1)-1
            replace=any(rows==i) .and. any(cols==a%col_ind(k))
            if (replace) cycle
            p=p+1; rr(p)=i; cc(p)=a%col_ind(k); vv(p)=a%values(k)
         end do
      end do
      do j=1,size(cols)
         do i=1,size(rows)
            p=p+1; rr(p)=rows(i); cc(p)=cols(j); vv(p)=value(i,j)
         end do
      end do
      call csr_from_triplet(nr,nc,rr(:p),cc(:p),vv(:p),a,info)
   end subroutine csr_assign_dense_block

end module matrixextra_slice
