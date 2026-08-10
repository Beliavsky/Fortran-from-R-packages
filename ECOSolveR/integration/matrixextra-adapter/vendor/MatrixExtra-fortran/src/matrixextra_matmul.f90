! SPDX-License-Identifier: GPL-3.0-only
module matrixextra_matmul
   use matrix_kinds, only : dp
   use matrix_sparse, only : csr_matrix, csc_matrix, csr_from_triplet
   implicit none
   private
   public :: csr_matvec_extra, csr_dense_matmul, dense_csc_matmul
   public :: csr_csr_matmul, csr_crossprod, csr_tcrossprod
   public :: csr_sparse_vector_matmul

contains

   function csr_matvec_extra(a,x) result(y)
      type(csr_matrix), intent(in) :: a
      real(dp), intent(in) :: x(:)
      real(dp), allocatable :: y(:)
      integer :: i,k
      allocate(y(a%nrow),source=0.0_dp)
      if (size(x)/=a%ncol) return
      do i=1,a%nrow
         do k=a%row_ptr(i),a%row_ptr(i+1)-1
            y(i)=y(i)+a%values(k)*x(a%col_ind(k))
         end do
      end do
   end function csr_matvec_extra

   function csr_sparse_vector_matmul(a,idx,val,n) result(y)
      type(csr_matrix), intent(in) :: a
      integer, intent(in) :: idx(:), n
      real(dp), intent(in) :: val(:)
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: x(:)
      integer :: k
      allocate(x(n),source=0.0_dp)
      if (size(idx)==size(val)) then
         do k=1,size(idx)
            if (idx(k)>=1 .and. idx(k)<=n) x(idx(k))=x(idx(k))+val(k)
         end do
      end if
      y=csr_matvec_extra(a,x)
   end function csr_sparse_vector_matmul

   function csr_dense_matmul(a,b) result(c)
      type(csr_matrix), intent(in) :: a
      real(dp), intent(in) :: b(:,:)
      real(dp), allocatable :: c(:,:)
      integer :: i,j,k
      allocate(c(a%nrow,size(b,2)),source=0.0_dp)
      if (a%ncol/=size(b,1)) return
      do i=1,a%nrow
         do k=a%row_ptr(i),a%row_ptr(i+1)-1
            do j=1,size(b,2)
               c(i,j)=c(i,j)+a%values(k)*b(a%col_ind(k),j)
            end do
         end do
      end do
   end function csr_dense_matmul

   function dense_csc_matmul(a,b) result(c)
      real(dp), intent(in) :: a(:,:)
      type(csc_matrix), intent(in) :: b
      real(dp), allocatable :: c(:,:)
      integer :: i,j,k
      allocate(c(size(a,1),b%ncol),source=0.0_dp)
      if (size(a,2)/=b%nrow) return
      do j=1,b%ncol
         do k=b%col_ptr(j),b%col_ptr(j+1)-1
            do i=1,size(a,1)
               c(i,j)=c(i,j)+a(i,b%row_ind(k))*b%values(k)
            end do
         end do
      end do
   end function dense_csc_matmul

   subroutine csr_csr_matmul(a,b,c,info,tol)
      type(csr_matrix), intent(in) :: a,b
      type(csr_matrix), intent(out) :: c
      integer, intent(out), optional :: info
      real(dp), intent(in), optional :: tol
      real(dp), allocatable :: acc(:), vals(:)
      integer, allocatable :: mark(:), rows(:), cols(:), touched(:)
      integer :: i,ka,kb,j,p,nz,cap,istat,nt,q
      real(dp) :: eps
      eps=0.0_dp
      if (present(tol)) eps=max(0.0_dp,tol)
      if (a%ncol/=b%nrow) then
         call csr_from_triplet(0,0,[integer ::],[integer ::],[real(dp) ::],c,istat)
         if (present(info)) info=1
         return
      end if
      allocate(acc(b%ncol),source=0.0_dp)
      allocate(mark(b%ncol),source=0)
      allocate(touched(max(1,b%ncol)))
      cap=max(16,a%nnz()+b%nnz())
      allocate(rows(cap),cols(cap),vals(cap))
      nz=0
      do i=1,a%nrow
         nt=0
         do ka=a%row_ptr(i),a%row_ptr(i+1)-1
            q=a%col_ind(ka)
            do kb=b%row_ptr(q),b%row_ptr(q+1)-1
               j=b%col_ind(kb)
               if (mark(j)==0) then
                  nt=nt+1; touched(nt)=j; mark(j)=1
               end if
               acc(j)=acc(j)+a%values(ka)*b%values(kb)
            end do
         end do
         call sort_int(touched(:nt))
         do p=1,nt
            j=touched(p)
            if (abs(acc(j))>eps) then
               if (nz==cap) call grow_triplets(rows,cols,vals,cap)
               nz=nz+1; rows(nz)=i; cols(nz)=j; vals(nz)=acc(j)
            end if
            acc(j)=0.0_dp; mark(j)=0
         end do
      end do
      call csr_from_triplet(a%nrow,b%ncol,rows(:nz),cols(:nz),vals(:nz),c,istat)
      if (present(info)) info=istat
   end subroutine csr_csr_matmul

   function csr_crossprod(a,b) result(c)
      type(csr_matrix), intent(in) :: a,b
      real(dp), allocatable :: c(:,:)
      integer :: i,j,k
      allocate(c(a%ncol,b%ncol),source=0.0_dp)
      if (a%nrow/=b%nrow) return
      ! Sparse accumulation over common rows.
      do i=1,a%nrow
         do j=a%row_ptr(i),a%row_ptr(i+1)-1
            do k=b%row_ptr(i),b%row_ptr(i+1)-1
               c(a%col_ind(j),b%col_ind(k)) = c(a%col_ind(j),b%col_ind(k)) + &
                  a%values(j)*b%values(k)
            end do
         end do
      end do
   end function csr_crossprod

   function csr_tcrossprod(a,b) result(c)
      type(csr_matrix), intent(in) :: a,b
      real(dp), allocatable :: c(:,:)
      integer :: i,j,ka,kb,ca,cb
      allocate(c(a%nrow,b%nrow),source=0.0_dp)
      if (a%ncol/=b%ncol) return
      do i=1,a%nrow
         do j=1,b%nrow
            ka=a%row_ptr(i); kb=b%row_ptr(j)
            do while (ka<a%row_ptr(i+1) .and. kb<b%row_ptr(j+1))
               ca=a%col_ind(ka); cb=b%col_ind(kb)
               if (ca==cb) then
                  c(i,j)=c(i,j)+a%values(ka)*b%values(kb)
                  ka=ka+1; kb=kb+1
               else if (ca<cb) then
                  ka=ka+1
               else
                  kb=kb+1
               end if
            end do
         end do
      end do
   end function csr_tcrossprod

   subroutine sort_int(x)
      integer, intent(inout) :: x(:)
      integer :: i,j,t
      do i=2,size(x)
         t=x(i); j=i-1
         do while (j>=1)
            if (x(j)<=t) exit
            x(j+1)=x(j); j=j-1
         end do
         x(j+1)=t
      end do
   end subroutine sort_int

   subroutine grow_triplets(r,c,v,cap)
      integer, allocatable, intent(inout) :: r(:),c(:)
      real(dp), allocatable, intent(inout) :: v(:)
      integer, intent(inout) :: cap
      integer, allocatable :: rr(:),cc(:)
      real(dp), allocatable :: vv(:)
      integer :: old
      old=cap; cap=2*cap
      allocate(rr(cap),cc(cap),vv(cap))
      rr(:old)=r; cc(:old)=c; vv(:old)=v
      call move_alloc(rr,r); call move_alloc(cc,c); call move_alloc(vv,v)
   end subroutine grow_triplets

end module matrixextra_matmul
