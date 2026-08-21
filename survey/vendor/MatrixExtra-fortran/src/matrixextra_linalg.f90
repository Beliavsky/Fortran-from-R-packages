! SPDX-License-Identifier: GPL-3.0-only
module matrixextra_linalg
   use matrix_kinds, only : dp
   use matrix_sparse, only : csr_matrix, csr_to_dense, csr_from_dense
   use matrixextra_matmul, only : csr_matvec_extra
   implicit none
   private
   public :: csr_norm, csr_diag, csr_set_diag

contains

   real(dp) function csr_norm(a,norm_type) result(v)
      type(csr_matrix), intent(in) :: a
      character(len=*), intent(in), optional :: norm_type
      character(len=1) :: t
      real(dp), allocatable :: sums(:),x(:),y(:),z(:)
      real(dp) :: old
      integer :: i,k,it
      t='O'; if (present(norm_type)) t=norm_type(1:1)
      select case(t)
      case('O','o','1')
         allocate(sums(a%ncol),source=0.0_dp)
         do k=1,a%nnz(); sums(a%col_ind(k))=sums(a%col_ind(k))+abs(a%values(k)); end do
         if (size(sums)>0) then; v=maxval(sums); else; v=0.0_dp; end if
      case('I','i')
         allocate(sums(a%nrow),source=0.0_dp)
         do i=1,a%nrow
            do k=a%row_ptr(i),a%row_ptr(i+1)-1; sums(i)=sums(i)+abs(a%values(k)); end do
         end do
         if (size(sums)>0) then; v=maxval(sums); else; v=0.0_dp; end if
      case('F','f')
         v=sqrt(sum(a%values*a%values))
      case('M','m')
         if (a%nnz()>0) then; v=maxval(abs(a%values)); else; v=0.0_dp; end if
      case('2')
         if (a%ncol==0 .or. a%nrow==0) then
            v=0.0_dp; return
         end if
         allocate(x(a%ncol),source=1.0_dp/sqrt(real(a%ncol,dp)))
         allocate(y(a%nrow),z(a%ncol))
         old=0.0_dp
         do it=1,200
            y=csr_matvec_extra(a,x)
            z=transpose_matvec(a,y)
            v=sqrt(sum(z*z))
            if (v<=tiny(1.0_dp)) then; v=0.0_dp; return; end if
            x=z/v
            v=sqrt(sum(y*y))
            if (abs(v-old)<=1.0e-12_dp*max(1.0_dp,v)) exit
            old=v
         end do
      case default
         v=-1.0_dp
      end select
   end function csr_norm

   function transpose_matvec(a,y) result(z)
      type(csr_matrix), intent(in) :: a
      real(dp), intent(in) :: y(:)
      real(dp), allocatable :: z(:)
      integer :: i,k
      allocate(z(a%ncol),source=0.0_dp)
      do i=1,a%nrow
         do k=a%row_ptr(i),a%row_ptr(i+1)-1
            z(a%col_ind(k))=z(a%col_ind(k))+a%values(k)*y(i)
         end do
      end do
   end function transpose_matvec

   function csr_diag(a) result(d)
      type(csr_matrix), intent(in) :: a
      real(dp), allocatable :: d(:)
      integer :: i,k,n
      n=min(a%nrow,a%ncol); allocate(d(n),source=0.0_dp)
      do i=1,n
         do k=a%row_ptr(i),a%row_ptr(i+1)-1
            if (a%col_ind(k)==i) then; d(i)=a%values(k); exit; end if
         end do
      end do
   end function csr_diag

   subroutine csr_set_diag(a,d)
      type(csr_matrix), intent(inout) :: a
      real(dp), intent(in) :: d(:)
      real(dp), allocatable :: x(:,:)
      integer :: i,n
      x=csr_to_dense(a); n=min(size(d),min(a%nrow,a%ncol))
      do i=1,n; x(i,i)=d(i); end do
      call csr_from_dense(x,a)
   end subroutine csr_set_diag

end module matrixextra_linalg
