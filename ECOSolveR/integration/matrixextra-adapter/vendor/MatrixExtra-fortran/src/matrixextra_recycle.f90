! SPDX-License-Identifier: GPL-3.0-only
module matrixextra_recycle
   use matrix_kinds, only : dp
   use matrix_sparse, only : csr_matrix
   implicit none
   private
   public :: csr_multiply_vector, csr_divide_vector, csr_power_vector
   public :: csr_mod_vector, csr_intdiv_vector

contains

   subroutine csr_multiply_vector(a,v)
      type(csr_matrix), intent(inout) :: a
      real(dp), intent(in) :: v(:)
      integer :: i,k,idx
      if (size(v)==0) return
      do i=1,a%nrow
         do k=a%row_ptr(i),a%row_ptr(i+1)-1
            idx=1+mod((i-1)+(a%col_ind(k)-1)*a%nrow,size(v))
            a%values(k)=a%values(k)*v(idx)
         end do
      end do
   end subroutine csr_multiply_vector

   subroutine csr_divide_vector(a,v)
      type(csr_matrix), intent(inout) :: a
      real(dp), intent(in) :: v(:)
      integer :: i,k,idx
      if (size(v)==0) return
      do i=1,a%nrow
         do k=a%row_ptr(i),a%row_ptr(i+1)-1
            idx=1+mod((i-1)+(a%col_ind(k)-1)*a%nrow,size(v))
            a%values(k)=a%values(k)/v(idx)
         end do
      end do
   end subroutine csr_divide_vector

   subroutine csr_power_vector(a,v)
      type(csr_matrix), intent(inout) :: a
      real(dp), intent(in) :: v(:)
      integer :: i,k,idx
      if (size(v)==0) return
      do i=1,a%nrow
         do k=a%row_ptr(i),a%row_ptr(i+1)-1
            idx=1+mod((i-1)+(a%col_ind(k)-1)*a%nrow,size(v))
            a%values(k)=a%values(k)**v(idx)
         end do
      end do
   end subroutine csr_power_vector

   subroutine csr_mod_vector(a,v)
      type(csr_matrix), intent(inout) :: a
      real(dp), intent(in) :: v(:)
      integer :: i,k,idx
      if (size(v)==0) return
      do i=1,a%nrow
         do k=a%row_ptr(i),a%row_ptr(i+1)-1
            idx=1+mod((i-1)+(a%col_ind(k)-1)*a%nrow,size(v))
            a%values(k)=modulo(a%values(k),v(idx))
         end do
      end do
   end subroutine csr_mod_vector

   subroutine csr_intdiv_vector(a,v)
      type(csr_matrix), intent(inout) :: a
      real(dp), intent(in) :: v(:)
      integer :: i,k,idx
      if (size(v)==0) return
      do i=1,a%nrow
         do k=a%row_ptr(i),a%row_ptr(i+1)-1
            idx=1+mod((i-1)+(a%col_ind(k)-1)*a%nrow,size(v))
            a%values(k)=floor(a%values(k)/v(idx))
         end do
      end do
   end subroutine csr_intdiv_vector

end module matrixextra_recycle
