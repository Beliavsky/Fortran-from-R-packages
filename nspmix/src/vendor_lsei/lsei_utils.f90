! SPDX-License-Identifier: GPL-2.0-or-later
module lsei_utils
   use lsei_kinds, only : dp
   implicit none
   private
   public :: indx, mat_maxs
contains
   subroutine indx(x,v,ind)
      real(dp), intent(in) :: x(:),v(:)
      integer, intent(out) :: ind(:)
      integer :: i,left,right,mid,n
      n=size(v)
      if(size(ind)/=size(x)) return
      do i=1,size(x)
         if(n==0) then; ind(i)=0; cycle; end if
         if(x(i)<v(1)) then; ind(i)=0; cycle; end if
         if(x(i)>=v(n)) then; ind(i)=n; cycle; end if
         left=1; right=n
         do while(left<right-1)
            mid=nint(0.5_dp*real(left+right,dp))
            if(x(i)>=v(mid)) then; left=mid; else; right=mid; end if
         end do
         ind(i)=left
      end do
   end subroutine indx

   function mat_maxs(x,dim) result(v)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in), optional :: dim
      real(dp), allocatable :: v(:)
      integer :: d
      d=1; if(present(dim)) d=dim
      if(d==1) then
         allocate(v(size(x,1))); if(size(x,2)>0) v=maxval(x,dim=2)
      else
         allocate(v(size(x,2))); if(size(x,1)>0) v=maxval(x,dim=1)
      end if
   end function mat_maxs
end module lsei_utils
