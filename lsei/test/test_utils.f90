program test_utils
   use lsei
   implicit none
   real(dp) :: x(7),v(6),m(3,2)
   integer :: ind(7)
   real(dp), allocatable :: r(:),c(:)
   x=[0d0,1d0,2d0,3d0,4d0,5d0,6d0]; v=[1d0,2d0,3d0,4d0,5d0,5d0]
   call indx(x,v,ind)
   if (any(ind/=[0,1,2,3,4,6,6])) error stop 1
   m=reshape([1d0,5d0,-2d0,4d0,2d0,8d0],[3,2])
   r=mat_maxs(m,1); c=mat_maxs(m,2)
   if (maxval(abs(r-[4d0,5d0,8d0]))>0d0) error stop 2
   if (maxval(abs(c-[5d0,8d0]))>0d0) error stop 3
   print *, 'PASS test_utils'
end program
