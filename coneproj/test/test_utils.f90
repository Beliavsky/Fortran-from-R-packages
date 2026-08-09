program test_utils
   use coneproj
   implicit none
   real(dp) :: x(3,2), edges(2,3)
   type(qr_result) :: qr
   integer, allocatable :: keep(:), red(:), eq(:)
   integer :: status
   x = reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp,1.0_dp,1.0_dp],[3,2])
   call qr_decomp(x,qr)
   if (qr%rank /= 2) error stop 'qr rank'
   edges(:,1)=[1.0_dp,0.0_dp]
   edges(:,2)=[2.0_dp,0.0_dp]
   edges(:,3)=[0.0_dp,1.0_dp]
   call check_irreducible(edges,keep,red,eq,status)
   if (status /= coneproj_success) error stop 'check_irreducible status'
   if (size(red) < 1) error stop 'expected reducible edge'
   print *, 'test_utils: PASS'
end program test_utils
