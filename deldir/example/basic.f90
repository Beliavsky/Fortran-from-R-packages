program basic
use deldir
implicit none
real(dp) :: x(5)=[0.0_dp,1.0_dp,1.0_dp,0.0_dp,0.5_dp]
real(dp) :: y(5)=[0.0_dp,0.0_dp,1.0_dp,1.0_dp,0.5_dp]
type(deldir_result) :: fit
type(voronoi_tile), allocatable :: tiles(:)
integer :: i
call deldir_compute(x,y,fit,rw=[0.0_dp,1.0_dp,0.0_dp,1.0_dp])
call deldir_tiles(fit,tiles)
print '(a,i0)','Delaunay edges: ',size(fit%delsgs)
print '(a,i0)','Voronoi shared edges: ',size(fit%dirsgs)
do i=1,size(tiles)
 print '(i3,2f10.4,f12.6)',i,fit%summary(i)%x,fit%summary(i)%y,tiles(i)%area
end do
end program
