program test_grid_stats_variogram
use fields, only: dp,bilinear_interpolate,describe_values,variogram_grid_result,variogram_grid,points_in_polygon
implicit none
real(dp)::xg(3),yg(3),z(3,3),loc(2,2),s(9),dat(4,4),poly(4,2),pts(3,2)
real(dp),allocatable::v(:)
logical,allocatable::inside(:)
type(variogram_grid_result)::vg
integer::i,j
xg=[0._dp,1._dp,2._dp];yg=xg
do j=1,3;do i=1,3;z(i,j)=xg(i)+2*yg(j);end do;end do
loc=reshape([0.5_dp,1.5_dp,0.5_dp,1.25_dp],[2,2]);v=bilinear_interpolate(xg,yg,z,loc)
call check(maxval(abs(v-(loc(:,1)+2*loc(:,2))))<1e-12_dp,'bilinear linear surface')
s=describe_values([1._dp,2._dp,3._dp,4._dp]);call check(abs(s(2)-2.5_dp)<1e-12_dp,'describe mean')
do j=1,4;do i=1,4;dat(i,j)=real(i+j,dp);end do;end do
vg=variogram_grid(dat,1._dp,1._dp,2._dp)
call check(size(vg%distance)>0 .and. all(vg%variogram>=0._dp),'variogram')
poly=reshape([0._dp,1._dp,1._dp,0._dp,0._dp,0._dp,1._dp,1._dp],[4,2])
pts=reshape([0.5_dp,1.5_dp,0.0_dp,0.5_dp,0.5_dp,0.2_dp],[3,2]);inside=points_in_polygon(pts,poly)
call check(inside(1) .and. .not.inside(2),'polygon')
print *,'test_grid_stats_variogram: PASS'
contains
subroutine check(ok,msg)
logical,intent(in)::ok;character(len=*),intent(in)::msg
if(.not.ok)then;print *,'FAIL: ',msg;error stop 1;end if
end subroutine
end program
