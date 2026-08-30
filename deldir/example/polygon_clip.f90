program polygon_clip
use deldir
implicit none
real(dp) :: x(5)=[0.0_dp,1.0_dp,1.0_dp,0.0_dp,0.5_dp]
real(dp) :: y(5)=[0.0_dp,0.0_dp,1.0_dp,1.0_dp,0.5_dp]
type(deldir_result) :: fit
type(voronoi_tile), allocatable :: tiles(:)
type(poly_set) :: window
integer :: i, status

allocate(window%path(1))
window%path(1)=make_path([0.1_dp,0.9_dp,0.7_dp,0.3_dp], &
                         [0.1_dp,0.1_dp,0.9_dp,0.9_dp])

call deldir_compute(x,y,fit,rw=[0.0_dp,1.0_dp,0.0_dp,1.0_dp],status=status)
if(status/=0)error stop 'deldir_compute failed'
call deldir_tiles(fit,tiles,clipp=window,status=status)
if(status/=0)error stop 'polygon clipping failed'

print '(a,i0)','clipped tiles: ',size(tiles)
do i=1,size(tiles)
    print '(i3,2x,i2,2x,f10.6)',tiles(i)%point_number,tiles(i)%n_components(),tiles(i)%area
end do
end program polygon_clip
