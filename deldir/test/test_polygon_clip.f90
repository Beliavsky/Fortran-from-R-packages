program test_polygon_clip
use deldir
implicit none
real(dp) :: x(5), y(5), area_sum
real(dp), allocatable :: cents(:,:)
type(deldir_result) :: res, cvres
type(voronoi_tile), allocatable :: tiles(:)
type(poly_set) :: clip
type(tile_info_result) :: info
integer :: i, kcenter, stat, it

x=[0.0_dp,1.0_dp,1.0_dp,0.0_dp,0.5_dp]
y=[0.0_dp,0.0_dp,1.0_dp,1.0_dp,0.5_dp]
call deldir_compute(x,y,res,rw=[0.0_dp,1.0_dp,0.0_dp,1.0_dp],sort_points=.true.,status=stat)
call check(stat==0,'base fit')

! Clipping by the rectangular window itself should preserve the partition.
allocate(clip%path(1))
clip%path(1)=make_path([0.0_dp,1.0_dp,1.0_dp,0.0_dp],[0.0_dp,0.0_dp,1.0_dp,1.0_dp])
call deldir_tiles(res,tiles,clipp=clip,clip_eps=1.0e-10_dp,status=stat)
call check(stat==0,'rectangular clip status')
call check(size(tiles)==5,'rectangular clip tile count')
area_sum=sum([(tiles(i)%area,i=1,size(tiles))])
call check(abs(area_sum-1.0_dp)<5.0e-9_dp,'rectangular clip area')

! A disconnected clipping set can split one Voronoi cell into two components.
deallocate(clip%path);allocate(clip%path(2))
clip%path(1)=make_path([0.2_dp,0.4_dp,0.4_dp,0.2_dp],[0.4_dp,0.4_dp,0.6_dp,0.6_dp])
clip%path(2)=make_path([0.6_dp,0.8_dp,0.8_dp,0.6_dp],[0.4_dp,0.4_dp,0.6_dp,0.6_dp])
call deldir_tiles(res,tiles,clipp=clip,clip_eps=1.0e-10_dp,status=stat)
call check(stat==0,'disconnected clip status')
call check(size(tiles)==1,'empty clipped tiles dropped')
call check(tiles(1)%point_number==5,'retained point number')
call check(tiles(1)%n_components()==2,'multi-component tile')
call check(.not.allocated(tiles(1)%x),'multi-component top-level coordinates')
call check(abs(tiles(1)%area-0.08_dp)<5.0e-9_dp,'multi-component area')
call check(all(abs(tile_centroid(tiles(1))-[0.5_dp,0.5_dp])<5.0e-9_dp),'multi-component centroid')
call deldir_tile_info(res,include_boundary=.false.,info=info,clipp=clip,clip_eps=1.0e-10_dp,status=stat)
call check(stat==0,'clipped tileInfo status')
call check(size(info%tiles)==1,'clipped tileInfo count')
call check(info%tiles(1)%point==5,'clipped tileInfo point')
call check(info%tiles(1)%n_edges==8,'clipped tileInfo edge count')
call check(abs(info%areas(1)-0.08_dp)<5.0e-9_dp,'clipped tileInfo area')

! Even-odd polygon sets support holes; signed component areas preserve hole area.
deallocate(clip%path);allocate(clip%path(2))
clip%path(1)=make_path([0.1_dp,0.9_dp,0.9_dp,0.1_dp],[0.1_dp,0.1_dp,0.9_dp,0.9_dp])
clip%path(2)=make_path([0.4_dp,0.4_dp,0.6_dp,0.6_dp],[0.4_dp,0.6_dp,0.6_dp,0.4_dp])
call deldir_tiles(res,tiles,clipp=clip,clip_eps=1.0e-10_dp,status=stat)
call check(stat==0,'hole clip status')
area_sum=sum([(tiles(i)%area,i=1,size(tiles))])
call check(abs(area_sum-0.60_dp)<5.0e-9_dp,'hole clip partition area')
kcenter=0
do i=1,size(tiles)
    if(tiles(i)%point_number==5)then
        kcenter=i
        exit
    end if
end do
call check(kcenter>0,'hole clip center retained')
call check(tiles(kcenter)%n_components()==2,'hole represented by signed ring')
call check(abs(tiles(kcenter)%area-0.42_dp)<5.0e-9_dp,'center area with hole')
call check(all(abs(tile_centroid(tiles(kcenter))-[0.5_dp,0.5_dp])<5.0e-9_dp),'hole centroid')

! A clipping set disjoint from rw produces an empty tile list rather than placeholders.
deallocate(clip%path);allocate(clip%path(1))
clip%path(1)=make_path([2.0_dp,3.0_dp,3.0_dp,2.0_dp],[2.0_dp,2.0_dp,3.0_dp,3.0_dp])
call deldir_tiles(res,tiles,clipp=clip,clip_eps=1.0e-10_dp,status=stat)
call check(stat==0 .and. size(tiles)==0,'empty clipped tessellation')
call centroidal_voronoi(res,cvres,cents,maxit=2,status=stat,clipp=clip,clip_eps=1.0e-10_dp)
call check(stat==8 .and. size(cents,1)==0,'empty polygon CVT status')

! Polygon-clipped Lloyd iteration is a native extension using the same clip set each iteration.
deallocate(clip%path);allocate(clip%path(1))
clip%path(1)=make_path([0.1_dp,0.9_dp,0.9_dp,0.1_dp],[0.1_dp,0.1_dp,0.9_dp,0.9_dp])
call centroidal_voronoi(res,cvres,cents,maxit=2,tol=1.0e-14_dp,iterations=it,status=stat, &
    clipp=clip,clip_eps=1.0e-10_dp)
call check(stat==0,'polygon CVT status')
call check(size(cents,2)==2 .and. size(cents,1)==cvres%n_data,'polygon CVT dimensions')
call check(all(cents(:,1)>=0.1_dp-1.0e-8_dp .and. cents(:,1)<=0.9_dp+1.0e-8_dp),'polygon CVT x bounds')
call check(all(cents(:,2)>=0.1_dp-1.0e-8_dp .and. cents(:,2)<=0.9_dp+1.0e-8_dp),'polygon CVT y bounds')

print '(a)','test_polygon_clip: PASS'
contains
subroutine check(ok,msg)
logical,intent(in)::ok
character(len=*),intent(in)::msg
if(.not.ok)then
    print '(a)','FAIL: '//trim(msg)
    error stop 1
end if
end subroutine check
end program test_polygon_clip
