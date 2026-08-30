program test_core
use deldir
implicit none
real(dp) :: x(5), y(5), expected_area, d
real(dp), allocatable :: cents(:,:)
type(deldir_result) :: res, cvres
type(voronoi_tile), allocatable :: tiles(:)
type(triangle3), allocatable :: tris(:)
type(dividing_segment), allocatable :: chain(:)
type(tile_info_result) :: ti
integer, allocatable :: nbr(:), tm(:,:)
integer :: i, stat, it
integer :: tags(5)
logical :: bnd
x=[0.0_dp,1.0_dp,1.0_dp,0.0_dp,0.5_dp]
y=[0.0_dp,0.0_dp,1.0_dp,1.0_dp,0.5_dp]
call deldir_compute(x,y,res,rw=[-0.2_dp,1.2_dp,-0.2_dp,1.2_dp],sort_points=.true.,status=stat)
call check(stat==0,'fit status')
call check(res%n_data==5,'n_data')
call check(size(res%delsgs)==8,'Delaunay edge count')
call check(size(res%dirsgs)==8,'Dirichlet shared edge count')
call check(abs(res%del_area-1.0_dp)<1.0e-12_dp,'convex hull area')
expected_area=(1.2_dp-(-0.2_dp))**2
call check(abs(res%dir_area-expected_area)<1.0e-11_dp,'window area')
call check(abs(sum([(res%summary(i)%dir_wt,i=1,5)])-1.0_dp)<1.0e-12_dp,'Dirichlet weights')
call deldir_tiles(res,tiles)
call check(size(tiles)==5,'tile count')
call check(abs(sum([(tiles(i)%area,i=1,5)])-expected_area)<1.0e-10_dp,'tile areas partition window')
call check(abs(tile_area(tiles(5))-0.5_dp)<1.0e-12_dp,'center tile area')
call check(all(abs(tile_centroid(tiles(5))-[0.5_dp,0.5_dp])<1.0e-12_dp),'center centroid')
call deldir_get_neighbors(res,5,nbr)
call check(size(nbr)==4,'center neighbors')
call deldir_triangle_matrix(res,tm)
call check(size(tm,1)==4,'triangle count')
call deldir_triangles(res,tris)
call check(size(tris)==4,'triangle-list count')
d=mean_nearest_neighbor_distance(x,y)
call check(abs(d-sqrt(0.5_dp))<1.0e-12_dp,'mean nearest neighbor distance')
call check(which_tile(0.49_dp,0.51_dp,tiles)==5,'which tile')
call check(inside_rect(0.5_dp,0.5_dp,[0.0_dp,1.0_dp,0.0_dp,1.0_dp]),'inside rectangle')
call check(inside_polygon(0.5_dp,0.5_dp,[0.0_dp,1.0_dp,1.0_dp,0.0_dp], &
    [0.0_dp,0.0_dp,1.0_dp,1.0_dp],on_boundary=bnd),'inside polygon')
call check(.not.bnd,'polygon interior boundary flag')
call check(inside_polygon(0.0_dp,0.5_dp,[0.0_dp,1.0_dp,1.0_dp,0.0_dp], &
    [0.0_dp,0.0_dp,1.0_dp,1.0_dp],on_boundary=bnd).and.bnd,'polygon boundary')
call centroidal_voronoi(res,cvres,cents,maxit=3,tol=1.0e-12_dp,iterations=it,status=stat)
call check(stat==0,'CVT status')
call check(size(cents,1)==5,'CVT centroids')
tags=[1,2,2,1,1]
call deldir_dividing_chain(res,tags,chain)
call check(size(chain)>0,'dividing chain')
call deldir_tile_info(res,include_boundary=.true.,info=ti)
call check(size(ti%tiles)==5,'tile info count')
call check(size(ti%all_edge_lengths)>0,'tile info edges')
call check(ti%total_perimeter>0.0_dp,'tile info perimeter')
call deldir_tile_info(res,include_boundary=.false.,info=ti)
call check(size(ti%tiles)==1 .and. ti%tiles(1)%point==5,'interior tile info')
print '(a)','test_core: PASS'
contains
subroutine check(ok,msg)
logical,intent(in)::ok
character(len=*),intent(in)::msg
if(.not.ok)then
 print '(a)','FAIL: '//trim(msg);error stop 1
end if
end subroutine
end program
