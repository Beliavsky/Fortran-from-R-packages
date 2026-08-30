program test_edge_cases
use deldir
implicit none
real(dp)::x(7),y(7)
type(deldir_result)::r1,r2
logical,allocatable::dup(:)
real(dp),allocatable::xs(:),ys(:)
integer,allocatable::ind(:),rind(:)
integer::stat
x=[0.0_dp,1.0_dp,1.0_dp,0.0_dp,0.5_dp,0.5_dp,9.0_dp]
y=[0.0_dp,0.0_dp,1.0_dp,1.0_dp,0.5_dp,0.5_dp,9.0_dp]
call duplicated_xy(x,y,dup)
call check(count(dup)==1 .and. dup(6),'duplicate mask')
call deldir_compute(x,y,r1,rw=[-0.2_dp,1.2_dp,-0.2_dp,1.2_dp],sort_points=.true.,status=stat)
call check(stat==0.and.r1%n_data==5,'duplicate elimination/window clipping')
call deldir_compute(x,y,r2,rw=[-0.2_dp,1.2_dp,-0.2_dp,1.2_dp],sort_points=.false.,status=stat)
call check(stat==0.and.r2%n_data==5,'unsorted fit')
call check(abs(r1%del_area-r2%del_area)<1.0e-12_dp,'sorted/unsorted hull area')
call check(abs(r1%dir_area-r2%dir_area)<1.0e-12_dp,'sorted/unsorted tile area')
call deldir_bin_sort(x(1:5),y(1:5),[-0.2_dp,1.2_dp,-0.2_dp,1.2_dp],xs,ys,ind,rind)
call check(all(ind(rind)==[1,2,3,4,5]),'bin-sort inverse indices')
print '(a)','test_edge_cases: PASS'
contains
subroutine check(ok,msg)
logical,intent(in)::ok
character(len=*),intent(in)::msg
if(.not.ok)then;print '(a)','FAIL: '//trim(msg);error stop 1;end if
end subroutine
end program
