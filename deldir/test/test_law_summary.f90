program test_law_summary
use deldir
implicit none
integer,parameter::m=10,n=m*m
real(dp)::x(n),y(n)
type(deldir_result)::fit
type(law_summary_result)::law
integer::i,j,k,stat
k=0
do j=0,m-1;do i=0,m-1
 k=k+1;x(k)=real(i,dp)+0.021_dp*sin(real(3*k+1,dp));y(k)=real(j,dp)+0.019_dp*cos(real(5*k+2,dp))
end do;end do
call deldir_compute(x,y,fit,rw=[-0.5_dp,9.5_dp,-0.5_dp,9.5_dp],status=stat)
call check(stat==0,'grid fit')
call deldir_law_summary(fit,law)
call check(size(law%layer1)>0,'law layer 1')
call check(size(law%layer2)>0,'law layer 2')
call check(size(law%layer3)>0,'law layer 3')
call check(size(law%kept)>0,'law interior kept')
call check(all(law%tile_areas>0.0_dp),'law areas')
print '(a)','test_law_summary: PASS'
contains
subroutine check(ok,msg)
logical,intent(in)::ok
character(len=*),intent(in)::msg
if(.not.ok)then;print '(a)','FAIL: '//trim(msg);error stop 1;end if
end subroutine
end program
