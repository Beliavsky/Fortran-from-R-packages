program test_covariance_distance
use fields, only: dp,fields_rdist,fields_rdist_earth,exponential,matern,wendland,wendland_covariance
implicit none
real(dp) :: x(2,2),y(1,2),d1
real(dp),allocatable :: d(:,:),k(:,:)
x=reshape([0.0_dp,3.0_dp,0.0_dp,4.0_dp],[2,2])
y=reshape([0.0_dp,0.0_dp],[1,2])
d=fields_rdist(x,y)
call check(abs(d(1,1))<1e-12_dp,'rdist origin')
call check(abs(d(2,1)-5.0_dp)<1e-12_dp,'rdist 3-4-5')
call check(abs(exponential(2.0_dp,2.0_dp)-exp(-1.0_dp))<1e-13_dp,'exponential')
call check(abs(matern(2.0_dp,2.0_dp,0.5_dp)-exp(-1.0_dp))<1e-13_dp,'matern nu=.5')
d1=0.3_dp
call check(abs(wendland(d1,1.0_dp,2,2)-(1-d1)**6*(35*d1*d1+18*d1+3)/3)<2e-12_dp,'wendland 2,2')
k=wendland_covariance(x,x,10.0_dp,2)
call check(maxval(abs(k-transpose(k)))<1e-13_dp,'wendland symmetry')
call check(maxval(abs([(k(1,1)-1.0_dp),(k(2,2)-1.0_dp)]))<1e-13_dp,'wendland diagonal')
print *,'test_covariance_distance: PASS'
contains
subroutine check(ok,msg)
logical,intent(in)::ok;character(len=*),intent(in)::msg
if(.not.ok) then; print *,'FAIL: ',msg; error stop 1; end if
end subroutine
end program
