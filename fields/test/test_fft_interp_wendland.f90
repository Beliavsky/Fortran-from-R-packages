program test_fft_interp_wendland
use fields
implicit none
real(dp)::x(5),y(5),z(5,5),center(1,2),coef(1)
real(dp),allocatable::h(:,:)
type(fft_interp_result)::q
integer::i,j
x=[(real(i-1,dp),i=1,5)];y=x
do i=1,5;do j=1,5;z(i,j)=2.0_dp+0.5_dp*cos(2.0_dp*acos(-1.0_dp)*real(i-1,dp)/5.0_dp);end do;end do
q=fft_interp_surface(x,y,z,3)
call check(size(q%z,1)==13 .and. size(q%z,2)==13,'fft interp shape')
do i=1,5
 call check(maxval(abs(q%z(1+(i-1)*3,1:13:3)-z(i,:)))<1e-10_dp,'fft interp knots')
end do
center(1,:)=[2.0_dp,2.0_dp];coef=1.0_dp
h=mult_wendland_grid(x,y,center,1.5_dp,coef)
call check(abs(h(3,3)-1.0_dp)<1e-14_dp,'wendland center')
call check(h(1,1)==0.0_dp,'wendland support')
print *,'test_fft_interp_wendland: PASS'
contains
subroutine check(ok,msg);logical,intent(in)::ok;character(*),intent(in)::msg
if(.not.ok)then;print *,'FAIL ',trim(msg);error stop;end if;end subroutine
end program
