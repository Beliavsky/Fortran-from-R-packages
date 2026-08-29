program test_random
use truncnorm, only: dp, rtruncnorm, etruncnorm, vtruncnorm, set_seed_int
implicit none
integer, parameter :: n = 30000
real(dp), allocatable :: x(:)
real(dp) :: m, s2, e, v
integer :: i

call set_seed_int(20260826)
allocate(x(n))
do i=1,n
   x(i)=rtruncnorm(-1.0_dp,2.0_dp,0.3_dp,1.4_dp)
   if (x(i) <= -1.0_dp .or. x(i) >= 2.0_dp) error stop 'random bound failure'
end do
m=sum(x)/real(n,dp)
s2=sum((x-m)**2)/real(n-1,dp)
e=etruncnorm(-1.0_dp,2.0_dp,0.3_dp,1.4_dp)
v=vtruncnorm(-1.0_dp,2.0_dp,0.3_dp,1.4_dp)
if (abs(m-e) > 0.025_dp) then
   write(*,*) 'mean mismatch',m,e
   error stop 1
end if
if (abs(sqrt(s2)-sqrt(v)) > 0.025_dp) then
   write(*,*) 'sd mismatch',sqrt(s2),sqrt(v)
   error stop 1
end if

print *, 'test_random: PASS'
end program test_random
