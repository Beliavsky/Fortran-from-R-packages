program test_sampler_branches
use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf, ieee_negative_inf, ieee_is_finite
use truncnorm, only: dp, rtruncnorm, etruncnorm, vtruncnorm, set_seed_int
implicit none
real(dp) :: pinf, ninf
pinf = ieee_value(0.0_dp, ieee_positive_inf)
ninf = ieee_value(0.0_dp, ieee_negative_inf)
call set_seed_int(90731)
call check_case(ninf, pinf, 0.0_dp, 1.0_dp, 'ordinary normal')
call check_case(-1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 'central two-sided')
call check_case(0.5_dp, 3.0_dp, 0.0_dp, 1.0_dp, 'positive half-normal branch')
call check_case(1.5_dp, 3.0_dp, 0.0_dp, 1.0_dp, 'positive exponential branch')
call check_case(-3.0_dp, -0.5_dp, 0.0_dp, 1.0_dp, 'negative symmetry branch')
call check_case(0.2_dp, pinf, 0.0_dp, 1.0_dp, 'left normal-rejection branch')
call check_case(1.5_dp, pinf, 0.0_dp, 1.0_dp, 'left exponential branch')
call check_case(ninf, -1.5_dp, 0.0_dp, 1.0_dp, 'right exponential branch')
print *, 'test_sampler_branches: PASS'
contains
subroutine check_case(a,b,mu,sd,label)
real(dp), intent(in) :: a,b,mu,sd
character(*), intent(in) :: label
integer, parameter :: n=8000
real(dp), allocatable :: x(:)
real(dp) :: m,s,e,v
x = rtruncnorm(n,a,b,mu,sd)
if (size(x) /= n) error stop 'rtruncnorm(n) size mismatch'
if (ieee_is_finite(a)) then
   if (any(x <= a)) then
      write(*,*) 'lower bound failure: ',trim(label)
      error stop 1
   end if
end if
if (ieee_is_finite(b)) then
   if (any(x >= b)) then
      write(*,*) 'upper bound failure: ',trim(label)
      error stop 1
   end if
end if
m=sum(x)/real(n,dp)
s=sqrt(sum((x-m)**2)/real(n-1,dp))
e=etruncnorm(a,b,mu,sd)
v=vtruncnorm(a,b,mu,sd)
if (abs(m-e) > 0.055_dp*max(sd,0.25_dp)) then
   write(*,*) 'mean failure: ',trim(label),m,e
   error stop 1
end if
if (abs(s-sqrt(v)) > 0.055_dp*max(sd,0.25_dp)) then
   write(*,*) 'sd failure: ',trim(label),s,sqrt(v)
   error stop 1
end if
end subroutine check_case
end program test_sampler_branches
