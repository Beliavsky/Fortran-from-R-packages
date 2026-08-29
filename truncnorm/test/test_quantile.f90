program test_quantile
use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf, ieee_negative_inf
use truncnorm, only: dp, ptruncnorm, qtruncnorm
implicit none
real(dp), parameter :: probs(*) = [0.001_dp, 0.01_dp, 0.1_dp, 0.25_dp, 0.5_dp, 0.9_dp, 0.99_dp, 0.999_dp]
real(dp) :: p, q, pp, pinf, ninf
integer :: i
pinf = ieee_value(0.0_dp, ieee_positive_inf)
ninf = ieee_value(0.0_dp, ieee_negative_inf)

call check_case(-1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp)
call check_case(-1.0_dp, 1.0_dp, 4.0_dp, 1.0_dp)
call check_case(1.0_dp, 2.0_dp, 0.0_dp, 3.0_dp)
call check_case(-1.0_dp, pinf, 4.0_dp, 1.0_dp)
call check_case(ninf, 1.0_dp, 0.0_dp, 3.0_dp)

q = qtruncnorm(0.0_dp, -1.0_dp, 2.0_dp, 0.0_dp, 1.0_dp)
if (q /= -1.0_dp) error stop 'q(0) endpoint failed'
q = qtruncnorm(1.0_dp, -1.0_dp, 2.0_dp, 0.0_dp, 1.0_dp)
if (q /= 2.0_dp) error stop 'q(1) endpoint failed'

print *, 'test_quantile: PASS'
contains
subroutine check_case(a,b,mu,sd)
real(dp), intent(in) :: a,b,mu,sd
integer :: j
do j=1,size(probs)
   p = probs(j)
   q = qtruncnorm(p,a,b,mu,sd)
   pp = ptruncnorm(q,a,b,mu,sd)
   if (abs(pp-p) > 2.0e-10_dp) then
      write(*,'(a,5(1x,es20.12))') 'FAIL p/q', a,b,p,q,pp
      error stop 1
   end if
end do
end subroutine check_case
end program test_quantile
