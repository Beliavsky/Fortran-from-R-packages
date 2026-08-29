program test_recycle
use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf
use truncnorm, only: dp, dtruncnorm, dtruncnorm_recycle, ptruncnorm_recycle, qtruncnorm_recycle, &
   etruncnorm_recycle, vtruncnorm_recycle
implicit none
real(dp), allocatable :: y(:), y2(:), p(:), q(:), e(:), v(:)
real(dp) :: means(3), pinf
integer :: i
pinf=ieee_value(0.0_dp,ieee_positive_inf)
means=[-1.0_dp,0.0_dp,1.0_dp]
y=dtruncnorm_recycle([1.0_dp],[0.0_dp],[pinf],means,[1.0_dp])
if (size(y)/=3) error stop 'recycle length'
do i=1,3
   if (abs(y(i)-dtruncnorm(1.0_dp,0.0_dp,pinf,means(i),1.0_dp)) > 1e-14_dp) error stop 'density recycle'
end do
y2=dtruncnorm([0.0_dp,0.5_dp,1.0_dp],0.0_dp,2.0_dp,0.0_dp,1.0_dp)
if (size(y2)/=3) error stop 'generic vector density length'
p=ptruncnorm_recycle([0.5_dp],[0.0_dp],[pinf],means,[1.0_dp])
q=qtruncnorm_recycle([0.5_dp],[0.0_dp],[pinf],means,[1.0_dp])
e=etruncnorm_recycle([0.0_dp],[pinf],means,[1.0_dp])
v=vtruncnorm_recycle([0.0_dp],[pinf],means,[1.0_dp])
if (any(p < 0.0_dp) .or. any(p > 1.0_dp)) error stop 'cdf recycle'
if (any(q < 0.0_dp)) error stop 'quantile recycle'
if (any(v <= 0.0_dp)) error stop 'variance recycle'
if (size(e)/=3) error stop 'moment recycle'
print *, 'test_recycle: PASS'
end program test_recycle
