program test_utils
  use quantreg, only : dp, rq_result, rq_fit_pfn, recursive_least_squares, combinations
  implicit none
  real(dp) :: x(120,2), y(120), bp(2,120), cov(2,2)
  integer, allocatable :: c(:,:)
  type(rq_result) :: fit
  integer :: i,info,nc
  do i=1,120
    x(i,1)=1.0_dp
    x(i,2)=real(i-60,dp)/20.0_dp
    y(i)=2.0_dp+0.5_dp*x(i,2)+0.1_dp*sin(real(i,dp))
  end do
  call rq_fit_pfn(x,y,0.5_dp,fit,123)
  if (fit%info/=0) error stop 'pfn info'
  if (maxval(abs(fit%coefficients-[2.0_dp,0.5_dp]))>0.05_dp) error stop 'pfn coef'
  call recursive_least_squares(x,y,bp,cov,info)
  if (info/=0) error stop 'rls info'
  if (maxval(abs(bp(:,120)-[2.0_dp,0.5_dp]))>0.05_dp) error stop 'rls coef'
  call combinations(5,2,c,nc)
  if (nc/=10 .or. size(c,2)/=10) error stop 'combinations'
  print *, 'test_utils: PASS'
end program
