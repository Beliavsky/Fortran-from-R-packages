program test_algorithm478
  use l1pack
  implicit none
  integer, parameter :: n=12,p=3
  real(dp) :: x(n,p),y(n)
  real(dp), parameter :: expected_coef(3)=[ &
    7.482405544798950e-1_dp,-1.264420015657083_dp,2.974850629551447e-1_dp]
  real(dp), parameter :: expected_sad=4.742805814748581_dp
  integer :: i
  type(l1fit_result) :: fit

  do i=1,n
    x(i,1)=1.0_dp
    x(i,2)=real(i-6,dp)/3.0_dp
    x(i,3)=sin(real(i,dp))
    y(i)=0.7_dp-1.2_dp*x(i,2)+0.4_dp*x(i,3)+0.2_dp*cos(real(3*i,dp))
    if(mod(i,5)==0)y(i)=y(i)+2.0_dp
  end do
  call l1fit(x,y,fit,intercept=.false.)
  if(maxval(abs(fit%coefficients-expected_coef))>1.0e-11_dp)error stop 1
  if(abs(fit%minimum-expected_sad)>1.0e-11_dp)error stop 2
  if(fit%rank/=3.or.fit%info/=1.or.fit%iterations/=6)error stop 3
  print '(a)','test_algorithm478: PASS'
end program test_algorithm478
