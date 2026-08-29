program test_moments
  use tmvtnorm
  implicit none
  real(dp)::mu(2),s(2,2),lo(2),up(2)
  type(tmvnorm_moments_t)::m
  integer::i
  mu=0.0_dp
  s=0.0_dp
  do i=1,2
  s(i,i)=1.0_dp
  end do
  lo=-1.0_dp
  up=1.0_dp
  call mtmvnorm(mu,s,lo,up,m)
  if(.not.m%ok) error stop 'moments failed'
  call close(m%mean(1),0.0_dp,2e-5_dp,'mean1')
  call close(m%mean(2),0.0_dp,2e-5_dp,'mean2')
  call close(m%covariance(1,1),0.29112509477279314_dp,5e-5_dp,'var1')
  call close(m%covariance(2,2),0.29112509477279314_dp,5e-5_dp,'var2')
  call close(m%covariance(1,2),0.0_dp,5e-5_dp,'cov12')
  ! Partial truncation exercises the Johnson-Kotz reduction.
  s=reshape([1.0_dp,0.5_dp,0.5_dp,1.0_dp],[2,2])
  lo=[-1.0_dp,-huge(1.0_dp)]
  up=[1.0_dp,huge(1.0_dp)]
  call mtmvnorm(mu,s,lo,up,m)
  if(.not.m%ok) error stop 'partial moments failed'
  call close(m%mean(1),0.0_dp,5e-5_dp,'partial mean1')
  call close(m%mean(2),0.0_dp,5e-5_dp,'partial mean2')
  call close(m%covariance(1,1),0.29112509477279314_dp,8e-5_dp,'partial var1')
  call close(m%covariance(1,2),0.5_dp*0.29112509477279314_dp,8e-5_dp,'partial cov')
  print *, 'test_moments: ok'
contains
  subroutine close(x,y,tol,msg)
    real(dp),intent(in)::x,y,tol
    character(len=*),intent(in)::msg
    if(abs(x-y)>tol*max(1.0_dp,abs(y))) then
    print *,'FAIL ',trim(msg),x,y
    error stop 1
    end if
  end subroutine
end program
