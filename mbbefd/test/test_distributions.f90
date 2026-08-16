program test_distributions
  use mbbefd, only : dp, dmbbefd, pmbbefd, qmbbefd, dmbbefd_gb, pmbbefd_gb, &
    dgbeta, pgbeta, qgbeta, dstpareto, pstpareto, g2a, ecmbbefd, mmbbefd, tlmbbefd
  implicit none
  real(dp) :: a,b,g,x,tol
  tol=2.0e-10_dp;a=0.5_dp;b=0.3_dp;g=(a+b)/((a+1.0_dp)*b);x=0.4_dp
  call check(dmbbefd(x,a,b),0.4464754368263065_dp,tol,'dmbbefd')
  call check(pmbbefd(x,a,b),0.17096030533461248_dp,tol,'pmbbefd')
  call check(qmbbefd(0.25_dp,a,b),0.5757166424934449_dp,tol,'qmbbefd')
  call check(dmbbefd(1.0_dp,a,b),0.5625_dp,tol,'atom')
  call check(g2a(g,b),a,5.0e-14_dp,'g2a')
  call check(dmbbefd_gb(x,g,b),dmbbefd(x,a,b),tol,'gb density equivalence')
  call check(pmbbefd_gb(x,g,b),pmbbefd(x,a,b),tol,'gb cdf equivalence')
  call check(dgbeta(0.4_dp,1.5_dp,2.2_dp,3.1_dp),1.5539665867396075_dp,5.0e-10_dp,'dgbeta')
  call check(pgbeta(0.4_dp,1.5_dp,2.2_dp,3.1_dp),0.23395724424111222_dp,5.0e-10_dp,'pgbeta')
  call check(qgbeta(0.65_dp,1.5_dp,2.2_dp,3.1_dp),0.6223941041723167_dp,2.0e-9_dp,'qgbeta')
  call check(dstpareto(0.3_dp,2.3_dp),1.2142017005561896_dp,5.0e-12_dp,'dstpareto')
  call check(pstpareto(0.3_dp,2.3_dp),0.5685165703999067_dp,5.0e-12_dp,'pstpareto')
  if(abs(ecmbbefd(1.0_dp,a,b)-1.0_dp)>1.0e-12_dp) error stop 'exposure endpoint'
  if(mmbbefd(1.0_dp,a,b)<=0.0_dp .or. mmbbefd(1.0_dp,a,b)>=1.0_dp) error stop 'moment range'
  call check(tlmbbefd(a,b),0.5625_dp,tol,'total loss')
  print '(a)', 'test_distributions: PASS'
contains
  subroutine check(got,want,eps,label)
    real(dp),intent(in)::got,want,eps;character(len=*),intent(in)::label
    if(abs(got-want)>eps*(1.0_dp+abs(want)))then
      print *,trim(label),got,want;error stop 1
    end if
  end subroutine check
end program test_distributions
