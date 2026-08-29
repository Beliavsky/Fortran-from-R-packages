program test_distributions
  use tmvtnorm
  use mvtnorm_distributions, only : dmvt_one
  implicit none
  real(dp) :: mu2(2),s2(2,2),lo2(2),up2(2),x2(2),p,d,marg,q
  real(dp) :: mu3(3),s3(3,3),lo3(3),up3(3)
  integer :: i
  mu2=0.0_dp
  s2=0.0_dp
  do i=1,2
  s2(i,i)=1.0_dp
  end do
  lo2=-1.0_dp
  up2=1.0_dp
  x2=0.0_dp
  d=dtmvnorm_one(x2,mu2,s2,lo2,up2)
  call assert_close(d,0.3414866223977847_dp,2.0e-6_dp,'joint normal density')
  p=ptmvnorm(lo2,[0.0_dp,0.0_dp],mu2,s2,lo2,up2)
  call assert_close(p,0.25_dp,2.0e-5_dp,'truncated normal rectangle probability')
  marg=dtmvnorm_marginal(0.0_dp,1,mu2,s2,lo2,up2)
  call assert_close(marg,0.5843685672568167_dp,2.0e-6_dp,'univariate marginal density')
  q=qtmvnorm_marginal(0.5_dp,1,mu2,s2,lo2,up2)
  call assert_close(q,0.0_dp,2.0e-4_dp,'marginal median')

  mu3=0.0_dp
  s3=0.0_dp
  do i=1,3
  s3(i,i)=1.0_dp
  end do
  lo3=-1.0_dp
  up3=1.0_dp
  d=dtmvnorm_marginal2(0.0_dp,0.0_dp,1,2,mu3,s3,lo3,up3)
  call assert_close(d,0.5843685672568167_dp**2,4.0e-6_dp,'bivariate marginal density')

  lo2=-huge(1.0_dp)
  up2=huge(1.0_dp)
  x2=[0.2_dp,-0.4_dp]
  d=dtmvt_one(x2,mu2,s2,5.0_dp,lo2,up2)
  call assert_close(d,dmvt_one(x2,mu2,s2,5.0_dp,.false.),2.0e-8_dp,'untruncated t density')
  print *, 'test_distributions: ok'
contains
  subroutine assert_close(x,y,tol,msg)
    real(dp),intent(in)::x,y,tol
    character(len=*),intent(in)::msg
    if(abs(x-y)>tol*max(1.0_dp,abs(y))) then
      print *, 'FAIL ',trim(msg),x,y
      error stop 1
    end if
  end subroutine
end program
