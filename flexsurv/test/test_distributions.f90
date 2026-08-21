program test_distributions
  use flexsurv, only : dp, dexp_fs, pexp_fs, qexp_fs, dweibull_fs, pweibull_fs, &
    qweibull_fs, dgompertz, pgompertz, qgompertz, pgengamma, qgengamma, pgenf, qgenf, &
    dist_mean, dist_rmst, dist_exponential, dist_weibull
  implicit none
  integer::fails
  real(dp)::x,p,q
  fails=0
  call chk(dexp_fs(2.0_dp,0.5_dp),0.5_dp*exp(-1.0_dp),1e-12_dp)
  call chk(pexp_fs(2.0_dp,0.5_dp),1.0_dp-exp(-1.0_dp),1e-12_dp)
  call chk(qexp_fs(0.7_dp,0.5_dp),-log(0.3_dp)/0.5_dp,1e-12_dp)
  x=1.7_dp;p=pweibull_fs(x,1.4_dp,2.3_dp);q=qweibull_fs(p,1.4_dp,2.3_dp);call chk(q,x,2e-10_dp)
  x=1.2_dp;p=pgompertz(x,0.15_dp,0.7_dp);q=qgompertz(p,0.15_dp,0.7_dp);call chk(q,x,2e-9_dp)
  x=2.0_dp;p=pgengamma(x,0.2_dp,0.8_dp,0.7_dp);q=qgengamma(p,0.2_dp,0.8_dp,0.7_dp);call chk(q,x,2e-7_dp)
  x=1.4_dp;p=pgenf(x,0.1_dp,0.8_dp,0.3_dp,0.7_dp);q=qgenf(p,0.1_dp,0.8_dp,0.3_dp,0.7_dp);call chk(q,x,5e-6_dp)
  call chk(dist_mean(dist_exponential,[0.5_dp]),2.0_dp,1e-12_dp)
  call chk(dist_rmst(dist_exponential,3.0_dp,[0.5_dp]),2.0_dp*(1.0_dp-exp(-1.5_dp)),2e-9_dp)
  call chk(dist_mean(dist_weibull,[1.0_dp,2.0_dp]),2.0_dp,1e-10_dp)
  if(fails>0)error stop 1
  print *,'test_distributions: PASS'
contains
  subroutine chk(a,b,tol)
    real(dp),intent(in)::a,b,tol
    if(abs(a-b)>tol*(1.0_dp+abs(b)))then
      print *,'FAIL ',a,b;fails=fails+1
    end if
  end subroutine chk
end program test_distributions
