program test_probability
  use matrixNormal
  use mvtnorm_special, only : normal_cdf
  use test_support
  implicit none
  real(dp)::m(2,2),u(2,2),v(2,2),lo(2,2),up(2,2),p,pleg,pref,plegref
  real(dp)::lov(4),upv(4),varc(4),varl(4)
  type(probability_result)::res
  integer::i
  m=0.0_dp
  u=reshape([1.0_dp,0.0_dp,0.0_dp,4.0_dp],[2,2])
  v=reshape([9.0_dp,0.0_dp,0.0_dp,16.0_dp],[2,2])
  lo=reshape([-1.0_dp,-2.0_dp,-3.0_dp,-4.0_dp],[2,2])
  up=reshape([0.5_dp,1.0_dp,2.0_dp,3.0_dp],[2,2])
  lov=vec(lo)
  upv=vec(up)
  varc=[9.0_dp,36.0_dp,16.0_dp,64.0_dp]
  varl=[9.0_dp,16.0_dp,36.0_dp,64.0_dp]
  pref=1.0_dp
  plegref=1.0_dp
  do i=1,4
    pref=pref*(normal_cdf(upv(i)/sqrt(varc(i)))-normal_cdf(lov(i)/sqrt(varc(i))))
    plegref=plegref*(normal_cdf(upv(i)/sqrt(varl(i)))-normal_cdf(lov(i)/sqrt(varl(i))))
  end do
  call assert_close(pref,0.006074900306399315_dp,5e-15_dp,'stored independent probability reference')
  call assert_close(plegref,0.0062069338267269954_dp,5e-15_dp,'stored legacy probability reference')
  res=pmatnorm(lo,up,m,u,v)
  p=res%value
  call assert_true(res%inform==0,'pmatnorm status: '//trim(res%message))
  call assert_close(p,pref,2e-5_dp,'pmatnorm correct V kron U ordering')
  res=pmatnorm(lo,up,m,u,v,legacy_covariance_order=.true.)
  pleg=res%value
  call assert_true(res%inform==0,'legacy pmatnorm status')
  call assert_close(pleg,plegref,2e-5_dp,'legacy U kron V ordering')
  res=pmatnorm(m,u,v)
  call assert_close(res%value,1.0_dp,1e-12_dp,'full support probability')
  res=pmatnorm(up,m,u,v)
  call assert_true(res%value>0.0_dp .and. res%value<1.0_dp,'upper-only CDF')
  write(*,'(a)') 'test_probability: PASS'
end program test_probability
