program test_paircopula
  use tscopula
  use test_utils
  implicit none
  type(pair_copula)::c
  real(dp)::u,v,p,x,tau
  u=0.31_dp;v=0.72_dp
  c=paircop('indep');call assert_close(pair_cdf(c,u,v),u*v,1.0e-12_dp,'independence CDF')
  call assert_close(pair_density(c,u,v),1.0_dp,1.0e-12_dp,'independence density')
  call assert_close(pair_h2(c,u,v),u,1.0e-12_dp,'independence h')
  c=paircop('gauss',0.65_dp);p=0.37_dp;x=pair_hinv2(c,p,v);call assert_close(pair_h2(c,x,v),p,2.0e-8_dp,'Gaussian inverse h')
  tau=pair_kendall(c);call assert_close(kendall_to_parameter(cop_gauss,tau),0.65_dp,1.0e-10_dp,'Gaussian tau conversion')
  c=paircop('clayton',2.0_dp);call assert_true(pair_density(c,u,v)>0.0_dp,'Clayton positive density')
  x=pair_hinv2(c,p,v);call assert_close(pair_h2(c,x,v),p,3.0e-7_dp,'Clayton inverse h')
  c=paircop('frank',3.0_dp);x=pair_cdf(c,u,v);call assert_true(x>=0.0_dp.and.x<=min(u,v),'Frank CDF bounds')
  c=paircop('gumbel',1.8_dp,rotation=180);call assert_true(pair_density(c,u,v)>0.0_dp,'rotated Gumbel density')
  c=paircop('t',0.4_dp,6.0_dp);call assert_true(pair_density(c,u,v)>0.0_dp,'Student copula density')
  call pass('pair copulas')
end program test_paircopula
