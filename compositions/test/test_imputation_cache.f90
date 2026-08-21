program test_imputation_cache
  use compositions
  implicit none
  real(dp) :: comp(4,3),pred(4,3),cov(3,3),norms(9,2),dl(4,3),h(3,3)
  integer :: mt(4,3),i
  type(clr_expectation_result) :: r
  comp=reshape([ &
    0.2_dp,0.0_dp,0.3_dp,0.25_dp, &
    0.3_dp,0.4_dp,0.0_dp,0.25_dp, &
    0.5_dp,0.6_dp,0.7_dp,0.50_dp],[4,3])
  mt=mt_observed; mt(2,1)=mt_bdl; mt(3,2)=mt_mar
  dl=0.0_dp; dl(2,1)=0.12_dp
  pred=0.0_dp
  do i=1,4
    pred(i,:)=log([0.25_dp,0.30_dp,0.45_dp])
    pred(i,:)=pred(i,:)-sum(pred(i,:))/3.0_dp
  end do
  h=-1.0_dp/3.0_dp; do i=1,3; h(i,i)=h(i,i)+1.0_dp; end do
  cov=0.12_dp*h
  norms=reshape([ &
    -2.0_dp,-1.5_dp,-1.0_dp,-0.5_dp,0.0_dp,0.5_dp,1.0_dp,1.5_dp,2.0_dp, &
     0.7_dp,-0.7_dp,0.3_dp,-0.3_dp,0.0_dp,1.0_dp,-1.0_dp,1.5_dp,-1.5_dp],[9,2])
  r=acomp_clr_expectation(comp,pred,mt,cov,norms,dl)
  if(.not.r%ok) error stop 'expectation not ok'
  if(abs(sum(r%mean_clr))>1.0e-11_dp) error stop 'mean not CLR centered'
  if(maxval(abs(matmul(r%covariance_clr,[1.0_dp,1.0_dp,1.0_dp])))>1.0e-10_dp) &
    error stop 'CLR covariance not centered'
  if(size(r%cache)<2) error stop 'missing-pattern cache not created'
  if(r%used(2)<0) error stop 'BDL Monte Carlo not exercised'
  if(abs(r%xlr_imputed(1,3))>1.0e-12_dp) error stop 'ALR reference not zero'
  print *, 'test_imputation_cache: PASS'
end program test_imputation_cache
