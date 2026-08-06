program test_dvine
  use tscopula
  use test_utils
  implicit none
  type(pair_copula)::pairs(2)
  type(dvine_copula)::d,ind
  real(dp),allocatable::u(:),res(:),innov(:)
  real(dp)::p,x,dens
  integer::i
  pairs(1)=paircop('gauss',0.55_dp);pairs(2)=paircop('clayton',1.0_dp)
  d=dvinecopula(pairs);call set_seed(9876);u=sim_dvine(d,120)
  call assert_true(all(u>0.0_dp).and.all(u<1.0_dp),'D-vine simulation support')
  dens=rblattdens(d,u(1:100),0.4_dp);call assert_true(dens>0.0_dp,'conditional density positive')
  p=rblatt(d,u(1:100),0.4_dp);x=irblatt(d,u(1:100),p);call assert_close(x,0.4_dp,8.0e-4_dp,'Rosenblatt inversion')
  res=resid_dvine(d,u);call assert_true(all(res>0.0_dp).and.all(res<1.0_dp),'D-vine residual support')
  ind=dvinecopula([paircop('indep'),paircop('indep')]);allocate(innov(50));do i=1,50;innov(i)=(real(i,dp)-0.5_dp)/50.0_dp;end do
  u=sim_dvine(ind,50,innov=innov);call assert_true(maxval(abs(u-innov))<5.0e-5_dp,'independent D-vine innovations')
  call assert_close(dvine_loglik(ind,u),0.0_dp,1.0e-9_dp,'independent D-vine likelihood')
  call pass('D-vine copulas')
end program test_dvine
