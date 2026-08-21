program test_pairwise_full
  use bzinb
  use test_support
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  integer,parameter::n=30
  integer::data(3,n),x(n),y(n),fail
  real(dp)::p(4)
  type(pairwise_bzinb_result)::r
  fail=0;p=[0.55_dp,0.15_dp,0.20_dp,0.10_dp]
  call set_bzinb_seed(777)
  call rbzinb_sample(n,1.2_dp,0.8_dp,1.0_dp,0.8_dp,1.0_dp,p,x,y)
  data(1,:)=x;data(2,:)=y
  call rbzinb_sample(n,0.9_dp,1.1_dp,0.7_dp,1.0_dp,0.9_dp,p,x,y)
  data(3,:)=y
  call pairwise_bzinb_full(data,r,full_param=.true.,maxiter=450)
  call assert_true('three pairs',r%npairs==3,fail)
  call assert_true('pair order',all(r%first==[1,1,2]).and.all(r%second==[2,3,3]),fail)
  call assert_true('rho finite',all(ieee_is_finite(r%rho)).and.all(r%rho>0.0_dp).and.all(r%rho<1.0_dp),fail)
  call assert_true('full params positive',all(r%param(1:5,:)>0.0_dp),fail)
  call assert_true('full p sums',maxval(abs(sum(r%param(6:9,:),dim=1)-1.0_dp))<1.0e-10_dp,fail)
  call assert_close('nonzero first',r%nonzero_first(1),real(count(data(1,:)/=0),dp)/real(n,dp),1.0e-14_dp,fail)
  call set_bzinb_seed(991)
  call pairwise_bzinb_full(data,r,full_param=.false.,nsample=2,maxiter=250)
  call assert_true('sample two pairs',r%npairs==2,fail)
  call assert_true('sample distinct',r%first(1)/=r%first(2).or.r%second(1)/=r%second(2),fail)
  if(fail==0)then;print *,'test_pairwise_full: PASS';else;error stop 1;end if
end program test_pairwise_full
