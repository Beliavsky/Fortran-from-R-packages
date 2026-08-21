program test_weighted_pairwise
  use bzinb
  use test_support
  implicit none
  integer,parameter::n=100
  integer::x(n),y(n),data(2,n),fail
  real(dp)::p(4),par(9),w,rho(2,2),ser(2,2)
  logical::conv(2,2)
  fail=0;p=[0.6_dp,0.15_dp,0.15_dp,0.10_dp];par=[1.2_dp,0.8_dp,1.0_dp,0.7_dp,0.9_dp,p]
  call set_bzinb_seed(314);call rbzinb_sample(n,par(1),par(2),par(3),par(4),par(5),p,x,y)
  w=weighted_pearson_correlation(x,y,par)
  call assert_true('weighted correlation finite',abs(w)<=1.000001_dp,fail)
  data(1,:)=x;data(2,:)=y
  call pairwise_bzinb(data,rho,ser,conv,maxiter=500)
  call assert_close('pair diag',rho(1,1),1.0_dp,0.0_dp,fail)
  call assert_close('pair symmetry',rho(1,2),rho(2,1),1e-14_dp,fail)
  if(fail==0)then;print *,'test_weighted_pairwise: PASS';else;error stop 1;end if
end program
