program test_bnb_bzinb
  use bzinb
  use test_support
  implicit none
  integer,parameter::n=180
  integer::x(n),y(n),fail
  real(dp)::p4(4),init9(9),ll0
  type(bnb_fit_result)::fb
  type(bzinb_fit_result)::fz
  fail=0;call set_bzinb_seed(90210)
  call rbnb_sample(n,1.4_dp,0.9_dp,1.1_dp,0.8_dp,1.0_dp,x,y)
  fb=fit_bnb(x,y,maxiter=900,initial=[1.3_dp,1.0_dp,1.0_dp,0.9_dp,0.9_dp])
  call assert_true('bnb finite',all(fb%param>0.0_dp).and.fb%loglik<0.0_dp,fail)
  call assert_close('bnb rho',fb%rho,true_correlation(1.4_dp,0.9_dp,1.1_dp,0.8_dp,1.0_dp),0.13_dp,fail)
  p4=[0.55_dp,0.15_dp,0.20_dp,0.10_dp]
  call rbzinb_sample(n,1.5_dp,0.8_dp,1.1_dp,0.7_dp,1.0_dp,p4,x,y)
  init9=[1.4_dp,0.9_dp,1.0_dp,0.8_dp,0.9_dp,0.50_dp,0.18_dp,0.20_dp,0.12_dp]
  ll0=loglik_bzinb(x,y,init9)
  fz=fit_bzinb(x,y,maxiter=1200,initial=init9,tol=2e-6_dp)
  call assert_true('bzinb likelihood improves',fz%loglik>=ll0-1e-6_dp,fail)
  call assert_close('bzinb psum',sum(fz%param(6:9)),1.0_dp,1e-12_dp,fail)
  call assert_true('bzinb rho valid',fz%rho>0.0_dp.and.fz%rho<1.0_dp,fail)
  if(fail==0)then;print *,'test_bnb_bzinb: PASS';else;error stop 1;end if
end program
