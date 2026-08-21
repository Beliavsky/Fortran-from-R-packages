program test_poisson_fits
  use bzinb
  use test_support
  implicit none
  integer,parameter::n=4000
  integer::x(n),y(n),fail
  type(bp_fit_result)::fbp
  type(bzip_a_fit_result)::fa
  type(bzip_b_fit_result)::fb
  real(dp)::p(4)
  fail=0;call set_bzinb_seed(77)
  call rbp_sample(n,0.8_dp,1.1_dp,1.5_dp,x,y);fbp=fit_bp(x,y)
  call assert_close('bp m0',fbp%param(1),0.8_dp,0.10_dp,fail)
  call assert_close('bp m1',fbp%param(2),1.1_dp,0.10_dp,fail)
  call assert_close('bp m2',fbp%param(3),1.5_dp,0.10_dp,fail)
  call rbzip_a_sample(n,0.6_dp,1.2_dp,0.9_dp,0.25_dp,x,y);fa=fit_bzip_a(x,y,maxiter=2000)
  call assert_close('bzipA p',fa%param(4),0.25_dp,0.05_dp,fail)
  p=[0.55_dp,0.20_dp,0.15_dp,0.10_dp]
  call rbzip_b_sample(n,0.5_dp,1.0_dp,0.8_dp,p,x,y);fb=fit_bzip_b(x,y,maxiter=3000)
  call assert_close('bzipB p1',fb%param(4),p(1),0.07_dp,fail)
  call assert_close('bzipB p4',fb%param(7),p(4),0.07_dp,fail)
  if(fail==0)then;print *,'test_poisson_fits: PASS';else;error stop 1;end if
end program
