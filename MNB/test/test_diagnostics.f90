program test_diagnostics
  use mnb, only : dp,mnb_residual_result,mnb_global_result,mnb_local_result,residuals_mnb,&
    global_influence_mnb,local_influence_mnb,local_weight,local_covariate
  implicit none
  integer,parameter::n=5,mi=2,nn=n*mi
  real(dp)::x(nn,2),y(nn),start(3)
  integer::i
  type(mnb_residual_result)::r
  type(mnb_global_result)::g
  type(mnb_local_result)::l,lc
  x(:,1)=1.0_dp;do i=1,nn;x(i,2)=real(mod((i-1)/mi,2),dp);end do
  y=[1.0_dp,0.0_dp,2.0_dp,1.0_dp,0.0_dp,1.0_dp,3.0_dp,2.0_dp,1.0_dp,1.0_dp]
  start=[2.0_dp,0.0_dp,0.1_dp]
  r=residuals_mnb(start,y,x,n,mi)
  if(size(r%weighted)/=nn .or. size(r%deviance)/=n)error stop 1
  if(any(r%leverage>=1.0_dp))error stop 2
  g=global_influence_mnb(start,y,x,n,mi);if(size(g%cook_distance)/=n)error stop 3
  l=local_influence_mnb(start,y,x,n,mi,local_weight);if(size(l%direction)/=n)error stop 4
  lc=local_influence_mnb(start,y,x,n,mi,local_covariate,covariate=2)
  if(maxval(abs(lc%total_curvature))>1.0e-5_dp)error stop 5
  print *, 'test_diagnostics: PASS'
end program
