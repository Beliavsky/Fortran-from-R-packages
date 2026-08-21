program test_local_schemes
  use mnb, only : dp,mnb_local_result,local_influence_mnb,local_weight_obs,local_dispersion
  implicit none
  integer,parameter::n=4,mi=2,nn=n*mi
  real(dp)::x(nn,2),y(nn),start(3)
  type(mnb_local_result)::lo,ld
  integer::i
  x(:,1)=1.0_dp;do i=1,nn;x(i,2)=real(mod((i-1)/mi,2),dp);end do
  y=[1.0_dp,0.0_dp,2.0_dp,1.0_dp,0.0_dp,1.0_dp,3.0_dp,2.0_dp];start=[2.0_dp,0.0_dp,0.1_dp]
  lo=local_influence_mnb(start,y,x,n,mi,local_weight_obs)
  if(size(lo%direction)/=nn .or. size(lo%total_curvature)/=nn)error stop 1
  ld=local_influence_mnb(start,y,x,n,mi,local_dispersion)
  if(size(ld%direction)/=n .or. size(ld%total_curvature)/=n)error stop 2
  if(any(abs(ld%total_curvature)>huge(1.0_dp)/10.0_dp))error stop 3
  print *, 'test_local_schemes: PASS'
end program
