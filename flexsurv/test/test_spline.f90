program test_spline
  use flexsurv
  implicit none
  type(survspline_model)::m
  real(dp)::kn(3),b(3),db(3),s,h,c
  real(dp)::xmat(4,2),dx(3,2),y(4)
  real(dp),allocatable::coef(:)
  integer::fails,st
  fails=0;kn=[-2.0_dp,0.0_dp,2.0_dp]
  allocate(m%knots(3),m%gamma(3));m%knots=kn;m%gamma=[0.0_dp,1.0_dp,0.0_dp]
  m%scale=spline_scale_hazard;m%timescale=spline_time_log
  call rp_basis(kn,0.5_dp,b);call rp_dbasis(kn,0.5_dp,db)
  if(abs(b(1)-1.0_dp)>1e-14_dp.or.abs(db(2)-1.0_dp)>1e-14_dp)fails=fails+1
  s=survspline_survival(m,2.0_dp);h=survspline_hazard(m,2.0_dp);c=survspline_cdf(m,2.0_dp)
  if(abs(s-exp(-2.0_dp))>1e-10_dp)then;print *,'spline S ',s;fails=fails+1;end if
  if(abs(h-1.0_dp)>1e-10_dp)then;print *,'spline h ',h;fails=fails+1;end if
  if(abs(c-(1.0_dp-exp(-2.0_dp)))>1e-10_dp)fails=fails+1
  ! QP initialization: fit y = 1 + 2x with derivative >= 0.
  xmat=reshape([1.0_dp,1.0_dp,1.0_dp,1.0_dp, 0.0_dp,1.0_dp,2.0_dp,3.0_dp],[4,2])
  y=[1.0_dp,3.0_dp,5.0_dp,7.0_dp]
  dx(:,1)=0.0_dp;dx(:,2)=1.0_dp
  call spline_qp_initial(y,xmat,dx,coef,st)
  if(st/=0.or.abs(coef(1)-1.0_dp)>1e-5_dp.or.abs(coef(2)-2.0_dp)>1e-5_dp)then
    print *,'qp ',st,coef;fails=fails+1
  end if
  if(fails>0)error stop 1
  print *,'test_spline: PASS'
end program test_spline
