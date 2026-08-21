program test_splines2ns
  use flexsurv_kinds, only : dp
  use flexsurv_splines2ns, only : splines2ns_basis, splines2ns_dbasis
  use flexsurv_spline, only : survspline_model,survspline_eta,survspline_deta_dt, &
    spline_basis_splines2ns,spline_time_identity
  implicit none
  real(dp)::k2(2),b2(2),d2(2),b3(3),d3(3),bm(3),bp(3),h
  real(dp)::k3(3)
  type(survspline_model)::m
  k2=[0.0_dp,1.0_dp]
  call splines2ns_basis(k2,0.3_dp,b2);call splines2ns_dbasis(k2,0.3_dp,d2)
  if(maxval(abs(b2-[1.0_dp,0.15_dp]))>2.0e-12_dp)error stop 'ns two-knot basis'
  if(maxval(abs(d2-[0.0_dp,0.5_dp]))>2.0e-11_dp)error stop 'ns two-knot derivative'
  call splines2ns_basis(k2,-1.0_dp,b2)
  if(abs(b2(2)+0.5_dp)>2.0e-10_dp)error stop 'ns left extrapolation'
  call splines2ns_basis(k2,2.0_dp,b2)
  if(abs(b2(2)-1.0_dp)>2.0e-10_dp)error stop 'ns right extrapolation'
  k3=[0.0_dp,0.4_dp,1.0_dp];h=1.0e-6_dp
  call splines2ns_basis(k3,0.37_dp-h,bm);call splines2ns_basis(k3,0.37_dp+h,bp)
  call splines2ns_basis(k3,0.37_dp,b3);call splines2ns_dbasis(k3,0.37_dp,d3)
  if(maxval(abs((bp-bm)/(2.0_dp*h)-d3))>2.0e-6_dp)error stop 'ns derivative finite diff'
  call splines2ns_dbasis(k3,-0.5_dp,bm);call splines2ns_dbasis(k3,-1.0_dp,bp)
  if(maxval(abs(bm-bp))>2.0e-10_dp)error stop 'ns natural left derivative'
  m%knots=k2;m%gamma=[0.2_dp,2.0_dp];m%basis=spline_basis_splines2ns
  m%timescale=spline_time_identity
  if(abs(survspline_eta(m,0.3_dp)-0.5_dp)>2.0e-11_dp)error stop 'ns model eta'
  if(abs(survspline_deta_dt(m,0.3_dp)-1.0_dp)>2.0e-10_dp)error stop 'ns model derivative'
  print *,'test_splines2ns: PASS'
end program
