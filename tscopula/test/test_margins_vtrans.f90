program test_margins_vtrans
  use tscopula
  use test_utils
  implicit none
  type(margin_spec)::m
  type(vtransform_spec)::vt
  real(dp)::p,x,v,u0
  m=margin('gauss',[1.0_dp,2.0_dp]);p=0.13_dp;x=qmarg(m,p)
  call assert_close(pmarg(m,x),p,2.0e-7_dp,'Gaussian margin inverse')
  call assert_close(dmarg(m,1.0_dp),1.0_dp/(2.0_dp*sqrt(2.0_dp*pi)),1.0e-12_dp,'Gaussian density')
  m=margin('laplace',[0.0_dp,1.5_dp]);x=qmarg(m,0.8_dp);call assert_close(pmarg(m,x),0.8_dp,1.0e-12_dp,'Laplace inverse')
  m=margin('st',[7.0_dp,0.0_dp,1.0_dp]);x=qmarg(m,0.2_dp);call assert_close(pmarg(m,x),0.2_dp,2.0e-7_dp,'Student inverse')
  vt=vlinear(0.3_dp);u0=0.18_dp;v=vtrans(vt,u0);call assert_close(vinverse(vt,v),u0,1.0e-10_dp,'linear lower inverse')
  call assert_close(vdownprob(vt,0.4_dp),0.3_dp,1.0e-12_dp,'linear down probability')
  vt=vsymmetric();call assert_close(vtrans(vt,0.2_dp),0.6_dp,1.0e-12_dp,'symmetric transform')
  call assert_close(pcoincide(vt),0.5_dp,1.0e-12_dp,'symmetric coincidence')
  vt=v2p(0.4_dp,1.3_dp);v=vtrans(vt,0.2_dp);call assert_close(vtrans(vt,vinverse(vt,v)),v,2.0e-7_dp,'V2p inverse')
  call pass('margins and V-transforms')
end program test_margins_vtrans
