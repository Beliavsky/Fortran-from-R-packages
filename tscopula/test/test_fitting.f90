program test_fitting
  use tscopula
  use test_utils
  implicit none
  type(arma_copula)::a
  type(arma_fit_result)::af
  type(margin_spec)::m
  type(margin_fit_result)::mf
  type(dvine_copula)::d
  type(dvine_fit_result)::df
  type(tscopula_spec)::base
  type(vtscopula_spec)::vc
  type(tscm_spec)::tm
  type(tscm_fit_result)::tf
  real(dp),allocatable::u(:),x(:)
  real(dp)::p,q
  call set_seed(112233)
  a=armacopula(ar=[0.5_dp],ma=[real(dp)::]);u=sim(a,350)
  af=fit(a,u,500);call assert_true(af%convergence==0.or.af%convergence==1,'AR fit completed')
  call assert_true(abs(af%model%ar(1)-0.5_dp)<0.25_dp,'AR fit recovery')
  m=margin('gauss',[0.0_dp,1.0_dp]);x=1.5_dp*normal_quantile(u)+0.7_dp
  mf=fit(m,x,500);call assert_true(abs(mf%margin%mu-0.7_dp)<0.2_dp,'margin mean recovery')
  call assert_true(abs(mf%margin%sigma-1.5_dp)<0.2_dp,'margin scale recovery')
  d=dvinecopula([paircop('gauss',0.5_dp)]);u=sim(d,250);df=fit(d,u)
  call assert_true(df%model%pairs(1)%par1>0.1_dp,'D-vine sequential fit')
  base=tscopula_from_arma(armacopula(ar=[0.2_dp],ma=[real(dp)::]));vc=vtscopula(base,vlinear(0.5_dp));tm=tscm(vc,m)
  x=sim(tm,180);tf=fit(tm,x,300);call assert_true(tf%log_likelihood>-huge(1.0_dp)/1000.0_dp,'full-model fit likelihood')
  p=pcondvtarma(a,vlinear(0.4_dp),u,0.3_dp);q=qcondvtarma(a,vlinear(0.4_dp),u,p)
  call assert_close(q,0.3_dp,3.0e-3_dp,'conditional VT-ARMA inversion')
  call pass('fitting and compatibility APIs')
end program test_fitting
