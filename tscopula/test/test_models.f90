program test_models
  use tscopula
  use test_utils
  implicit none
  type(arma_copula)::a
  type(tscopula_spec)::base
  type(vtscopula_spec)::vc,vcw
  type(tscm_spec)::model
  type(margin_spec)::m
  type(empirical_distribution)::e
  real(dp),allocatable::x(:),u(:),prof(:)
  real(dp)::p,q
  call set_seed(54321);a=armacopula(ar=[0.45_dp],ma=[real(dp)::]);base=tscopula_from_arma(a)
  vc=vtscopula(base,vlinear(0.35_dp));vcw=vtscopula(base,vlinear(0.35_dp),paircop('frank',2.0_dp));m=margin('gauss',[0.5_dp,1.2_dp]);model=tscm(vc,m)
  x=sim_tscm(model,220);call assert_true(size(x)==220,'full-model simulation size')
  p=predict_tscm_cdf(model,x(1:200),0.3_dp);call assert_true(p>=0.0_dp.and.p<=1.0_dp,'full-model CDF')
  q=predict_tscm_quantile(model,x(1:200),p);call assert_close(predict_tscm_cdf(model,x(1:200),q),p,2.0e-3_dp,'full-model forecast inversion')
  u=pmarg(m,x);call assert_true(tscopula_loglik(base,vtrans(vc%transform,u))<huge(1.0_dp),'full-model finite likelihood')
  prof=profilefulcrum(vc,u,[0.25_dp,0.35_dp,0.45_dp]);call assert_true(size(prof)==3,'fulcrum profile')
  call assert_true(wobjective(vcw%wcopula,vcw%transform,u)>=0.0_dp,'W-copula objective')
  e=fit_edf(x);call assert_true(pedf(e,huge(1.0_dp))>0.999_dp,'empirical CDF upper tail')
  call assert_true(predict_empirical(e,0.0_dp,'density')>=0.0_dp,'empirical KDE density')
  call assert_close(aicc(-100.0_dp,4,200),208.0_dp+40.0_dp/195.0_dp,1.0e-12_dp,'AICc')
  call pass('V-copula and full models')
end program test_models
