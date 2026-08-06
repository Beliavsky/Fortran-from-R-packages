program test_timeseries
  use tscopula
  use test_utils
  implicit none
  type(arma_copula)::a
  type(sarma_copula)::s
  type(arma_filter_result)::f
  real(dp),allocatable::u(:),g(:),pac(:),r(:)
  call set_seed(12345)
  a=armacopula(ar=[0.55_dp],ma=[-0.2_dp]);call assert_true(.not.non_stat(a%ar),'AR stationarity')
  call assert_true(.not.non_invert(a%ma),'MA invertibility')
  g=arma_autocovariance(a,3,sigmastarma(a)**2);call assert_close(g(1),1.0_dp,2.0e-5_dp,'standardized variance')
  u=sim_arma_copula(a,500);call assert_true(all(u>0.0_dp).and.all(u<1.0_dp),'ARMA uniform simulation')
  f=kfilter(a,u);call assert_true(f%objective<huge(1.0_dp)/1000.0_dp,'finite ARMA likelihood')
  r=resid_arma_copula(a,u);call assert_true(all(r>0.0_dp).and.all(r<1.0_dp),'ARMA residuals')
  pac=kpacf_arma([0.5_dp],[real(dp)::],5);call assert_close(pac(1),2.0_dp*asin(0.5_dp)/pi,3.0e-3_dp,'AR(1) Kendall PACF')
  s=sarmacopula(ar=[0.2_dp],sar=[0.3_dp],period=4);a=sarma2arma(s);call assert_true(size(a%ar)==5,'seasonal AR expansion')
  call assert_close(predict_arma_cdf(a,u,0.5_dp),min(max(predict_arma_cdf(a,u,0.5_dp),0.0_dp),1.0_dp),1.0e-14_dp,'forecast CDF bounds')
  call pass('ARMA and SARMA copulas')
end program test_timeseries
