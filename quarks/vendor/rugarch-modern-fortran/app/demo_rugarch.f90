program demo_rugarch
   use rugarch
   implicit none

   integer, parameter :: n = 1200
   type(garch_spec) :: spec
   type(garch_fit_result) :: fit
   type(arfima_spec) :: aspec
   type(arfima_fit_result) :: afit
   type(var_test_result) :: vtest
   real(dp), allocatable :: y(:), sigma(:), eps(:), forecast(:)
   real(dp), allocatable :: ax(:), var_series(:)
   integer :: i

   call seed_rng(20260722)

   spec = make_garch_spec(1,1,model_gjrgarch,dist_std)
   spec%mean = 0.0002_dp
   spec%omega = 2.0e-6_dp
   spec%alpha = 0.06_dp
   spec%gamma = 0.08_dp
   spec%beta = 0.88_dp
   spec%shape = 7.0_dp
   allocate(y(n),sigma(n),eps(n))
   call simulate_garch(spec,n,y,sigma,eps)

   fit = fit_gjrgarch11(y,dist_std,fit_mean=.true.,fit_shape=.true.,max_iterations=1200)
   allocate(forecast(5))
   call forecast_volatility(fit%spec,fit%residuals,fit%sigma,5,forecast)

   print '(a,i0)', 'observations: ', n
   print '(a,f10.6)', 'true omega: ', spec%omega
   print '(a,f10.6)', 'fit omega:  ', fit%spec%omega
   print '(a,f10.6)', 'fit alpha:  ', fit%spec%alpha(1)
   print '(a,f10.6)', 'fit gamma:  ', fit%spec%gamma(1)
   print '(a,f10.6)', 'fit beta:   ', fit%spec%beta(1)
   print '(a,f10.4)', 'fit shape:  ', fit%spec%shape
   print '(a,f12.3)', 'log likelihood: ', fit%log_likelihood
   print '(a,5(1x,f9.6))', 'sigma forecast:', forecast

   aspec = make_arfima_spec(1,0)
   aspec%d = 0.25_dp
   aspec%ar(1) = 0.35_dp
   aspec%innovation_sd = 0.8_dp
   allocate(ax(1500))
   call simulate_arfima(aspec,size(ax),ax)
   afit = fit_arfima(ax,p=1,q=0)
   print '(a,f9.4)', 'true ARFIMA d: ', aspec%d
   print '(a,f9.4)', 'GPH ARFIMA d:  ', afit%spec%d
   print '(a,f9.4)', 'fitted AR(1):  ', afit%spec%ar(1)

   allocate(var_series(n))
   do i=1,n
      var_series(i)=value_at_risk(0.05_dp,fit%spec%mean,fit%sigma(i), &
         fit%spec%cond_dist,fit%spec%shape,fit%spec%skew,fit%spec%lambda)
   end do
   vtest=var_test(0.05_dp,y,var_series)
   print '(a,i0)', 'VaR exceedances: ', vtest%actual_exceedances
   print '(a,f9.5)', 'VaR UC p-value:   ', vtest%uc_p_value
   print '(a,f10.6)', 'standardized NIG density at zero: ', &
      dsnig(0.0_dp,0.0_dp,1.0_dp,0.2_dp,1.0_dp)
   print '(a,f10.6)', 'standardized GH skew-t density at zero: ', &
      dsghst(0.0_dp,0.0_dp,1.0_dp,0.5_dp,8.0_dp)
   spec=make_garch_spec(1,1,model_fgarch,dist_norm)
   call configure_fgarch_submodel(spec,fgarch_allgarch)
   print '(a,a)', 'available fGARCH example submodel: ', trim(fgarch_submodel_name(spec%fgarch_submodel))
end program demo_rugarch
