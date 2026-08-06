program test_backtest
  use ghyp_kinds, only : dp, i8
  use tsgarch
  use test_support
  implicit none
  real(dp)::dummy(30)
  type(garch_spec)::spec
  type(garch_parameters)::par
  type(garch_simulation)::sim
  type(backtest_result)::bt
  type(fit_options)::opt
  integer::i
  dummy=[(0.01_dp*sin(real(i,dp)),i=1,30)]
  spec=standard_spec('garch','norm')
  par=standard_parameters(dummy,spec)
  sim=simulate_garch(spec,par,105,burn=20,seed=991_i8)
  opt%max_iterations=350
  opt%tolerance=1.0e-5_dp
  opt%compute_inference=.false.
  bt=backtest_var(sim%series(:,1),spec,65,probability=0.05_dp,refit_every=20,window=60,options=opt)
  call assert_true(bt%status==tsg_success,'backtest failed: '//trim(bt%message))
  call assert_true(size(bt%actual)==40.and.bt%refits==2,'backtest dimensions/refits')
  call assert_true(bt%coverage>=0.0_dp.and.bt%coverage<=1.0_dp,'coverage range')
  call assert_true(bt%kupiec_pvalue>=0.0_dp.and.bt%kupiec_pvalue<=1.0_dp,'Kupiec p-value range')
  write(*,'(a,f8.4)')'test_backtest: PASS coverage ',bt%coverage
end program test_backtest
