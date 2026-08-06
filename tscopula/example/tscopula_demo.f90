program tscopula_demo
  use tscopula
  implicit none
  type(arma_copula) :: arma
  type(dvine2_copula) :: vine2
  type(tscopula_spec) :: core
  type(vtscopula_spec) :: vcop
  type(tscm_spec) :: model
  type(margin_spec) :: marg
  real(dp), allocatable :: x(:), u(:), tau(:)
  real(dp) :: q01, cdf0

  call set_seed(20260804)
  arma = armacopula(ar=[0.55_dp], ma=[-0.20_dp])
  vine2 = arma2dvine('gauss', arma%ar, arma%ma, 4)
  core = tscopula_from_dvine(vine2%vine)
  vcop = vtscopula(core, v2p(delta=0.42_dp,kappa=1.25_dp))
  marg = margin('sst',[7.0_dp,1.20_dp,0.0_dp,1.0_dp])
  model = tscm(vcop,marg)

  x = sim(model,500)
  u = pmarg(marg,x)
  tau = kendall(core,4)
  q01 = predict_tscm_quantile(model,x(1:499),0.01_dp)
  cdf0 = predict_tscm_cdf(model,x(1:499),0.0_dp)

  write(*,'(a,f10.5)') 'sample mean: ',sum(x)/real(size(x),dp)
  write(*,'(a,4f10.5)') 'lag Kendall values: ',tau
  write(*,'(a,f10.5)') 'one-step 1% quantile: ',q01
  write(*,'(a,f10.5)') 'one-step CDF at zero: ',cdf0
  write(*,'(a,f10.5)') 'copula log likelihood: ',tscopula_loglik(core,vtrans(vcop%transform,u))
end program tscopula_demo
