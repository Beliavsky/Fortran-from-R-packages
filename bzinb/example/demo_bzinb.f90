program demo_bzinb
  use bzinb
  implicit none
  integer,parameter::n=120
  integer::x(n),y(n)
  real(dp)::p(4)
  type(bzinb_fit_result)::fit
  p=[0.55_dp,0.15_dp,0.20_dp,0.10_dp]
  call set_bzinb_seed(42)
  call rbzinb_sample(n,1.5_dp,0.8_dp,1.1_dp,0.7_dp,1.0_dp,p,x,y)
  fit=fit_bzinb(x,y,maxiter=2500,tol=1.0e-5_dp,initial=[1.4_dp,0.9_dp,1.0_dp,0.8_dp,0.9_dp,0.5_dp,0.18_dp,0.2_dp,0.12_dp])
  print '(a,l1)','converged: ',fit%converged
  print '(a,9(1x,f8.4))','params:',fit%param
  print '(a,f10.5)','rho: ',fit%rho
  print '(a,f12.4)','logLik: ',fit%loglik
end program
