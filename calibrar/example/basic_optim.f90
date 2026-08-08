program basic_optim
  use calibrar, only : dp, optim_options, optim_result, optim2, rosenbrock_objective
  implicit none
  type(optim_options) :: op
  type(optim_result) :: res
  real(dp) :: x0(2),lo(2),hi(2)
  x0=[-1.2_dp,1.0_dp];lo=-5.0_dp;hi=5.0_dp
  op%maxit=1000;op%maxfeval=30000;op%gradient_method="central"
  call optim2(x0,rosenbrock_objective,res,"BFGS",lo,hi,op)
  print '(a,2f12.6)', 'par = ',res%par
  print '(a,f14.8)', 'value = ',res%value
end program basic_optim
