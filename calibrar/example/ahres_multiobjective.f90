program ahres_multiobjective
  use calibrar, only : dp, ahr_options, ahr_result, ahres, two_component_objective
  implicit none
  type(ahr_options) :: op
  type(ahr_result) :: res
  real(dp) :: x0(2),lo(2),hi(2)
  x0=[1.5_dp,-1.0_dp];lo=-3.0_dp;hi=3.0_dp
  op%maxgen=150;op%popsize=16;op%seed=44;op%termination=0
  call ahres(x0,two_component_objective,2,res,lo,hi,op,[1.0_dp,1.0_dp])
  print '(a,2f12.6)', 'par = ',res%par
  print '(a,2f12.6)', 'partial = ',res%partial
  print '(a,f14.8)', 'weighted value = ',res%value
end program ahres_multiobjective
