program test_ahres
  use calibrar, only : dp, ahr_options, ahr_result, ahres, two_component_objective
  implicit none
  type(ahr_options) :: op
  type(ahr_result) :: res
  real(dp) :: x0(3),lo(3),hi(3),w(2)
  x0=[2.0_dp,-1.0_dp,1.5_dp];lo=-3.0_dp;hi=3.0_dp;w=[1.0_dp,1.0_dp]
  op%maxgen=220;op%popsize=18;op%seed=8123;op%selection=0.5_dp;op%alpha=0.1_dp;op%step=0.5_dp
  op%termination=0
  call ahres(x0,two_component_objective,2,res,lo,hi,op,w)
  if(res%value>0.8_dp) then
    print *, res%par,res%partial,res%value
    error stop "AHR-ES failed"
  end if
  if(any(res%par<lo) .or. any(res%par>hi)) error stop "AHR-ES bounds failed"
  print *, "PASS test_ahres"
end program test_ahres
