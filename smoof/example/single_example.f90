program single_example
  use smoof_kinds, only : dp
  use smoof_single, only : ackley, rosenbrock
  implicit none
  print '(a,es16.8)','Ackley(0) = ',ackley([0.0_dp,0.0_dp,0.0_dp])
  print '(a,es16.8)','Rosenbrock(1) = ',rosenbrock([1.0_dp,1.0_dp,1.0_dp])
end program single_example
