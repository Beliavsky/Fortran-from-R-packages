program test_cmaes
  use cmaes, only : dp, cma_control, cma_result, cma_es
  use cmaes_functions, only : f_sphere
  implicit none
  type(cma_control) :: ctrl
  type(cma_result) :: res
  real(dp) :: par(2), lo(2), hi(2)

  par = [4.0_dp, -3.0_dp]
  lo = -10.0_dp
  hi = 10.0_dp
  ctrl%seed = 20260810
  ctrl%maxit = 400
  ctrl%stopfitness = 1.0e-10_dp
  res = cma_es(par, f_sphere, lo, hi, ctrl)
  if (res%value > 1.0e-10_dp) error stop "sphere did not reach stopfitness"
  if (res%convergence /= 0) error stop "sphere convergence code"
  if (res%function_evaluations /= res%iterations * (4 + floor(3.0_dp * log(2.0_dp)))) &
    error stop "evaluation count mismatch"
  print '(a,es12.4)', 'sphere value: ', res%value
end program test_cmaes
