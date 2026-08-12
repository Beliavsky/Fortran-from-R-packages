program sphere_example
  use cmaes, only : dp, cma_control, cma_result, cma_es
  use cmaes_functions, only : f_sphere
  implicit none
  type(cma_control) :: control
  type(cma_result) :: result
  real(dp) :: par(5), lower(5), upper(5)

  par = [4.0_dp, -3.0_dp, 2.0_dp, 1.0_dp, -4.0_dp]
  lower = -10.0_dp
  upper = 10.0_dp
  control%seed = 12345
  control%maxit = 1000
  control%stopfitness = 1.0e-12_dp

  result = cma_es(par, f_sphere, lower, upper, control)
  print '(a,es16.8)', 'value = ', result%value
  print '(a,*(f12.6,1x))', 'par   = ', result%par
  print '(a,i0)', 'evaluations = ', result%function_evaluations
end program sphere_example
