program test_maximize_bounds
  use cmaes, only : dp, cma_control, cma_result, cma_es
  implicit none
  type(cma_control) :: ctrl
  type(cma_result) :: res
  real(dp) :: par(2), lo(2), hi(2)

  par = [5.0_dp, 5.0_dp]
  lo = -10.0_dp
  hi = 10.0_dp
  ctrl%seed = 77
  ctrl%maxit = 500
  ctrl%fnscale = -1.0_dp
  ctrl%stopfitness = -1.0e-8_dp
  res = cma_es(par, negsphere, lo, hi, ctrl)
  if (res%value <= -1.0e-8_dp) error stop "maximization did not approach zero"

  par = [2.0_dp, 2.0_dp]
  lo = [0.5_dp, -10.0_dp]
  hi = 10.0_dp
  ctrl = cma_control()
  ctrl%seed = 1234
  ctrl%maxit = 800
  res = cma_es(par, sphere, lo, hi, ctrl)
  if (any(res%par < lo - 1.0e-6_dp) .or. any(res%par > hi + 1.0e-6_dp)) error stop "best point violates bounds"
  if (abs(res%par(1) - 0.5_dp) > 3.0e-3_dp) error stop "boundary optimum not reached"
  print '(a,2f12.6)', 'bounded par: ', res%par
contains
  function negsphere(x) result(v)
    real(dp), intent(in) :: x(:)
    real(dp) :: v
    v = -dot_product(x, x)
  end function negsphere

  function sphere(x) result(v)
    real(dp), intent(in) :: x(:)
    real(dp) :: v
    v = dot_product(x, x)
  end function sphere
end program test_maximize_bounds
