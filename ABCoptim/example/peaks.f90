program peaks_example
  use abcoptim, only : dp, abc_control, abc_result, abc_cpp
  implicit none

  type(abc_control) :: control
  type(abc_result) :: result
  real(dp) :: par(2), lb(1), ub(1)
  real(dp), parameter :: pi = acos(-1.0_dp)

  par = 0.0_dp
  lb = -10.0_dp
  ub = 10.0_dp
  control%seed = 213
  control%criter = 80

  call abc_cpp(par, objective, lb, ub, result, control)

  print '(a,es16.8)', 'value = ', result%value
  print '(a,2f16.8)', 'par   = ', result%par
  print '(a,i0)', 'cycles = ', result%counts
  print '(a,i0)', 'evaluations = ', result%objective_evaluations

contains

  function objective(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value

    value = -cos(x(1)) * cos(x(2)) * &
      exp(-((x(1) - pi)**2 + (x(2) - pi)**2))
  end function objective

end program peaks_example
