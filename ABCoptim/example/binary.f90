program binary_example
  use abcoptim, only : dp, abc_control, abc_result, abc_optim
  implicit none

  type(abc_control) :: control
  type(abc_result) :: result
  real(dp) :: par(4), lb(1), ub(1)
  real(dp), parameter :: target(4) = [1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp]

  par = 0.0_dp
  lb = 0.0_dp
  ub = 1.0_dp
  control%optiinteger = .true.
  control%seed = 1234
  control%criter = 50

  call abc_optim(par, objective, lb, ub, result, control)

  print '(a,es16.8)', 'value = ', result%value
  print '(a,4f8.2)', 'par   = ', result%par

contains

  function objective(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value

    value = sum(abs(x - target))
  end function objective

end program binary_example
