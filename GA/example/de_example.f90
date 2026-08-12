program de_example
  use ga
  implicit none
  type(de_control_type) :: control
  type(ga_real_result) :: result
  real(dp) :: lower(3), upper(3)

  lower = -5.0_dp
  upper =  5.0_dp
  control%pop_size = 40
  control%max_iter = 250
  control%run = 70
  control%stepsize = 0.8_dp
  control%pcrossover = 0.7_dp
  control%seed = 2468

  call de_real(fitness, lower, upper, control, result)
  print '(a,f16.10)', 'best fitness = ', result%fitness_value
  print '(a,*(f12.7,1x))', 'solution     = ', result%solution
contains
  function fitness(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f
    f = -sum((x-[1.0_dp,-2.0_dp,0.5_dp])**2)
  end function fitness
end program de_example
