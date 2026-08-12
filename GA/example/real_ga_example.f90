program real_ga_example
  use ga
  implicit none
  type(ga_control_type) :: control
  type(ga_real_result) :: result
  real(dp) :: lower(2), upper(2)

  lower = -5.12_dp
  upper =  5.12_dp
  control%pop_size = 80
  control%max_iter = 300
  control%run = 80
  control%seed = 12345

  call ga_real(rastrigin_fitness, lower, upper, control, result)

  print '(a,f16.8)', 'best fitness = ', result%fitness_value
  print '(a,*(f12.6,1x))', 'solution     = ', result%solution
  print '(a,i0)', 'iterations   = ', result%iter
  print '(a,i0)', 'evaluations  = ', result%evaluations
contains
  function rastrigin_fitness(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f
    real(dp), parameter :: pi = acos(-1.0_dp)
    f = -(20.0_dp + sum(x*x - 10.0_dp*cos(2.0_dp*pi*x)))
  end function rastrigin_fitness
end program real_ga_example
