program permutation_example
  use ga
  implicit none
  type(ga_control_type) :: control
  type(ga_int_result) :: result
  integer, parameter :: target(8) = [5,2,8,1,7,3,6,4]

  control%pop_size = 100
  control%max_iter = 300
  control%run = 90
  control%pmutation = 0.25_dp
  control%seed = 13579

  call ga_permutation(fitness, 1, 8, control, result)
  print '(a,f10.3)', 'best fitness = ', result%fitness_value
  print '(a,*(i0,1x))', 'solution     = ', result%solution
contains
  function fitness(x) result(f)
    integer, intent(in) :: x(:)
    real(dp) :: f
    f = real(count(x==target),dp)
  end function fitness
end program permutation_example
