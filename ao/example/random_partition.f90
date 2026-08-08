program random_partition_example
  use ao
  implicit none
  type(ao_result) :: result
  type(ao_options) :: options
  real(dp) :: initial(4)
  initial = 0.0_dp
  options%partition = AO_PARTITION_RANDOM
  options%new_block_probability = 0.35_dp
  options%minimum_block_number = 2
  options%iteration_limit = 30
  call ao_seed(2026)
  call ao_optimize(objective, initial, result, options, gradient=gradient)
  print '(a,4f12.6)', 'estimate: ', result%estimate
  print '(a,es14.6)', 'value:    ', result%value
contains
  function objective(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value
    value = sum((x-1.0_dp)**2)
  end function objective
  subroutine gradient(x,g)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    g = 2.0_dp*(x-1.0_dp)
  end subroutine gradient
end program random_partition_example
