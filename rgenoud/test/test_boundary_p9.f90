program test_boundary_p9
  use rgenoud, only : dp, genoud_options, genoud_result, genoud_optimize
  implicit none
  type(genoud_options) :: opt
  type(genoud_result) :: res
  real(dp) :: lower(1), upper(1)

  lower = 0.0_dp
  upper = 1.0_dp
  opt%pop_size = 40
  opt%max_generations = 30
  opt%wait_generations = 6
  opt%boundary_enforcement = 2
  opt%operator_weights = [10.0_dp, 10.0_dp, 10.0_dp, 10.0_dp, 10.0_dp, &
    10.0_dp, 10.0_dp, 10.0_dp, 20.0_dp]
  opt%p9_mix = 1.0_dp
  opt%seed = 90210
  call genoud_optimize(outside_optimum, lower, upper, opt, res)
  if (res%par(1) < lower(1) .or. res%par(1) > upper(1)) then
    error stop "boundary enforcement test failed"
  end if
  if (abs(res%par(1) - 1.0_dp) > 1.0e-5_dp) error stop "bounded optimum test failed"
  if (mod(res%operators(6), 2) /= 0 .or. mod(res%operators(8), 2) /= 0) then
    error stop "crossover operator count parity failed"
  end if
contains
  real(dp) function outside_optimum(x) result(f)
    real(dp), intent(in) :: x(:)
    f = (x(1) - 3.0_dp)**2
  end function outside_optimum
end program test_boundary_p9
