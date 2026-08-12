program test_initial_population
  use iso_fortran_env, only : int64
  use rmalschains, only : dp, mals_control, mals_result, malschains_control, malschains_optimize
  implicit none
  type(mals_control) :: ctrl
  type(mals_result) :: res
  real(dp) :: lo(3), hi(3), initial(3, 2)
  lo = -2.0_dp; hi = 2.0_dp
  initial(:, 1) = 0.0_dp
  initial(:, 2) = 1.0_dp
  ctrl = malschains_control(popsize=20, ls='sw', istep=50, effort=0.5_dp, alpha=0.5_dp, &
    optimum=0.0_dp, threshold=1.0e-12_dp, seed=int(9, int64))
  res = malschains_optimize(sphere, lo, hi, 500, ctrl, initial)
  if (res%fitness > 1.0e-14_dp) error stop 'test_initial_population: supplied optimum was not retained'
  if (res%actual_nfe /= res%num_eval_ea + res%num_eval_ls + 2) &
    error stop 'test_initial_population: initial-population accounting mismatch'
  print *, 'PASS test_initial_population'
contains
  function sphere(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f
    f = sum(x*x)
  end function sphere
end program test_initial_population
