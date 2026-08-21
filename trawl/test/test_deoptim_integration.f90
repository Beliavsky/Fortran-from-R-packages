program test_deoptim_integration
  use trawl_kinds, only : dp
  use trawl_optimize, only : differential_evolution
  use trawl_rng, only : set_trawl_seed
  use deoptim, only : de_control, de_result, deoptim_solve, de_success, i8
  implicit none

  real(dp) :: lower(2), upper(2), best(2), bestval, best2(2), bestval2
  type(de_control) :: control
  type(de_result) :: direct
  integer :: status, failures
  integer(i8), parameter :: seed = 24681357_i8

  failures = 0
  lower = [-3.0_dp, -2.0_dp]
  upper = [ 3.0_dp,  4.0_dp]

  call differential_evolution(rosenbrock, lower, upper, best, bestval, &
    itermax=250, npop=40, status=status, seed=seed)

  control%itermax = 250
  control%strategy = 2
  control%cr = 0.5_dp
  control%f = 0.8_dp
  control%bs = .false.
  control%trace = 0
  control%np = 40
  control%p = 0.2_dp
  control%c = 0.0_dp
  control%reltol = sqrt(epsilon(1.0_dp))
  control%steptol = 250
  control%storepopfrom = 251
  control%storepopfreq = 1
  control%seed = seed

  call deoptim_solve(rosenbrock, lower, upper, direct, control=control)

  if (status /= 0 .or. direct%status /= de_success) failures = failures + 1
  if (.not. allocated(direct%bestmem)) then
    failures = failures + 1
  else
    if (maxval(abs(best - direct%bestmem)) > 0.0_dp) failures = failures + 1
  end if
  if (abs(bestval - direct%bestval) > 0.0_dp) failures = failures + 1
  if (bestval > 1.0e-8_dp) failures = failures + 1
  if (maxval(abs(best - [1.0_dp, 1.0_dp])) > 5.0e-4_dp) failures = failures + 1

  ! The trawl package has one public RNG seeding entry point.  The DEoptim
  ! wrapper consumes that RNG when no explicit optimizer seed is supplied, so
  ! repeated fits after the same set_trawl_seed() call are reproducible.
  call set_trawl_seed(97531)
  call differential_evolution(rosenbrock, lower, upper, best, bestval, &
    itermax=80, npop=20, status=status)
  call set_trawl_seed(97531)
  call differential_evolution(rosenbrock, lower, upper, best2, bestval2, &
    itermax=80, npop=20, status=status)
  if (maxval(abs(best - best2)) > 0.0_dp) failures = failures + 1
  if (abs(bestval - bestval2) > 0.0_dp) failures = failures + 1

  if (failures == 0) then
    print '(a)', 'test_deoptim_integration: PASS'
  else
    print '(a,i0)', 'test_deoptim_integration: FAIL ', failures
    error stop 1
  end if

contains

  function rosenbrock(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f
    f = 100.0_dp * (x(2) - x(1) * x(1))**2 + (1.0_dp - x(1))**2
  end function rosenbrock

end program test_deoptim_integration
