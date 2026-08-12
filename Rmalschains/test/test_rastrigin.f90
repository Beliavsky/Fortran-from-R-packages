program test_rastrigin
  use iso_fortran_env, only : int64
  use rmalschains, only : dp, mals_control, mals_result, malschains_control, malschains_optimize
  implicit none
  type(mals_control) :: ctrl
  type(mals_result) :: res
  real(dp) :: lo(10), hi(10)
  lo = -1.0_dp; hi = 1.0_dp
  ctrl = malschains_control(popsize=20, ls='simplex', istep=100, effort=0.8_dp, alpha=0.3_dp, seed=int(777, int64))
  res = malschains_optimize(rastrigin, lo, hi, 12000, ctrl)
  if (res%fitness > 1.0_dp) then
    print *, res%fitness
    error stop 'test_rastrigin: insufficient convergence'
  end if
  print '(a,es12.4)', 'PASS test_rastrigin fitness=', res%fitness
contains
  function rastrigin(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f
    real(dp), parameter :: twopi = 6.2831853071795864769_dp
    f = sum(x*x - 10.0_dp*cos(twopi*x) + 10.0_dp)
  end function rastrigin
end program test_rastrigin
