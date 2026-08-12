program example_rastrigin
  use iso_fortran_env, only : int64
  use rmalschains, only : dp, mals_control, mals_result, malschains_control, malschains_optimize
  implicit none
  type(mals_control) :: ctrl
  type(mals_result) :: res
  real(dp) :: lo(30), hi(30)
  lo = -1.0_dp; hi = 1.0_dp
  ctrl = malschains_control(popsize=20, ls='simplex', istep=100, effort=0.8_dp, alpha=0.3_dp, seed=int(2026, int64))
  res = malschains_optimize(rastrigin, lo, hi, 30000, ctrl)
  print '(a,es14.6)', 'Rastrigin fitness: ', res%fitness
  print '(a,i0)', 'actual objective evaluations: ', res%actual_nfe
contains
  function rastrigin(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f
    real(dp), parameter :: twopi = 6.2831853071795864769_dp
    f = sum(x*x - 10.0_dp*cos(twopi*x) + 10.0_dp)
  end function rastrigin
end program example_rastrigin
