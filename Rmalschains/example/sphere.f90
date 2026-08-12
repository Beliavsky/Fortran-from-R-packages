program example_sphere
  use iso_fortran_env, only : int64
  use rmalschains, only : dp, mals_control, mals_result, malschains_control, malschains_optimize
  implicit none
  type(mals_control) :: ctrl
  type(mals_result) :: res
  real(dp) :: lo(10), hi(10)
  lo = -5.0_dp; hi = 5.0_dp
  ctrl = malschains_control(popsize=20, ls='cmaes', istep=150, effort=0.5_dp, alpha=0.5_dp, seed=int(1234, int64))
  res = malschains_optimize(sphere, lo, hi, 8000, ctrl)
  print '(a,es14.6)', 'fitness: ', res%fitness
  print '(a,*(f10.5,1x))', 'solution: ', res%sol
  print '(a,i0,a,i0)', 'EA evaluations: ', res%num_eval_ea, ', LS evaluations: ', res%num_eval_ls
contains
  function sphere(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f
    f = sum(x*x)
  end function sphere
end program example_sphere
