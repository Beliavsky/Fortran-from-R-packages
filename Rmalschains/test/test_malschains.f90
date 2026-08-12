program test_malschains
  use iso_fortran_env, only : int64
  use rmalschains, only : dp, mals_control, mals_result, malschains_control, malschains_optimize
  implicit none
  type(mals_control) :: ctrl
  type(mals_result) :: res
  real(dp) :: lo(5), hi(5)
  lo = -5.0_dp; hi = 5.0_dp
  ctrl = malschains_control(popsize=20, ls='cmaes', istep=100, effort=0.5_dp, alpha=0.5_dp, seed=int(12345, int64))
  res = malschains_optimize(sphere, lo, hi, 5000, ctrl)
  if (res%fitness > 1.0e-5_dp) error stop 'test_malschains: CMA-ES chain did not converge'
  if (res%num_eval_ea <= 0 .or. res%num_eval_ls <= 0) error stop 'test_malschains: missing evaluation accounting'
  if (res%actual_nfe < res%num_eval_ea + res%num_eval_ls) error stop 'test_malschains: actual count too small'
  print '(a,es12.4,a,i0,a,i0)', 'PASS test_malschains fitness=', res%fitness, &
    ' EA=', res%num_eval_ea, ' LS=', res%num_eval_ls
contains
  function sphere(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f
    f = sum(x*x)
  end function sphere
end program test_malschains
