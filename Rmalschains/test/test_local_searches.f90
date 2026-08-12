program test_local_searches
  use iso_fortran_env, only : int64
  use rmalschains, only : dp, mals_control, mals_result, malschains_control, malschains_optimize
  implicit none
  character(len=8), parameter :: methods(6) = [character(len=8) :: 'sw', 'ssw', 'simplex', 'mts1', 'mts2', 'cmaes']
  type(mals_control) :: ctrl
  type(mals_result) :: res
  real(dp) :: lo(4), hi(4), init(4, 1)
  integer :: k
  lo = -4.0_dp; hi = 4.0_dp; init(:, 1) = 2.5_dp
  do k = 1, size(methods)
    ctrl = malschains_control(popsize=20, ls=trim(methods(k)), istep=100, effort=0.5_dp, &
      threshold=0.0_dp, seed=int(100+k, int64))
    ctrl%ls_only = .true.
    ctrl%legacy_ls_only_zero_start = .false.
    res = malschains_optimize(shifted_sphere, lo, hi, 3500, ctrl, init)
    if (res%fitness > 5.0e-3_dp) then
      print *, 'method=', trim(methods(k)), ' fitness=', res%fitness
      error stop 'test_local_searches: local search did not converge'
    end if
  end do
  print *, 'PASS test_local_searches'
contains
  function shifted_sphere(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f
    f = sum((x - 0.75_dp)**2)
  end function shifted_sphere
end program test_local_searches
