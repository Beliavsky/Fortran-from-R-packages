program test_compatibility
  use iso_fortran_env, only : int64
  use rmalschains, only : dp, mals_control, mals_result, malschains_control, malschains_optimize
  implicit none
  type(mals_control) :: legacy, fixed
  type(mals_result) :: a, b
  real(dp) :: lo(2), hi(2)
  lo = -3.0_dp; hi = 3.0_dp
  legacy = malschains_control(popsize=23, ls='sw', istep=100, seed=int(19, int64))
  if (legacy%popsize /= 20) error stop 'test_compatibility: population rounding mismatch'
  legacy%ls_only = .true.
  a = malschains_optimize(shifted, lo, hi, 500, legacy)
  if (abs(a%fitness) > 1.0e-15_dp .or. maxval(abs(a%sol)) > 1.0e-15_dp) &
    error stop 'test_compatibility: legacy zero-start behavior not reproduced'
  fixed = legacy
  fixed%legacy_ls_only_zero_start = .false.
  b = malschains_optimize(shifted, lo, hi, 1500, fixed)
  if (b%fitness > 1.0e-3_dp) error stop 'test_compatibility: corrected LS-only mode failed'
  if (b%actual_nfe <= b%num_eval_ls) error stop 'test_compatibility: setup evaluations should be visible in actual_nfe'
  print *, 'PASS test_compatibility'
contains
  function shifted(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f
    f = sum((x - 1.0_dp)**2)
  end function shifted
end program test_compatibility
