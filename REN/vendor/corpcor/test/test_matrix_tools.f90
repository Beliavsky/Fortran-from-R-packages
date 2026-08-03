program test_matrix_tools
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use corpcor, only : dp, covariance_decomposition, precision_decomposition, sm2vec, &
    sm_index, vec2sm, rebuild_cov, decompose_cov, rebuild_invcov, decompose_invcov
  implicit none
  real(dp), parameter :: tol = 1.0e-12_dp
  real(dp) :: a(3, 3), r(3, 3), v(3), pr(3, 3), pv(3)
  real(dp), allocatable :: vec(:), back(:, :), cov(:, :), precision(:, :)
  integer, allocatable :: idx(:, :)
  type(covariance_decomposition) :: dc
  type(precision_decomposition) :: di

  a = reshape([1.0_dp, 2.0_dp, 3.0_dp, 2.0_dp, 4.0_dp, 5.0_dp, 3.0_dp, 5.0_dp, 6.0_dp], [3, 3])
  vec = sm2vec(a)
  if (maxval(abs(vec - [2.0_dp, 3.0_dp, 5.0_dp])) > tol) error stop 'sm2vec ordering failed'
  back = vec2sm(vec)
  if (.not. ieee_is_nan(back(1, 1))) error stop 'vec2sm diagonal should be NaN'
  if (maxval(abs([back(2,1), back(3,1), back(3,2)] - vec)) > tol) error stop 'vec2sm failed'
  idx = sm_index(3)
  if (any(idx(:, 1) /= [1, 1, 2]) .or. any(idx(:, 2) /= [2, 3, 3])) error stop 'sm_index failed'

  r = reshape([1.0_dp, 0.2_dp, -0.1_dp, 0.2_dp, 1.0_dp, 0.3_dp, -0.1_dp, 0.3_dp, 1.0_dp], [3, 3])
  v = [4.0_dp, 9.0_dp, 16.0_dp]
  cov = rebuild_cov(r, v)
  dc = decompose_cov(cov)
  if (maxval(abs(dc%correlation - r)) > tol .or. maxval(abs(dc%variance-v)) > tol) &
    error stop 'covariance decomposition failed'

  pr = r
  pv = [2.0_dp, 3.0_dp, 5.0_dp]
  precision = rebuild_invcov(pr, pv)
  di = decompose_invcov(precision)
  if (maxval(abs(di%partial_correlation-pr)) > tol .or. maxval(abs(di%partial_variance-pv)) > tol) &
    error stop 'precision decomposition failed'
  print '(a)', 'test_matrix_tools: PASS'
end program test_matrix_tools
