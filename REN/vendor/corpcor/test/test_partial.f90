program test_partial
  use corpcor, only : dp, matrix_shrinkage_result, vector_shrinkage_result, &
    cor2pcor, pcor2cor, pcor_shrink, pvar_shrink
  implicit none
  real(dp), parameter :: tol = 5.0e-8_dp
  real(dp) :: r(4, 4), x(6, 4), w(6)
  real(dp), allocatable :: p(:, :), rr(:, :)
  type(matrix_shrinkage_result) :: ps
  type(vector_shrinkage_result) :: pv

  r = reshape([ &
    1.0_dp, 0.69213620_dp, 0.63920343_dp, 0.71889083_dp, &
    0.69213620_dp, 1.0_dp, 0.49694209_dp, 0.71609211_dp, &
    0.63920343_dp, 0.49694209_dp, 1.0_dp, 0.54282977_dp, &
    0.71889083_dp, 0.71609211_dp, 0.54282977_dp, 1.0_dp], [4, 4])
  p = cor2pcor(r)
  rr = pcor2cor(p)
  if (maxval(abs(rr-r)) > tol) error stop 'partial-correlation roundtrip failed'

  call make_data(x, w)
  ps = pcor_shrink(x, w=w)
  if (maxval(abs(ps%value-p)) > tol) error stop 'pcor shrink mismatch'
  pv = pvar_shrink(x, w=w)
  if (any(pv%value <= 0.0_dp)) error stop 'partial variances must be positive'
  if (size(ps%standardized_partial_variance) /= 4) error stop 'missing standardized partial variances'
  print '(a)', 'test_partial: PASS'
contains
  subroutine make_data(a, weights)
    real(dp), intent(out) :: a(6, 4), weights(6)
    a = reshape([ &
      1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, &
      2.0_dp, 1.0_dp, 5.0_dp, 4.0_dp, 7.0_dp, 8.0_dp, &
      3.0_dp, 4.0_dp, 2.0_dp, 5.0_dp, 6.0_dp, 7.0_dp, &
      4.0_dp, 3.0_dp, 6.0_dp, 7.0_dp, 8.0_dp, 9.0_dp], [6, 4])
    weights = [1.0_dp, 2.0_dp, 1.0_dp, 3.0_dp, 2.0_dp, 1.0_dp]
    weights = weights / sum(weights)
  end subroutine make_data
end program test_partial
