program test_shrinkage
  use corpcor, only : dp, matrix_shrinkage_result, vector_shrinkage_result, &
    estimate_lambda, estimate_lambda_var, cor_shrink, var_shrink, cov_shrink, invcov_shrink
  implicit none
  real(dp), parameter :: tol = 2.0e-8_dp
  real(dp) :: x(6, 4), w(6), lambda, lambda_var
  real(dp) :: r_ref(4, 4), cov_ref(4, 4)
  type(matrix_shrinkage_result) :: cr, cv, prec
  type(vector_shrinkage_result) :: vr

  call make_data(x, w)
  lambda = estimate_lambda(x, w)
  lambda_var = estimate_lambda_var(x, w)
  if (abs(lambda - 0.2380242847868838_dp) > tol) error stop 'lambda reference failed'
  if (abs(lambda_var - 1.0_dp) > tol) error stop 'lambda.var reference failed'

  r_ref = reshape([ &
    1.0_dp, 0.69213620_dp, 0.63920343_dp, 0.71889083_dp, &
    0.69213620_dp, 1.0_dp, 0.49694209_dp, 0.71609211_dp, &
    0.63920343_dp, 0.49694209_dp, 1.0_dp, 0.54282977_dp, &
    0.71889083_dp, 0.71609211_dp, 0.54282977_dp, 1.0_dp], [4, 4])
  cr = cor_shrink(x, w=w)
  if (maxval(abs(cr%value-r_ref)) > tol) error stop 'correlation shrinkage reference failed'
  vr = var_shrink(x, w=w)
  if (maxval(abs(vr%value - 4.0_dp)) > tol) error stop 'variance shrinkage reference failed'

  cov_ref = 4.0_dp * r_ref
  cv = cov_shrink(x, w=w)
  if (maxval(abs(cv%value-cov_ref)) > 5.0e-8_dp) error stop 'covariance shrinkage reference failed'
  prec = invcov_shrink(x, w=w)
  if (maxval(abs(matmul(cv%value, prec%value) - identity(4))) > 5.0e-7_dp) &
    error stop 'precision inverse identity failed'
  print '(a)', 'test_shrinkage: PASS'
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
  pure function identity(n) result(a)
    integer, intent(in) :: n
    real(dp) :: a(n, n)
    integer :: i
    a = 0.0_dp
    do i = 1, n
      a(i, i) = 1.0_dp
    end do
  end function identity
end program test_shrinkage
