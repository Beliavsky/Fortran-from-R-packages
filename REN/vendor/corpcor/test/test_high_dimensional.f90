program test_high_dimensional
  use corpcor, only : dp, matrix_shrinkage_result, cov_shrink, invcov_shrink, &
    is_positive_definite, powcor_shrink, crossprod_powcor_shrink
  implicit none
  real(dp) :: x(5, 12), y(12, 3)
  type(matrix_shrinkage_result) :: cv, iv, p2, cp
  integer :: i, j

  do j = 1, 12
    do i = 1, 5
      x(i, j) = sin(0.17_dp * real(i*j, dp)) + 0.05_dp * real(j, dp) + 0.1_dp * real(i, dp)
    end do
  end do
  do j = 1, 3
    do i = 1, 12
      y(i, j) = cos(0.11_dp * real(i*j, dp))
    end do
  end do
  cv = cov_shrink(x)
  iv = invcov_shrink(x)
  if (.not. is_positive_definite(cv%value)) error stop 'high-dimensional covariance not positive definite'
  if (maxval(abs(matmul(cv%value, iv%value)-identity(12))) > 1.0e-5_dp) &
    error stop 'high-dimensional inverse failed'
  p2 = powcor_shrink(x, 2.0_dp)
  cp = crossprod_powcor_shrink(x, y, 2.0_dp)
  if (maxval(abs(cp%value - matmul(p2%value, y))) > 1.0e-8_dp) error stop 'crossprod power failed'
  print '(a)', 'test_high_dimensional: PASS'
contains
  pure function identity(n) result(a)
    integer, intent(in) :: n
    real(dp) :: a(n, n)
    integer :: k
    a = 0.0_dp
    do k = 1, n
      a(k, k) = 1.0_dp
    end do
  end function identity
end program test_high_dimensional
