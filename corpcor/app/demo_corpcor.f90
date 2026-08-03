program demo_corpcor
  use corpcor, only : dp, matrix_shrinkage_result, cov_shrink, invcov_shrink, pcor_shrink
  implicit none
  real(dp) :: x(20, 8)
  type(matrix_shrinkage_result) :: covariance, precision, partial
  integer :: i, j

  do j = 1, 8
    do i = 1, 20
      x(i, j) = 0.5_dp * sin(0.09_dp * real(i*j, dp)) + &
        0.3_dp * cos(0.14_dp * real(i, dp)) + 0.04_dp * real(j, dp)
    end do
  end do
  covariance = cov_shrink(x)
  precision = invcov_shrink(x)
  partial = pcor_shrink(x)

  print '(a)', 'corpcor modern Fortran demonstration'
  print '(a,i0,a,i0)', 'observations: ', size(x, 1), ', variables: ', size(x, 2)
  print '(a,f9.6)', 'estimated correlation lambda: ', covariance%lambda
  print '(a,f9.6)', 'estimated variance lambda:    ', covariance%lambda_var
  print '(a,es12.4)', 'inverse consistency error:    ', &
    maxval(abs(matmul(covariance%value, precision%value) - identity(8)))
  print '(a,f9.6)', 'largest absolute partial correlation: ', &
    maxval(abs(partial%value - identity(8)))
contains
  pure function identity(n) result(out)
    integer, intent(in) :: n
    real(dp) :: out(n, n)
    integer :: k
    out = 0.0_dp
    do k = 1, n
      out(k, k) = 1.0_dp
    end do
  end function identity
end program demo_corpcor
