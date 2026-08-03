program example_covariance_shrinkage
  use corpcor, only : dp, matrix_shrinkage_result, cov_shrink
  implicit none
  real(dp) :: x(8, 5)
  type(matrix_shrinkage_result) :: fit
  integer :: i, j

  do j = 1, 5
    do i = 1, 8
      x(i, j) = sin(0.2_dp * real(i*j, dp)) + 0.15_dp * real(j, dp)
    end do
  end do
  fit = cov_shrink(x)
  print '(a,f8.5)', 'correlation shrinkage intensity: ', fit%lambda
  print '(a,f8.5)', 'variance shrinkage intensity:    ', fit%lambda_var
  print '(a)'; print '(a)', 'shrinkage covariance:'
  do i = 1, 5
    print '(*(f11.6,1x))', fit%value(i, :)
  end do
end program example_covariance_shrinkage
