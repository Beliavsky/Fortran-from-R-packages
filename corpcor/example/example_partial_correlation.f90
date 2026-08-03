program example_partial_correlation
  use corpcor, only : dp, matrix_shrinkage_result, pcor_shrink
  implicit none
  real(dp) :: x(10, 4)
  type(matrix_shrinkage_result) :: fit
  integer :: i, j

  do j = 1, 4
    do i = 1, 10
      x(i, j) = cos(0.13_dp * real(i*(j+1), dp)) + 0.03_dp * real(i*j, dp)
    end do
  end do
  fit = pcor_shrink(x)
  print '(a)'; print '(a)', 'shrinkage partial correlation:'
  do i = 1, 4
    print '(*(f10.5,1x))', fit%value(i, :)
  end do
  print '(a)'; print '(a)', 'standardized partial variances:'
  print '(*(f10.5,1x))', fit%standardized_partial_variance
end program example_partial_correlation
