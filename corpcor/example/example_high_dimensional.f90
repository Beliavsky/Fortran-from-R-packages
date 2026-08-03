program example_high_dimensional
  use corpcor, only : dp, matrix_shrinkage_result, invcov_shrink, rank_condition
  use corpcor, only : rank_condition_result
  implicit none
  real(dp) :: x(12, 40)
  type(matrix_shrinkage_result) :: precision
  type(rank_condition_result) :: rc
  integer :: i, j

  do j = 1, 40
    do i = 1, 12
      x(i, j) = sin(0.031_dp * real(i*j, dp)) + cos(0.07_dp * real(i+j, dp))
    end do
  end do
  precision = invcov_shrink(x)
  rc = rank_condition(precision%value)
  print '(a,i0,a,i0)', 'data dimensions: ', size(x, 1), ' x ', size(x, 2)
  print '(a,f9.6)', 'estimated lambda: ', precision%lambda
  print '(a,i0)', 'precision rank: ', rc%rank
  print '(a,es12.4)', 'precision condition number: ', rc%condition
end program example_high_dimensional
