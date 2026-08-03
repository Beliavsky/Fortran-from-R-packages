! SPDX-License-Identifier: BSD-2-Clause
program test_linkage
  use fastcluster
  implicit none

  real(dp), allocatable :: d(:, :)
  real(dp) :: x(4, 2)
  type(hclust_result) :: result

  x = reshape([0.0_dp, 1.0_dp, 0.0_dp, 5.0_dp, &
               0.0_dp, 0.0_dp, 2.0_dp, 5.0_dp], shape(x))
  call pairwise_distances(x, 'euclidean', d)

  call hclust_matrix(d, 'single', result)
  call check(result%ok(), 'single status')
  call check_heights(result%height, [1.0_dp, 2.0_dp, sqrt(34.0_dp)], 'single')
  call check_merge(result)

  call hclust_matrix(d, 'complete', result)
  call check(result%ok(), 'complete status')
  call check_heights(result%height, [1.0_dp, sqrt(5.0_dp), sqrt(50.0_dp)], 'complete')

  call hclust_matrix(d, 'average', result)
  call check(result%ok(), 'average status')
  call check_heights(result%height, [1.0_dp, (2.0_dp + sqrt(5.0_dp)) / 2.0_dp, &
    (sqrt(50.0_dp) + sqrt(41.0_dp) + sqrt(34.0_dp)) / 3.0_dp], 'average')

  call hclust_matrix(d, 'mcquitty', result)
  call check(result%ok(), 'mcquitty status')
  call check_heights(result%height, [1.0_dp, (2.0_dp + sqrt(5.0_dp)) / 2.0_dp, &
    (sqrt(34.0_dp) + (sqrt(50.0_dp) + sqrt(41.0_dp)) / 2.0_dp) / 2.0_dp], &
    'mcquitty')

  call hclust_matrix(d, 'ward.D2', result)
  call check(result%ok(), 'ward.D2 status')
  call check_heights(result%height, [1.0_dp, sqrt(17.0_dp / 3.0_dp), &
    sqrt(1825.0_dp / 30.0_dp)], 'ward.D2')

  call hclust_matrix(d(1:3, 1:3), 'ward.D', result)
  call check(result%ok(), 'ward.D status')
  call check_heights(result%height, [1.0_dp, (3.0_dp + 2.0_dp * sqrt(5.0_dp)) / 3.0_dp], &
    'ward.D')

  call hclust_matrix(d(1:3, 1:3), 'centroid', result)
  call check(result%ok(), 'centroid matrix status')
  call check_heights(result%height, [1.0_dp, 0.75_dp + 0.5_dp * sqrt(5.0_dp)], &
    'centroid matrix')

  call hclust_matrix(d(1:3, 1:3), 'median', result)
  call check(result%ok(), 'median matrix status')
  call check_heights(result%height, [1.0_dp, 0.75_dp + 0.5_dp * sqrt(5.0_dp)], &
    'median matrix')

  print '(a)', 'test_linkage: PASS'

contains

  subroutine check_merge(value)
    type(hclust_result), intent(in) :: value
    integer :: expected_merge(3, 2)

    expected_merge(:, 1) = [-1, -3, -4]
    expected_merge(:, 2) = [-2, 1, 2]
    call check(all(value%merge == expected_merge), 'merge matrix')
    call check(all(value%order == [4, 3, 1, 2]), 'dendrogram order')
  end subroutine check_merge

  subroutine check_heights(actual, expected, label)
    real(dp), intent(in) :: actual(:), expected(:)
    character(len=*), intent(in) :: label

    call check(size(actual) == size(expected), label//' height size')
    call check(maxval(abs(actual - expected)) <= 2.0e-12_dp * &
      max(1.0_dp, maxval(abs(expected))), label//' heights')
  end subroutine check_heights

  subroutine check(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label

    if (.not. condition) then
      write (*, '(a)') 'FAIL: '//label
      error stop 1
    end if
  end subroutine check

end program test_linkage
