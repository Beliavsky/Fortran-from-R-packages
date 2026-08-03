! SPDX-License-Identifier: BSD-2-Clause
program test_vector
  use fastcluster
  implicit none

  real(dp), allocatable :: d(:, :)
  real(dp) :: x(4, 2)
  type(hclust_result) :: from_matrix, from_vector

  x = reshape([0.0_dp, 1.0_dp, 0.0_dp, 5.0_dp, &
               0.0_dp, 0.0_dp, 2.0_dp, 5.0_dp], shape(x))

  call hclust_vector(x, 'ward', from_vector)
  call check(from_vector%ok(), 'vector Ward status')
  call check_heights(from_vector%height, [1.0_dp, sqrt(17.0_dp / 3.0_dp), &
    sqrt(1825.0_dp / 30.0_dp)], 'vector Ward')

  call hclust_vector(x, 'centroid', from_vector)
  call check(from_vector%ok(), 'vector centroid status')
  call check_heights(from_vector%height, [1.0_dp, sqrt(17.0_dp / 4.0_dp), &
    sqrt(365.0_dp / 9.0_dp)], 'vector centroid')

  call hclust_vector(x, 'median', from_vector)
  call check(from_vector%ok(), 'vector median status')
  call check_heights(from_vector%height, [1.0_dp, sqrt(17.0_dp / 4.0_dp), &
    sqrt(617.0_dp / 16.0_dp)], 'vector median')

  call pairwise_distances(x, 'manhattan', d)
  call hclust_matrix(d, 'single', from_matrix)
  call hclust_vector(x, 'single', from_vector, metric='manhattan')
  call check(from_vector%ok(), 'single Manhattan status')
  call check(all(from_vector%merge == from_matrix%merge), 'single Manhattan merge')
  call check(maxval(abs(from_vector%height - from_matrix%height)) < 1.0e-13_dp, &
    'single Manhattan heights')

  call pairwise_distances(x, 'maximum', d)
  call hclust_matrix(d, 'single', from_matrix)
  call hclust_vector(x, 'single', from_vector, metric='maximum')
  call check(all(from_vector%merge == from_matrix%merge), 'single maximum merge')

  call hclust_vector(x, 'single', from_vector, metric='minkowski', p=3.0_dp)
  call check(from_vector%ok(), 'single Minkowski status')
  call check(from_vector%metric == 'minkowski', 'metric name')

  print '(a)', 'test_vector: PASS'

contains

  subroutine check_heights(actual, expected, label)
    real(dp), intent(in) :: actual(:), expected(:)
    character(len=*), intent(in) :: label

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

end program test_vector
