! SPDX-License-Identifier: BSD-2-Clause
program example_vector_metrics
  use fastcluster
  implicit none

  real(dp) :: x(5, 3)
  type(hclust_result) :: result
  character(len=10) :: metrics(6)
  integer :: i

  x = reshape([0.0_dp, 1.0_dp, 2.0_dp, 4.0_dp, 5.0_dp, &
               0.0_dp, 1.0_dp, 0.0_dp, 4.0_dp, 3.0_dp, &
               1.0_dp, 1.0_dp, 0.0_dp, 2.0_dp, 2.0_dp], shape(x))
  metrics = [character(len=10) :: 'euclidean', 'maximum', 'manhattan', &
    'canberra', 'binary', 'minkowski']

  do i = 1, size(metrics)
    if (trim(metrics(i)) == 'minkowski') then
      call hclust_vector(x, 'single', result, metric=trim(metrics(i)), p=3.0_dp)
    else
      call hclust_vector(x, 'single', result, metric=trim(metrics(i)))
    end if
    if (.not. result%ok()) error stop result%message
    write (*, '(a10,2x,*(f10.5,1x))') trim(metrics(i)), result%height
  end do
end program example_vector_metrics
