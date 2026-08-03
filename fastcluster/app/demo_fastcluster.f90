! SPDX-License-Identifier: BSD-2-Clause
program demo_fastcluster
  use fastcluster
  implicit none

  real(dp) :: x(6, 2)
  type(hclust_result) :: result
  integer :: i

  x = reshape([0.0_dp, 0.2_dp, 1.0_dp, 5.0_dp, 5.2_dp, 8.0_dp, &
               0.0_dp, 0.1_dp, 1.2_dp, 5.0_dp, 4.8_dp, 1.0_dp], shape(x))

  call hclust_vector(x, 'ward', result)
  if (.not. result%ok()) error stop result%message

  write (*, '(a)') 'Ward hierarchical clustering'
  write (*, '(a)') ' step  left right       height'
  do i = 1, result%n - 1
    write (*, '(i5,2i6,f13.6)') i, result%merge(i, 1), result%merge(i, 2), &
      result%height(i)
  end do
  write (*, '(a,*(i0,1x))') 'leaf order: ', result%order
end program demo_fastcluster
