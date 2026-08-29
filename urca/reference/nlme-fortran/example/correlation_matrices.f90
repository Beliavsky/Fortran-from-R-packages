! SPDX-License-Identifier: GPL-2.0-or-later
program correlation_matrices
  use nlme
  implicit none
  real(dp) :: time(5)=[0.0_dp,1.0_dp,2.0_dp,3.0_dp,4.0_dp]
  real(dp),allocatable :: matrix(:,:)
  type(correlation_spec) :: corr
  integer :: status,i

  corr%kind=COR_CAR1
  allocate(corr%par(1))
  corr%par=0.65_dp
  call correlation_matrix(corr,time,matrix,status)
  write(*,'(a)')'continuous AR(1) correlation matrix:'
  do i=1,size(matrix,1)
    write(*,'(*(f9.5,1x))')matrix(i,:)
  end do
end program correlation_matrices
