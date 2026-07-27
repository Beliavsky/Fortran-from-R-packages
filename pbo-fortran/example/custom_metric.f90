! SPDX-License-Identifier: MIT
! Copyright (c) 2014 Matthew R. Barry
module custom_metric_mod
  use pbo, only : dp
  implicit none
  private
  public :: downside_score
contains
  subroutine downside_score(data, values)
    real(dp), intent(in) :: data(:,:)
    real(dp), intent(out) :: values(:)
    integer :: j

    do j = 1, size(data,2)
      values(j) = sum(data(:,j)) - 2.0_dp * sum(max(-data(:,j),0.0_dp))
    end do
  end subroutine downside_score
end module custom_metric_mod

program custom_metric
  use pbo, only : dp, pbo_result, compute_pbo
  use custom_metric_mod, only : downside_score
  implicit none
  type(pbo_result) :: result
  real(dp) :: x(12,3)
  integer :: i

  do i = 1, 12
    x(i,1) = 0.01_dp + 0.002_dp * sin(real(i,dp))
    x(i,2) = 0.005_dp + 0.004_dp * cos(real(i,dp))
    x(i,3) = 0.008_dp - 0.0005_dp * real(i,dp)
  end do
  call compute_pbo(x, 4, downside_score, result)
  if (.not. result%success) error stop result%message
  print '(a,f8.4)', 'PBO for custom downside score = ', result%phi
end program custom_metric
