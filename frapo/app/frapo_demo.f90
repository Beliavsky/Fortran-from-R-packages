! SPDX-License-Identifier: GPL-3.0-or-later
program frapo_demo
  use frapo
  implicit none

  real(dp) :: returns(60, 4)
  real(dp), allocatable :: covariance(:, :)
  type(portfolio_result) :: gmv, diversified, erc
  integer :: i

  do i = 1, size(returns, 1)
    returns(i, 1) = 0.010_dp * sin(0.17_dp * real(i, dp))
    returns(i, 2) = 0.007_dp * cos(0.11_dp * real(i, dp))
    returns(i, 3) = 0.009_dp * sin(0.07_dp * real(i, dp) + 0.4_dp)
    returns(i, 4) = 0.006_dp * cos(0.19_dp * real(i, dp) - 0.2_dp)
  end do

  call sample_covariance(returns, covariance)
  gmv = pgmv(returns, percentage=.false.)
  diversified = pmd(returns, percentage=.false.)
  erc = perc(covariance, percentage=.false.)

  write(*, '(a)') 'FRAPO modern Fortran demonstration'
  call print_weights('Global minimum variance', gmv%weights)
  call print_weights('Most diversified', diversified%weights)
  call print_weights('Equal risk contribution', erc%weights)

contains

  subroutine print_weights(label, weights)
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: weights(:)
    write(*, '(a)') trim(label)
    write(*, '(*(f11.6,1x))') weights
  end subroutine print_weights
end program frapo_demo
