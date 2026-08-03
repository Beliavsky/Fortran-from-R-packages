! SPDX-License-Identifier: Apache-2.0
module test_support
  use intraday_model
  implicit none
  private
  public :: check, make_parameters
contains
  subroutine check(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write(*, '(a)') 'FAIL: ' // trim(label)
      error stop 1
    end if
  end subroutine check

  function make_parameters(n_bin) result(par)
    integer, intent(in) :: n_bin
    type(volume_parameters) :: par
    integer :: i
    real(dp) :: x

    allocate(par%phi(n_bin))
    do i = 1, n_bin
      x = 2.0_dp * acos(-1.0_dp) * real(i - 1, dp) / real(n_bin, dp)
      par%phi(i) = -0.18_dp * cos(x) + 0.08_dp * cos(2.0_dp * x)
    end do
    par%phi = par%phi - sum(par%phi) / real(n_bin, dp)
    par%a_eta = 0.85_dp
    par%a_mu = 0.45_dp
    par%var_eta = 0.015_dp
    par%var_mu = 0.010_dp
    par%r = 0.005_dp
    par%x0 = [3.0_dp, 0.0_dp]
    par%v0 = 0.0_dp
    par%v0(1, 1) = 0.02_dp
    par%v0(2, 2) = 0.01_dp
  end function make_parameters
end module test_support
