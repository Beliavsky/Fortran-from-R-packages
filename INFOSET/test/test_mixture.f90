! SPDX-License-Identifier: GPL-2.0-or-later
program test_mixture
  use infoset, only : dp, tail_mixture, infoset_estimate
  use infoset, only : tail_mixture_result, information_set_result
  use infoset, only : infoset_success, infoset_not_converged
  implicit none
  integer, parameter :: n = 400
  real(dp) :: y(n)
  type(tail_mixture_result) :: fit
  type(information_set_result) :: information
  integer :: i

  do i = 1, 170
    y(i) = exp(-0.085_dp + 0.018_dp*sin(0.37_dp*real(i,dp)) &
      + 0.006_dp*cos(0.11_dp*real(i,dp)))
  end do
  do i = 171, n
    y(i) = exp(0.035_dp + 0.022_dp*sin(0.29_dp*real(i,dp)) &
      + 0.005_dp*cos(0.17_dp*real(i,dp)))
  end do

  call tail_mixture(y, 0.0_dp, 1, fit)
  call check(fit%status == infoset_success .or. &
    fit%status == infoset_not_converged, 'tail mixture status')
  call check(fit%left_mean < fit%right_mean, 'ordered mixture means')
  call check(fit%left_probability > 0.2_dp .and. &
    fit%left_probability < 0.8_dp, 'mixture probability')
  call check(fit%change_point > minval(y) .and. &
    fit%change_point < maxval(y), 'change point range')
  call check(fit%first_type_error >= 0.0_dp .and. &
    fit%first_type_error <= 1.0_dp, 'first error range')
  call check(fit%second_type_error >= 0.0_dp .and. &
    fit%second_type_error <= 1.0_dp, 'second error range')

  call infoset_estimate(y, information)
  call check(information%n_change_points >= 1, 'information set split')
  call check(information%change_points(1) < 1.0_dp, 'left change point')
  print '(a)', 'test_mixture: PASS'
contains
  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*,'(a)') 'FAIL: '//message
      error stop 1
    end if
  end subroutine check
end program test_mixture
