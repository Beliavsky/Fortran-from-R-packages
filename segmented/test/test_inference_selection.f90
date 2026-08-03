! SPDX-License-Identifier: GPL-2.0-or-later
program test_inference_selection
  use segmented
  implicit none
  integer, parameter :: n = 140
  real(dp) :: y(n), x(n, 2), z(n), value, power
  real(dp), allocatable :: bic(:)
  type(test_result) :: davies, score
  type(segmented_result) :: selected
  integer :: i, status

  do i = 1, n
    value = 10.0_dp * real(i - 1, dp) / real(n - 1, dp)
    x(i, :) = [1.0_dp, value]
    z(i) = value
    y(i) = 0.7_dp + 0.4_dp * value + 1.1_dp * max(0.0_dp, value - 4.2_dp) + &
        0.08_dp * sin(1.7_dp * real(i, dp))
  end do
  call davies_test(y, x, z, 35, davies)
  call check(davies%status == SEG_SUCCESS, 'Davies status')
  call check(davies%p_value < 1.0e-6_dp, 'Davies significance')
  call check(abs(davies%breakpoint - 4.2_dp) < 0.5_dp, 'Davies location')
  call pscore_test(y, x, z, 4.2_dp, score)
  call check(score%status == SEG_SUCCESS .and. abs(score%statistic) > 20.0_dp, &
      'score test')
  power = pwr_seg(z, x, 4.2_dp, 1.1_dp, 0.08_dp, status=status)
  call check(status == SEG_SUCCESS .and. power > 0.999_dp, 'power calculation')
  call select_breakpoints_bic(y, x, z, 3, FAMILY_GAUSSIAN, selected, bic)
  call check(selected%status == SEG_SUCCESS, 'BIC selection status')
  call check(selected%n_break == 1, 'BIC selects one breakpoint')
  call check(abs(selected%breakpoints(1) - 4.2_dp) < 0.05_dp, &
      'selected breakpoint')
  print *, 'test_inference_selection: PASS'
contains
  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') 'FAIL: ' // trim(message)
      error stop 1
    end if
  end subroutine check
end program test_inference_selection
