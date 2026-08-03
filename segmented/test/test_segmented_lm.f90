! SPDX-License-Identifier: GPL-2.0-or-later
program test_segmented_lm
  use segmented
  implicit none
  integer, parameter :: n = 180
  real(dp) :: y(n), x(n, 2), z(n, 2), psi(2), value
  real(dp), allocatable :: prediction(:), slopes(:), slope_se(:), ci(:,:)
  type(segmented_result) :: fit
  integer :: i, status

  do i = 1, n
    value = 12.0_dp * real(i - 1, dp) / real(n - 1, dp)
    x(i, :) = [1.0_dp, value]
    z(i, :) = value
    y(i) = 1.0_dp + 0.2_dp * value + max(0.0_dp, value - 3.0_dp) - &
        0.8_dp * max(0.0_dp, value - 8.0_dp) + 0.02_dp * sin(real(i, dp))
  end do
  psi = [4.5_dp, 6.5_dp]
  call fit_segmented_lm(y, x, z, psi, fit)
  call check(fit%status == SEG_SUCCESS, 'segmented LM status')
  call check(fit%converged, 'segmented LM convergence')
  call check(maxval(abs(fit%breakpoints - [3.0_dp, 8.0_dp])) < 0.03_dp, &
      'two breakpoints')
  call check(maxval(abs(fit%coefficients - [1.0_dp, 0.2_dp, 1.0_dp, -0.8_dp])) &
      < 0.01_dp, 'piecewise coefficients')
  call predict_segmented(fit, x, z, prediction, status)
  call check(status == SEG_SUCCESS, 'prediction status')
  call check(maxval(abs(prediction - fit%fitted)) < 1.0e-12_dp, 'fitted prediction')
  call segment_slopes(fit, 2, slopes, slope_se, status)
  call check(status == SEG_SUCCESS, 'slope status')
  call check(maxval(abs(slopes - [0.2_dp, 1.2_dp, 0.4_dp])) < 0.01_dp, &
      'segment slopes')
  call breakpoint_confint(fit, 0.95_dp, ci, status)
  call check(status == SEG_SUCCESS .and. size(ci, 1) == 2, 'breakpoint intervals')
  call check(all(fit%breakpoint_se > 0.0_dp), 'breakpoint standard errors')
  call check(abs(aapc(fit, 0.0_dp, 12.0_dp, 2, status=status) - (8.2_dp / 12.0_dp)) &
      < 0.01_dp, 'average slope')
  call check(abs(broken_line_values(10.0_dp, 1.0_dp, 0.2_dp, &
      [1.0_dp, -0.8_dp], [3.0_dp, 8.0_dp]) - 8.4_dp) < 1.0e-12_dp, &
      'broken line evaluation')
  print *, 'test_segmented_lm: PASS'
contains
  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') 'FAIL: ' // trim(message)
      error stop 1
    end if
  end subroutine check
end program test_segmented_lm
