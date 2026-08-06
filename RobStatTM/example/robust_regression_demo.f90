program robust_regression_demo
  use robstattm, only : dp, robstattm_control, regression_result, lmrobdet_control, &
    lmrobdetmm, lmrobdetdcml
  implicit none
  integer, parameter :: n = 40
  real(dp) :: x(n, 2), y(n), t
  type(robstattm_control) :: control
  type(regression_result) :: mm_fit, dcml_fit
  integer :: i

  do i = 1, n
    t = -1.0_dp + 2.0_dp * real(i - 1, dp) / real(n - 1, dp)
    x(i, :) = [1.0_dp, t]
    y(i) = 1.0_dp + 2.5_dp * t + 0.05_dp * sin(9.0_dp * t)
  end do
  y(4) = y(4) + 8.0_dp
  y(31) = y(31) - 7.0_dp

  control = lmrobdet_control(n_resample=150)
  call lmrobdetmm(x, y, mm_fit, control)
  call lmrobdetdcml(x, y, dcml_fit, control)

  write(*, '(a, *(f10.5, 1x))') 'MM coefficients:   ', mm_fit%coefficients
  write(*, '(a, *(f10.5, 1x))') 'DCML coefficients: ', dcml_fit%coefficients
  write(*, '(a, f10.5)') 'DCML LS mixing:    ', dcml_fit%mixing
end program robust_regression_demo
