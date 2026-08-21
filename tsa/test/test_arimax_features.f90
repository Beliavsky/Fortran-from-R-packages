program test_arimax_features
  use tsa
  use tseries_random, only : seed_random, random_normal
  implicit none

  integer :: failures
  failures = 0
  call test_xreg(failures)
  call test_integrated(failures)
  call test_seasonal(failures)

  if (failures /= 0) then
    write(*,'(a,i0)') 'test_arimax_features: FAIL ', failures
    error stop 1
  end if
  write(*,'(a)') 'test_arimax_features: PASS'

contains

  subroutine check(condition, name, failures)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: name
    integer, intent(inout) :: failures
    if (.not. condition) then
      failures = failures + 1
      write(*,'(a,a)') 'FAIL: ', trim(name)
    end if
  end subroutine check

  subroutine test_xreg(failures)
    integer, intent(inout) :: failures
    integer, parameter :: n = 320
    real(dp) :: y(n), e(n), xr(n,1)
    integer :: i
    type(arimax_result) :: fit

    call seed_random(73191)
    e(1) = 0.1_dp*random_normal()
    xr(1,1) = sin(0.031_dp)
    y(1) = 1.25_dp + 2.0_dp*xr(1,1) + e(1)
    do i = 2, n
      xr(i,1) = sin(0.031_dp*real(i,dp)) + &
        0.3_dp*cos(0.071_dp*real(i,dp))
      e(i) = 0.55_dp*e(i-1) + 0.12_dp*random_normal()
      y(i) = 1.25_dp + 2.0_dp*xr(i,1) + e(i)
    end do
    fit = arimax_fit(y, 1, 0, 0, xr)
    call check(fit%status == 0, 'ARIMAX xreg status', failures)
    call check(size(fit%regression) == 2, 'ARIMAX regression size', failures)
    if (size(fit%regression) == 2) then
      call check(abs(fit%regression(1)-1.25_dp) < 0.15_dp, &
        'ARIMAX intercept', failures)
      call check(abs(fit%regression(2)-2.0_dp) < 0.15_dp, &
        'ARIMAX slope', failures)
    end if
    if (size(fit%ar) == 1) then
      call check(abs(fit%ar(1)-0.55_dp) < 0.15_dp, 'ARIMAX AR(1)', failures)
    end if
  end subroutine test_xreg

  subroutine test_integrated(failures)
    integer, intent(inout) :: failures
    integer, parameter :: n = 220
    real(dp) :: y(n)
    integer :: i
    type(arimax_result) :: fit

    call seed_random(1193)
    y(1) = 0.0_dp
    do i = 2, n
      y(i) = y(i-1) + 0.3_dp*random_normal()
    end do
    fit = arima_fit(y, 0, 1, 0)
    call check(fit%status == 0, 'ARIMA(0,1,0) status', failures)
    call check(size(fit%coefficients) == 0, &
      'ARIMA(0,1,0) zero parameters', failures)
    call check(fit%sigma2 > 0.04_dp .and. fit%sigma2 < 0.16_dp, &
      'ARIMA(0,1,0) variance', failures)
  end subroutine test_integrated

  subroutine test_seasonal(failures)
    integer, intent(inout) :: failures
    integer, parameter :: n = 480, s = 12
    real(dp) :: y(n)
    integer :: i
    type(arimax_result) :: fit

    call seed_random(81021)
    y(:s) = 0.15_dp
    do i = s+1, n
      y(i) = 0.62_dp*y(i-s) + 0.16_dp*random_normal()
    end do
    fit = arima_fit(y, 0, 0, 0, seasonal_p=1, period=s)
    call check(fit%status == 0, 'seasonal AR status', failures)
    if (size(fit%sar) == 1) then
      call check(abs(fit%sar(1)-0.62_dp) < 0.15_dp, &
        'seasonal AR coefficient', failures)
    else
      call check(.false., 'seasonal AR coefficient size', failures)
    end if
  end subroutine test_seasonal
end program test_arimax_features
