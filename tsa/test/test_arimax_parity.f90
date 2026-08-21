program test_arimax_parity
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
  use tsa
  use tseries_random, only : seed_random, random_normal
  implicit none

  integer :: failures
  failures = 0
  call test_transfer(failures)
  call test_io(failures)
  call test_fixed(failures)
  call test_ml_missing(failures)
  call test_ma_ml(failures)
  call test_ma_invert(failures)

  if (failures /= 0) then
    write(*,'(a,i0)') 'test_arimax_parity: FAIL ', failures
    error stop 1
  end if
  write(*,'(a)') 'test_arimax_parity: PASS'

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

  subroutine test_transfer(failures)
    integer, intent(inout) :: failures
    integer, parameter :: n = 360
    real(dp) :: y(n), xt(n,1), empty(n,0)
    real(dp), allocatable :: effect(:)
    real(dp) :: phi(1), theta(2)
    type(transfer_spec) :: tr(1)
    type(arimax_result) :: fit
    integer :: i

    call seed_random(41931)
    do i = 1, n
      xt(i,1) = sin(0.047_dp*real(i,dp)) + 0.5_dp*cos(0.113_dp*real(i,dp))
    end do
    phi = [0.42_dp]
    theta = [1.10_dp,-0.36_dp]
    call transfer_filter(xt(:,1),phi,theta,effect)
    do i = 1, n
      y(i) = effect(i) + 0.035_dp*random_normal()
    end do
    tr(1)%ar_order = 1
    tr(1)%ma_order = 1
    fit = arimax_fit(y,0,0,0,empty,include_mean=.false.,xtransf=xt, &
      transfer=tr,method='CSS')
    call check(fit%status == 0,'transfer status',failures)
    call check(size(fit%transfer) == 3,'transfer parameter count',failures)
    if (size(fit%transfer) == 3) then
      call check(abs(fit%transfer(1)-phi(1)) < 0.12_dp,'transfer AR',failures)
      call check(abs(fit%transfer(2)-theta(1)) < 0.12_dp,'transfer MA0',failures)
      call check(abs(fit%transfer(3)-theta(2)) < 0.12_dp,'transfer MA1',failures)
    end if
  end subroutine test_transfer

  subroutine test_io(failures)
    integer, intent(inout) :: failures
    integer, parameter :: n = 420, point = 180
    real(dp) :: y(n), e(n), empty(n,0)
    real(dp), allocatable :: pulse_effect(:)
    type(arimax_result) :: fit
    integer :: i, st

    call seed_random(58123)
    call io_regressor(n,point,[0.55_dp],[real(dp)::],0,0,1,pulse_effect,st)
    e(1) = 0.10_dp*random_normal()
    y(1) = e(1) + 2.2_dp*pulse_effect(1)
    do i = 2, n
      e(i) = 0.55_dp*e(i-1) + 0.10_dp*random_normal()
      y(i) = e(i) + 2.2_dp*pulse_effect(i)
    end do
    fit = arimax_fit(y,1,0,0,empty,include_mean=.false.,io=[point],method='CSS')
    call check(fit%status == 0,'IO status',failures)
    if (size(fit%ar) == 1) call check(abs(fit%ar(1)-0.55_dp) < 0.12_dp,'IO AR',failures)
    if (size(fit%regression) == 1) then
      call check(abs(fit%regression(1)-2.2_dp) < 0.30_dp,'IO coefficient',failures)
    else
      call check(.false.,'IO coefficient size',failures)
    end if
  end subroutine test_io

  subroutine test_fixed(failures)
    integer, intent(inout) :: failures
    integer, parameter :: n = 260
    real(dp) :: y(n), fixed(2), nanv
    type(arimax_result) :: fit
    integer :: i

    nanv = ieee_value(0.0_dp,ieee_quiet_nan)
    call seed_random(8401)
    y(1) = 1.4_dp + 0.1_dp*random_normal()
    do i = 2, n
      y(i) = 1.4_dp + 0.50_dp*(y(i-1)-1.4_dp) + 0.1_dp*random_normal()
    end do
    fixed = [0.50_dp,nanv]
    fit = arima_fit(y,1,0,0,fixed=fixed,method='CSS')
    call check(fit%status == 0,'fixed status',failures)
    call check(abs(fit%ar(1)-0.50_dp) < 1.0e-12_dp,'fixed AR exact',failures)
    call check(.not. fit%estimated(1) .and. fit%estimated(2),'fixed mask',failures)
    call check(abs(fit%regression(1)-1.4_dp) < 0.15_dp,'fixed intercept',failures)
  end subroutine test_fixed

  subroutine test_ml_missing(failures)
    integer, intent(inout) :: failures
    integer, parameter :: n = 360
    real(dp) :: y(n), nanv
    type(arimax_result) :: fit
    integer :: i

    nanv = ieee_value(0.0_dp,ieee_quiet_nan)
    call seed_random(15517)
    y(1) = 0.12_dp*random_normal()
    do i = 2, n
      y(i) = 0.68_dp*y(i-1) + 0.12_dp*random_normal()
    end do
    y(91) = nanv
    y(217) = nanv
    fit = arima_fit(y,1,0,0,include_mean=.false.,method='ML')
    call check(fit%status == 0,'ML missing status',failures)
    call check(abs(fit%ar(1)-0.68_dp) < 0.15_dp,'ML missing AR',failures)
    call check(ieee_is_finite(fit%loglik) .and. ieee_is_finite(fit%aic), &
      'ML likelihood/AIC',failures)
  end subroutine test_ml_missing

  subroutine test_ma_ml(failures)
    integer, intent(inout) :: failures
    integer, parameter :: n = 420
    real(dp) :: y(n), e(n)
    type(arimax_result) :: fit
    integer :: i

    call seed_random(99127)
    e = 0.0_dp
    e(1) = 0.11_dp*random_normal()
    y(1) = e(1)
    do i = 2, n
      e(i) = 0.11_dp*random_normal()
      y(i) = e(i) + 0.45_dp*e(i-1)
    end do
    fit = arima_fit(y,0,0,1,include_mean=.false.,method='ML')
    call check(fit%status == 0,'MA ML status',failures)
    call check(abs(fit%ma(1)-0.45_dp) < 0.15_dp,'MA direct coefficient',failures)
  end subroutine test_ma_ml

  subroutine test_ma_invert(failures)
    integer, intent(inout) :: failures
    integer, parameter :: n = 300
    real(dp) :: y(n), e(n), init(1)
    type(arimax_result) :: fit
    integer :: i

    call seed_random(67433)
    e(1) = 0.10_dp*random_normal()
    y(1) = e(1)
    do i = 2, n
      e(i) = 0.10_dp*random_normal()
      y(i) = e(i)+0.50_dp*e(i-1)
    end do
    init = [2.0_dp]
    fit = arima_fit(y,0,0,1,include_mean=.false.,method='ML',init=init)
    call check(fit%status == 0,'MA invertibility status',failures)
    call check(abs(fit%ma(1)) < 1.0_dp,'MA normalized to invertible root',failures)
  end subroutine test_ma_invert

end program test_arimax_parity
