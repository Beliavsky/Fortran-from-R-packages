program test_tsa
  use tsa
  use leaps, only : regsubsets_result
  use tseries_random, only : seed_random, random_normal
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none

  integer :: failures
  failures = 0

  call test_moments(failures)
  call test_acf_spectrum(failures)
  call test_runs_and_tests(failures)
  call test_simulation(failures)
  call test_subset_and_eacf(failures)
  call test_tar(failures)
  call test_arima(failures)
  call test_boxcox(failures)

  if (failures /= 0) then
    write(*,'(a,i0)') 'test_tsa: FAIL ', failures
    error stop 1
  end if
  write(*,'(a)') 'test_tsa: PASS'

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

  subroutine test_moments(failures)
    integer, intent(inout) :: failures
    real(dp) :: x(5)
    x = [-2.0_dp,-1.0_dp,0.0_dp,1.0_dp,2.0_dp]
    call check(abs(skewness(x)) < 1.0e-14_dp, 'skewness symmetric', failures)
    call check(abs(kurtosis(x) + 1.3_dp) < 1.0e-12_dp, &
      'TSA population excess kurtosis', failures)
  end subroutine test_moments

  subroutine test_acf_spectrum(failures)
    integer, intent(inout) :: failures
    real(dp), allocatable :: ac(:)
    real(dp) :: x(64), pi
    type(spectrum_result) :: sp
    integer :: i, imax

    pi = acos(-1.0_dp)
    do i = 1, size(x)
      x(i) = sin(2.0_dp*pi*0.125_dp*real(i,dp))
    end do
    call autocorrelation(x, 4, ac)
    call check(size(ac) == 4, 'acf size', failures)
    call check(abs(ac(1)-cos(2.0_dp*pi*0.125_dp)) < 0.08_dp, &
      'acf sinusoid lag 1', failures)
    sp = periodogram(x)
    imax = maxloc(sp%spectrum, dim=1)
    call check(abs(sp%frequency(imax)-0.125_dp) < 1.0e-12_dp, &
      'periodogram peak', failures)
  end subroutine test_acf_spectrum

  subroutine test_runs_and_tests(failures)
    integer, intent(inout) :: failures
    real(dp) :: x(20), e(80)
    type(runs_result) :: rr
    type(tsa_test_result) :: lb, kt, tt
    integer :: i

    do i = 1, size(x)
      x(i) = merge(-1.0_dp,1.0_dp,mod(i,2)==0)
    end do
    rr = runs(x)
    call check(rr%observed_runs == 20, 'runs alternating count', failures)
    call check(rr%p_value >= 0.0_dp .and. rr%p_value <= 1.0_dp, &
      'runs p-value range', failures)

    do i = 1, size(e)
      e(i) = sin(0.7_dp*real(i,dp)) + 0.2_dp*cos(1.3_dp*real(i,dp))
    end do
    lb = lb_test(e, 8, 1)
    call check(lb%status == 0, 'LB test status', failures)
    kt = keenan_test(e, 2)
    call check(kt%status == 0, 'Keenan test status', failures)
    tt = tsay_test(e, 2)
    call check(tt%status == 0, 'Tsay test status', failures)
  end subroutine test_runs_and_tests

  subroutine test_simulation(failures)
    integer, intent(inout) :: failures
    real(dp), allocatable :: x(:)
    real(dp) :: q(50)
    integer :: st

    call seed_random(12345)
    call qar_sim(q, phi0=0.1_dp, phi1=0.05_dp, sigma=0.1_dp)
    call check(all(ieee_is_finite(q)), 'QAR finite', failures)

    call garch_sim([0.1_dp,0.1_dp], [0.8_dp], 200, 50, x, st)
    call check(st == 0 .and. size(x) == 200, 'GARCH simulation', failures)
    call check(sum(x*x) > 0.0_dp, 'GARCH nonzero', failures)
  end subroutine test_simulation

  subroutine test_subset_and_eacf(failures)
    integer, intent(inout) :: failures
    real(dp) :: y(120)
    type(regsubsets_result) :: subs
    real(dp), allocatable :: em(:,:)
    character(len=1), allocatable :: sym(:,:)
    integer :: i, st

    call seed_random(22011)
    y(1:2) = [0.2_dp,-0.1_dp]
    do i = 3, size(y)
      y(i) = 0.65_dp*y(i-1) - 0.2_dp*y(i-2) + &
        0.08_dp*random_normal()
    end do
    call armasubsets_fit(y, 3, 2, subs, nvmax=4, nbest=2, status=st)
    call check(st == 0, 'armasubsets status', failures)
    call check(subs%nvar == 5, 'armasubsets predictors', failures)

    call eacf(y, 3, 4, em, sym, st)
    call check(st == 0, 'eacf status', failures)
    call check(all(shape(em) == [4,5]), 'eacf dimensions', failures)
  end subroutine test_subset_and_eacf

  subroutine test_tar(failures)
    integer, intent(inout) :: failures
    real(dp) :: y(300)
    type(tar_result) :: fit
    type(tsa_test_result) :: lr
    integer :: i, st

    call seed_random(92831)
    y(1) = 0.2_dp
    do i = 2, size(y)
      if (y(i-1) <= 0.0_dp) then
        y(i) = -0.15_dp + 0.45_dp*y(i-1)
      else
        y(i) = 0.20_dp - 0.35_dp*y(i-1)
      end if
      y(i) = y(i) + 0.05_dp*random_normal()
    end do
    call tar_fit(y, 1, 1, 1, fit, a=0.1_dp, b=0.9_dp, status=st)
    call check(st == 0, 'TAR fit status', failures)
    call check(abs(fit%threshold) < 0.15_dp, 'TAR threshold vicinity', failures)
    call check(fit%n1 > 20 .and. fit%n2 > 20, 'TAR regime sizes', failures)

    lr = tlrt_test(y, 1, 1)
    call check(lr%status == 0, 'TLRT status', failures)
    call check(lr%p_value >= 0.0_dp .and. lr%p_value <= 1.0_dp, &
      'TLRT p-value range', failures)
  end subroutine test_tar

  subroutine test_arima(failures)
    integer, intent(inout) :: failures
    real(dp), allocatable :: y(:)
    type(arimax_result) :: fit

    call seed_random(90210)
    call arima_sim([0.65_dp], [real(dp) ::], 0, 260, y, sigma=0.35_dp, &
      ntrans=150)
    fit = arima_fit(y, 1, 0, 0, include_mean=.false.)
    call check(fit%status == 0, 'ARIMA fit status', failures)
    call check(size(fit%ar) == 1, 'ARIMA ar size', failures)
    if (size(fit%ar) == 1) then
      call check(abs(fit%ar(1)-0.65_dp) < 0.18_dp, 'ARIMA AR(1) recovery', failures)
    end if
    call check(fit%sigma2 > 0.0_dp, 'ARIMA variance positive', failures)
  end subroutine test_arima

  subroutine test_boxcox(failures)
    integer, intent(inout) :: failures
    real(dp) :: y(60), lambda(5), mle, ci(2)
    real(dp), allocatable :: ll(:)
    integer :: i, st

    do i = 1, size(y)
      y(i) = exp(0.01_dp*real(i,dp) + 0.08_dp*sin(0.4_dp*real(i,dp)))
    end do
    lambda = [-1.0_dp,-0.5_dp,0.0_dp,0.5_dp,1.0_dp]
    call boxcox_ar(y, lambda, 1, ll, mle, ci, st)
    call check(st == 0 .and. size(ll) == 5, 'BoxCox AR status', failures)
    call check(mle >= -1.0_dp .and. mle <= 1.0_dp, 'BoxCox mle grid', failures)
    call check(ci(1) <= ci(2), 'BoxCox CI ordering', failures)
  end subroutine test_boxcox
end program test_tsa
