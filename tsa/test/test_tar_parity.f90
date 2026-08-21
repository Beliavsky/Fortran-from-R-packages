program test_tar_parity
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
  use tsa, only : dp, tar_result, tar_multi_result, tar_fit, tar_fit_multi
  use tseries_random, only : seed_random, random_normal
  implicit none

  integer :: failures
  failures = 0
  call test_univariate_missing(failures)
  call test_multiseries_cls(failures)

  if (failures == 0) then
    print '(a)', 'test_tar_parity: PASS'
  else
    print '(a,i0)', 'test_tar_parity: FAIL ', failures
    error stop 1
  end if

contains

  subroutine test_univariate_missing(failures)
    integer, intent(inout) :: failures
    real(dp) :: y(240), nanv
    integer :: i, st
    type(tar_result) :: fit

    nanv = ieee_value(0.0_dp,ieee_quiet_nan)
    y(1) = 1.0_dp
    do i = 2, size(y)
      if (y(i-1) <= 1.0_dp) then
        y(i) = 0.8_dp + 0.15_dp*y(i-1) + 0.015_dp*sin(0.31_dp*real(i,dp))
      else
        y(i) = 1.2_dp - 0.20_dp*y(i-1) + 0.015_dp*cos(0.23_dp*real(i,dp))
      end if
      y(i) = max(y(i),0.05_dp)
    end do
    y(80) = nanv
    y(151) = nanv

    call tar_fit(y,1,1,1,fit,estimate_threshold=.false.,threshold=0.0_dp, &
      method='CLS',transform='sqrt',standard=.true.,status=st)
    call check(st == 0,'univariate TAR accepts missing windows',failures)
    call check(ieee_is_finite(fit%transform_mean),'finite transform mean',failures)
    call check(fit%n1+fit%n2 < size(y)-1,'missing windows omitted',failures)
    call check(.not. ieee_is_finite(fit%residuals(80)),'missing residual retained as NaN',failures)
  end subroutine test_univariate_missing

  subroutine test_multiseries_cls(failures)
    integer, intent(inout) :: failures
    real(dp) :: y(360,2), nanv, e
    real(dp), parameter :: b10=0.45_dp,b11=0.65_dp,b20=-0.45_dp,b21=0.65_dp
    real(dp), parameter :: c10=0.55_dp,c11=0.50_dp,c20=-0.55_dp,c21=0.50_dp
    integer :: i, s, st
    type(tar_multi_result) :: fit

    nanv = ieee_value(0.0_dp,ieee_quiet_nan)
    call seed_random(7319)
    y(1,:) = [-0.4_dp,0.35_dp]
    do i = 2, size(y,1)
      do s = 1, 2
        e = 0.08_dp*random_normal()
        if (y(i-1,s) <= 0.0_dp) then
          if (s == 1) then
            y(i,s) = b10+b11*y(i-1,s)+e
          else
            y(i,s) = c10+c11*y(i-1,s)+e
          end if
        else
          if (s == 1) then
            y(i,s) = b20+b21*y(i-1,s)+e
          else
            y(i,s) = c20+c21*y(i-1,s)+e
          end if
        end if
      end do
    end do
    y(121,1)=nanv
    y(247,2)=nanv

    call tar_fit_multi(y,1,1,1,fit,estimate_threshold=.false.,threshold=0.0_dp, &
      method='CLS',order_select=.false.,status=st)
    call check(st == 0,'multiseries TAR status',failures)
    call check(fit%nseries == 2,'multiseries TAR nseries',failures)
    call check(size(fit%coefficients1) == 4,'lower coefficient block size',failures)
    call check(size(fit%coefficients2) == 4,'upper coefficient block size',failures)
    if (st == 0 .and. size(fit%coefficients1) == 4) then
      call check(abs(fit%coefficients1(1)-b10) < 0.06_dp,'lower baseline intercept',failures)
      call check(abs(fit%coefficients1(2)-b11) < 0.10_dp,'lower baseline AR',failures)
      call check(abs(fit%coefficients1(3)-(c10-b10)) < 0.08_dp,'lower series delta intercept',failures)
      call check(abs(fit%coefficients1(4)-(c11-b11)) < 0.13_dp,'lower series delta AR',failures)
      call check(abs(fit%coefficients2(1)-b20) < 0.06_dp,'upper baseline intercept',failures)
      call check(abs(fit%coefficients2(2)-b21) < 0.10_dp,'upper baseline AR',failures)
      call check(abs(fit%coefficients2(3)-(c20-b20)) < 0.08_dp,'upper series delta intercept',failures)
      call check(abs(fit%coefficients2(4)-(c21-b21)) < 0.13_dp,'upper series delta AR',failures)
    end if
    call check(all(fit%n1 > 20) .and. all(fit%n2 > 20),'each series represented in both regimes',failures)
    call check(.not. ieee_is_finite(fit%residuals(121,1)),'matrix missing residual NaN',failures)
  end subroutine test_multiseries_cls

  subroutine check(ok,label,failures)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: label
    integer, intent(inout) :: failures
    if (.not. ok) then
      print '(a)', 'FAIL: '//trim(label)
      failures = failures+1
    end if
  end subroutine check

end program test_tar_parity
