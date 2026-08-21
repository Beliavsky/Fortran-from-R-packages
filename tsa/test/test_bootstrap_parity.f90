program test_bootstrap_parity
  use tsa, only : dp, arimax_result, arima_bootstrap_sample
  implicit none

  integer :: failures
  failures = 0
  call test_conditional_prefix(failures)
  call test_tsa_ma_convolution(failures)

  if (failures == 0) then
    print '(a)', 'test_bootstrap_parity: PASS'
  else
    print '(a,i0)', 'test_bootstrap_parity: FAIL ', failures
    error stop 1
  end if

contains

  subroutine test_conditional_prefix(failures)
    integer, intent(inout) :: failures
    type(arimax_result) :: fit
    real(dp), allocatable :: x(:)
    real(dp) :: initial(2)
    integer :: st, i

    fit%p = 1
    fit%d = 1
    fit%q = 0
    fit%include_mean = .false.
    fit%sigma2 = 0.04_dp
    allocate(fit%ar(1),fit%ma(0),fit%residuals(80),fit%series(80))
    fit%ar = [0.25_dp]
    fit%residuals = 0.0_dp
    fit%series = [(real(i,dp),i=1,80)]
    initial = [3.25_dp,-1.50_dp]

    call arima_bootstrap_sample(fit,x,normal=.true.,cond_boot=.true., &
      init=initial,status=st)
    call check(st == 0,'conditional sample status',failures)
    call check(size(x) == 80,'conditional sample length',failures)
    call check(maxval(abs(x(:2)-initial)) < 1.0e-12_dp, &
      'conditional bootstrap preserves p+d prefix',failures)
  end subroutine test_conditional_prefix

  subroutine test_tsa_ma_convolution(failures)
    integer, intent(inout) :: failures
    type(arimax_result) :: fit
    real(dp), allocatable :: x(:)
    integer :: st

    fit%p = 0
    fit%d = 0
    fit%q = 1
    fit%include_mean = .false.
    fit%sigma2 = 1.0_dp
    allocate(fit%ar(0),fit%ma(1),fit%residuals(30),fit%series(30))
    fit%ma = [0.5_dp]
    fit%residuals = 2.0_dp
    fit%series = 0.0_dp

    call arima_bootstrap_sample(fit,x,normal=.false.,cond_boot=.false., &
      ntrans=0,status=st)
    call check(st == 0,'MA bootstrap sample status',failures)
    call check(maxval(abs(x-1.0_dp)) < 1.0e-12_dp, &
      'TSA source MA convolution semantics',failures)
  end subroutine test_tsa_ma_convolution

  subroutine check(ok,label,failures)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: label
    integer, intent(inout) :: failures
    if (.not. ok) then
      print '(a)', 'FAIL: '//trim(label)
      failures = failures+1
    end if
  end subroutine check

end program test_bootstrap_parity
