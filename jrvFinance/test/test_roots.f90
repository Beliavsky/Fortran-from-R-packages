! SPDX-License-Identifier: GPL-2.0-or-later
module root_test_callbacks
  use jrvfinance, only: dp
  implicit none
  type :: empty_context
    integer :: unused = 0
  end type empty_context
contains
  subroutine trig_callback(x, context, value, gradient)
    real(dp), intent(in) :: x
    class(*), intent(in) :: context
    real(dp), intent(out) :: value, gradient
    select type(context)
    type is(empty_context)
      value = sin(x)-cos(x)
      gradient = cos(x)+sin(x)
    class default
      value = 0.0_dp
      gradient = 0.0_dp
    end select
  end subroutine trig_callback
  subroutine sin_callback(x, context, value, gradient)
    real(dp), intent(in) :: x
    class(*), intent(in) :: context
    real(dp), intent(out) :: value, gradient
    select type(context)
    type is(empty_context)
      value = sin(x)
      gradient = cos(x)
    class default
      value = 0.0_dp
      gradient = 0.0_dp
    end select
  end subroutine sin_callback
end module root_test_callbacks

program test_roots
  use jrvfinance, only: dp, root_result, newton_raphson_root, bisection_root, JRV_OK
  use root_test_callbacks
  implicit none
  type(empty_context) :: ctx
  type(root_result) :: res
  real(dp), parameter :: pi = acos(-1.0_dp)

  res = newton_raphson_root(trig_callback,ctx,0.0_dp,-2.0_dp,2.0_dp)
  call check(res%status == JRV_OK .and. abs(res%root-pi/4.0_dp) < 1.0e-7_dp, &
    'Newton root')
  res = bisection_root(sin_callback,ctx,7.0_dp,1.0_dp,13.0_dp)
  call check(res%status == JRV_OK .and. abs(sin(res%root)) < 1.0e-5_dp, &
    'bisection root')

  print '(a)', 'test_roots: PASS'
contains
  subroutine check(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write(*,'(a)') 'FAIL: '//trim(label)
      error stop 1
    end if
  end subroutine check
end program test_roots
