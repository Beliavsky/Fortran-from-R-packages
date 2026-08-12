program test_cpp_continuous
  use abcoptim, only : dp, abc_control, abc_result, abc_cpp
  implicit none

  type(abc_control) :: ctl
  type(abc_result) :: ans
  real(dp) :: par(2), lb(1), ub(1)
  real(dp), parameter :: pi = acos(-1.0_dp)

  par = 0.0_dp
  lb = -10.0_dp
  ub = 10.0_dp
  ctl%food_number = 30
  ctl%limit = 80
  ctl%max_cycle = 1200
  ctl%criter = 120
  ctl%seed = 213

  call abc_cpp(par, peaks, lb, ub, ans, ctl)
  call check(ans%value < -0.9999_dp, "peak objective")
  call check(maxval(abs(ans%par - pi)) < 2.0e-2_dp, "peak parameters")
  call check(size(ans%hist, 2) == ans%counts + 1, "C++ history includes initial best")
  print *, "test_cpp_continuous: PASS", ans%value, ans%par

contains

  function peaks(x) result(y)
    real(dp), intent(in) :: x(:)
    real(dp) :: y
    y = -cos(x(1)) * cos(x(2)) * exp(-((x(1) - pi)**2 + (x(2) - pi)**2))
  end function peaks

  subroutine check(ok, message)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: message
    if (.not. ok) then
      print *, "FAIL: ", trim(message)
      error stop 1
    end if
  end subroutine check

end program test_cpp_continuous
