program test_r_wild
  use abcoptim, only : dp, abc_control, abc_result, abc_optim
  implicit none

  type(abc_control) :: ctl
  type(abc_result) :: ans
  real(dp) :: par(1), lb(1), ub(1)

  par = 50.0_dp
  lb = -100.0_dp
  ub = 100.0_dp
  ctl%food_number = 20
  ctl%limit = 100
  ctl%max_cycle = 1600
  ctl%criter = 120
  ctl%seed = 213

  call abc_optim(par, wild, lb, ub, ans, ctl)
  call check(abs(ans%par(1) + 15.81515_dp) < 2.0e-3_dp, "wild-function minimizer")
  call check(size(ans%hist, 2) == ans%counts, "R history excludes initial best")
  print *, "test_r_wild: PASS", ans%value, ans%par

contains

  function wild(x) result(y)
    real(dp), intent(in) :: x(:)
    real(dp) :: y
    y = 10.0_dp * sin(0.3_dp*x(1)) * sin(1.3_dp*x(1)**2) + &
      1.0e-5_dp*x(1)**4 + 0.2_dp*x(1) + 80.0_dp
  end function wild

  subroutine check(ok, message)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: message
    if (.not. ok) then
      print *, "FAIL: ", trim(message)
      error stop 1
    end if
  end subroutine check

end program test_r_wild
