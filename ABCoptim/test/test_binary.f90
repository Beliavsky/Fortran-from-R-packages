program test_binary
  use abcoptim, only : dp, abc_control, abc_result, abc_optim
  implicit none

  type(abc_control) :: ctl
  type(abc_result) :: ans
  real(dp) :: par(6), lb(1), ub(1)
  real(dp), parameter :: target(6) = [1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp]

  par = 0.0_dp
  lb = 0.0_dp
  ub = 1.0_dp
  ctl%optiinteger = .true.
  ctl%food_number = 24
  ctl%limit = 30
  ctl%max_cycle = 600
  ctl%criter = 80
  ctl%seed = 9183

  call abc_optim(par, hamming, lb, ub, ans, ctl)
  call check(ans%value < 1.0e-15_dp, "binary optimum")
  call check(maxval(abs(ans%par - target)) < 1.0e-15_dp, "binary parameters")
  call check(all((abs(ans%par) < 1.0e-15_dp) .or. &
    (abs(ans%par - 1.0_dp) < 1.0e-15_dp)), "binary result values")
  print *, "test_binary: PASS", ans%par

contains

  function hamming(x) result(y)
    real(dp), intent(in) :: x(:)
    real(dp) :: y
    y = sum(abs(x - target))
  end function hamming

  subroutine check(ok, message)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: message
    if (.not. ok) then
      print *, "FAIL: ", trim(message)
      error stop 1
    end if
  end subroutine check

end program test_binary
