program test_scaling_and_initialization
  use abcoptim, only : dp, abc_control, abc_result, abc_cpp
  implicit none

  type(abc_control) :: ctl
  type(abc_result) :: ans
  real(dp) :: par(2), lb(2), ub(2), expected(2,4)
  real(dp) :: fnscale
  integer :: i

  par = 0.0_dp
  lb = [-2.0_dp, 10.0_dp]
  ub = [4.0_dp, 16.0_dp]
  ctl%food_number = 4
  ctl%limit = 100
  ctl%max_cycle = 5
  ctl%criter = 0
  ctl%seed = 42

  do i = 1, 4
    expected(:,i) = lb + (ub - lb) * real(i - 1, dp) / 3.0_dp
  end do
  call abc_cpp(par, constant, lb, ub, ans, ctl)
  call check(maxval(abs(ans%foods - expected)) < 1.0e-14_dp, "equally-spaced initial foods")

  par = 0.0_dp
  lb = -5.0_dp
  ub = 5.0_dp
  ctl%food_number = 25
  ctl%max_cycle = 800
  ctl%criter = 100
  ctl%seed = 777
  fnscale = -1.0_dp
  call abc_cpp(par, gaussian_peak, lb, ub, ans, ctl, fnscale=fnscale)
  call check(gaussian_peak(ans%par) > 0.9999_dp, "fnscale maximization")
  call check(abs(ans%value + gaussian_peak(ans%par)) < 1.0e-12_dp, "scaled return value")
  print *, "test_scaling_and_initialization: PASS"

contains

  function constant(x) result(y)
    real(dp), intent(in) :: x(:)
    real(dp) :: y
    y = 0.0_dp * sum(x)
  end function constant

  function gaussian_peak(x) result(y)
    real(dp), intent(in) :: x(:)
    real(dp) :: y
    y = exp(-sum((x - 1.25_dp)**2))
  end function gaussian_peak

  subroutine check(ok, message)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: message
    if (.not. ok) then
      print *, "FAIL: ", trim(message)
      error stop 1
    end if
  end subroutine check

end program test_scaling_and_initialization
