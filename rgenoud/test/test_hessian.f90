program test_hessian
  use rgenoud, only : dp, genoud_options, genoud_result, genoud_optimize
  implicit none
  type(genoud_options) :: opt
  type(genoud_result) :: res
  real(dp) :: lower(2), upper(2)

  lower = -5.0_dp
  upper = 5.0_dp
  opt%pop_size = 30
  opt%max_generations = 30
  opt%wait_generations = 5
  opt%seed = 17
  call genoud_optimize(quadratic, lower, upper, opt, res, gradient=quad_grad, hessian=.true.)
  if (.not. allocated(res%hessian)) error stop "hessian missing"
  if (maxval(abs(res%hessian - reshape([2.0_dp, 0.0_dp, 0.0_dp, 4.0_dp], [2, 2]))) > 1.0e-4_dp) then
    error stop "hessian test failed"
  end if
contains
  real(dp) function quadratic(x) result(f)
    real(dp), intent(in) :: x(:)
    f = x(1)**2 + 2.0_dp * x(2)**2
  end function quadratic
  subroutine quad_grad(x, g)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    g = [2.0_dp * x(1), 4.0_dp * x(2)]
  end subroutine quad_grad
end program test_hessian
