program test_spg
  use bb, only: dp, spg, spg_control, spg_result
  implicit none
  type(spg_control) :: ctrl
  type(spg_result) :: fit
  real(dp) :: x0(2)

  x0 = [-1.2_dp, 1.0_dp]
  ctrl = spg_control(method=3, m=10, maxit=5000, gtol=1.0e-7_dp, &
    ftol=1.0e-14_dp, trace=.false.)
  fit = spg(x0, rosenbrock, ctrl, rosen_grad)
  if (.not. fit%succeeded()) error stop 'SPG Rosenbrock did not converge'
  if (maxval(abs(fit%par - [1.0_dp, 1.0_dp])) > 2.0e-4_dp) error stop 'SPG Rosenbrock parameters wrong'
  if (fit%value > 1.0e-8_dp) error stop 'SPG Rosenbrock value wrong'

  ctrl%maximize = .true.
  fit = spg([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], maxfn, ctrl)
  if (abs(fit%value - 10.0_dp) > 1.0e-10_dp) error stop 'SPG maximize exact-start value wrong'
  if (maxval(abs(fit%par - [1.0_dp,2.0_dp,3.0_dp,4.0_dp])) > 1.0e-10_dp) &
    error stop 'SPG maximize exact-start parameters wrong'

  print *, 'test_spg: PASS'
contains
  function rosenbrock(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f
    f = 100.0_dp * (x(2) - x(1)**2)**2 + (1.0_dp - x(1))**2
  end function rosenbrock

  subroutine rosen_grad(x, g)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    g(1) = -400.0_dp*x(1)*(x(2)-x(1)**2) - 2.0_dp*(1.0_dp-x(1))
    g(2) = 200.0_dp*(x(2)-x(1)**2)
  end subroutine rosen_grad

  function maxfn(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f
    integer :: i
    f = 0.0_dp
    do i = 1, size(x)
      f = f + (x(i) - real(i,dp))**2
    end do
    f = 10.0_dp - f*f
  end function maxfn
end program test_spg
