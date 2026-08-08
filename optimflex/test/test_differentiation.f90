program test_differentiation
  use optimflex
  implicit none
  real(dp) :: x(2), g(2), h(2,2), j(2,2)
  x = [1.5_dp, -0.5_dp]
  call fast_grad(quad,x,g,diff_central)
  if (maxval(abs(g-[3.0_dp,-3.0_dp])) > 1.0e-7_dp) error stop 'central gradient failed'
  call fast_grad(quad,x,g,diff_richardson)
  if (maxval(abs(g-[3.0_dp,-3.0_dp])) > 1.0e-7_dp) error stop 'Richardson gradient failed'
  call fast_hess(quad,x,h,diff_central)
  if (maxval(abs(h-reshape([2.0_dp,0.0_dp,0.0_dp,6.0_dp],[2,2]))) > 2.0e-5_dp) error stop 'Hessian failed'
  call fast_jac(resid,x,j,diff_central)
  if (maxval(abs(j-reshape([1.0_dp,2.0_dp,1.0_dp,-1.0_dp],[2,2]))) > 1.0e-7_dp) error stop 'Jacobian failed'
  print *, 'test_differentiation: PASS'
contains
  real(dp) function quad(z) result(f)
    real(dp), intent(in) :: z(:)
    f = z(1)*z(1) + 3.0_dp*z(2)*z(2)
  end function quad
  function resid(z) result(r)
    real(dp), intent(in) :: z(:)
    real(dp), allocatable :: r(:)
    allocate(r(2))
    r = [z(1)+z(2), 2.0_dp*z(1)-z(2)]
  end function resid
end program test_differentiation
