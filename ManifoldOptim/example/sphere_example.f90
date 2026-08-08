! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of the ManifoldOptim/ROPTLIB computational interface.
! See LICENSE and UPSTREAM_PROVENANCE.md for upstream copyright and provenance.
program sphere_example
  use manifoldoptim
  implicit none

  type(manifold_domain) :: domain
  type(solver_options) :: opt
  type(solver_result) :: res
  real(dp) :: x0(5), a(5)

  allocate(domain%component(1))
  domain%component(1) = make_component(MANI_SPHERE, 5)

  a = [1.0_dp, 2.0_dp, -3.0_dp, 0.5_dp, 4.0_dp]
  x0 = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
  opt%tolerance = 1.0e-8_dp

  call manifold_optimize(domain, x0, objective, gradient, 'LRBFGS', res, opt)

  write(*,'(a,1x,es14.6)') 'f(x*) =', res%fval
  write(*,'(a,*(1x,f10.6))') 'x* =', res%xopt
  write(*,'(a,1x,es14.6)') '||grad f|| =', res%normgf

contains

  subroutine objective(x, f)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f
    f = -dot_product(a, x)
  end subroutine objective

  subroutine gradient(x, g)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    if (size(x) /= size(g)) error stop 'gradient: size mismatch'
    g = -a
  end subroutine gradient

end program sphere_example
