! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of the ManifoldOptim/ROPTLIB computational interface.
! See LICENSE and UPSTREAM_PROVENANCE.md for upstream copyright and provenance.
program test_sphere
  use manifoldoptim
  implicit none
  type(manifold_domain)::domain
  type(solver_options)::opt
  type(solver_result)::res
  real(dp)::x0(5),a(5),truth(5)
  allocate(domain%component(1))
  domain%component(1)=make_component(MANI_SPHERE,5)
  a=[1.0_dp,2.0_dp,-3.0_dp,0.5_dp,4.0_dp]
  truth=a/sqrt(dot_product(a,a))
  x0=[1.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp]
  opt%tolerance=1.0e-8_dp
  opt%max_iteration=300
  call manifold_optimize(domain,x0,obj,grad,'LRBFGS',res,opt)
  if(.not.point_is_valid(domain,res%xopt,1.0e-8_dp))error stop 'sphere validity'
  if(vecnorm(res%xopt-truth)>2.0e-4_dp)then
  write(*,*)res%xopt,truth
  error stop 'sphere optimum'
  end if
  write(*,*)'PASS test_sphere'
contains
  subroutine obj(x,f)
    real(dp),intent(in)::x(:)
    real(dp),intent(out)::f
    f=-dot_product(a,x)
  end subroutine
  subroutine grad(x,g)
    real(dp),intent(in)::x(:)
    real(dp),intent(out)::g(:)
    if (size(x) /= size(g)) error stop 'grad: size mismatch'
    g=-a
  end subroutine
end program
