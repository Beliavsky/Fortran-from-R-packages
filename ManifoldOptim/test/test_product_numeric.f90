! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of the ManifoldOptim/ROPTLIB computational interface.
! See LICENSE and UPSTREAM_PROVENANCE.md for upstream copyright and provenance.
program test_product_numeric
  use manifoldoptim
  implicit none
  type(manifold_domain)::domain
  type(solver_options)::opt
  type(solver_result)::res
  real(dp)::x0(5),a(3),t(2)
  allocate(domain%component(2))
  domain%component(1)=make_component(MANI_SPHERE,3)
  domain%component(2)=make_component(MANI_EUCLIDEAN,2,m=1)
  a=[1.0_dp,-2.0_dp,0.5_dp]
  t=[3.0_dp,-1.0_dp]
  x0=[1.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp]
  opt%tolerance=2.0e-6_dp
  opt%max_iteration=400
  opt%eps_numerical_grad=1.0e-6_dp
  call manifold_optimize(domain,x0,obj,'LRBFGS',res,opt)
  if(.not.point_is_valid(domain,res%xopt,1.0e-7_dp))error stop 'product validity'
  if(vecnorm(res%xopt(1:3)-a/vecnorm(a))>2.0e-3_dp.or.maxval(abs(res%xopt(4:5)-t))>2.0e-3_dp)then
    write(*,*)res%xopt
    error stop 'numeric gradient product optimum'
  end if
  write(*,*)'PASS test_product_numeric'
contains
  subroutine obj(x,f)
    real(dp),intent(in)::x(:)
    real(dp),intent(out)::f
    f=-dot_product(a,x(1:3))+0.5_dp*sum((x(4:5)-t)**2)
  end subroutine
end program
