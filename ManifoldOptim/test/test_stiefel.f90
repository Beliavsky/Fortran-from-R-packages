! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of the ManifoldOptim/ROPTLIB computational interface.
! See LICENSE and UPSTREAM_PROVENANCE.md for upstream copyright and provenance.
program test_stiefel
  use manifoldoptim
  implicit none
  integer,parameter::n=5,p=2
  type(manifold_domain)::domain
  type(solver_options)::opt
  type(solver_result)::res
  real(dp)::x0(n*p),avec(n),xmat(n,p),v
  allocate(domain%component(1))
  domain%component(1)=make_component(MANI_STIEFEL,n,p=p)
  xmat=0.0_dp
  xmat(1,1)=1.0_dp
  xmat(2,2)=1.0_dp
  x0=reshape(xmat,[n*p])
  avec=[0.5_dp,-1.0_dp,3.0_dp,2.0_dp,-0.25_dp]
  opt%tolerance=1.0e-7_dp
  opt%max_iteration=400
  call manifold_optimize(domain,x0,obj,grad,'RCG',res,opt)
  if(.not.point_is_valid(domain,res%xopt,2.0e-7_dp))error stop 'stiefel validity'
  xmat=reshape(res%xopt,[n,p])
  v=abs(dot_product(xmat(:,1),avec))/sqrt(dot_product(avec,avec))
  if(v<0.999_dp)then
  write(*,*)v,res%fval
  error stop 'stiefel optimum'
  end if
  write(*,*)'PASS test_stiefel'
contains
  subroutine obj(x,f)
    real(dp),intent(in)::x(:)
    real(dp),intent(out)::f
    real(dp)::xm(n,p)
    xm=reshape(x,[n,p])
    f=-dot_product(avec,xm(:,1))
  end subroutine
  subroutine grad(x,g)
    real(dp),intent(in)::x(:)
    real(dp),intent(out)::g(:)
    real(dp)::gm(n,p)
    if (size(x) /= size(g)) error stop 'grad: size mismatch'
    gm=0.0_dp
    gm(:,1)=-avec
    g=reshape(gm,[n*p])
  end subroutine
end program
