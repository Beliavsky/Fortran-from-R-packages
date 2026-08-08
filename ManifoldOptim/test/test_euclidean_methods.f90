! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of the ManifoldOptim/ROPTLIB computational interface.
! See LICENSE and UPSTREAM_PROVENANCE.md for upstream copyright and provenance.
program test_euclidean_methods
  use manifoldoptim
  implicit none
  type(manifold_domain) :: domain
  type(solver_options) :: opt
  type(solver_result) :: res
  real(dp) :: x0(4), target(4)
  character(len=16), parameter :: methods(11) = [character(len=16) :: &
    'RSD','RCG','RNewton','RBFGS','LRBFGS','RBroydenFamily','RWRBFGS', &
    'RTRSD','RTRNewton','RTRSR1','LRTRSR1']
  integer :: i
  allocate(domain%component(1))
  domain%component(1)=make_component(MANI_EUCLIDEAN,4,m=1)
  x0=[4.0_dp,-3.0_dp,2.0_dp,5.0_dp]
  target=[1.0_dp,2.0_dp,-1.0_dp,0.5_dp]
  opt%tolerance=1.0e-7_dp
  opt%max_iteration=300
  opt%trust_radius=1.0_dp
  do i=1,size(methods)
    call manifold_optimize(domain,x0,obj,grad,hess,trim(methods(i)),res,opt)
    if(maxval(abs(res%xopt-target))>3.0e-4_dp)then
      write(*,*)'FAIL ',trim(methods(i)),res%status,res%fval,res%normgf,res%xopt
      error stop 1
    end if
  end do
  write(*,*)'PASS test_euclidean_methods'
contains
  subroutine obj(x,f)
    real(dp),intent(in)::x(:)
    real(dp),intent(out)::f
    f=0.5_dp*sum((x-target)**2)
  end subroutine
  subroutine grad(x,g)
    real(dp),intent(in)::x(:)
    real(dp),intent(out)::g(:)
    g=x-target
  end subroutine
  subroutine hess(x,v,hv)
    real(dp),intent(in)::x(:),v(:)
    real(dp),intent(out)::hv(:)
    if (size(x) /= size(v)) error stop 'hess: size mismatch'
    hv=v
  end subroutine
end program
