! SPDX-License-Identifier: GPL-3.0-or-later
module cccp_solver
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use cccp_kinds, only : dp
   use cccp_types
   use cccp_linalg, only : solve_system, equality_particular, vector_norm2, symmetrize
   use cccp_cones, only : cone_barrier, cones_interior, phase_shift_needed, barrier_parameter, &
      cone_slacks, cone_duals, cone_dual_gradient, minimum_cone_slack, nnoc
   implicit none
   private

   integer, parameter :: objective_linear = 1
   integer, parameter :: objective_quadratic = 2
   integer, parameter :: objective_user = 3

   public :: ctrl, solve_lp, solve_qp, solve_dnl, solve_dcp, dnl, dcp
   public :: cccp_lp, cccp_qp, cccp, cccp_solve
   public :: getx, gety, gets, getz, getstate, getstatus, getniter, getparams

   procedure(constraints_callback), pointer, save :: phase_nl_target => null()
   integer, save :: phase_nl_n = 0, phase_nl_m = 0

   interface cccp
      module procedure cccp_lp
      module procedure cccp_qp
   end interface cccp
   interface cccp_solve
      module procedure cccp_lp
      module procedure cccp_qp
   end interface cccp_solve
   interface dnl
      module procedure solve_dnl
   end interface dnl
   interface dcp
      module procedure solve_dcp
   end interface dcp

contains

   function ctrl(maxiters, abstol, reltol, feastol, stepadj, beta, trace) result(c)
      integer, intent(in), optional :: maxiters
      real(dp), intent(in), optional :: abstol, reltol, feastol, stepadj, beta
      logical, intent(in), optional :: trace
      type(cccp_control) :: c
      if (present(maxiters)) c%maxiters = maxiters
      if (present(abstol)) c%abstol = abstol
      if (present(reltol)) c%reltol = reltol
      if (present(feastol)) c%feastol = feastol
      if (present(stepadj)) c%stepadj = stepadj
      if (present(beta)) c%beta = beta
      if (present(trace)) c%trace = trace
   end function ctrl

   function cccp_lp(q, a, b, cones, control) result(sol)
      real(dp), intent(in) :: q(:)
      real(dp), intent(in), optional :: a(:,:), b(:)
      type(cone_constraint), intent(in), optional :: cones(:)
      type(cccp_control), intent(in), optional :: control
      type(cccp_solution) :: sol
      call solve_lp(q, a, b, cones, control, sol)
   end function cccp_lp

   function cccp_qp(p, q, a, b, cones, control) result(sol)
      real(dp), intent(in) :: p(:,:), q(:)
      real(dp), intent(in), optional :: a(:,:), b(:)
      type(cone_constraint), intent(in), optional :: cones(:)
      type(cccp_control), intent(in), optional :: control
      type(cccp_solution) :: sol
      call solve_qp(p, q, a, b, cones, control, sol)
   end function cccp_qp

   subroutine solve_lp(q, a, b, cones, control, sol)
      real(dp), intent(in) :: q(:)
      real(dp), intent(in), optional :: a(:,:), b(:)
      type(cone_constraint), intent(in), optional :: cones(:)
      type(cccp_control), intent(in), optional :: control
      type(cccp_solution), intent(out) :: sol
      type(objective_spec) :: obj
      type(cccp_control) :: ctl
      type(cone_constraint), allocatable :: cc(:)
      real(dp), allocatable :: aa(:,:), bb(:)
      integer :: n

      n = size(q)
      obj%mode = objective_linear
      obj%n = n
      allocate(obj%q(n)); obj%q = q
      call normalize_inputs(n, a, b, cones, control, aa, bb, cc, ctl, sol%info)
      if (sol%info /= cccp_success) then
         sol%status = 'invalid input'
         return
      end if
      call solve_core(obj, aa, bb, cc, ctl, sol=sol)
   end subroutine solve_lp

   subroutine solve_qp(p, q, a, b, cones, control, sol)
      real(dp), intent(in) :: p(:,:), q(:)
      real(dp), intent(in), optional :: a(:,:), b(:)
      type(cone_constraint), intent(in), optional :: cones(:)
      type(cccp_control), intent(in), optional :: control
      type(cccp_solution), intent(out) :: sol
      type(objective_spec) :: obj
      type(cccp_control) :: ctl
      type(cone_constraint), allocatable :: cc(:)
      real(dp), allocatable :: aa(:,:), bb(:)
      integer :: n

      n = size(q)
      if (any(shape(p) /= [n,n])) then
         sol%info = cccp_invalid_input
         sol%status = 'invalid input'
         return
      end if
      obj%mode = objective_quadratic
      obj%n = n
      allocate(obj%q(n), obj%p(n,n)); obj%q = q; obj%p = symmetrize(p)
      call normalize_inputs(n, a, b, cones, control, aa, bb, cc, ctl, sol%info)
      if (sol%info /= cccp_success) then
         sol%status = 'invalid input'
         return
      end if
      call solve_core(obj, aa, bb, cc, ctl, sol=sol)
   end subroutine solve_qp

   subroutine solve_dnl(q, x0, nl_eval, mnl, a, b, cones, control, sol)
      real(dp), intent(in) :: q(:), x0(:)
      procedure(constraints_callback) :: nl_eval
      integer, intent(in) :: mnl
      real(dp), intent(in), optional :: a(:,:), b(:)
      type(cone_constraint), intent(in), optional :: cones(:)
      type(cccp_control), intent(in), optional :: control
      type(cccp_solution), intent(out) :: sol
      type(objective_spec) :: obj
      type(nonlinear_spec) :: nl
      type(cccp_control) :: ctl
      type(cone_constraint), allocatable :: cc(:)
      real(dp), allocatable :: aa(:,:), bb(:)
      integer :: n

      n = size(q)
      if (size(x0) /= n .or. mnl < 1) then
         sol%info = cccp_invalid_input; sol%status = 'invalid input'; return
      end if
      obj%mode = objective_linear; obj%n = n
      allocate(obj%q(n)); obj%q = q
      nl%m = mnl; nl%eval => nl_eval
      call normalize_inputs(n,a,b,cones,control,aa,bb,cc,ctl,sol%info)
      if(sol%info/=cccp_success) then; sol%status='invalid input'; return; end if
      call solve_core(obj,aa,bb,cc,ctl,x0=x0,nl=nl,sol=sol)
   end subroutine solve_dnl

   subroutine solve_dcp(x0, obj_eval, a, b, cones, control, sol, nl_eval, mnl)
      real(dp), intent(in) :: x0(:)
      procedure(objective_callback) :: obj_eval
      real(dp), intent(in), optional :: a(:,:), b(:)
      type(cone_constraint), intent(in), optional :: cones(:)
      type(cccp_control), intent(in), optional :: control
      type(cccp_solution), intent(out) :: sol
      procedure(constraints_callback), optional :: nl_eval
      integer, intent(in), optional :: mnl
      type(objective_spec) :: obj
      type(nonlinear_spec) :: nl
      type(cccp_control) :: ctl
      type(cone_constraint), allocatable :: cc(:)
      real(dp), allocatable :: aa(:,:), bb(:)
      integer :: n

      n=size(x0); obj%mode=objective_user; obj%n=n; obj%eval=>obj_eval
      call normalize_inputs(n,a,b,cones,control,aa,bb,cc,ctl,sol%info)
      if(sol%info/=cccp_success) then; sol%status='invalid input'; return; end if
      if(present(nl_eval)) then
         if(.not.present(mnl).or.mnl<1) then; sol%info=cccp_invalid_input; sol%status='invalid input'; return; end if
         nl%m=mnl; nl%eval=>nl_eval
         call solve_core(obj,aa,bb,cc,ctl,x0=x0,nl=nl,sol=sol)
      else
         call solve_core(obj,aa,bb,cc,ctl,x0=x0,sol=sol)
      end if
   end subroutine solve_dcp

   subroutine normalize_inputs(n,a,b,cones,control,aa,bb,cc,ctl,info)
      integer,intent(in)::n
      real(dp),intent(in),optional::a(:,:),b(:)
      type(cone_constraint),intent(in),optional::cones(:)
      type(cccp_control),intent(in),optional::control
      real(dp),allocatable,intent(out)::aa(:,:),bb(:)
      type(cone_constraint),allocatable,intent(out)::cc(:)
      type(cccp_control),intent(out)::ctl
      integer,intent(out)::info
      integer::k
      info=cccp_success; ctl=cccp_control(); if(present(control))ctl=control
      if(present(a))then
         if(size(a,2)/=n .or. .not.present(b) .or. size(b)/=size(a,1))then; info=cccp_invalid_input; return; end if
         allocate(aa(size(a,1),n),bb(size(a,1))); aa=a; bb=b
      else
         if(present(b))then; info=cccp_invalid_input; return; end if
         allocate(aa(0,n),bb(0))
      end if
      if(present(cones))then
         allocate(cc(size(cones))); cc=cones
         do k=1,size(cc)
            if(size(cc(k)%g,2)/=n .or. size(cc(k)%g,1)/=size(cc(k)%h))then; info=cccp_invalid_input; return; end if
            if(cc(k)%kind==cone_psdc .and. size(cc(k)%h)/=cc(k)%dim*cc(k)%dim)then; info=cccp_invalid_input; return; end if
         end do
      else
         allocate(cc(0))
      end if
      if(ctl%maxiters<1 .or. ctl%abstol<=0.0_dp .or. ctl%feastol<=0.0_dp .or. &
         ctl%stepadj<=0.0_dp .or. ctl%stepadj>1.0_dp .or. ctl%beta<=0.0_dp .or. ctl%beta>=1.0_dp)then
         info=cccp_invalid_input
      end if
   end subroutine normalize_inputs

   subroutine evaluate_objective(obj,x,f,g,h,info)
      type(objective_spec),intent(in)::obj
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::f,g(:),h(:,:)
      integer,intent(out)::info
      select case(obj%mode)
      case(objective_linear)
         f=dot_product(obj%q,x); g=obj%q; h=0.0_dp; info=0
      case(objective_quadratic)
         g=matmul(obj%p,x)+obj%q; h=obj%p
         f=0.5_dp*dot_product(x,matmul(obj%p,x))+dot_product(obj%q,x); info=0
      case(objective_user)
         if(.not.associated(obj%eval))then; info=1; return; end if
         call obj%eval(x,f,g,h,info)
      case default
         info=1
      end select
   end subroutine evaluate_objective

   subroutine nonlinear_barrier(nl,x,phi,gphi,hphi,ok)
      type(nonlinear_spec),intent(in)::nl
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::phi,gphi(:),hphi(:,:)
      logical,intent(out)::ok
      real(dp),allocatable::f(:),g(:,:),h(:,:,:)
      integer::i,info,n
      n=size(x); phi=0.0_dp; gphi=0.0_dp; hphi=0.0_dp; ok=.false.
      if(nl%m<=0 .or. .not.associated(nl%eval))then; ok=.true.; return; end if
      allocate(f(nl%m),g(nl%m,n),h(n,n,nl%m))
      call nl%eval(x,f,g,h,info)
      if(info/=0 .or. any(.not.ieee_is_finite(f)) .or. any(f>=0.0_dp))return
      do i=1,nl%m
         phi=phi-log(-f(i))
         gphi=gphi-g(i,:)/f(i)
         hphi=hphi-h(:,:,i)/f(i)+matmul(reshape(g(i,:),[n,1]),reshape(g(i,:),[1,n]))/(f(i)*f(i))
      end do
      ok=.true.
   end subroutine nonlinear_barrier

   logical function nonlinear_interior(nl,x) result(ok)
      type(nonlinear_spec),intent(in)::nl
      real(dp),intent(in)::x(:)
      real(dp),allocatable::f(:),g(:,:),h(:,:,:)
      integer::info,n
      n=size(x); allocate(f(nl%m),g(nl%m,n),h(n,n,nl%m))
      call nl%eval(x,f,g,h,info)
      ok=info==0 .and. all(ieee_is_finite(f)) .and. all(f<0.0_dp)
   end function nonlinear_interior

   subroutine solve_core(obj,a,b,cones,ctl,x0,nl,sol)
      type(objective_spec),intent(in)::obj
      real(dp),intent(in)::a(:,:),b(:)
      type(cone_constraint),intent(in)::cones(:)
      type(cccp_control),intent(in)::ctl
      real(dp),intent(in),optional::x0(:)
      type(nonlinear_spec),intent(in),optional::nl
      type(cccp_solution),intent(out)::sol
      real(dp),allocatable::x(:),xeq(:)
      integer::info,n
      logical::ok,oknl
      n=obj%n; allocate(x(n),xeq(n))
      if(present(x0))then
         x=x0
         if(size(x0)/=n)then; sol%info=cccp_invalid_input; sol%status='invalid input'; return; end if
         if(size(a,1)>0 .and. vector_norm2(matmul(a,x)-b)>10.0_dp*ctl%feastol)then
            sol%info=cccp_infeasible_start; sol%status='infeasible start'; return
         end if
         ok=cones_interior(cones,x)
         if(present(nl))then
            oknl=nonlinear_interior(nl,x)
            if((.not.ok).or.(.not.oknl))then
               call find_combined_interior(a,b,cones,nl,ctl,x,info)
               if(info/=0)then; sol%info=cccp_infeasible_start; sol%status='infeasible'; return; end if
            end if
         else if(.not.ok)then
            if(size(cones)>0)then
               call find_cone_interior(a,b,cones,ctl,x,info)
               if(info/=0)then; sol%info=cccp_infeasible_start; sol%status='infeasible'; return; end if
            end if
         end if
      else
         call equality_particular(a,b,xeq,info)
         if(info/=0)then; sol%info=cccp_singular_system; sol%status='singular equality'; return; end if
         x=xeq
         if(size(cones)>0)then
            ok=cones_interior(cones,x)
            if(.not.ok)then
               call find_cone_interior(a,b,cones,ctl,x,info)
               if(info/=0)then; sol%info=cccp_infeasible_start; sol%status='infeasible'; return; end if
            end if
         end if
         if(present(nl))then; sol%info=cccp_infeasible_start; sol%status='x0 required'; return; end if
      end if
      if(present(nl))then
         call barrier_optimize(obj,a,b,cones,ctl,x,sol, nl=nl)
      else
         call barrier_optimize(obj,a,b,cones,ctl,x,sol)
      end if
   end subroutine solve_core

   subroutine find_combined_interior(a,b,cones,nl,ctl,x,info)
      real(dp),intent(in)::a(:,:),b(:)
      type(cone_constraint),intent(in)::cones(:)
      type(nonlinear_spec),intent(in)::nl
      type(cccp_control),intent(in)::ctl
      real(dp),intent(inout)::x(:)
      integer,intent(out)::info
      type(objective_spec)::phase_obj
      type(nonlinear_spec)::phase_nl
      type(cccp_solution),allocatable::phase_sol
      type(cccp_control)::pctl
      real(dp),allocatable::y(:),aa(:,:),bb(:),f(:),g(:,:),h(:,:,:)
      integer::n,p,eval_info
      real(dp),parameter::delta=1.0e-6_dp
      allocate(phase_sol)
      n=size(x);p=size(a,1)
      allocate(y(n+1),aa(p,n+1),bb(p),f(nl%m),g(nl%m,n),h(n,n,nl%m))
      call nl%eval(x,f,g,h,eval_info)
      if(eval_info/=0 .or. any(.not.ieee_is_finite(f)))then;info=1;return;end if
      y(1:n)=x
      y(n+1)=max(phase_shift_needed(cones,x),maxval(f)+1.0_dp)
      aa=0.0_dp;if(p>0)aa(:,1:n)=a;bb=b
      phase_obj%mode=objective_linear;phase_obj%n=n+1
      allocate(phase_obj%q(n+1));phase_obj%q=0.0_dp;phase_obj%q(n+1)=1.0_dp
      phase_nl_target=>nl%eval;phase_nl_n=n;phase_nl_m=nl%m
      phase_nl%m=nl%m;phase_nl%eval=>phase_nonlinear_callback
      pctl=ctl;pctl%abstol=min(ctl%abstol,1.0e-9_dp);pctl%reltol=min(ctl%reltol,1.0e-9_dp)
      call barrier_optimize(phase_obj,aa,bb,cones,pctl,y,phase_sol,nl=phase_nl,phase=.true.,phase_delta=delta)
      nullify(phase_nl_target);phase_nl_n=0;phase_nl_m=0
      if((phase_sol%info==cccp_success.or.phase_sol%info==cccp_max_iterations).and.allocated(phase_sol%x))then
         if(phase_sol%x(n+1)<-0.05_dp*delta)then;x=phase_sol%x(1:n);info=0;return;end if
      end if
      info=1
   end subroutine find_combined_interior

   subroutine phase_nonlinear_callback(y,f,g,h,info)
      real(dp),intent(in)::y(:)
      real(dp),intent(out)::f(:),g(:,:),h(:,:,:)
      integer,intent(out)::info
      real(dp),allocatable::f0(:),g0(:,:),h0(:,:,:)
      if(.not.associated(phase_nl_target).or.size(y)/=phase_nl_n+1)then;info=1;return;end if
      allocate(f0(phase_nl_m),g0(phase_nl_m,phase_nl_n),h0(phase_nl_n,phase_nl_n,phase_nl_m))
      call phase_nl_target(y(1:phase_nl_n),f0,g0,h0,info)
      if(info/=0)return
      f=f0-y(size(y));g=0.0_dp;g(:,1:phase_nl_n)=g0;g(:,size(y))=-1.0_dp
      h=0.0_dp;h(1:phase_nl_n,1:phase_nl_n,:)=h0
   end subroutine phase_nonlinear_callback

   subroutine find_cone_interior(a,b,cones,ctl,x,info)
      real(dp),intent(in)::a(:,:),b(:)
      type(cone_constraint),intent(in)::cones(:)
      type(cccp_control),intent(in)::ctl
      real(dp),intent(inout)::x(:)
      integer,intent(out)::info
      type(objective_spec)::phase_obj
      type(cccp_solution),allocatable::phase_sol
      type(cccp_control)::pctl
      real(dp),allocatable::y(:),aa(:,:),bb(:)
      integer::n,p
      real(dp),parameter::delta=1.0e-6_dp
      allocate(phase_sol)
      n=size(x);p=size(a,1)
      allocate(y(n+1),aa(p,n+1),bb(p))
      y(1:n)=x; y(n+1)=phase_shift_needed(cones,x)
      aa=0.0_dp; if(p>0)aa(:,1:n)=a; bb=b
      phase_obj%mode=objective_linear; phase_obj%n=n+1
      allocate(phase_obj%q(n+1)); phase_obj%q=0.0_dp; phase_obj%q(n+1)=1.0_dp
      pctl=ctl; pctl%abstol=min(ctl%abstol,1.0e-9_dp); pctl%reltol=min(ctl%reltol,1.0e-9_dp)
      call barrier_optimize(phase_obj,aa,bb,cones,pctl,y,phase_sol,phase=.true.,phase_delta=delta)
      if((phase_sol%info==cccp_success .or. phase_sol%info==cccp_max_iterations) .and. &
         allocated(phase_sol%x) .and. phase_sol%x(n+1)<-0.05_dp*delta)then
         x=phase_sol%x(1:n); info=0
      else
         info=1
      end if
   end subroutine find_cone_interior

   subroutine total_value(obj,cones,x,tbar,value,grad,hess,ok,nl,phase,phase_delta)
      type(objective_spec),intent(in)::obj
      type(cone_constraint),intent(in)::cones(:)
      real(dp),intent(in)::x(:),tbar
      real(dp),intent(out)::value,grad(:),hess(:,:)
      logical,intent(out)::ok
      type(nonlinear_spec),intent(in),optional::nl
      logical,intent(in),optional::phase
      real(dp),intent(in),optional::phase_delta
      real(dp)::f,phi,phinl,delta
      real(dp),allocatable::g(:),hh(:,:),gb(:),hb(:,:),gn(:),hn(:,:)
      integer::n,info
      logical::okc,okn,ph
      n=size(x); allocate(g(n),hh(n,n),gb(n),hb(n,n),gn(n),hn(n,n))
      call evaluate_objective(obj,x,f,g,hh,info)
      if(info/=0 .or. .not.ieee_is_finite(f))then; ok=.false.; return; end if
      ph=.false.; if(present(phase))ph=phase
      if(ph)then
         call cone_barrier(cones,x,phi,gb,hb,okc,phase=.true.,tau=x(n))
         delta=phase_delta+x(n)
         if(delta<=0.0_dp)okc=.false.
         if(okc)then
            phi=phi-log(delta); gb(n)=gb(n)-1.0_dp/delta; hb(n,n)=hb(n,n)+1.0_dp/(delta*delta)
         end if
      else
         call cone_barrier(cones,x,phi,gb,hb,okc)
      end if
      if(.not.okc)then; ok=.false.; return; end if
      phinl=0.0_dp;gn=0.0_dp;hn=0.0_dp;okn=.true.
      if(present(nl))call nonlinear_barrier(nl,x,phinl,gn,hn,okn)
      if(.not.okn)then; ok=.false.; return; end if
      value=tbar*f+phi+phinl
      grad=tbar*g+gb+gn
      hess=tbar*hh+hb+hn
      ok=ieee_is_finite(value).and.all(ieee_is_finite(grad)).and.all(ieee_is_finite(hess))
   end subroutine total_value

   subroutine barrier_optimize(obj,a,b,cones,ctl,x,sol,nl,phase,phase_delta)
      type(objective_spec),intent(in)::obj
      real(dp),intent(in)::a(:,:),b(:)
      type(cone_constraint),intent(in)::cones(:)
      type(cccp_control),intent(in)::ctl
      real(dp),intent(inout)::x(:)
      type(cccp_solution),intent(out)::sol
      type(nonlinear_spec),intent(in),optional::nl
      logical,intent(in),optional::phase
      real(dp),intent(in),optional::phase_delta
      real(dp),allocatable::grad(:),hess(:,:),kkt(:,:),rhs(:),answer(:),stepv(:),trial(:),lambda(:),gobj(:),hobj(:,:),dualgrad(:)
      real(dp)::value,newvalue,tbar,alpha,gd,dec,fobj,eqres,nu_gap
      integer::n,p,outer,inner,info,total_iter,nu,iobj
      logical::ok,ok2,ph
      n=size(x);p=size(a,1); allocate(grad(n),hess(n,n),stepv(n),trial(n),gobj(n),hobj(n,n),dualgrad(n))
      allocate(kkt(n+p,n+p),rhs(n+p),answer(n+p),lambda(p))
      ph=.false.;if(present(phase))ph=phase
      tbar=1.0_dp;total_iter=0;nu=barrier_parameter(cones)
      if(present(nl))nu=nu+nl%m
      if(ph)nu=nu+1
      do outer=1,ctl%max_outer
         do inner=1,ctl%maxiters
            total_iter=total_iter+1
            if(present(nl))then
               call total_value(obj,cones,x,tbar,value,grad,hess,ok,nl=nl,phase=ph,phase_delta=phase_delta)
            else
               call total_value(obj,cones,x,tbar,value,grad,hess,ok,phase=ph,phase_delta=phase_delta)
            end if
            if(.not.ok)then; sol%info=cccp_domain_error; sol%status='domain error'; return; end if
            do iobj=1,n
               hess(iobj,iobj)=hess(iobj,iobj)+1.0e-10_dp*max(1.0_dp,abs(hess(iobj,iobj)))
            end do
            kkt=0.0_dp;kkt(1:n,1:n)=hess;rhs(1:n)=-grad
            if(p>0)then
               kkt(1:n,n+1:n+p)=transpose(a);kkt(n+1:n+p,1:n)=a;rhs(n+1:n+p)=0.0_dp
            end if
            call solve_system(kkt,rhs,answer,info)
            if(info/=0)then; sol%info=cccp_singular_system; sol%status='singular KKT'; return; end if
            stepv=answer(1:n); if(p>0)lambda=answer(n+1:n+p)
            gd=dot_product(grad,stepv);dec=-gd
            if(dec/2.0_dp<=ctl%abstol)exit
            alpha=1.0_dp
            do
               trial=x+alpha*stepv
               if(present(nl))then
                  call total_value(obj,cones,trial,tbar,newvalue,gobj,hobj,ok2,nl=nl,phase=ph,phase_delta=phase_delta)
               else
                  call total_value(obj,cones,trial,tbar,newvalue,gobj,hobj,ok2,phase=ph,phase_delta=phase_delta)
               end if
               if(ok2)then
                  if(newvalue<=value+0.01_dp*alpha*gd)exit
               end if
               alpha=alpha*ctl%beta
               if(alpha<1.0e-14_dp)then; sol%info=cccp_domain_error; sol%status='line search failed'; return; end if
            end do
            x=trial
            if(ctl%trace)write(*,'(a,i0,a,i0,a,es12.4,a,es12.4)')'outer=',outer,' inner=',inner,' value=',newvalue,' decrement=',dec
         end do
         if(nu==0)exit
         nu_gap=real(nu,dp)/tbar
         if(nu_gap<=max(ctl%abstol,ctl%reltol))exit
         tbar=tbar*ctl%barrier_growth
      end do
      call evaluate_objective(obj,x,fobj,gobj,hobj,info)
      allocate(sol%x(n));sol%x=x
      allocate(sol%y(p)); if(p>0)then;sol%y=lambda/tbar;else;sol%y=0.0_dp;end if
      if(ph)then
         allocate(sol%s(0),sol%z(0),sol%cone_offsets(0,2))
      else
         call cone_slacks(cones,x,sol%s,sol%cone_offsets)
         call cone_duals(cones,x,tbar,sol%z)
      end if
      eqres=0.0_dp;if(p>0)eqres=vector_norm2(matmul(a,x)-b)
      sol%state%pobj=fobj;sol%state%dobj=fobj;sol%state%dgap=real(nu,dp)/tbar
      sol%state%rdgap=sol%state%dgap/max(1.0_dp,abs(fobj));sol%state%certp=eqres
      dualgrad=gobj
      if(p>0)dualgrad=dualgrad+matmul(transpose(a),sol%y)
      if((.not.ph).and.size(sol%z)>0)then
         call cone_dual_gradient(cones,sol%z,grad)
         dualgrad=dualgrad+grad
      end if
      sol%state%certd=vector_norm2(dualgrad)
      if(ph)then
         sol%state%pslack=minimum_cone_slack(cones,x(1:n-1))
      else
         sol%state%pslack=minimum_cone_slack(cones,x)
      end if
      if(size(sol%z)>0)sol%state%dslack=minval(sol%z)
      sol%niter=total_iter
      if(nu==0 .or. real(nu,dp)/tbar<=10.0_dp*max(ctl%abstol,ctl%reltol))then
         sol%info=cccp_success;sol%status='optimal'
      else
         sol%info=cccp_max_iterations;sol%status='unknown'
      end if
   end subroutine barrier_optimize

   function getx(sol) result(x)
      type(cccp_solution),intent(in)::sol
      real(dp),allocatable::x(:)
      x=sol%x
   end function getx
   function gety(sol) result(y)
      type(cccp_solution),intent(in)::sol
      real(dp),allocatable::y(:)
      y=sol%y
   end function gety
   function gets(sol) result(s)
      type(cccp_solution),intent(in)::sol
      real(dp),allocatable::s(:)
      s=sol%s
   end function gets
   function getz(sol) result(z)
      type(cccp_solution),intent(in)::sol
      real(dp),allocatable::z(:)
      z=sol%z
   end function getz
   function getstate(sol) result(state)
      type(cccp_solution),intent(in)::sol
      type(cccp_state)::state
      state=sol%state
   end function getstate
   function getstatus(sol) result(status)
      type(cccp_solution),intent(in)::sol
      character(len=:),allocatable::status
      status=trim(sol%status)
   end function getstatus
   integer function getniter(sol) result(n)
      type(cccp_solution),intent(in)::sol
      n=sol%niter
   end function getniter
   function getparams(control) result(params)
      type(cccp_control),intent(in)::control
      type(cccp_control)::params
      params=control
   end function getparams

end module cccp_solver
