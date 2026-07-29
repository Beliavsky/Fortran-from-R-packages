! SPDX-License-Identifier: GPL-3.0-or-later
module cccp_special
   use cccp_kinds, only : dp
   use cccp_types
   use cccp_cones, only : nnoc, cones_interior
   use cccp_solver, only : solve_lp, solve_qp, solve_dcp
   use cccp_linalg, only : solve_system, vector_norm2
   implicit none
   private

   type, public :: gp_function
      real(dp), allocatable :: f(:,:)
      real(dp), allocatable :: g(:)
   end type gp_function

   public :: l1, rp, gp, logsumexp_value_gradient_hessian

   real(dp), allocatable, save :: rp_p(:,:), rp_mrc(:)
   real(dp), allocatable, save :: gp_f0(:,:), gp_g0(:)
   type(gp_function), allocatable, save :: gp_constraints(:)

contains

   function l1(p, q, control) result(sol)
      real(dp), intent(in) :: p(:,:)
      real(dp), intent(in), optional :: q(:)
      type(cccp_control), intent(in), optional :: control
      type(cccp_solution) :: sol
      real(dp), allocatable :: target(:), gg(:,:), hh(:), qq(:)
      type(cone_constraint) :: c(1)
      integer :: m,n
      m=size(p,1);n=size(p,2)
      allocate(qq(m));qq=0.0_dp;if(present(q))then
         if(size(q)/=m)then;sol%info=cccp_invalid_input;sol%status='invalid input';return;end if
         qq=q
      end if
      allocate(target(n+m),gg(2*m,n+m),hh(2*m));target=0.0_dp;target(n+1:)=1.0_dp;gg=0.0_dp
      gg(1:m,1:n)=p;gg(1:m,n+1:n+m)=-identity_matrix(m)
      gg(m+1:2*m,1:n)=-p;gg(m+1:2*m,n+1:n+m)=-identity_matrix(m)
      hh(1:m)=qq;hh(m+1:)=-qq;c(1)=nnoc(gg,hh)
      call solve_lp(target,cones=c,control=control,sol=sol)
   end function l1

   function rp(x0, p, mrc, control) result(sol)
      real(dp), intent(in) :: x0(:), p(:,:), mrc(:)
      type(cccp_control), intent(in), optional :: control
      type(cccp_solution) :: sol
      type(cone_constraint) :: c(1)
      real(dp), allocatable :: g(:,:),h(:)
      integer::n,i
      n=size(x0)
      if(any(shape(p)/=[n,n]).or.size(mrc)/=n.or.any(x0<=0.0_dp).or.sum(mrc)<=0.0_dp)then
         sol%info=cccp_invalid_input;sol%status='invalid input';return
      end if
      if(allocated(rp_p))deallocate(rp_p,rp_mrc)
      allocate(rp_p(n,n),rp_mrc(n));rp_p=2.0_dp*0.5_dp*(p+transpose(p));rp_mrc=mrc/sum(mrc)
      allocate(g(n,n),h(n));g=0.0_dp
      do i=1,n;g(i,i)=-1.0_dp;end do
      h=0.0_dp;c(1)=nnoc(g,h)
      call solve_dcp(x0,rp_objective,cones=c,control=control,sol=sol)
      if(allocated(sol%x).and.sum(sol%x)>0.0_dp)sol%x=sol%x/sum(sol%x)
   end function rp

   subroutine rp_objective(x,f,g,h,info)
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::f,g(:),h(:,:)
      integer,intent(out)::info
      integer::i,n
      n=size(x)
      if(any(x<=0.0_dp))then;info=1;f=huge(1.0_dp);g=0.0_dp;h=0.0_dp;return;end if
      f=0.5_dp*dot_product(x,matmul(rp_p,x))-dot_product(rp_mrc,log(x))
      g=matmul(rp_p,x)-rp_mrc/x
      h=rp_p
      do i=1,n;h(i,i)=h(i,i)+rp_mrc(i)/(x(i)*x(i));end do
      info=0
   end subroutine rp_objective

   function gp(f0,g0,constraints,nno,a,b,control) result(sol)
      real(dp),intent(in)::f0(:,:),g0(:)
      type(gp_function),intent(in),optional::constraints(:)
      type(cone_constraint),intent(in),optional::nno
      real(dp),intent(in),optional::a(:,:),b(:)
      type(cccp_control),intent(in),optional::control
      type(cccp_solution)::sol
      type(cone_constraint),allocatable::cones(:)
      real(dp),allocatable::x0(:)
      integer::n,m
      n=size(f0,2)
      if(size(g0)/=size(f0,1))then;sol%info=cccp_invalid_input;sol%status='invalid input';return;end if
      if(allocated(gp_f0))deallocate(gp_f0,gp_g0)
      allocate(gp_f0(size(f0,1),n),gp_g0(size(g0)));gp_f0=f0;gp_g0=g0
      m=0;if(present(constraints))m=size(constraints)
      if(allocated(gp_constraints))deallocate(gp_constraints)
      allocate(gp_constraints(m));if(m>0)gp_constraints=constraints
      if(present(nno))then;allocate(cones(1));cones(1)=nno;else;allocate(cones(0));end if
      allocate(x0(n));x0=0.0_dp
      call find_gp_start(x0,cones,a,b,control,sol%info)
      if(sol%info/=cccp_success)then;sol%status='infeasible start';return;end if
      if(m>0)then
         call solve_dcp(x0,gp_objective,a,b,cones,control,sol,gp_nonlinear,m)
      else
         call solve_dcp(x0,gp_objective,a,b,cones,control,sol)
      end if
      if(allocated(sol%x))sol%x=exp(sol%x)
   end function gp



   subroutine find_gp_start(x,cones,a,b,control,info)
      real(dp),intent(inout)::x(:)
      type(cone_constraint),intent(in)::cones(:)
      real(dp),intent(in),optional::a(:,:),b(:)
      type(cccp_control),intent(in),optional::control
      integer,intent(out)::info
      type(cccp_solution),allocatable::qsol
      type(cccp_control)::ctl
      real(dp),allocatable::eye(:,:),zero(:),f(:),gg(:,:),hh(:,:,:),grad(:),trial(:), &
         aat(:,:),rhs(:),lam(:)
      real(dp)::pen,newpen,alpha,v,margin
      integer::n,m,k,iter,linfo,p
      logical::ok
      allocate(qsol)
      n=size(x);m=size(gp_constraints);ctl=cccp_control();if(present(control))ctl=control
      allocate(eye(n,n),zero(n));eye=identity_matrix(n);zero=0.0_dp
      if(size(cones)>0)then
         if(present(a))then
            call solve_qp(1.0e-6_dp*eye,zero,a,b,cones,ctl,qsol)
         else
            call solve_qp(1.0e-6_dp*eye,zero,cones=cones,control=ctl,sol=qsol)
         end if
         if(.not.allocated(qsol%x))then;info=1;return;end if
         x=qsol%x
      else if(present(a))then
         p=size(a,1);allocate(aat(p,p),rhs(p),lam(p))
         aat=matmul(a,transpose(a));call solve_system(aat,b,lam,linfo)
         if(linfo/=0)then;info=1;return;end if
         x=matmul(transpose(a),lam)
         deallocate(aat,rhs,lam)
      end if
      if(m==0)then;info=0;return;end if
      allocate(f(m),gg(m,n),hh(n,n,m),grad(n),trial(n))
      margin=1.0e-5_dp
      do iter=1,2000
         call gp_nonlinear(x,f,gg,hh,linfo)
         if(linfo/=0)then;info=1;return;end if
         ok=all(f < -margin)
         if(ok)then
            if(size(cones)==0)then
               info=0;return
            else
               ok=cones_interior(cones,x)
               if(ok)then;info=0;return;end if
            end if
         end if
         pen=0.0_dp;grad=0.0_dp
         do k=1,m
            v=max(0.0_dp,f(k)+margin)
            pen=pen+0.5_dp*v*v
            grad=grad+v*gg(k,:)
         end do
         if(present(a))then
            p=size(a,1)
            allocate(aat(p,p),rhs(p),lam(p))
            aat=matmul(a,transpose(a));rhs=matmul(a,grad)
            call solve_system(aat,rhs,lam,linfo)
            if(linfo==0)grad=grad-matmul(transpose(a),lam)
            deallocate(aat,rhs,lam)
         end if
         if(vector_norm2(grad)<1.0e-12_dp)exit
         alpha=1.0_dp
         do
            trial=x-alpha*grad
            if(size(cones)>0)then
               if(.not.cones_interior(cones,trial))then
                  alpha=alpha*0.5_dp
                  if(alpha<1e-14_dp)exit
                  cycle
               end if
            end if
            call gp_nonlinear(trial,f,gg,hh,linfo)
            newpen=0.0_dp
            do k=1,m
               v=max(0.0_dp,f(k)+margin);newpen=newpen+0.5_dp*v*v
            end do
            if(newpen<pen)exit
            alpha=alpha*0.5_dp
            if(alpha<1e-14_dp)exit
         end do
         if(alpha<1e-14_dp)exit
         x=trial
      end do
      info=1
   end subroutine find_gp_start

   subroutine gp_objective(x,f,g,h,info)
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::f,g(:),h(:,:)
      integer,intent(out)::info
      call logsumexp_value_gradient_hessian(x,gp_f0,gp_g0,f,g,h)
      info=0
   end subroutine gp_objective

   subroutine gp_nonlinear(x,f,g,h,info)
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::f(:),g(:,:),h(:,:,:)
      integer,intent(out)::info
      integer::k
      do k=1,size(gp_constraints)
         call logsumexp_value_gradient_hessian(x,gp_constraints(k)%f,gp_constraints(k)%g, &
            f(k),g(k,:),h(:,:,k))
      end do
      info=0
   end subroutine gp_nonlinear

   subroutine logsumexp_value_gradient_hessian(x,fmat,gvec,value,grad,hess)
      real(dp),intent(in)::x(:),fmat(:,:),gvec(:)
      real(dp),intent(out)::value,grad(:),hess(:,:)
      real(dp),allocatable::y(:),prob(:),centered(:,:)
      real(dp)::ymax,ysum
      integer::m,n,i
      m=size(fmat,1);n=size(fmat,2)
      allocate(y(m),prob(m),centered(m,n))
      y=matmul(fmat,x)+gvec;ymax=maxval(y);prob=exp(y-ymax);ysum=sum(prob);prob=prob/ysum
      value=ymax+log(ysum);grad=matmul(transpose(fmat),prob)
      do i=1,m;centered(i,:)=sqrt(prob(i))*(fmat(i,:)-grad);end do
      hess=matmul(transpose(centered),centered)
   end subroutine logsumexp_value_gradient_hessian

   pure function identity_matrix(n) result(a)
      integer,intent(in)::n
      real(dp)::a(n,n)
      integer::i
      a=0.0_dp;do i=1,n;a(i,i)=1.0_dp;end do
   end function identity_matrix

end module cccp_special
