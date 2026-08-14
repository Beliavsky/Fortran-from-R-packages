module rebayes_rlr
   use rebayes_kinds, only : dp
   use rebayes_math, only : solve_spd, soft_threshold, log1pexp
   implicit none
   private
   public :: rlr_control, rlr_result, rlr_fit
   type :: rlr_control
      integer :: max_admm = 20000
      integer :: max_newton = 50
      real(dp) :: rho = 1.0_dp
      real(dp) :: tol = 1.0e-9_dp
      real(dp) :: ridge = 1.0e-10_dp
   end type rlr_control
   type :: rlr_result
      real(dp),allocatable :: coef(:)
      real(dp) :: loglik = -huge(1.0_dp)
      integer :: iterations = 0
      integer :: status = 1
   end type rlr_result
contains
   subroutine rlr_fit(x,y,dmat,lambda,result,trials,control)
      real(dp),intent(in)::x(:,:),dmat(:,:),lambda
      real(dp),intent(in)::y(:)
      type(rlr_result),intent(out)::result
      real(dp),intent(in),optional::trials(:)
      type(rlr_control),intent(in),optional::control
      type(rlr_control)::ctl
      real(dp),allocatable::theta(:),z(:),u(:),zold(:),dtz(:),ntr(:)
      real(dp)::rnorm,snorm
      integer::it,p,m,n
      n=size(x,1);p=size(x,2);m=size(dmat,1)
      if(size(y)/=n.or.size(dmat,2)/=p)error stop "rlr_fit: dimensions"
      ctl=rlr_control();if(present(control))ctl=control
      allocate(theta(p),z(m),u(m),zold(m),dtz(m),ntr(n));theta=0;z=0;u=0
      if(present(trials))then;if(size(trials)/=n)error stop "rlr_fit: trials";ntr=trials;else;ntr=1.0_dp;end if
      do it=1,ctl%max_admm
         call theta_step(x,y,ntr,dmat,z,u,ctl,theta)
         dtz=matmul(dmat,theta)+u;zold=z;z=soft_threshold(dtz,lambda/ctl%rho)
         u=u+matmul(dmat,theta)-z
         rnorm=sqrt(sum((matmul(dmat,theta)-z)**2))
         snorm=ctl%rho*sqrt(sum(matmul(transpose(dmat),z-zold)**2))
         if(max(rnorm,snorm)<ctl%tol*sqrt(real(max(1,p),dp)))exit
      end do
      allocate(result%coef(p));result%coef=theta
      result%loglik=-logistic_loss(x,y,ntr,theta)
      result%iterations=min(it,ctl%max_admm);result%status=merge(0,1,it<=ctl%max_admm)
   end subroutine rlr_fit

   subroutine theta_step(x,y,ntr,dmat,z,u,ctl,theta)
      real(dp),intent(in)::x(:,:),y(:),ntr(:),dmat(:,:),z(:),u(:)
      type(rlr_control),intent(in)::ctl
      real(dp),intent(inout)::theta(:)
      real(dp),allocatable::eta(:),prob(:),grad(:),hess(:,:),step(:),r(:),cand(:)
      real(dp)::obj,objnew,alpha,wi
      logical::ok
      integer::i,j,k,p,it,n
      n=size(x,1);p=size(x,2)
      allocate(eta(n),prob(n),grad(p),hess(p,p),step(p),r(size(z)),cand(p))
      do it=1,ctl%max_newton
         eta=matmul(x,theta);prob=1.0_dp/(1.0_dp+exp(-max(-40.0_dp,min(40.0_dp,eta))))
         r=matmul(dmat,theta)-z+u
         grad=matmul(transpose(x),ntr*prob-y)+ctl%rho*matmul(transpose(dmat),r)
         hess=ctl%rho*matmul(transpose(dmat),dmat)
         do j=1,p;hess(j,j)=hess(j,j)+ctl%ridge;end do
         do i=1,n
            wi=ntr(i)*prob(i)*(1.0_dp-prob(i))
            do j=1,p;do k=1,p;hess(j,k)=hess(j,k)+wi*x(i,j)*x(i,k);end do;end do
         end do
         step=-grad;call solve_spd(hess,step,ok);if(.not.ok)exit
         if(sqrt(sum(step*step))<ctl%tol*(1.0_dp+sqrt(sum(theta*theta))))exit
         obj=logistic_loss(x,y,ntr,theta)+0.5_dp*ctl%rho*sum(r*r);alpha=1.0_dp
         do
            cand=theta+alpha*step;r=matmul(dmat,cand)-z+u
            objnew=logistic_loss(x,y,ntr,cand)+0.5_dp*ctl%rho*sum(r*r)
            if(objnew<=obj+1.0e-4_dp*alpha*sum(grad*step).or.alpha<1.0e-8_dp)exit
            alpha=0.5_dp*alpha
         end do
         theta=cand
      end do
   end subroutine theta_step

   real(dp) function logistic_loss(x,y,ntr,theta) result(v)
      real(dp),intent(in)::x(:,:),y(:),ntr(:),theta(:)
      real(dp)::eta
      integer::i
      v=0.0_dp
      do i=1,size(y)
         eta=dot_product(x(i,:),theta)
         v=v+ntr(i)*log1pexp(eta)-y(i)*eta
      end do
   end function logistic_loss
end module rebayes_rlr
