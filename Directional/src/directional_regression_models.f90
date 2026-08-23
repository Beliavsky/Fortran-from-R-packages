module directional_regression_models
   use directional_kinds, only : dp, pi
   use directional_linalg, only : solve_linear, svd3, det3
   use directional_special, only : normal_cdf, normal_pdf
   implicit none
   private

   type, public :: directional_regression_result
      real(dp), allocatable :: beta(:,:)
      real(dp), allocatable :: fitted(:,:)
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: fit = 0.0_dp
      real(dp) :: shape1 = 0.0_dp
      real(dp) :: shape2 = 0.0_dp
      integer :: iterations = 0
      logical :: converged = .false.
   end type

   public :: iag_reg, sipc_reg, cipc_reg, spcauchy_reg, pkbd_reg
   public :: spml_reg, gcpc_reg, esag_reg, sespc_reg, spher_reg

contains

   function iag_reg(y,x,tol,maxit) result(res)
      real(dp), intent(in) :: y(:,:), x(:,:)
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxit
      type(directional_regression_result) :: res
      real(dp), allocatable :: beta(:,:)
      real(dp) :: eps
      integer :: nit
      eps=1.0e-6_dp; if(present(tol)) eps=tol
      nit=3000; if(present(maxit)) nit=maxit
      call ols_beta(x,y,beta)
      call optimize_beta(beta,x,y,1,eps,nit,res%iterations,res%converged)
      call finish_directional_regression(beta,x,y,res)
      res%loglik=-iag_reg_nll(beta,x,y)-1.5_dp*real(size(y,1),dp)*log(2.0_dp*pi)
   end function iag_reg

   function sipc_reg(y,x,tol,maxit) result(res)
      real(dp), intent(in) :: y(:,:), x(:,:)
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxit
      type(directional_regression_result) :: res
      real(dp), allocatable :: beta(:,:)
      real(dp) :: eps
      integer :: nit
      eps=1.0e-6_dp; if(present(tol)) eps=tol
      nit=3000; if(present(maxit)) nit=maxit
      call ols_beta(x,y,beta)
      call optimize_beta(beta,x,y,2,eps,nit,res%iterations,res%converged)
      call finish_directional_regression(beta,x,y,res)
      res%loglik=-sipc_reg_nll(beta,x,y)-real(size(y,1),dp)*log(4.0_dp*pi*pi)
   end function sipc_reg

   function cipc_reg(y,x,tol,maxit) result(res)
      real(dp), intent(in) :: y(:,:), x(:,:)
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxit
      type(directional_regression_result) :: res
      if(size(y,2)/=2) error stop 'cipc_reg: y must have two columns'
      res=spcauchy_reg(y,x,tol,maxit)
   end function cipc_reg

   function spcauchy_reg(y,x,tol,maxit) result(res)
      real(dp), intent(in) :: y(:,:), x(:,:)
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxit
      type(directional_regression_result) :: res
      real(dp), allocatable :: beta(:,:)
      real(dp) :: eps, d
      integer :: nit
      eps=1.0e-6_dp; if(present(tol)) eps=tol
      nit=3000; if(present(maxit)) nit=maxit
      call ols_beta(x,y,beta)
      call optimize_beta(beta,x,y,3,eps,nit,res%iterations,res%converged)
      call finish_directional_regression(beta,x,y,res)
      d=real(size(y,2)-1,dp)
      res%loglik=-spcauchy_reg_nll(beta,x,y)
      res%loglik=res%loglik+real(size(y,1),dp)*log_gamma(0.5_dp*(d+1.0_dp))
      res%loglik=res%loglik-0.5_dp*real(size(y,1),dp)*(d+1.0_dp)*log(pi)
      res%loglik=res%loglik-real(size(y,1),dp)*log(2.0_dp)
   end function spcauchy_reg

   function pkbd_reg(y,x,tol,maxit) result(res)
      real(dp), intent(in) :: y(:,:), x(:,:)
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxit
      type(directional_regression_result) :: res
      real(dp), allocatable :: beta(:,:)
      real(dp) :: eps, d
      integer :: nit
      eps=1.0e-6_dp; if(present(tol)) eps=tol
      nit=3000; if(present(maxit)) nit=maxit
      call ols_beta(x,y,beta)
      call optimize_beta(beta,x,y,4,eps,nit,res%iterations,res%converged)
      call finish_directional_regression(beta,x,y,res)
      d=real(size(y,2)-1,dp)
      res%loglik=-pkbd_reg_nll(beta,x,y)
      res%loglik=res%loglik-0.5_dp*real(size(y,1),dp)*(d-1.0_dp)*log(2.0_dp)
      res%loglik=res%loglik+real(size(y,1),dp)*log_gamma(0.5_dp*(d+1.0_dp))
      res%loglik=res%loglik-0.5_dp*real(size(y,1),dp)*(d+1.0_dp)*log(pi)
      res%loglik=res%loglik-real(size(y,1),dp)*log(2.0_dp)
   end function pkbd_reg

   function spml_reg(y,x,tol,maxit) result(res)
      real(dp), intent(in) :: y(:,:), x(:,:)
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxit
      type(directional_regression_result) :: res
      real(dp), allocatable :: beta(:,:)
      real(dp) :: eps
      integer :: nit
      if(size(y,2)/=2) error stop 'spml_reg: y must have two columns'
      eps=1.0e-6_dp; if(present(tol)) eps=tol
      nit=3000; if(present(maxit)) nit=maxit
      call ols_beta(x,y,beta)
      call optimize_beta(beta,x,y,5,eps,nit,res%iterations,res%converged)
      call finish_directional_regression(beta,x,y,res)
      res%loglik=-spml_reg_nll(beta,x,y)-real(size(y,1),dp)*log(2.0_dp*pi)
   end function spml_reg

   function gcpc_reg(y,x,tol,maxit) result(res)
      real(dp), intent(in) :: y(:,:), x(:,:)
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxit
      type(directional_regression_result) :: res
      real(dp), allocatable :: beta(:,:), par(:)
      real(dp) :: eps
      integer :: nit
      if(size(y,2)/=2) error stop 'gcpc_reg: y must have two columns'
      eps=1.0e-6_dp; if(present(tol)) eps=tol
      nit=5000; if(present(maxit)) nit=maxit
      call ols_beta(x,y,beta)
      allocate(par(1+size(beta)))
      par(1)=log(0.5_dp)
      par(2:)=reshape(beta,[size(beta)])
      call optimize_vector_gcpc(par,x,y,eps,nit,res%iterations,res%converged)
      beta=reshape(par(2:),shape(beta))
      call finish_directional_regression(beta,x,y,res)
      res%shape1=exp(max(-20.0_dp,min(20.0_dp,par(1))))
      res%loglik=-gcpc_reg_nll(par,x,y)-real(size(y,1),dp)*log(2.0_dp*pi)
   end function gcpc_reg

   function esag_reg(y,x,tol,maxit) result(res)
      real(dp), intent(in) :: y(:,:), x(:,:)
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxit
      type(directional_regression_result) :: res
      real(dp), allocatable :: beta(:,:), par(:)
      real(dp) :: eps, rho, den
      integer :: nit
      if(size(y,2)/=3) error stop 'esag_reg: y must have three columns'
      eps=1.0e-6_dp; if(present(tol)) eps=tol
      nit=6000; if(present(maxit)) nit=maxit
      call ols_beta(x,y,beta)
      allocate(par(2+size(beta))); par=0.0_dp
      par(3:)=reshape(beta,[size(beta)])
      call optimize_vector_esag(par,x,y,.false.,eps,nit,res%iterations,res%converged)
      beta=reshape(par(3:),shape(beta))
      call finish_directional_regression(beta,x,y,res)
      rho=sqrt(par(1)**2+par(2)**2+1.0_dp)-sqrt(par(1)**2+par(2)**2)
      den=max(1.0e-14_dp,1.0_dp/rho-rho)
      res%shape1=rho
      res%shape2=0.5_dp*acos(max(-1.0_dp,min(1.0_dp,2.0_dp*par(1)/den)))
      res%loglik=-esag_reg_nll(par,x,y,.false.)-real(size(y,1),dp)*log(2.0_dp*pi)
   end function esag_reg

   function sespc_reg(y,x,tol,maxit) result(res)
      real(dp), intent(in) :: y(:,:), x(:,:)
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxit
      type(directional_regression_result) :: res
      real(dp), allocatable :: beta(:,:), par(:)
      real(dp) :: eps, rho, den
      integer :: nit
      if(size(y,2)/=3) error stop 'sespc_reg: y must have three columns'
      eps=1.0e-6_dp; if(present(tol)) eps=tol
      nit=6000; if(present(maxit)) nit=maxit
      call ols_beta(x,y,beta)
      allocate(par(2+size(beta))); par=0.0_dp
      par(3:)=reshape(beta,[size(beta)])
      call optimize_vector_esag(par,x,y,.true.,eps,nit,res%iterations,res%converged)
      beta=reshape(par(3:),shape(beta))
      call finish_directional_regression(beta,x,y,res)
      rho=sqrt(par(1)**2+par(2)**2+1.0_dp)-sqrt(par(1)**2+par(2)**2)
      den=max(1.0e-14_dp,1.0_dp/rho-rho)
      res%shape1=rho
      res%shape2=0.5_dp*acos(max(-1.0_dp,min(1.0_dp,2.0_dp*par(1)/den)))
      res%loglik=-esag_reg_nll(par,x,y,.true.)-real(size(y,1),dp)*log(2.0_dp*pi)
   end function sespc_reg

   function spher_reg(y,x) result(res)
      real(dp), intent(in) :: y(:,:), x(:,:)
      type(directional_regression_result) :: res
      real(dp) :: xy(3,3), u(3,3), s(3), v(3,3), a(3,3), du
      integer :: i
      if(size(x,2)/=3 .or. size(y,2)/=3 .or. size(x,1)/=size(y,1)) &
         error stop 'spher_reg: x and y must be n by 3'
      xy=matmul(transpose(x),y)
      call svd3(xy,u,s,v)
      a=matmul(v,transpose(u))
      du=det3(a)
      if(du<0.0_dp)then
         u(:,3)=-u(:,3)
         a=matmul(v,transpose(u))
      end if
      allocate(res%beta(3,3),res%fitted(size(x,1),3))
      res%beta=a
      res%fitted=matmul(x,transpose(a))
      do i=1,size(x,1)
         res%fit=res%fit+dot_product(y(i,:),res%fitted(i,:))
      end do
      res%converged=.true.; res%iterations=1
   end function spher_reg

   subroutine ols_beta(x,y,beta)
      real(dp), intent(in) :: x(:,:), y(:,:)
      real(dp), allocatable, intent(out) :: beta(:,:)
      real(dp) :: xtx(size(x,2),size(x,2)), rhs(size(x,2)), sol(size(x,2))
      integer :: j, info
      allocate(beta(size(x,2),size(y,2)))
      xtx=matmul(transpose(x),x)
      do j=1,size(y,2)
         rhs=matmul(transpose(x),y(:,j))
         call solve_linear(xtx,rhs,sol,info)
         if(info/=0)then
            beta(:,j)=0.0_dp
         else
            beta(:,j)=sol
         end if
      end do
   end subroutine ols_beta

   subroutine finish_directional_regression(beta,x,y,res)
      real(dp), intent(in) :: beta(:,:), x(:,:), y(:,:)
      type(directional_regression_result), intent(inout) :: res
      real(dp) :: nr
      integer :: i
      allocate(res%beta(size(beta,1),size(beta,2)))
      allocate(res%fitted(size(x,1),size(beta,2)))
      res%beta=beta
      res%fitted=matmul(x,beta)
      res%fit=0.0_dp
      do i=1,size(x,1)
         nr=sqrt(sum(res%fitted(i,:)**2))
         if(nr>tiny(1.0_dp)) res%fitted(i,:)=res%fitted(i,:)/nr
         res%fit=res%fit+dot_product(y(i,:),res%fitted(i,:))
      end do
   end subroutine finish_directional_regression

   subroutine optimize_beta(beta,x,y,kind,eps,maxit,iters,conv)
      real(dp), intent(inout) :: beta(:,:)
      real(dp), intent(in) :: x(:,:), y(:,:), eps
      integer, intent(in) :: kind, maxit
      integer, intent(out) :: iters
      logical, intent(out) :: conv
      real(dp) :: best, old, val, step
      integer :: i,j,it
      step=0.25_dp
      best=beta_nll(beta,x,y,kind)
      conv=.false.; iters=0
      do it=1,maxit
         old=best
         do j=1,size(beta,2)
            do i=1,size(beta,1)
               beta(i,j)=beta(i,j)+step
               val=beta_nll(beta,x,y,kind)
               if(val<best)then
                  best=val
               else
                  beta(i,j)=beta(i,j)-2.0_dp*step
                  val=beta_nll(beta,x,y,kind)
                  if(val<best)then
                     best=val
                  else
                     beta(i,j)=beta(i,j)+step
                  end if
               end if
            end do
         end do
         iters=it
         if(old-best<=eps*max(1.0_dp,abs(best)))then
            step=0.5_dp*step
            if(step<eps)then; conv=.true.; exit; end if
         end if
      end do
   end subroutine optimize_beta

   real(dp) function beta_nll(beta,x,y,kind) result(v)
      real(dp),intent(in)::beta(:,:),x(:,:),y(:,:)
      integer,intent(in)::kind
      select case(kind)
      case(1); v=iag_reg_nll(beta,x,y)
      case(2); v=sipc_reg_nll(beta,x,y)
      case(3); v=spcauchy_reg_nll(beta,x,y)
      case(4); v=pkbd_reg_nll(beta,x,y)
      case(5); v=spml_reg_nll(beta,x,y)
      case default; v=huge(1.0_dp)
      end select
   end function beta_nll

   real(dp) function iag_reg_nll(beta,x,y) result(v)
      real(dp),intent(in)::beta(:,:),x(:,:),y(:,:)
      real(dp)::mu(size(y,1),size(y,2)),a,m2,pdf,cdf
      integer::i
      if(size(y,2)/=3)then;v=huge(1.0_dp);return;end if
      mu=matmul(x,beta);v=0.5_dp*sum(mu*mu)
      do i=1,size(y,1)
         a=dot_product(y(i,:),mu(i,:));pdf=normal_pdf(a);cdf=normal_cdf(a)
         m2=a+(1.0_dp+a*a)*cdf/max(pdf,tiny(1.0_dp))
         v=v-log(max(m2,tiny(1.0_dp)))
      end do
   end function iag_reg_nll

   real(dp) function sipc_reg_nll(beta,x,y) result(v)
      real(dp),intent(in)::beta(:,:),x(:,:),y(:,:)
      real(dp)::mu(size(y,1),size(y,2)),a,rl,d,sqd,up
      integer::i
      if(size(y,2)/=3)then;v=huge(1.0_dp);return;end if
      mu=matmul(x,beta);v=0.0_dp
      do i=1,size(y,1)
         a=dot_product(y(i,:),mu(i,:));rl=sum(mu(i,:)**2)
         d=max(tiny(1.0_dp),rl+1.0_dp-a*a);sqd=sqrt(d)
         up=(rl+1.0_dp)*sqd*(atan2(sqd,-a)-atan2(sqd,a)+pi)+2.0_dp*a*d
         v=v-log(max(up,tiny(1.0_dp)))+2.0_dp*log(d)
      end do
   end function sipc_reg_nll

   real(dp) function spcauchy_reg_nll(beta,x,y) result(v)
      real(dp),intent(in)::beta(:,:),x(:,:),y(:,:)
      real(dp)::mu(size(y,1),size(y,2)),a,g2,com2,d
      integer::i
      mu=matmul(x,beta); d=real(size(y,2)-1,dp); v=0.0_dp
      do i=1,size(y,1)
         g2=sum(mu(i,:)**2);a=dot_product(y(i,:),mu(i,:))
         com2=sqrt(1.0_dp+g2)-a
         v=v+d*log(max(com2,tiny(1.0_dp)))
      end do
   end function spcauchy_reg_nll

   real(dp) function pkbd_reg_nll(beta,x,y) result(v)
      real(dp),intent(in)::beta(:,:),x(:,:),y(:,:)
      real(dp)::mu(size(y,1),size(y,2)),a,g2,com,com2,d
      integer::i
      mu=matmul(x,beta);d=real(size(y,2)-1,dp);v=0.0_dp
      do i=1,size(y,1)
         g2=max(sum(mu(i,:)**2),1.0e-20_dp);a=dot_product(y(i,:),mu(i,:))
         com=sqrt(g2+1.0_dp);com2=max(com-a,tiny(1.0_dp))
         v=v+0.5_dp*(d+1.0_dp)*log(com2)
         v=v+0.5_dp*(d-1.0_dp)*(log(max(com-1.0_dp,tiny(1.0_dp)))-log(g2))
      end do
   end function pkbd_reg_nll

   real(dp) function spml_reg_nll(beta,x,y) result(v)
      real(dp),intent(in)::beta(:,:),x(:,:),y(:,:)
      real(dp)::mu(size(y,1),2),tau,pt,term
      integer::i
      if(size(y,2)/=2)then;v=huge(1.0_dp);return;end if
      mu=matmul(x,beta);v=0.5_dp*sum(mu*mu)
      do i=1,size(y,1)
         tau=dot_product(y(i,:),mu(i,:));pt=normal_cdf(tau)
         term=1.0_dp+tau*pt*sqrt(2.0_dp*pi)*exp(0.5_dp*tau*tau)
         v=v-log(max(term,tiny(1.0_dp)))
      end do
   end function spml_reg_nll

   subroutine optimize_vector_gcpc(par,x,y,eps,maxit,iters,conv)
      real(dp),intent(inout)::par(:);real(dp),intent(in)::x(:,:),y(:,:),eps
      integer,intent(in)::maxit;integer,intent(out)::iters;logical,intent(out)::conv
      real(dp)::best,old,val,step;integer::i,it
      step=0.2_dp;best=gcpc_reg_nll(par,x,y);conv=.false.
      do it=1,maxit
         old=best
         do i=1,size(par)
            par(i)=par(i)+step;val=gcpc_reg_nll(par,x,y)
            if(val<best)then;best=val
            else
               par(i)=par(i)-2.0_dp*step;val=gcpc_reg_nll(par,x,y)
               if(val<best)then;best=val;else;par(i)=par(i)+step;end if
            end if
         end do
         if(old-best<=eps*max(1.0_dp,abs(best)))then
            step=0.5_dp*step;if(step<eps)then;conv=.true.;exit;end if
         end if
      end do
      iters=min(it,maxit)
   end subroutine optimize_vector_gcpc

   real(dp) function gcpc_reg_nll(par,x,y) result(v)
      real(dp),intent(in)::par(:),x(:,:),y(:,:)
      real(dp)::beta(size(x,2),2),mu(size(y,1),2),rho,g2,nr,ksi(2),s1,s12,s2,a,b,den
      integer::i
      rho=exp(max(-20.0_dp,min(20.0_dp,par(1))));beta=reshape(par(2:),shape(beta));mu=matmul(x,beta)
      v=0.5_dp*real(size(y,1),dp)*log(rho)
      do i=1,size(y,1)
         g2=max(sum(mu(i,:)**2),1.0e-20_dp);nr=sqrt(g2);ksi=mu(i,:)/nr
         s1=ksi(1)**2+ksi(2)**2/rho;s12=ksi(1)*ksi(2)*(1.0_dp-1.0_dp/rho)
         s2=ksi(2)**2+ksi(1)**2/rho;a=dot_product(y(i,:),mu(i,:))
         b=y(i,1)**2*s1+2.0_dp*y(i,1)*y(i,2)*s12+y(i,2)**2*s2
         den=b*sqrt(g2+1.0_dp)-a*sqrt(max(b,tiny(1.0_dp)))
         v=v+log(max(den,tiny(1.0_dp)))
      end do
   end function gcpc_reg_nll

   subroutine optimize_vector_esag(par,x,y,sespc,eps,maxit,iters,conv)
      real(dp),intent(inout)::par(:);real(dp),intent(in)::x(:,:),y(:,:),eps
      logical,intent(in)::sespc;integer,intent(in)::maxit;integer,intent(out)::iters
      logical,intent(out)::conv
      real(dp)::best,old,val,step;integer::i,it
      step=0.15_dp;best=esag_reg_nll(par,x,y,sespc);conv=.false.
      do it=1,maxit
         old=best
         do i=1,size(par)
            par(i)=par(i)+step;val=esag_reg_nll(par,x,y,sespc)
            if(val<best)then;best=val
            else
               par(i)=par(i)-2.0_dp*step;val=esag_reg_nll(par,x,y,sespc)
               if(val<best)then;best=val;else;par(i)=par(i)+step;end if
            end if
         end do
         if(old-best<=eps*max(1.0_dp,abs(best)))then
            step=0.5_dp*step;if(step<eps)then;conv=.true.;exit;end if
         end if
      end do
      iters=min(it,maxit)
   end subroutine optimize_vector_esag

   real(dp) function esag_reg_nll(par,x,y,sespc) result(v)
      real(dp),intent(in)::par(:),x(:,:),y(:,:);logical,intent(in)::sespc
      real(dp)::beta(size(x,2),3),m(size(y,1),3),g1,g2,a,m2,heta,m0,rl
      real(dp)::x1b(3),x2b(3),tx1(3,3),tx2(3,3),vinv(3,3),q(3,3),e,sqe,up
      integer::i
      beta=reshape(par(3:),shape(beta));m=matmul(x,beta);heta=sqrt(par(1)**2+par(2)**2+1.0_dp)-1.0_dp;v=0.0_dp
      do i=1,size(y,1)
         rl=max(sum(m(i,:)**2),1.0e-20_dp);m0=sqrt(max(m(i,2)**2+m(i,3)**2,1.0e-20_dp))
         x1b=[-m0*m0,m(i,1)*m(i,2),m(i,1)*m(i,3)]/(m0*sqrt(rl))
         x2b=[0.0_dp,-m(i,3),m(i,2)]/m0
         tx1=outer3(x1b,x1b);tx2=outer3(x2b,x2b)
         vinv=outer3(x1b,x2b)+outer3(x2b,x1b)
         q=par(2)*vinv+par(1)*(tx1-tx2)+heta*(tx1+tx2)
         g1=1.0_dp+dot_product(y(i,:),matmul(q,y(i,:)))
         g1=max(g1,1.0e-14_dp);g2=dot_product(y(i,:),m(i,:))
         if(.not.sespc)then
            a=g2/sqrt(g1);m2=(1.0_dp+a*a)*normal_cdf(a)+a*normal_pdf(a)
            v=v-0.5_dp*a*a+0.5_dp*rl+1.5_dp*log(g1)-log(max(m2,tiny(1.0_dp)))
         else
            e=max(g1*rl+g1-g2*g2,1.0e-20_dp);sqe=sqrt(e)
            up=g1*(rl+1.0_dp)*sqe*(atan2(sqe,-g2)-atan2(sqe,g2)+pi)+2.0_dp*g2*e
            v=v-log(max(up,tiny(1.0_dp)))+log(g1)+2.0_dp*log(e)
         end if
      end do
   end function esag_reg_nll

   pure function outer3(a,b) result(c)
      real(dp),intent(in)::a(3),b(3);real(dp)::c(3,3);integer::i,j
      do j=1,3;do i=1,3;c(i,j)=a(i)*b(j);end do;end do
   end function outer3

end module directional_regression_models
