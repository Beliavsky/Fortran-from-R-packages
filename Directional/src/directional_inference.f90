module directional_inference
   use directional_kinds, only : dp, pi
   use directional_special, only : log_bessel_i
   use directional_linalg, only : solve_linear
   implicit none
   private
   type, public :: circular_mle_result
      real(dp) :: mu=0.0_dp, rho=0.0_dp, lambda=0.0_dp, loglik=-huge(1.0_dp)
   end type
   type, public :: vmf_mle_result
      real(dp), allocatable :: mu(:)
      real(dp) :: kappa=0.0_dp, loglik=-huge(1.0_dp)
   end type
   public :: cardio_mle, circexp_mle, vmf_mle, spcauchy_mle, pkbd_mle
contains
   function cardio_mle(x,rads) result(res)
      real(dp),intent(in)::x(:); logical,intent(in),optional::rads
      type(circular_mle_result)::res
      real(dp)::z(size(x)),best,mu,rho,val,stepm,stepr
      logical::rr; integer::i,j,it
      rr=.false.;if(present(rads))rr=rads;z=x;if(.not.rr)z=z*pi/180
      mu=modulo(atan2(sum(sin(z)),sum(cos(z))),2*pi);rho=min(0.5_dp,sqrt(sum(cos(z))**2+sum(sin(z))**2)/size(z))
      stepm=pi/4; stepr=0.1_dp; best=cardio_ll(z,mu,rho)
      do it=1,80
         do i=-1,1; do j=-1,1
            val=cardio_ll(z,modulo(mu+i*stepm,2*pi),max(0.0_dp,min(0.5_dp,rho+j*stepr)))
            if(val>best)then;best=val;mu=modulo(mu+i*stepm,2*pi);rho=max(0.0_dp,min(0.5_dp,rho+j*stepr));end if
         end do;end do
         stepm=stepm*0.7_dp;stepr=stepr*0.7_dp
      end do
      res%mu=merge(mu,mu*180/pi,rr);res%rho=rho;res%loglik=best
   end function
   pure real(dp) function cardio_ll(x,mu,rho) result(v)
      real(dp),intent(in)::x(:),mu,rho
      v=-size(x)*log(2*pi)+sum(log(max(tiny(1.0_dp),1+2*rho*cos(x-mu))))
   end function

   function circexp_mle(x,rads,tol) result(res)
      real(dp),intent(in)::x(:);logical,intent(in),optional::rads;real(dp),intent(in),optional::tol
      type(circular_mle_result)::res
      real(dp)::z(size(x)),a,b,c,d,fc,fd,gr,eps
      logical::rr; integer::it
      rr=.false.;if(present(rads))rr=rads;z=x;if(.not.rr)z=z*pi/180
      eps=1e-8_dp;if(present(tol))eps=tol;a=1e-5_dp;b=1000.0_dp;gr=(sqrt(5.0_dp)-1)/2
      c=b-gr*(b-a);d=a+gr*(b-a);fc=circexp_ll(z,c);fd=circexp_ll(z,d)
      do it=1,300
         if(abs(b-a)<eps*max(1.0_dp,c+d))exit
         if(fc>fd)then;b=d;d=c;fd=fc;c=b-gr*(b-a);fc=circexp_ll(z,c)
         else;a=c;c=d;fc=fd;d=a+gr*(b-a);fd=circexp_ll(z,d);end if
      end do
      res%lambda=0.5_dp*(a+b);res%loglik=circexp_ll(z,res%lambda)
   end function
   pure real(dp) function circexp_ll(x,lam) result(v)
      real(dp),intent(in)::x(:),lam
      v=size(x)*log(lam)-lam*sum(x)-size(x)*log(max(tiny(1.0_dp),1-exp(-2*pi*lam)))
   end function

   function vmf_mle(x) result(res)
      real(dp),intent(in)::x(:,:)
      type(vmf_mle_result)::res
      real(dp)::m(size(x,2)),r,apk,k1,k2,nu
      integer::p,it
      p=size(x,2);allocate(res%mu(p));m=sum(x,dim=1)/size(x,1);r=sqrt(sum(m*m));res%mu=m/max(r,tiny(1.0_dp))
      if(r<1e-12_dp)then;res%kappa=0;res%loglik=-size(x,1)*0.5_dp*p*log(2*pi);return;end if
      k1=max(1e-8_dp,r*(p-r*r)/max(1e-10_dp,1-r*r));nu=0.5_dp*p-1
      do it=1,100
         apk=exp(log_bessel_i(nu+1,k1)-log_bessel_i(nu,k1))
         k2=max(1e-10_dp,k1-(apk-r)/max(1e-12_dp,1-apk*apk-real(p-1,dp)*apk/k1))
         if(abs(k2-k1)<1e-10_dp*max(1.0_dp,k1))exit;k1=k2
      end do
      res%kappa=k2
      res%loglik=size(x,1)*(nu*log(res%kappa)-0.5_dp*p*log(2*pi)-log_bessel_i(nu,res%kappa))+res%kappa*sum(matmul(x,res%mu))
   end function

   subroutine spcauchy_mle(x,mu,rho,loglik,tol,maxit)
      real(dp),intent(in)::x(:,:);real(dp),intent(out)::mu(size(x,2)),rho,loglik
      real(dp),intent(in),optional::tol;integer,intent(in),optional::maxit
      real(dp)::mes(size(x,2)),g2,a(size(x,1)),com,com2(size(x,1)),up(size(x,1),size(x,2)),der(size(x,2)),h(size(x,2),size(x,2)),step(size(x,2)),old,eps
      integer::i,j,k,it,nit,info,d,n;eps=1e-6_dp;if(present(tol))eps=tol;nit=100;if(present(maxit))nit=maxit
      n=size(x,1);d=size(x,2)-1;mes=sum(x,dim=1)/n;old=-huge(1.0_dp)
      do it=1,nit
         g2=sum(mes*mes);com=sqrt(g2+1);a=matmul(x,mes);com2=1/(com-a)
         do i=1,n;up(i,:)=x(i,:)-mes/com;end do
         der=d*sum(up*spread(com2,2,size(x,2)),dim=1);h=0
         do j=1,size(x,2);do k=1,size(x,2)
            h(j,k)=d*(sum(up(:,j)*up(:,k)*com2*com2)-sum(com2)*((merge(com,0.0_dp,j==k)-mes(j)*mes(k)/com)/com**2))
         end do;end do
         call solve_linear(h,der,step,info);if(info/=0)exit;mes=mes-step
         loglik=d*sum(log(1/(sqrt(sum(mes*mes)+1)-matmul(x,mes))))+n*log_gamma(.5_dp*(d+1))-0.5_dp*n*(d+1)*log(pi)-n*log(2.0_dp)
         if(abs(loglik-old)<eps)exit;old=loglik
      end do
      rho=(sqrt(sum(mes*mes)+1)-1)/max(sqrt(sum(mes*mes)),tiny(1.0_dp));mu=mes/max(sqrt(sum(mes*mes)),tiny(1.0_dp))
   end subroutine

   subroutine pkbd_mle(x,mu,rho,loglik,tol,maxit)
      real(dp),intent(in)::x(:,:);real(dp),intent(out)::mu(size(x,2)),rho,loglik
      real(dp),intent(in),optional::tol;integer,intent(in),optional::maxit
      ! Stable numerical maximization in the equivalent mesos parameterization.
      real(dp)::mes(size(x,2)),grad(size(x,2)),cand(size(x,2)),base,val,h,eps,step,gamma,com
      integer::j,it,nit,n,d
      eps=1e-6_dp;if(present(tol))eps=tol;nit=200;if(present(maxit))nit=maxit;n=size(x,1);d=size(x,2)-1;mes=sum(x,dim=1)/n;step=0.2_dp
      base=pkbd_mes_ll(x,mes)
      do it=1,nit
         do j=1,size(x,2)
            h=1e-5_dp*max(1.0_dp,abs(mes(j)));cand=mes;cand(j)=cand(j)+h;val=pkbd_mes_ll(x,cand);cand(j)=mes(j)-h
            grad(j)=(val-pkbd_mes_ll(x,cand))/(2*h)
         end do
         cand=mes+step*grad/max(sqrt(sum(grad*grad)),tiny(1.0_dp));val=pkbd_mes_ll(x,cand)
         if(val>base)then;if(abs(val-base)<eps)then;mes=cand;base=val;exit;end if;mes=cand;base=val;step=min(1.0_dp,step*1.1_dp)
         else;step=step*0.5_dp;if(step<1e-10_dp)exit;end if
      end do
      gamma=sqrt(sum(mes*mes));com=sqrt(gamma*gamma+1);rho=(com-1)/max(gamma,tiny(1.0_dp));mu=mes/max(gamma,tiny(1.0_dp))
      loglik=base-0.5_dp*n*(d-1)*log(2.0_dp)+n*log_gamma(.5_dp*(d+1))-0.5_dp*n*(d+1)*log(pi)-n*log(2.0_dp)
   end subroutine
   pure real(dp) function pkbd_mes_ll(x,mes) result(v)
      real(dp),intent(in)::x(:,:),mes(:);real(dp)::g2,com;integer::n,d
      n=size(x,1);d=size(x,2)-1;g2=max(sum(mes*mes),1e-16_dp);com=sqrt(g2+1)
      v=-0.5_dp*(d+1)*sum(log(max(tiny(1.0_dp),com-matmul(x,mes))))-0.5_dp*n*(d-1)*(log(com-1)-log(g2))
   end function
end module directional_inference
