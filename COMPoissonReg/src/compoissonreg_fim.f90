module compoissonreg_fim
   use compoissonreg_kinds, only : dp
   use compoissonreg_types, only : cmp_control_t
   use compoissonreg_normalizer, only : cmp_logz_hybrid
   use compoissonreg_distributions, only : ecmp, cmp_stats, rcmp, rzicmp
   use compoissonreg_numerics, only : logistic
   implicit none
   private
   public :: fim_cmp, fim_cmp_mc, fim_zicmp, fim_zicmp_mc, fim_zicmp_reg

contains

   subroutine logz_derivs(lambda,nu,control,g,h)
      real(dp),intent(in)::lambda,nu
      type(cmp_control_t),intent(in)::control
      real(dp),intent(out)::g(2),h(2,2)
      real(dp)::x(2),hi,hj,f0,fp,fm,fpp,fpm,fmp,fmm
      integer::i,j
      x=[lambda,nu];f0=cmp_logz_hybrid(x(1),x(2),control)
      do i=1,2
         hi=1.0e-5_dp*max(1.0_dp,abs(x(i)))
         if(x(i)-hi<=0.0_dp)then
            fp=eval_shift(x,i,hi,control);g(i)=(fp-f0)/hi
            h(i,i)=(eval_shift(x,i,2.0_dp*hi,control)-2.0_dp*fp+f0)/(hi*hi)
         else
            fp=eval_shift(x,i,hi,control);fm=eval_shift(x,i,-hi,control)
            g(i)=(fp-fm)/(2.0_dp*hi);h(i,i)=(fp-2.0_dp*f0+fm)/(hi*hi)
         end if
      end do
      i=1;j=2;hi=1.0e-4_dp*max(1.0_dp,abs(x(i)));hj=1.0e-4_dp*max(1.0_dp,abs(x(j)))
      if(x(1)-hi>0.0_dp.and.x(2)-hj>0.0_dp)then
         fpp=eval2(x,hi,hj,control);fpm=eval2(x,hi,-hj,control)
         fmp=eval2(x,-hi,hj,control);fmm=eval2(x,-hi,-hj,control)
         h(1,2)=(fpp-fpm-fmp+fmm)/(4.0_dp*hi*hj);h(2,1)=h(1,2)
      else
         fpp=eval2(x,hi,hj,control);fp=eval_shift(x,1,hi,control)
         fm=eval_shift(x,2,hj,control);h(1,2)=(fpp-fp-fm+f0)/(hi*hj);h(2,1)=h(1,2)
      end if
   contains
      function eval_shift(xx,k,delta,ctrl) result(v)
         real(dp),intent(in)::xx(2),delta;integer,intent(in)::k
         type(cmp_control_t),intent(in)::ctrl
         real(dp)::v,z(2);z=xx;z(k)=max(tiny(1.0_dp),z(k)+delta)
         v=cmp_logz_hybrid(z(1),z(2),ctrl)
      end function eval_shift
      function eval2(xx,d1,d2,ctrl) result(v)
         real(dp),intent(in)::xx(2),d1,d2;type(cmp_control_t),intent(in)::ctrl
         real(dp)::v,z(2);z=[max(tiny(1.0_dp),xx(1)+d1),max(tiny(1.0_dp),xx(2)+d2)]
         v=cmp_logz_hybrid(z(1),z(2),ctrl)
      end function eval2
   end subroutine logz_derivs

   subroutine fim_cmp(lambda,nu,fim,control)
      real(dp),intent(in)::lambda,nu
      real(dp),intent(out)::fim(2,2)
      type(cmp_control_t),intent(in),optional::control
      type(cmp_control_t)::ctrl;real(dp)::g(2),h(2,2)
      ctrl=cmp_control_t();if(present(control))ctrl=control
      call logz_derivs(lambda,nu,ctrl,g,h);fim=h
      fim(1,1)=fim(1,1)+ecmp(lambda,nu,ctrl)/(lambda*lambda)
   end subroutine fim_cmp

   subroutine fim_cmp_mc(lambda,nu,reps,fim,control)
      real(dp),intent(in)::lambda,nu
      integer,intent(in)::reps
      real(dp),intent(out)::fim(2,2)
      type(cmp_control_t),intent(in),optional::control
      type(cmp_control_t)::ctrl;integer,allocatable::x(:)
      real(dp)::logz,ey,vy,elogf,s(2);integer::r
      ctrl=cmp_control_t();if(present(control))ctrl=control
      allocate(x(reps));call rcmp(reps,lambda,nu,x,ctrl)
      call cmp_stats(lambda,nu,ctrl,logz,ey,vy,elogf);fim=0.0_dp
      do r=1,reps
         s(1)=(real(x(r),dp)-ey)/lambda
         s(2)=elogf-log_gamma(real(x(r)+1,dp))
         fim=fim+outer(s,s)
      end do
      fim=fim/real(reps,dp)
   end subroutine fim_cmp_mc

   subroutine fim_zicmp(lambda,nu,p,fim,control)
      real(dp),intent(in)::lambda,nu,p
      real(dp),intent(out)::fim(3,3)
      type(cmp_control_t),intent(in),optional::control
      type(cmp_control_t)::ctrl
      real(dp)::g(2),h(2,2),z,meany,den
      ctrl=cmp_control_t();if(present(control))ctrl=control
      z=exp(cmp_logz_hybrid(lambda,nu,ctrl));meany=(1.0_dp-p)*ecmp(lambda,nu,ctrl)
      call logz_derivs(lambda,nu,ctrl,g,h);den=p*(z-1.0_dp)+1.0_dp
      fim=0.0_dp
      fim(1,1)=(1.0_dp-p)*h(1,1)-p*(1.0_dp-p)*g(1)**2/den+meany/lambda**2
      fim(2,2)=(1.0_dp-p)*h(2,2)-p*(1.0_dp-p)*g(2)**2/den
      fim(3,3)=(1.0_dp/z)*(z-1.0_dp)**2/den+(1.0_dp-1.0_dp/z)/(1.0_dp-p)
      fim(1,2)=(1.0_dp-p)*h(1,2)-p*(1.0_dp-p)*g(1)*g(2)/den;fim(2,1)=fim(1,2)
      fim(1,3)=-g(1)/den;fim(3,1)=fim(1,3)
      fim(2,3)=-g(2)/den;fim(3,2)=fim(2,3)
   end subroutine fim_zicmp

   subroutine fim_zicmp_mc(lambda,nu,p,reps,fim,control)
      real(dp),intent(in)::lambda,nu,p
      integer,intent(in)::reps
      real(dp),intent(out)::fim(3,3)
      type(cmp_control_t),intent(in),optional::control
      type(cmp_control_t)::ctrl;integer,allocatable::x(:)
      real(dp)::logz,ey,vy,elogf,f0,den,tau,s(3);integer::r
      ctrl=cmp_control_t();if(present(control))ctrl=control
      allocate(x(reps));call rzicmp(reps,lambda,nu,p,x,ctrl)
      call cmp_stats(lambda,nu,ctrl,logz,ey,vy,elogf);f0=exp(-logz);fim=0.0_dp
      do r=1,reps
         if(x(r)==0)then
            den=p+(1.0_dp-p)*f0;tau=(1.0_dp-p)*f0/den
            s(1)=-tau*ey/lambda;s(2)=tau*elogf;s(3)=(1.0_dp-f0)/den
         else
            s(1)=(real(x(r),dp)-ey)/lambda
            s(2)=elogf-log_gamma(real(x(r)+1,dp));s(3)=-1.0_dp/(1.0_dp-p)
         end if
         fim=fim+outer(s,s)
      end do
      fim=fim/real(reps,dp)
   end subroutine fim_zicmp_mc

   subroutine fim_zicmp_reg(xmat,smat,wmat,beta,gamma,zeta,offx,offs,offw,fim,control)
      real(dp),intent(in)::xmat(:,:),smat(:,:),wmat(:,:),beta(:),gamma(:),zeta(:)
      real(dp),intent(in)::offx(:),offs(:),offw(:)
      real(dp),intent(out)::fim(:,:)
      type(cmp_control_t),intent(in),optional::control
      type(cmp_control_t)::ctrl
      real(dp)::lambda,nu,p,one(3,3),di,dj
      integer::i,a,b,d1,d2,d3,ia,ib,n
      ctrl=cmp_control_t();if(present(control))ctrl=control
      n=size(xmat,1);d1=size(xmat,2);d2=size(smat,2);d3=size(wmat,2)
      if(size(fim,1)/=d1+d2+d3.or.size(fim,2)/=d1+d2+d3)error stop 'fim_zicmp_reg: size mismatch'
      fim=0.0_dp
      do i=1,n
         lambda=exp(dot_product(xmat(i,:),beta)+offx(i))
         nu=exp(dot_product(smat(i,:),gamma)+offs(i))
         p=logistic(dot_product(wmat(i,:),zeta)+offw(i))
         call fim_zicmp(lambda,nu,p,one,ctrl)
         do a=1,d1+d2+d3
            call design_entry(a,i,d1,d2,xmat,smat,wmat,lambda,nu,p,ia,di)
            do b=1,d1+d2+d3
               call design_entry(b,i,d1,d2,xmat,smat,wmat,lambda,nu,p,ib,dj)
               fim(a,b)=fim(a,b)+one(ia,ib)*di*dj
            end do
         end do
      end do
   contains
      subroutine design_entry(k,row,dd1,dd2,xm,sm,wm,lam,nuv,pv,which,val)
         integer,intent(in)::k,row,dd1,dd2
         real(dp),intent(in)::xm(:,:),sm(:,:),wm(:,:),lam,nuv,pv
         integer,intent(out)::which;real(dp),intent(out)::val
         if(k<=dd1)then;which=1;val=lam*xm(row,k)
         else if(k<=dd1+dd2)then;which=2;val=nuv*sm(row,k-dd1)
         else;which=3;val=pv*(1.0_dp-pv)*wm(row,k-dd1-dd2);end if
      end subroutine design_entry
   end subroutine fim_zicmp_reg

   pure function outer(a,b) result(c)
      real(dp),intent(in)::a(:),b(:);real(dp)::c(size(a),size(b));integer::i,j
      do j=1,size(b);do i=1,size(a);c(i,j)=a(i)*b(j);end do;end do
   end function outer

end module compoissonreg_fim
