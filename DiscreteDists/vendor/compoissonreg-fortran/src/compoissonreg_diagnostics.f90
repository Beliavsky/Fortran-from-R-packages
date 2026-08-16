module compoissonreg_diagnostics
   use compoissonreg_kinds, only : dp
   use compoissonreg_types, only : cmp_fit_t, zicmp_fit_t, cmp_init_t, cmp_fixed_t, cmp_offset_t
   use compoissonreg_types, only : equitest_t, default_init
   use compoissonreg_distributions, only : ecmp, vcmp, ezicmp, pcmp, pzicmp, rcmp_vec, rzicmp_vec
   use compoissonreg_distributions, only : dcmp, dzicmp, ncmp, tcmp
   use compoissonreg_normalizer, only : z_prodlogj
   use compoissonreg_regression, only : fitted_cmp, fitted_zicmp, fit_cmp_raw, fit_zicmp_raw
   use compoissonreg_numerics, only : qnorm_std, pchisq_upper, invert_matrix
   implicit none
   private
   public :: aic_cmp, bic_cmp, aic_zicmp, bic_zicmp
   public :: sdev_cmp, sdev_zicmp, predict_cmp, predict_zicmp
   public :: residuals_cmp_raw, residuals_cmp_quantile
   public :: residuals_zicmp_raw, residuals_zicmp_quantile
   public :: equitest_cmp, equitest_zicmp, leverage_cmp, deviance_cmp
   public :: bootstrap_cmp, bootstrap_zicmp

contains

   pure real(dp) function aic_cmp(fit)
      type(cmp_fit_t),intent(in)::fit
      aic_cmp=-2.0_dp*fit%loglik+2.0_dp*real(count(.not.fit%fixed%beta)+count(.not.fit%fixed%gamma),dp)
   end function aic_cmp

   pure real(dp) function bic_cmp(fit)
      type(cmp_fit_t),intent(in)::fit
      integer::k
      k=count(.not.fit%fixed%beta)+count(.not.fit%fixed%gamma)
      bic_cmp=-2.0_dp*fit%loglik+log(real(size(fit%y),dp))*real(k,dp)
   end function bic_cmp

   pure real(dp) function aic_zicmp(fit)
      type(zicmp_fit_t),intent(in)::fit
      integer::k
      k=count(.not.fit%fixed%beta)+count(.not.fit%fixed%gamma)+count(.not.fit%fixed%zeta)
      aic_zicmp=-2.0_dp*fit%loglik+2.0_dp*real(k,dp)
   end function aic_zicmp

   pure real(dp) function bic_zicmp(fit)
      type(zicmp_fit_t),intent(in)::fit
      integer::k
      k=count(.not.fit%fixed%beta)+count(.not.fit%fixed%gamma)+count(.not.fit%fixed%zeta)
      bic_zicmp=-2.0_dp*fit%loglik+log(real(size(fit%y),dp))*real(k,dp)
   end function bic_zicmp

   subroutine sdev_cmp(fit,sd)
      type(cmp_fit_t),intent(in)::fit
      real(dp),intent(out)::sd(:)
      integer::i
      if(size(sd)/=size(fit%covariance,1))error stop 'sdev_cmp: size mismatch'
      do i=1,size(sd);sd(i)=sqrt(max(0.0_dp,fit%covariance(i,i)));end do
   end subroutine sdev_cmp

   subroutine sdev_zicmp(fit,sd)
      type(zicmp_fit_t),intent(in)::fit
      real(dp),intent(out)::sd(:)
      integer::i
      if(size(sd)/=size(fit%covariance,1))error stop 'sdev_zicmp: size mismatch'
      do i=1,size(sd);sd(i)=sqrt(max(0.0_dp,fit%covariance(i,i)));end do
   end subroutine sdev_zicmp

   subroutine predict_cmp(fit,xmat,smat,offx,offs,lambda,nu,mean)
      type(cmp_fit_t),intent(in)::fit
      real(dp),intent(in)::xmat(:,:),smat(:,:),offx(:),offs(:)
      real(dp),intent(out)::lambda(:),nu(:)
      real(dp),intent(out),optional::mean(:)
      integer::i
      call fitted_cmp(xmat,smat,fit%beta,fit%gamma,offx,offs,lambda,nu)
      if(present(mean))then
         do i=1,size(mean);mean(i)=ecmp(lambda(i),nu(i),fit%control);end do
      end if
   end subroutine predict_cmp

   subroutine predict_zicmp(fit,xmat,smat,wmat,offx,offs,offw,lambda,nu,p,mean)
      type(zicmp_fit_t),intent(in)::fit
      real(dp),intent(in)::xmat(:,:),smat(:,:),wmat(:,:),offx(:),offs(:),offw(:)
      real(dp),intent(out)::lambda(:),nu(:),p(:)
      real(dp),intent(out),optional::mean(:)
      integer::i
      call fitted_zicmp(xmat,smat,wmat,fit%beta,fit%gamma,fit%zeta,offx,offs,offw,lambda,nu,p)
      if(present(mean))then
         do i=1,size(mean);mean(i)=ezicmp(lambda(i),nu(i),p(i),fit%control);end do
      end if
   end subroutine predict_zicmp

   subroutine residuals_cmp_raw(fit,res)
      type(cmp_fit_t),intent(in)::fit
      real(dp),intent(out)::res(:)
      real(dp),allocatable::lambda(:),nu(:)
      integer::i,n
      n=size(fit%y);allocate(lambda(n),nu(n))
      call fitted_cmp(fit%xmat,fit%smat,fit%beta,fit%gamma,fit%offset%x,fit%offset%s,lambda,nu)
      do i=1,n;res(i)=real(fit%y(i),dp)-ecmp(lambda(i),nu(i),fit%control);end do
   end subroutine residuals_cmp_raw

   subroutine residuals_cmp_quantile(fit,res)
      type(cmp_fit_t),intent(in)::fit
      real(dp),intent(out)::res(:)
      real(dp),allocatable::lambda(:),nu(:)
      real(dp)::fl,fu,u
      integer::i,n
      n=size(fit%y);allocate(lambda(n),nu(n))
      call fitted_cmp(fit%xmat,fit%smat,fit%beta,fit%gamma,fit%offset%x,fit%offset%s,lambda,nu)
      do i=1,n
         fl=pcmp(fit%y(i)-1,lambda(i),nu(i),fit%control)
         fu=pcmp(fit%y(i),lambda(i),nu(i),fit%control)
         call random_number(u);u=fl+(fu-fl)*u
         res(i)=qnorm_std(min(1.0_dp-epsilon(1.0_dp),max(tiny(1.0_dp),u)))
      end do
   end subroutine residuals_cmp_quantile

   subroutine residuals_zicmp_raw(fit,res)
      type(zicmp_fit_t),intent(in)::fit
      real(dp),intent(out)::res(:)
      real(dp),allocatable::lambda(:),nu(:),p(:)
      integer::i,n
      n=size(fit%y);allocate(lambda(n),nu(n),p(n))
      call fitted_zicmp(fit%xmat,fit%smat,fit%wmat,fit%beta,fit%gamma,fit%zeta, &
         fit%offset%x,fit%offset%s,fit%offset%w,lambda,nu,p)
      do i=1,n;res(i)=real(fit%y(i),dp)-ezicmp(lambda(i),nu(i),p(i),fit%control);end do
   end subroutine residuals_zicmp_raw

   subroutine residuals_zicmp_quantile(fit,res)
      type(zicmp_fit_t),intent(in)::fit
      real(dp),intent(out)::res(:)
      real(dp),allocatable::lambda(:),nu(:),p(:)
      real(dp)::fl,fu,u
      integer::i,n
      n=size(fit%y);allocate(lambda(n),nu(n),p(n))
      call fitted_zicmp(fit%xmat,fit%smat,fit%wmat,fit%beta,fit%gamma,fit%zeta, &
         fit%offset%x,fit%offset%s,fit%offset%w,lambda,nu,p)
      do i=1,n
         fl=pzicmp(fit%y(i)-1,lambda(i),nu(i),p(i),fit%control)
         fu=pzicmp(fit%y(i),lambda(i),nu(i),p(i),fit%control)
         call random_number(u);u=fl+(fu-fl)*u
         res(i)=qnorm_std(min(1.0_dp-epsilon(1.0_dp),max(tiny(1.0_dp),u)))
      end do
   end subroutine residuals_zicmp_quantile

   function equitest_cmp(fit) result(test)
      type(cmp_fit_t),intent(in)::fit
      type(equitest_t)::test
      type(cmp_fit_t)::fit0
      type(cmp_init_t)::ini
      type(cmp_fixed_t)::fix
      integer::d1,d2
      d1=size(fit%beta);d2=size(fit%gamma)
      if(any(fit%fixed%gamma))error stop 'equitest_cmp: gamma already fixed'
      ini=default_init(d1,d2);ini%beta=fit%beta;ini%gamma=0.0_dp
      fix=fit%fixed;fix%gamma=.true.
      call fit_cmp_raw(fit%y,fit%xmat,fit%smat,fit0,ini,fix,fit%offset,fit%control)
      test%statistic=max(0.0_dp,-2.0_dp*(fit0%loglik-fit%loglik));test%df=d2
      test%p_value=pchisq_upper(test%statistic,d2)
   end function equitest_cmp

   function equitest_zicmp(fit) result(test)
      type(zicmp_fit_t),intent(in)::fit
      type(equitest_t)::test
      type(zicmp_fit_t)::fit0
      type(cmp_init_t)::ini
      type(cmp_fixed_t)::fix
      integer::d1,d2,d3
      d1=size(fit%beta);d2=size(fit%gamma);d3=size(fit%zeta)
      if(any(fit%fixed%gamma))error stop 'equitest_zicmp: gamma already fixed'
      ini=default_init(d1,d2,d3);ini%beta=fit%beta;ini%gamma=0.0_dp;ini%zeta=fit%zeta
      fix=fit%fixed;fix%gamma=.true.
      call fit_zicmp_raw(fit%y,fit%xmat,fit%smat,fit%wmat,fit0,ini,fix,fit%offset,fit%control)
      test%statistic=max(0.0_dp,-2.0_dp*(fit0%loglik-fit%loglik));test%df=d2
      test%p_value=pchisq_upper(test%statistic,d2)
   end function equitest_zicmp

   subroutine leverage_cmp(fit,lev)
      type(cmp_fit_t),intent(in)::fit
      real(dp),intent(out)::lev(:)
      real(dp),allocatable::lambda(:),nu(:),ey(:),vy(:),curly(:,:),a(:,:),ainv(:,:),v(:)
      real(dp)::z,elogf,den
      integer::i,n,d,mx
      logical::ok
      n=size(fit%y);d=size(fit%beta)+1
      allocate(lambda(n),nu(n),ey(n),vy(n),curly(n,d),a(d,d),ainv(d,d),v(d))
      call fitted_cmp(fit%xmat,fit%smat,fit%beta,fit%gamma,fit%offset%x,fit%offset%s,lambda,nu)
      mx=0
      do i=1,n;mx=max(mx,tcmp(lambda(i),nu(i),fit%control));end do
      do i=1,n
         ey(i)=ecmp(lambda(i),nu(i),fit%control);vy(i)=vcmp(lambda(i),nu(i),fit%control)
         z=ncmp(lambda(i),nu(i),control=fit%control)
         elogf=z_prodlogj(lambda(i),nu(i),mx)/z
         den=real(fit%y(i),dp)-ey(i)
         if(abs(den)<1.0e-10_dp)den=sign(1.0e-10_dp,den+1.0e-30_dp)
         curly(i,1:d-1)=fit%xmat(i,:)
         curly(i,d)=(elogf-log_gamma(real(fit%y(i)+1,dp)))/den
      end do
      a=0.0_dp
      do i=1,n
         a=a+vy(i)*outer(curly(i,:),curly(i,:))
      end do
      call invert_matrix(a,ainv,ok)
      if(.not.ok)then;lev=0.0_dp;return;end if
      do i=1,n
         v=sqrt(max(0.0_dp,vy(i)))*curly(i,:)
         lev(i)=dot_product(v,matmul(ainv,v))
      end do
   contains
      pure function outer(x,y) result(aout)
         real(dp),intent(in)::x(:),y(:);real(dp)::aout(size(x),size(y));integer::ii,jj
         do jj=1,size(y);do ii=1,size(x);aout(ii,jj)=x(ii)*y(jj);end do;end do
      end function outer
   end subroutine leverage_cmp

   subroutine deviance_cmp(fit,dev)
      type(cmp_fit_t),intent(in)::fit
      real(dp),intent(out)::dev(:)
      type(cmp_fit_t)::onefit
      type(cmp_init_t)::ini
      type(cmp_fixed_t)::fix
      real(dp),allocatable::lev(:),lambda(:),nu(:),xone(:,:),sone(:,:)
      type(cmp_offset_t)::offone
      real(dp)::ll,llstar
      integer::i,n,d1,d2
      n=size(fit%y);d1=size(fit%beta);d2=size(fit%gamma)
      if(size(dev)/=n)error stop 'deviance_cmp: size mismatch'
      allocate(lev(n),lambda(n),nu(n),xone(1,d1),sone(1,d2))
      call fitted_cmp(fit%xmat,fit%smat,fit%beta,fit%gamma,fit%offset%x,fit%offset%s,lambda,nu)
      call leverage_cmp(fit,lev)
      do i=1,n
         ini=default_init(d1,d2);ini%beta=fit%beta;ini%gamma=0.0_dp
         fix=fit%fixed;fix%gamma=.true.
         xone(1,:)=fit%xmat(i,:);sone(1,:)=fit%smat(i,:)
         allocate(offone%x(1),offone%s(1),offone%w(1))
         offone%x(1)=fit%offset%x(i);offone%s(1)=fit%offset%s(i);offone%w(1)=fit%offset%w(i)
         call fit_cmp_raw([fit%y(i)],xone,sone,onefit,ini,fix,offone,fit%control)
         llstar=onefit%loglik
         ll=dcmp(fit%y(i),lambda(i),nu(i),.true.,fit%control)
         dev(i)=max(0.0_dp,-2.0_dp*(ll-llstar))/sqrt(max(1.0e-12_dp,1.0_dp-lev(i)))
         deallocate(offone%x,offone%s,offone%w)
      end do
   end subroutine deviance_cmp

   subroutine bootstrap_cmp(fit,reps,samples,success)
      type(cmp_fit_t),intent(in)::fit
      integer,intent(in)::reps
      real(dp),intent(out)::samples(:,:)
      logical,intent(out),optional::success(:)
      real(dp),allocatable::lambda(:),nu(:)
      integer,allocatable::yb(:)
      type(cmp_fit_t)::fb
      type(cmp_init_t)::ini
      integer::r,n,d1,d2
      n=size(fit%y);d1=size(fit%beta);d2=size(fit%gamma)
      if(size(samples,1)/=reps.or.size(samples,2)/=d1+d2)error stop 'bootstrap_cmp: size mismatch'
      allocate(lambda(n),nu(n),yb(n));call fitted_cmp(fit%xmat,fit%smat,fit%beta,fit%gamma, &
         fit%offset%x,fit%offset%s,lambda,nu)
      ini=default_init(d1,d2);ini%beta=fit%beta;ini%gamma=fit%gamma
      do r=1,reps
         call rcmp_vec(lambda,nu,yb,fit%control)
         call fit_cmp_raw(yb,fit%xmat,fit%smat,fb,ini,fit%fixed,fit%offset,fit%control)
         samples(r,1:d1)=fb%beta;samples(r,d1+1:)=fb%gamma
         if(present(success))success(r)=fb%converged
      end do
   end subroutine bootstrap_cmp

   subroutine bootstrap_zicmp(fit,reps,samples,success)
      type(zicmp_fit_t),intent(in)::fit
      integer,intent(in)::reps
      real(dp),intent(out)::samples(:,:)
      logical,intent(out),optional::success(:)
      real(dp),allocatable::lambda(:),nu(:),p(:)
      integer,allocatable::yb(:)
      type(zicmp_fit_t)::fb
      type(cmp_init_t)::ini
      integer::r,n,d1,d2,d3
      n=size(fit%y);d1=size(fit%beta);d2=size(fit%gamma);d3=size(fit%zeta)
      if(size(samples,1)/=reps.or.size(samples,2)/=d1+d2+d3)error stop 'bootstrap_zicmp: size mismatch'
      allocate(lambda(n),nu(n),p(n),yb(n));call fitted_zicmp(fit%xmat,fit%smat,fit%wmat, &
         fit%beta,fit%gamma,fit%zeta,fit%offset%x,fit%offset%s,fit%offset%w,lambda,nu,p)
      ini=default_init(d1,d2,d3);ini%beta=fit%beta;ini%gamma=fit%gamma;ini%zeta=fit%zeta
      do r=1,reps
         call rzicmp_vec(lambda,nu,p,yb,fit%control)
         call fit_zicmp_raw(yb,fit%xmat,fit%smat,fit%wmat,fb,ini,fit%fixed,fit%offset,fit%control)
         samples(r,1:d1)=fb%beta;samples(r,d1+1:d1+d2)=fb%gamma;samples(r,d1+d2+1:)=fb%zeta
         if(present(success))success(r)=fb%converged
      end do
   end subroutine bootstrap_zicmp

end module compoissonreg_diagnostics
