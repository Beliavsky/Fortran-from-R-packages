! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Derived from OptionPricing 0.1.2 by Wolfgang Hormann and Kemal Dingec.
module optionpricing_asian_mc
   use optionpricing_kinds, only : dp
   use optionpricing_math, only : normal_cdf, normal_pdf, mean_value, sample_sd
   use optionpricing_random, only : seed_rng, fill_normal
   use optionpricing_linalg, only : least_squares
   use optionpricing_types, only : greeks_result, conditional_result, moments_result
   use optionpricing_asian_analytic, only : eval_ecv, eval_lb, eval_eqcv, find_bcv
   implicit none
   private
   public :: asian_call_naive_mc, asian_call_ncv_lr_mc
   public :: asian_call_cmc_cv, asian_call_best_mc
   public :: simulate_asian_call_z, asian_call_naive_greeks_z
   public :: conditional_estimates_z, build_conditional_samples
contains

   function asian_call_naive_mc(n,t,d,k,r,sigma,s0,seed) result(res)
      integer, intent(in) :: n,d
      real(dp), intent(in) :: t,k,r,sigma,s0
      integer, intent(in), optional :: seed
      type(greeks_result) :: res
      real(dp), allocatable :: z(:), price(:), delta(:), gamma(:), logst(:), sumst(:), z1(:)
      real(dp) :: dt,drift,vol,a,disc,l1
      integer :: i
      if (.not.valid_inputs(n,t,d,k,sigma,s0,res)) return
      if (present(seed)) call seed_rng(seed)
      allocate(z(n),price(n),delta(n),gamma(n),logst(n),sumst(n),z1(n))
      dt=t/real(d,dp); drift=(r-0.5_dp*sigma*sigma)*dt
      vol=sigma*sqrt(dt); disc=exp(-r*t); l1=vol
      logst=log(s0); sumst=0.0_dp
      do i=1,d
         call fill_normal(z)
         if (i==1) z1=z
         logst=logst+drift+vol*z
         sumst=sumst+exp(logst)
      end do
      do i=1,n
         a=sumst(i)/real(d,dp)
         price(i)=disc*max(a-k,0.0_dp)
         delta(i)=disc*merge(a/s0,0.0_dp,a>k)
         gamma(i)=disc*merge(k*z1(i)/(s0*s0*l1),0.0_dp,a>k)
      end do
      call summarize_three(price,delta,gamma,res)
   end function asian_call_naive_mc

   function asian_call_ncv_lr_mc(n,t,d,k,r,sigma,s0,seed) result(res)
      integer, intent(in) :: n,d
      real(dp), intent(in) :: t,k,r,sigma,s0
      integer, intent(in), optional :: seed
      type(greeks_result) :: res
      real(dp), allocatable :: z(:), y(:), price(:), delta(:), gamma(:), logst(:), sumst(:), sumlog(:), z1(:)
      real(dp) :: dt,drift,vol,a,g,disc,l1
      type(moments_result) :: ew
      integer :: i
      if (.not.valid_inputs(n,t,d,k,sigma,s0,res)) return
      if (present(seed)) call seed_rng(seed)
      allocate(z(n),y(n),price(n),delta(n),gamma(n),logst(n),sumst(n),sumlog(n),z1(n))
      dt=t/real(d,dp); drift=(r-0.5_dp*sigma*sigma)*dt
      vol=sigma*sqrt(dt); disc=exp(-r*t); l1=vol
      logst=log(s0); sumst=0.0_dp; sumlog=0.0_dp
      do i=1,d
         call fill_normal(z)
         if (i==1) z1=z
         logst=logst+drift+vol*z
         sumst=sumst+exp(logst); sumlog=sumlog+logst
      end do
      ew=eval_ecv(t,d,k,r,sigma,s0)
      do i=1,n
         a=sumst(i)/real(d,dp); g=exp(sumlog(i)/real(d,dp))
         y(i)=disc*max(a-k,0.0_dp)*merge(1.0_dp,0.0_dp,g<k)
         price(i)=y(i)+ew%price
         delta(i)=y(i)*z1(i)/(s0*l1)+ew%delta
         gamma(i)=y(i)*(z1(i)**2-z1(i)*l1-1.0_dp)/(s0*l1)**2+ew%gamma
      end do
      call summarize_three(price,delta,gamma,res)
   end function asian_call_ncv_lr_mc

   function asian_call_cmc_cv(n,t,d,k,r,sigma,s0,np,maxiter,tol,seed) result(res)
      integer, intent(in) :: n,d,np,maxiter
      real(dp), intent(in) :: t,k,r,sigma,s0,tol
      integer, intent(in), optional :: seed
      type(greeks_result) :: res
      real(dp), allocatable :: z(:,:), y(:,:), x(:,:), beta(:,:), adjusted(:)
      real(dp) :: expected(2)
      type(moments_result) :: lb,eq,ew
      integer :: status,total
      if (n<=np .or. np<3) then
         res%status=1; res%message='n must exceed a pilot sample of at least three'; return
      end if
      if (.not.valid_inputs(n,t,d,k,sigma,s0,res)) return
      if (present(seed)) call seed_rng(seed)
      total=n
      allocate(z(d,total)); call fill_normal(z)
      call build_conditional_samples(z,t,d,k,r,sigma,s0,'std',y,x,status=status,maxiter=maxiter,tol=tol)
      if (status/=0) then
         res%status=status; res%message='conditional simulation failed'; return
      end if
      allocate(beta(3,1),adjusted(total-np))
      call least_squares(x(1:np,1:2),reshape(y(1:np,1),[np,1]),beta,status,.true.)
      if (status==1) then
         res%status=2; res%message='pilot regression failed'; return
      end if
      lb=eval_lb(t,d,k,r,sigma,s0,.false.); eq=eval_eqcv(t,d,k,r,sigma,s0); ew=eval_ecv(t,d,k,r,sigma,s0)
      expected=[lb%price,eq%price]
      adjusted=y(np+1:total,1)-matmul(x(np+1:total,1:2)-spread(expected,1,total-np),beta(2:3,1))+ew%price
      res%estimate=0.0_dp; res%error95=0.0_dp
      res%estimate(1)=mean_value(adjusted)
      res%error95(1)=1.96_dp*sample_sd(adjusted)/sqrt(real(size(adjusted),dp))
      res%message='price-only CMC/CV estimate in estimate(1)'
   end function asian_call_cmc_cv

   function asian_call_best_mc(n,t,d,k,r,sigma,s0,np,maxiter,tol,seed) result(res)
      integer, intent(in) :: n,d,np,maxiter
      real(dp), intent(in) :: t,k,r,sigma,s0,tol
      integer, intent(in), optional :: seed
      type(greeks_result) :: res
      real(dp), allocatable :: z(:,:), y(:,:), x(:,:), beta(:,:), adjusted(:,:)
      real(dp) :: expected(6), ewvec(3)
      type(moments_result) :: lb,eq,ew
      integer :: status,total,j
      if (n<=np .or. np<8) then
         res%status=1; res%message='n must exceed a pilot sample of at least eight'; return
      end if
      if (.not.valid_inputs(n,t,d,k,sigma,s0,res)) return
      if (present(seed)) call seed_rng(seed)
      total=n
      allocate(z(d,total)); call fill_normal(z)
      call build_conditional_samples(z,t,d,k,r,sigma,s0,'std',y,x,status=status,maxiter=maxiter,tol=tol)
      if (status/=0) then
         res%status=status; res%message='conditional simulation failed'; return
      end if
      allocate(beta(7,3),adjusted(total-np,3))
      call least_squares(x(1:np,:),y(1:np,:),beta,status,.true.)
      if (status==1) then
         res%status=2; res%message='pilot regression failed'; return
      end if
      lb=eval_lb(t,d,k,r,sigma,s0,.false.); eq=eval_eqcv(t,d,k,r,sigma,s0); ew=eval_ecv(t,d,k,r,sigma,s0)
      expected=[lb%price,eq%price,exp(r*t)*lb%delta,exp(r*t)*lb%gamma,eq%delta,eq%gamma]
      ewvec=[ew%price,ew%delta,ew%gamma]
      adjusted=y(np+1:total,:)-matmul(x(np+1:total,:)-spread(expected,1,total-np),beta(2:7,:))
      adjusted=adjusted+spread(ewvec,1,total-np)
      do j=1,3
         res%estimate(j)=mean_value(adjusted(:,j))
         res%error95(j)=1.96_dp*sample_sd(adjusted(:,j))/sqrt(real(total-np,dp))
      end do
   end function asian_call_best_mc

   pure function simulate_asian_call_z(zmat,t,k,r,sigma,s0) result(value)
      real(dp), intent(in) :: zmat(:,:),t,k,r,sigma,s0
      real(dp) :: value,dt,drift,vol,logst,sumst
      integer :: i,j,n,d
      d=size(zmat,2); n=size(zmat,1)
      dt=t/real(d,dp); drift=(r-0.5_dp*sigma*sigma)*dt; vol=sigma*sqrt(dt)
      value=0.0_dp
      do j=1,n
         logst=log(s0); sumst=0.0_dp
         do i=1,d
            logst=logst+drift+vol*zmat(j,i)
            sumst=sumst+exp(logst)
         end do
         value=value+exp(-r*t)*max(sumst/real(d,dp)-k,0.0_dp)
      end do
      value=value/real(n,dp)
   end function simulate_asian_call_z

   subroutine asian_call_naive_greeks_z(z,t,k,r,sigma,s0,estimate,qmat,mode,status)
      real(dp), intent(in) :: z(:,:),t,k,r,sigma,s0
      real(dp), intent(out) :: estimate(3)
      real(dp), intent(in), optional :: qmat(:,:)
      character(len=*), intent(in), optional :: mode
      integer, intent(out), optional :: status
      real(dp), allocatable :: transformed(:,:), st(:), logst(:), sumst(:), z1(:)
      real(dp) :: dt,drift,vol,a,disc
      integer :: d,n,i,j,istat
      character(len=16) :: method
      d=size(z,1); n=size(z,2); istat=0; method='std'; if (present(mode)) method=adjustl(mode)
      allocate(sumst(n),z1(n)); sumst=0.0_dp
      dt=t/real(d,dp); drift=(r-0.5_dp*sigma*sigma)*dt; vol=sigma*sqrt(dt); disc=exp(-r*t)
      if (trim(method)=='pca') then
         if (.not.present(qmat)) then
            istat=1; estimate=0.0_dp; if(present(status))status=istat; return
         end if
         allocate(transformed(d,n),st(n)); transformed=matmul(qmat,z)
         do j=1,n
            sumst(j)=sum(s0*exp([(drift*real(i,dp)+vol*transformed(i,j),i=1,d)]))
         end do
         z1=transformed(1,:)
      else
         allocate(logst(n)); logst=log(s0)
         do i=1,d
            logst=logst+drift+vol*z(i,:)
            sumst=sumst+exp(logst)
            if(i==1) z1=z(i,:)
         end do
      end if
      estimate=0.0_dp
      do j=1,n
         a=sumst(j)/real(d,dp)
         estimate(1)=estimate(1)+disc*max(a-k,0.0_dp)
         estimate(2)=estimate(2)+disc*merge(a/s0,0.0_dp,a>k)
         estimate(3)=estimate(3)+disc*merge(k*z1(j)/(s0*s0*vol),0.0_dp,a>k)
      end do
      estimate=estimate/real(n,dp)
      if(present(status))status=istat
   end subroutine asian_call_naive_greeks_z

   function conditional_estimates_z(z,t,d,k,r,sigma,s0,mode,qmat,maxiter,tol) result(res)
      real(dp), intent(in) :: z(:,:),t,k,r,sigma,s0
      integer, intent(in) :: d,maxiter
      character(len=*), intent(in) :: mode
      real(dp), intent(in), optional :: qmat(:,:)
      real(dp), intent(in), optional :: tol
      type(conditional_result) :: res
      real(dp), allocatable :: y(:,:),x(:,:)
      real(dp) :: eps
      integer :: status
      eps=1.0e-14_dp; if(present(tol))eps=tol
      if (present(qmat)) then
         call build_conditional_samples(z,t,d,k,r,sigma,s0,mode,y,x,qmat,status,maxiter,eps)
      else
         call build_conditional_samples(z,t,d,k,r,sigma,s0,mode,y,x,status=status,maxiter=maxiter,tol=eps)
      end if
      res%status=status
      if(status/=0) then
         res%message='conditional sample construction failed'; return
      end if
      res%y=sum(y,dim=1)/real(size(y,1),dp)
      res%controls=sum(x,dim=1)/real(size(x,1),dp)
   end function conditional_estimates_z

   subroutine build_conditional_samples(z,t,d,k,r,sigma,s0,mode,y,x,qmat,status,maxiter,tol)
      real(dp), intent(in) :: z(:,:),t,k,r,sigma,s0
      integer, intent(in) :: d
      character(len=*), intent(in) :: mode
      real(dp), allocatable, intent(out) :: y(:,:),x(:,:)
      real(dp), intent(in), optional :: qmat(:,:)
      integer, intent(out) :: status
      integer, intent(in), optional :: maxiter
      real(dp), intent(in), optional :: tol
      real(dp), allocatable :: v(:),a(:),zp(:,:),s(:,:)
      real(dp) :: dt,varx,drift,vol
      integer :: n,i,j,iters
      real(dp) :: eps
      n=size(z,2); status=0; iters=100; if(present(maxiter))iters=maxiter
      eps=1.0e-14_dp; if(present(tol))eps=tol
      if (d<1 .or. n<1 .or. size(z,1)<d-1) then
         status=1; allocate(y(0,3),x(0,6)); return
      end if
      allocate(v(d),a(d),s(d,n))
      dt=t/real(d,dp); varx=real(d*(d+1)*(2*d+1),dp)/6.0_dp
      drift=(r-0.5_dp*sigma*sigma)*dt; vol=sigma*sqrt(dt)
      do i=1,d
         v(i)=real(d-i+1,dp)/sqrt(varx)
         a(i)=vol*sum(v(1:i))
      end do
      select case(trim(adjustl(mode)))
      case('std')
         if(size(z,1)/=d) then
            status=2; allocate(y(0,3),x(0,6)); return
         end if
         allocate(zp(d,n)); zp=z-spread(matmul(v,z),1,d)*spread(v,2,n)
         do j=1,n
            s(1,j)=s0*exp(drift+vol*zp(1,j))
            do i=2,d
               s(i,j)=s(i-1,j)*exp(drift+vol*zp(i,j))
            end do
         end do
      case default
         if(.not.present(qmat)) then
            status=3; allocate(y(0,3),x(0,6)); return
         end if
         if(size(qmat,1)/=d .or. size(qmat,2)>size(z,1)) then
            status=4; allocate(y(0,3),x(0,6)); return
         end if
         allocate(zp(d,n)); zp=matmul(qmat,z(1:size(qmat,2),:))
         do j=1,n
            do i=1,d
               s(i,j)=s0*exp(drift*real(i,dp)+vol*zp(i,j))
            end do
         end do
      end select
      call samples_from_conditional_paths(s,a,t,d,k,r,sigma,s0,y,x,status,iters,eps)
   end subroutine build_conditional_samples

   subroutine samples_from_conditional_paths(s,a,t,d,k,r,sigma,s0,y,x,status,maxiter,tol)
      real(dp), intent(in) :: s(:,:),a(:),t,k,r,sigma,s0,tol
      integer, intent(in) :: d,maxiter
      real(dp), allocatable, intent(out) :: y(:,:),x(:,:)
      integer, intent(out) :: status
      real(dp) :: dt,varx,mus,sigmas,zcut,bcv,disc
      real(dp) :: h0bcv,hpbcv,hppbcv,bcv0,bcv00
      real(dp) :: fk,fdk,alpha,beta,b,func,deriv,intv,ak,apb,apk
      real(dp) :: intcv,cecv,abcv,apbcv,c1,c2,c3,deltacv,gammacv
      real(dp), allocatable :: m(:),mcv(:),expak(:),expab(:),expabcv(:)
      integer :: n,i,j,it,status_b
      n=size(s,2); status=0
      allocate(y(n,3),x(n,6),m(d),mcv(d),expak(d),expab(d),expabcv(d))
      dt=t/real(d,dp); varx=real(d*(d+1)*(2*d+1),dp)/6.0_dp
      mus=log(s0)+(r-0.5_dp*sigma*sigma)*dt*real(d+1,dp)/2.0_dp
      sigmas=sigma/real(d,dp)*sqrt(dt*varx); zcut=(log(k)-mus)/sigmas
      bcv=find_bcv(t,d,k,r,sigma,s0,status_b)
      if(status_b/=0) then
         status=10+status_b; return
      end if
      disc=exp(-r*t)
      mcv=exp(0.5_dp*a*a)*(normal_cdf(zcut-a)-normal_cdf(bcv-a))
      h0bcv=0.0_dp; hpbcv=0.0_dp; hppbcv=0.0_dp
      do i=1,d
         h0bcv=h0bcv+exp(a(i)*bcv+r*real(i,dp)*dt-0.5_dp*a(i)**2)
         hpbcv=hpbcv+a(i)*exp(a(i)*bcv+r*real(i,dp)*dt-0.5_dp*a(i)**2)
         hppbcv=hppbcv+a(i)**2*exp(a(i)*bcv+r*real(i,dp)*dt-0.5_dp*a(i)**2)
      end do
      h0bcv=h0bcv/real(d,dp); hpbcv=s0*hpbcv/real(d,dp); hppbcv=s0*hppbcv/real(d,dp)
      bcv0=-h0bcv/hpbcv
      bcv00=-bcv0*(2.0_dp/s0+hppbcv*bcv0/hpbcv)
      expak=exp(a*zcut)
      expabcv=exp(a*bcv)
      do j=1,n
         fk=sum(expak*s(:,j)); fdk=sum(a*expak*s(:,j))
         alpha=fk/real(d,dp); beta=fdk/real(d,dp)
         if(alpha<=0.0_dp .or. beta<=tiny(1.0_dp)) then
            b=zcut
         else
            b=zcut+alpha/beta*log(k/alpha)
         end if
         do it=1,maxiter
            expab=exp(a*b)
            func=sum(expab*s(:,j))/real(d,dp)-k
            deriv=sum(a*expab*s(:,j))/real(d,dp)
            if(abs(deriv)<=tiny(1.0_dp)) exit
            alpha=b-func/deriv
            if(abs(alpha-b)<tol) then
               b=alpha; exit
            end if
            b=alpha
         end do
         expab=exp(a*b)
         m=exp(0.5_dp*a*a)*(normal_cdf(zcut-a)-normal_cdf(b-a))
         intv=sum(m*s(:,j))/real(d,dp)
         y(j,1)=disc*(intv-k*(normal_cdf(zcut)-normal_cdf(b)))
         ak=sum(expak*s(:,j))/real(d,dp)
         y(j,2)=disc*(intv-(ak-k)*normal_pdf(zcut)/sigmas)/s0
         apb=sum(a*expab*s(:,j))/real(d,dp)
         apk=sum(a*expak*s(:,j))/real(d,dp)
         c1=(k*k*normal_pdf(b)/apb-ak*normal_pdf(zcut)/sigmas)/(s0*s0)
         c2=(apk-ak*sigmas+(ak-k)*(sigmas-zcut))*normal_pdf(zcut)/(s0*sigmas)**2
         y(j,3)=disc*(c1+c2)

         intcv=sum(mcv*s(:,j))/real(d,dp)
         cecv=intcv-k*(normal_cdf(zcut)-normal_cdf(bcv))
         x(j,1)=disc*cecv; x(j,2)=cecv*cecv
         abcv=sum(expabcv*s(:,j))/real(d,dp)
         deltacv=intcv/s0-(ak-k)*normal_pdf(zcut)/(s0*sigmas)-(abcv-k)*normal_pdf(bcv)*bcv0
         apbcv=sum(a*expabcv*s(:,j))/real(d,dp)
         c1=-ak*normal_pdf(zcut)/(s0*s0*sigmas)-abcv/s0*normal_pdf(bcv)*bcv0
         c2=normal_pdf(zcut)/(s0*sigmas)**2*(apk-ak*sigmas+(ak-k)*(sigmas-zcut))
         c3=normal_pdf(bcv)*(abcv/s0*bcv0+apbcv*bcv0*bcv0+(abcv-k)*(bcv00-bcv*bcv0*bcv0))
         gammacv=c1+c2-c3
         x(j,3)=deltacv; x(j,4)=gammacv
         x(j,5)=2.0_dp*deltacv*cecv
         x(j,6)=2.0_dp*(deltacv*deltacv+cecv*gammacv)
      end do
   end subroutine samples_from_conditional_paths

   subroutine summarize_three(price,delta,gamma,res)
      real(dp), intent(in) :: price(:),delta(:),gamma(:)
      type(greeks_result), intent(inout) :: res
      integer :: n
      n=size(price)
      res%estimate=[mean_value(price),mean_value(delta),mean_value(gamma)]
      res%error95=1.96_dp*[sample_sd(price),sample_sd(delta),sample_sd(gamma)]/sqrt(real(n,dp))
   end subroutine summarize_three

   logical function valid_inputs(n,t,d,k,sigma,s0,res) result(ok)
      integer, intent(in) :: n,d
      real(dp), intent(in) :: t,k,sigma,s0
      type(greeks_result), intent(inout) :: res
      ok=n>1 .and. d>0 .and. t>0.0_dp .and. k>0.0_dp .and. sigma>0.0_dp .and. s0>0.0_dp
      if(.not.ok) then
         res%status=1; res%message='invalid option or sample parameters'
      end if
   end function valid_inputs
end module optionpricing_asian_mc
