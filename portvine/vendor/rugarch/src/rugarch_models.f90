module rugarch_models
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use rugarch_kinds, only : dp
   use rugarch_distributions, only : distribution_pdf, random_innovation, absolute_moment
   use rugarch_types, only : garch_spec, model_sgarch, model_gjrgarch, model_egarch, &
      model_aparch, model_igarch, model_figarch, model_csgarch, model_realgarch, model_fgarch
   implicit none
   private

   public :: garch_filter, garch_log_likelihood, simulate_garch
   public :: forecast_volatility, garch_kappa, fgarch_kappa, true_persistence
   public :: news_impact, unconditional_variance, figarch_weights
   public :: realgarch_filter, realgarch_log_likelihood, simulate_realgarch

   interface figarch_weights
      module procedure figarch_weights_scalar
      module procedure figarch_weights_vector
   end interface figarch_weights

contains

   subroutine mean_residuals(y, spec, residuals)
      real(dp), intent(in) :: y(:)
      type(garch_spec), intent(in) :: spec
      real(dp), intent(out) :: residuals(size(y))
      integer :: i, j, ar_order, ma_order

      ar_order=size_or_zero(spec%ar)
      ma_order=size_or_zero(spec%ma)
      residuals=0.0_dp
      do i=1,size(y)
         residuals(i)=y(i)-spec%mean
         do j=1,min(ar_order,i-1)
            residuals(i)=residuals(i)-spec%ar(j)*(y(i-j)-spec%mean)
         end do
         do j=1,min(ma_order,i-1)
            residuals(i)=residuals(i)-spec%ma(j)*residuals(i-j)
         end do
      end do
   end subroutine mean_residuals

   subroutine garch_filter(y, spec, residuals, sigma, valid)
      real(dp), intent(in) :: y(:)
      type(garch_spec), intent(in) :: spec
      real(dp), intent(out) :: residuals(size(y)), sigma(size(y))
      logical, intent(out), optional :: valid
      real(dp), allocatable :: h(:), hp(:), logh(:), z(:), psi(:), qcomp(:)
      real(dp) :: v0, mabs, shock, intercept, zp, base, kdelta
      integer :: n, i, j, p, q, ar_order, ma_order, start, kmax
      logical :: ok

      n=size(y)
      p=size_or_zero(spec%alpha)
      q=size_or_zero(spec%beta)
      ar_order=size_or_zero(spec%ar)
      ma_order=size_or_zero(spec%ma)
      start=max(2,max(max(p,q),max(ar_order,ma_order))+1)
      ok=n>=start .and. spec%omega>=0.0_dp
      if(spec%model==model_aparch .and. spec%delta<=0.0_dp)ok=.false.
      if(size_or_zero(spec%gamma)/=p)ok=.false.
      if(spec%model==model_fgarch) then
         if(size_or_zero(spec%eta1)/=p .or. size_or_zero(spec%eta2)/=p)ok=.false.
         if(spec%fgarch_lambda<=0.0_dp .or. &
            spec%delta+spec%fgarch_fk*spec%fgarch_lambda<=0.0_dp)ok=.false.
      end if
      if(spec%model==model_realgarch)ok=.false.
      if(.not.ok)then
         residuals=0.0_dp;sigma=huge(1.0_dp)
         if(present(valid))valid=.false.
         return
      end if

      call mean_residuals(y,spec,residuals)
      v0=max(sample_variance(residuals),1.0e-12_dp)
      allocate(h(n),z(n))
      h=v0;sigma=sqrt(v0);z=residuals/sigma

      select case(spec%model)
      case(model_egarch)
         allocate(logh(n));logh=log(v0)
         mabs=absolute_moment(1.0_dp,spec%cond_dist,spec%shape,spec%skew,spec%lambda)
         if(.not.finite_value(mabs))mabs=sqrt(2.0_dp/acos(-1.0_dp))
         do i=start,n
            logh(i)=spec%omega
            do j=1,p
               logh(i)=logh(i)+spec%alpha(j)*(abs(z(i-j))-mabs)+spec%gamma(j)*z(i-j)
            end do
            do j=1,q
               logh(i)=logh(i)+spec%beta(j)*logh(i-j)
            end do
            logh(i)=max(-100.0_dp,min(100.0_dp,logh(i)))
            h(i)=exp(logh(i));sigma(i)=sqrt(h(i));z(i)=residuals(i)/sigma(i)
         end do

      case(model_aparch)
         allocate(hp(n));hp=v0**(0.5_dp*spec%delta)
         do i=start,n
            hp(i)=spec%omega
            do j=1,p
               shock=abs(residuals(i-j))-spec%gamma(j)*residuals(i-j)
               hp(i)=hp(i)+spec%alpha(j)*max(shock,0.0_dp)**spec%delta
            end do
            do j=1,q
               hp(i)=hp(i)+spec%beta(j)*hp(i-j)
            end do
            if(hp(i)<=0.0_dp)then;ok=.false.;exit;end if
            sigma(i)=hp(i)**(1.0_dp/spec%delta);h(i)=sigma(i)**2;z(i)=residuals(i)/sigma(i)
         end do

      case(model_fgarch)
         allocate(hp(n));hp=sqrt(v0)**spec%fgarch_lambda
         kdelta=spec%delta+spec%fgarch_fk*spec%fgarch_lambda
         do i=start,n
            hp(i)=spec%omega
            do j=1,p
               zp=residuals(i-j)/max(sigma(i-j),1.0e-12_dp)-spec%eta2(j)
               base=sqrt(1.0e-6_dp+zp*zp)-spec%eta1(j)*zp
               if(base<=0.0_dp)then;ok=.false.;exit;end if
               hp(i)=hp(i)+spec%alpha(j)*base**kdelta*sigma(i-j)**spec%fgarch_lambda
            end do
            if(.not.ok)exit
            do j=1,q
               hp(i)=hp(i)+spec%beta(j)*sigma(i-j)**spec%fgarch_lambda
            end do
            if(hp(i)<=0.0_dp)then;ok=.false.;exit;end if
            sigma(i)=hp(i)**(1.0_dp/spec%fgarch_lambda);h(i)=sigma(i)**2
            z(i)=residuals(i)/sigma(i)
         end do

      case(model_figarch)
         kmax=min(max(10,spec%figarch_truncation),n-1)
         allocate(psi(kmax))
         call figarch_weights(spec%frac_d,spec%alpha,spec%beta,psi)
         if(any(psi < -1.0e-10_dp))ok=.false.
         intercept=spec%omega/max(1.0e-8_dp,1.0_dp-first_or_zero(spec%beta))
         if(ok)then
            do i=2,n
               h(i)=intercept
               do j=1,min(kmax,i-1)
                  h(i)=h(i)+psi(j)*residuals(i-j)**2
               end do
               if(h(i)<=0.0_dp)then;ok=.false.;exit;end if
               sigma(i)=sqrt(h(i));z(i)=residuals(i)/sigma(i)
            end do
         end if

      case(model_csgarch)
         allocate(qcomp(n));qcomp=v0
         do i=start,n
            qcomp(i)=spec%omega+spec%rho*qcomp(i-1)+spec%phi*(residuals(i-1)**2-h(i-1))
            h(i)=qcomp(i)
            do j=1,p
               h(i)=h(i)+spec%alpha(j)*(residuals(i-j)**2-qcomp(i-j))
            end do
            do j=1,q
               h(i)=h(i)+spec%beta(j)*(h(i-j)-qcomp(i-j))
            end do
            if(qcomp(i)<=0.0_dp .or. h(i)<=0.0_dp)then;ok=.false.;exit;end if
            sigma(i)=sqrt(h(i));z(i)=residuals(i)/sigma(i)
         end do

      case default
         do i=start,n
            h(i)=spec%omega
            do j=1,p
               shock=residuals(i-j)**2
               h(i)=h(i)+spec%alpha(j)*shock
               if(spec%model==model_gjrgarch .and. residuals(i-j)<0.0_dp) &
                  h(i)=h(i)+spec%gamma(j)*shock
            end do
            do j=1,q
               h(i)=h(i)+spec%beta(j)*h(i-j)
            end do
            if(h(i)<=0.0_dp)then;ok=.false.;exit;end if
            sigma(i)=sqrt(h(i));z(i)=residuals(i)/sigma(i)
         end do
      end select

      if(any(.not.ieee_is_finite(sigma)) .or. any(sigma<=0.0_dp))ok=.false.
      if(present(valid))valid=ok
   end subroutine garch_filter

   function garch_log_likelihood(y, spec, residuals, sigma) result(loglik)
      real(dp), intent(in) :: y(:)
      type(garch_spec), intent(in) :: spec
      real(dp), intent(out), optional :: residuals(size(y)),sigma(size(y))
      real(dp) :: loglik,density
      real(dp),allocatable::eps(:),vol(:)
      integer::i,start
      logical::valid
      allocate(eps(size(y)),vol(size(y)))
      call garch_filter(y,spec,eps,vol,valid)
      if(.not.valid)then
         loglik=-huge(1.0_dp)
      else
         start=max(2,max(max(size_or_zero(spec%ar),size_or_zero(spec%ma)), &
            max(size_or_zero(spec%alpha),size_or_zero(spec%beta)))+1)
         loglik=0.0_dp
         do i=start,size(y)
            density=distribution_pdf(eps(i)/vol(i),spec%cond_dist,spec%shape,spec%skew,spec%lambda)/vol(i)
            if(density<=tiny(1.0_dp) .or. .not.finite_value(density))then
               loglik=-huge(1.0_dp);exit
            end if
            loglik=loglik+log(density)
         end do
      end if
      if(present(residuals))residuals=eps
      if(present(sigma))sigma=vol
   end function garch_log_likelihood

   subroutine realgarch_filter(returns,realized,spec,residuals,sigma,measurement_residuals,valid)
      real(dp),intent(in)::returns(:),realized(:)
      type(garch_spec),intent(in)::spec
      real(dp),intent(out)::residuals(size(returns)),sigma(size(returns)),measurement_residuals(size(returns))
      logical,intent(out),optional::valid
      real(dp),allocatable::h(:),z(:)
      real(dp)::logh,v0
      integer::i,j,n,p,q,start
      logical::ok

      n=size(returns);p=size_or_zero(spec%alpha);q=size_or_zero(spec%beta)
      start=max(2,max(p,q)+1)
      ok=size(realized)==n .and. n>=start .and. all(realized>0.0_dp) .and. &
         spec%measurement_sd>0.0_dp .and. p+q>0
      call mean_residuals(returns,spec,residuals)
      v0=max(sample_variance(residuals),1.0e-12_dp)
      allocate(h(n),z(n));h=v0;sigma=sqrt(v0);z=residuals/sigma
      measurement_residuals=0.0_dp
      if(ok)then
         do i=start,n
            logh=spec%omega
            do j=1,p
               logh=logh+spec%alpha(j)*log(realized(i-j))
            end do
            do j=1,q
               logh=logh+spec%beta(j)*log(h(i-j))
            end do
            h(i)=exp(max(-100.0_dp,min(100.0_dp,logh)))
            if(h(i)<=0.0_dp)then;ok=.false.;exit;end if
            sigma(i)=sqrt(h(i));z(i)=residuals(i)/sigma(i)
            measurement_residuals(i)=log(realized(i))-spec%xi-spec%phi*log(h(i))- &
               spec%tau1*z(i)-spec%tau2*(z(i)*z(i)-1.0_dp)
         end do
      end if
      if(any(.not.ieee_is_finite(sigma)))ok=.false.
      if(present(valid))valid=ok
   end subroutine realgarch_filter

   function realgarch_log_likelihood(returns,realized,spec,residuals,sigma,measurement_residuals) result(loglik)
      real(dp),intent(in)::returns(:),realized(:)
      type(garch_spec),intent(in)::spec
      real(dp),intent(out),optional::residuals(size(returns)),sigma(size(returns)), &
         measurement_residuals(size(returns))
      real(dp)::loglik,density,mdensity
      real(dp),allocatable::eps(:),vol(:),u(:)
      integer::i,start
      logical::ok
      allocate(eps(size(returns)),vol(size(returns)),u(size(returns)))
      call realgarch_filter(returns,realized,spec,eps,vol,u,ok)
      if(.not.ok)then
         loglik=-huge(1.0_dp)
      else
         loglik=0.0_dp
         start=max(2,max(size_or_zero(spec%alpha),size_or_zero(spec%beta))+1)
         do i=start,size(returns)
            density=distribution_pdf(eps(i)/vol(i),spec%cond_dist,spec%shape,spec%skew,spec%lambda)/vol(i)
            mdensity=exp(-0.5_dp*(u(i)/spec%measurement_sd)**2)/ &
               (sqrt(2.0_dp*acos(-1.0_dp))*spec%measurement_sd)
            if(density<=tiny(1.0_dp) .or. mdensity<=tiny(1.0_dp) .or. &
               .not.finite_value(density) .or. .not.finite_value(mdensity))then
               loglik=-huge(1.0_dp);exit
            end if
            loglik=loglik+log(density)+log(mdensity)
         end do
      end if
      if(present(residuals))residuals=eps
      if(present(sigma))sigma=vol
      if(present(measurement_residuals))measurement_residuals=u
   end function realgarch_log_likelihood

   subroutine simulate_garch(spec,n,y,sigma,residuals,burn_in)
      type(garch_spec),intent(in)::spec
      integer,intent(in)::n
      real(dp),intent(out)::y(n),sigma(n),residuals(n)
      integer,intent(in),optional::burn_in
      real(dp),allocatable::all_y(:),all_sigma(:),all_eps(:),qcomp(:)
      real(dp)::innovation,hnew
      integer::burn,nt,i,j,ar_order,ma_order,p,q
      burn=500;if(present(burn_in))burn=max(0,burn_in)
      nt=n+burn;p=size_or_zero(spec%alpha);q=size_or_zero(spec%beta)
      allocate(all_y(nt),all_sigma(nt),all_eps(nt),qcomp(nt))
      all_y=spec%mean;all_eps=0.0_dp;all_sigma=sqrt(max(unconditional_variance(spec),1.0e-8_dp));qcomp=all_sigma**2
      do i=1,nt
         if(i>1)then
            if(spec%model==model_csgarch)then
               qcomp(i)=spec%omega+spec%rho*qcomp(i-1)+spec%phi*(all_eps(i-1)**2-all_sigma(i-1)**2)
               hnew=qcomp(i)
               do j=1,min(p,i-1)
                  hnew=hnew+spec%alpha(j)*(all_eps(i-j)**2-qcomp(i-j))
               end do
               do j=1,min(q,i-1)
                  hnew=hnew+spec%beta(j)*(all_sigma(i-j)**2-qcomp(i-j))
               end do
               all_sigma(i)=sqrt(max(hnew,1.0e-20_dp))
            else
               call one_step_sigma(spec,all_eps,all_sigma,i)
            end if
         end if
         innovation=random_innovation(spec%cond_dist,spec%shape,spec%skew,spec%lambda)
         all_eps(i)=all_sigma(i)*innovation
         all_y(i)=spec%mean+all_eps(i)
         ar_order=size_or_zero(spec%ar);ma_order=size_or_zero(spec%ma)
         do j=1,min(ar_order,i-1);all_y(i)=all_y(i)+spec%ar(j)*(all_y(i-j)-spec%mean);end do
         do j=1,min(ma_order,i-1);all_y(i)=all_y(i)+spec%ma(j)*all_eps(i-j);end do
      end do
      y=all_y(burn+1:nt);sigma=all_sigma(burn+1:nt);residuals=all_eps(burn+1:nt)
   end subroutine simulate_garch

   subroutine simulate_realgarch(spec,n,returns,realized,sigma,residuals,burn_in)
      type(garch_spec),intent(in)::spec
      integer,intent(in)::n
      real(dp),intent(out)::returns(n),realized(n),sigma(n),residuals(n)
      integer,intent(in),optional::burn_in
      real(dp),allocatable::rall(:),xall(:),sall(:),eall(:)
      real(dp)::z,u,logh,denom
      integer::burn,nt,i,j,p,q,start
      burn=500;if(present(burn_in))burn=max(0,burn_in)
      p=size_or_zero(spec%alpha);q=size_or_zero(spec%beta);start=max(2,max(p,q)+1)
      nt=n+burn;allocate(rall(nt),xall(nt),sall(nt),eall(nt))
      denom=max(1.0e-6_dp,1.0_dp-sum(spec%beta)-spec%phi*sum(spec%alpha))
      sall=sqrt(max(exp((spec%omega+spec%xi*sum(spec%alpha))/denom),1.0e-8_dp))
      xall=sall*sall;rall=spec%mean;eall=0.0_dp
      do i=1,nt
         if(i>=start)then
            logh=spec%omega
            do j=1,p
               logh=logh+spec%alpha(j)*log(max(xall(i-j),1.0e-20_dp))
            end do
            do j=1,q
               logh=logh+spec%beta(j)*log(max(sall(i-j)**2,1.0e-20_dp))
            end do
            sall(i)=sqrt(exp(max(-100.0_dp,min(100.0_dp,logh))))
         end if
         z=random_innovation(spec%cond_dist,spec%shape,spec%skew,spec%lambda)
         eall(i)=sall(i)*z;rall(i)=spec%mean+eall(i)
         u=spec%measurement_sd*random_innovation(10,5.0_dp,1.0_dp)
         xall(i)=exp(spec%xi+spec%phi*log(sall(i)**2)+spec%tau1*z+spec%tau2*(z*z-1.0_dp)+u)
      end do
      returns=rall(burn+1:nt);realized=xall(burn+1:nt);sigma=sall(burn+1:nt);residuals=eall(burn+1:nt)
   end subroutine simulate_realgarch

   subroutine one_step_sigma(spec,eps,sigma,i)
      type(garch_spec),intent(in)::spec
      real(dp),intent(in)::eps(:)
      real(dp),intent(inout)::sigma(:)
      integer,intent(in)::i
      real(dp)::h,hp,z,mabs,shock,zp,base,kdelta,intercept
      real(dp),allocatable::psi(:)
      integer::j,p,q,kmax
      p=size_or_zero(spec%alpha);q=size_or_zero(spec%beta)
      select case(spec%model)
      case(model_egarch)
         h=spec%omega;mabs=absolute_moment(1.0_dp,spec%cond_dist,spec%shape,spec%skew,spec%lambda)
         do j=1,min(p,i-1)
            z=eps(i-j)/max(sigma(i-j),1.0e-12_dp)
            h=h+spec%alpha(j)*(abs(z)-mabs)+spec%gamma(j)*z
         end do
         do j=1,min(q,i-1);h=h+spec%beta(j)*log(max(sigma(i-j)**2,1.0e-20_dp));end do
         sigma(i)=exp(0.5_dp*max(-100.0_dp,min(100.0_dp,h)))
      case(model_aparch)
         hp=spec%omega
         do j=1,min(p,i-1)
            shock=abs(eps(i-j))-spec%gamma(j)*eps(i-j)
            hp=hp+spec%alpha(j)*max(shock,0.0_dp)**spec%delta
         end do
         do j=1,min(q,i-1);hp=hp+spec%beta(j)*sigma(i-j)**spec%delta;end do
         sigma(i)=max(hp,1.0e-20_dp)**(1.0_dp/spec%delta)
      case(model_fgarch)
         hp=spec%omega;kdelta=spec%delta+spec%fgarch_fk*spec%fgarch_lambda
         do j=1,min(p,i-1)
            zp=eps(i-j)/max(sigma(i-j),1.0e-12_dp)-spec%eta2(j)
            base=sqrt(1.0e-6_dp+zp*zp)-spec%eta1(j)*zp
            hp=hp+spec%alpha(j)*max(base,1.0e-12_dp)**kdelta*sigma(i-j)**spec%fgarch_lambda
         end do
         do j=1,min(q,i-1);hp=hp+spec%beta(j)*sigma(i-j)**spec%fgarch_lambda;end do
         sigma(i)=max(hp,1.0e-20_dp)**(1.0_dp/spec%fgarch_lambda)
      case(model_figarch)
         kmax=min(max(10,spec%figarch_truncation),i-1);allocate(psi(kmax))
         call figarch_weights(spec%frac_d,spec%alpha,spec%beta,psi)
         intercept=spec%omega/max(1.0e-8_dp,1.0_dp-first_or_zero(spec%beta));h=intercept
         do j=1,kmax;h=h+max(psi(j),0.0_dp)*eps(i-j)**2;end do
         sigma(i)=sqrt(max(h,1.0e-20_dp))
      case default
         h=spec%omega
         do j=1,min(p,i-1)
            h=h+spec%alpha(j)*eps(i-j)**2
            if(spec%model==model_gjrgarch .and. eps(i-j)<0.0_dp)h=h+spec%gamma(j)*eps(i-j)**2
         end do
         do j=1,min(q,i-1);h=h+spec%beta(j)*sigma(i-j)**2;end do
         sigma(i)=sqrt(max(h,1.0e-20_dp))
      end select
   end subroutine one_step_sigma

   subroutine forecast_volatility(spec,residuals,sigma,horizon,forecast)
      type(garch_spec),intent(in)::spec
      real(dp),intent(in)::residuals(:),sigma(:)
      integer,intent(in)::horizon
      real(dp),intent(out)::forecast(horizon)
      real(dp),allocatable::eps(:),vol(:),qcomp(:)
      real(dp)::hnew,shock_variance
      integer::n,h,i,j,p,q,lag_index
      n=size(residuals);p=size_or_zero(spec%alpha);q=size_or_zero(spec%beta)
      allocate(eps(n+horizon),vol(n+horizon))
      eps(1:n)=residuals;vol(1:n)=sigma
      if(spec%model==model_csgarch)then
         allocate(qcomp(n+horizon));qcomp=vol**2
         do i=2,n
            qcomp(i)=spec%omega+spec%rho*qcomp(i-1)+spec%phi*(eps(i-1)**2-vol(i-1)**2)
         end do
      end if
      do h=1,horizon
         eps(n+h)=0.0_dp
         if(spec%model==model_csgarch)then
            qcomp(n+h)=spec%omega+spec%rho*qcomp(n+h-1)
            hnew=qcomp(n+h)
            do j=1,min(p,n+h-1)
               lag_index=n+h-j
               if(lag_index<=n)then
                  shock_variance=eps(lag_index)**2
               else
                  shock_variance=vol(lag_index)**2
               end if
               hnew=hnew+spec%alpha(j)*(shock_variance-qcomp(lag_index))
            end do
            do j=1,min(q,n+h-1)
               lag_index=n+h-j
               hnew=hnew+spec%beta(j)*(vol(lag_index)**2-qcomp(lag_index))
            end do
            vol(n+h)=sqrt(max(hnew,1.0e-20_dp))
         else
            call one_step_sigma(spec,eps,vol,n+h)
         end if
         forecast(h)=vol(n+h)
      end do
   end subroutine forecast_volatility

   subroutine figarch_weights_scalar(d,alpha,beta,psi)
      real(dp),intent(in)::d,alpha,beta
      real(dp),intent(out)::psi(:)
      real(dp)::avec(1),bvec(1)
      avec(1)=alpha;bvec(1)=beta
      call figarch_weights_vector(d,avec,bvec,psi)
   end subroutine figarch_weights_scalar

   subroutine figarch_weights_vector(d,alpha,beta,psi)
      real(dp),intent(in)::d,alpha(:),beta(:)
      real(dp),intent(out)::psi(:)
      real(dp),allocatable::frac(:),product(:),quotient(:)
      integer::j,k,p,q,n

      n=size(psi);if(n==0)return
      p=size(alpha);q=size(beta)
      allocate(frac(0:n),product(0:n),quotient(0:n))
      frac=0.0_dp;product=0.0_dp;quotient=0.0_dp
      frac(0)=1.0_dp
      do j=1,n
         frac(j)=frac(j-1)*(real(j-1,dp)-d)/real(j,dp)
      end do
      do j=0,n
         product(j)=frac(j)
         do k=1,min(p,j)
            product(j)=product(j)-alpha(k)*frac(j-k)
         end do
      end do
      quotient(0)=1.0_dp
      do j=1,n
         quotient(j)=product(j)
         do k=1,min(q,j)
            quotient(j)=quotient(j)+beta(k)*quotient(j-k)
         end do
         psi(j)=-quotient(j)
      end do
   end subroutine figarch_weights_vector

   function unconditional_variance(spec) result(value)
      type(garch_spec),intent(in)::spec
      real(dp)::value,persistence
      persistence=true_persistence(spec)
      select case(spec%model)
      case(model_egarch)
         value=exp(spec%omega/max(1.0e-8_dp,1.0_dp-sum_if_allocated(spec%beta)))
      case(model_igarch,model_figarch)
         value=max(spec%omega,1.0e-8_dp)
      case(model_realgarch)
         value=exp((spec%omega+spec%xi*sum_if_allocated(spec%alpha))/ &
            max(1.0e-8_dp,1.0_dp-sum_if_allocated(spec%beta)-spec%phi*sum_if_allocated(spec%alpha)))
      case(model_fgarch)
         if(persistence<0.999999_dp)then
            value=(spec%omega/max(1.0e-8_dp,1.0_dp-persistence))**(2.0_dp/spec%fgarch_lambda)
         else;value=max(spec%omega,1.0e-8_dp);end if
      case default
         if(persistence<0.999999_dp)then;value=spec%omega/max(1.0e-8_dp,1.0_dp-persistence)
         else;value=max(spec%omega,1.0e-8_dp);end if
      end select
   end function unconditional_variance

   function garch_kappa(spec,leverage) result(value)
      type(garch_spec),intent(in)::spec
      real(dp),intent(in)::leverage
      real(dp)::value,dx,x,fx
      integer::i,weight
      integer,parameter::ngrid=4000
      real(dp),parameter::limit=30.0_dp
      dx=2.0_dp*limit/real(ngrid,dp);value=0.0_dp
      do i=0,ngrid
         x=-limit+real(i,dp)*dx
         fx=max(abs(x)-leverage*x,0.0_dp)**spec%delta* &
            distribution_pdf(x,spec%cond_dist,spec%shape,spec%skew,spec%lambda)
         if(i==0 .or. i==ngrid)then;weight=1;else if(mod(i,2)==0)then;weight=2;else;weight=4;end if
         value=value+real(weight,dp)*fx
      end do
      value=value*dx/3.0_dp
   end function garch_kappa

   function fgarch_kappa(spec,eta1,eta2) result(value)
      type(garch_spec),intent(in)::spec
      real(dp),intent(in)::eta1,eta2
      real(dp)::value,dx,x,fx,base
      integer::i,weight
      integer,parameter::ngrid=5000
      real(dp),parameter::limit=35.0_dp
      dx=2.0_dp*limit/real(ngrid,dp);value=0.0_dp
      do i=0,ngrid
         x=-limit+real(i,dp)*dx;base=sqrt(1.0e-6_dp+(x-eta2)**2)-eta1*(x-eta2)
         fx=max(base,0.0_dp)**(spec%delta+spec%fgarch_fk*spec%fgarch_lambda)* &
            distribution_pdf(x,spec%cond_dist,spec%shape,spec%skew,spec%lambda)
         if(i==0 .or. i==ngrid)then;weight=1;else if(mod(i,2)==0)then;weight=2;else;weight=4;end if
         value=value+real(weight,dp)*fx
      end do
      value=value*dx/3.0_dp
   end function fgarch_kappa

   function true_persistence(spec) result(value)
      type(garch_spec),intent(in)::spec
      real(dp)::value
      integer::j
      select case(spec%model)
      case(model_gjrgarch)
         value=sum_if_allocated(spec%beta)+sum_if_allocated(spec%alpha)+0.5_dp*sum_if_allocated(spec%gamma)
      case(model_aparch)
         value=sum_if_allocated(spec%beta)
         do j=1,size_or_zero(spec%alpha);value=value+spec%alpha(j)*garch_kappa(spec,spec%gamma(j));end do
      case(model_fgarch)
         value=sum_if_allocated(spec%beta)
         do j=1,size_or_zero(spec%alpha);value=value+spec%alpha(j)*fgarch_kappa(spec,spec%eta1(j),spec%eta2(j));end do
      case(model_egarch)
         value=sum_if_allocated(spec%beta)
      case(model_figarch,model_igarch)
         value=1.0_dp
      case(model_csgarch)
         value=max(spec%rho,sum_if_allocated(spec%alpha)+sum_if_allocated(spec%beta))
      case(model_realgarch)
         value=sum_if_allocated(spec%beta)+sum_if_allocated(spec%alpha)*spec%phi
      case default
         value=sum_if_allocated(spec%alpha)+sum_if_allocated(spec%beta)
      end select
   end function true_persistence

   function news_impact(shock,spec,previous_variance) result(value)
      real(dp),intent(in)::shock,previous_variance
      type(garch_spec),intent(in)::spec
      real(dp)::value,z,mabs,zp,base
      select case(spec%model)
      case(model_gjrgarch)
         value=spec%omega+first_or_zero(spec%alpha)*shock**2+first_or_zero(spec%beta)*previous_variance
         if(shock<0.0_dp)value=value+first_or_zero(spec%gamma)*shock**2
      case(model_egarch)
         z=shock/sqrt(max(previous_variance,1.0e-20_dp))
         mabs=absolute_moment(1.0_dp,spec%cond_dist,spec%shape,spec%skew,spec%lambda)
         value=exp(spec%omega+first_or_zero(spec%alpha)*(abs(z)-mabs)+first_or_zero(spec%gamma)*z+ &
            first_or_zero(spec%beta)*log(max(previous_variance,1.0e-20_dp)))
      case(model_aparch)
         value=(spec%omega+first_or_zero(spec%alpha)*(abs(shock)-first_or_zero(spec%gamma)*shock)**spec%delta+ &
            first_or_zero(spec%beta)*previous_variance**(0.5_dp*spec%delta))**(2.0_dp/spec%delta)
      case(model_fgarch)
         z=shock/sqrt(max(previous_variance,1.0e-20_dp));zp=z-first_or_zero(spec%eta2)
         base=sqrt(1.0e-6_dp+zp*zp)-first_or_zero(spec%eta1)*zp
         value=(spec%omega+first_or_zero(spec%alpha)* &
            base**(spec%delta+spec%fgarch_fk*spec%fgarch_lambda)* &
            previous_variance**(0.5_dp*spec%fgarch_lambda)+first_or_zero(spec%beta)* &
            previous_variance**(0.5_dp*spec%fgarch_lambda))**(2.0_dp/spec%fgarch_lambda)
      case default
         value=spec%omega+first_or_zero(spec%alpha)*shock**2+first_or_zero(spec%beta)*previous_variance
      end select
   end function news_impact

   pure integer function size_or_zero(x) result(n)
      real(dp),allocatable,intent(in)::x(:)
      if(allocated(x))then;n=size(x);else;n=0;end if
   end function size_or_zero

   pure function sum_if_allocated(x) result(value)
      real(dp),allocatable,intent(in)::x(:)
      real(dp)::value
      if(allocated(x))then;value=sum(x);else;value=0.0_dp;end if
   end function sum_if_allocated

   pure function first_or_zero(x) result(value)
      real(dp),allocatable,intent(in)::x(:)
      real(dp)::value
      if(allocated(x).and.size(x)>0)then;value=x(1);else;value=0.0_dp;end if
   end function first_or_zero

   pure function sample_variance(x) result(value)
      real(dp),intent(in)::x(:)
      real(dp)::value,meanx
      if(size(x)<2)then;value=0.0_dp;else;meanx=sum(x)/real(size(x),dp);value=sum((x-meanx)**2)/real(size(x)-1,dp);end if
   end function sample_variance

   pure elemental logical function finite_value(x) result(value)
      real(dp),intent(in)::x
      value=ieee_is_finite(x)
   end function finite_value

end module rugarch_models
