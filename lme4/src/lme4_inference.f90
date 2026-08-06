module lme4_inference
   use lme4_kinds, only : dp
   use lme4_types, only : random_term_t, lmm_control_t, glmm_control_t, &
      lmm_result_t, glmm_result_t, bootstrap_result_t, profile_result_t, &
      influence_result_t
   use lme4_lmm, only : fit_lmm
   use lme4_glmm, only : fit_glmm
   use lme4_simulation, only : simulate_lmm, simulate_glmm
   use lme4_linalg, only : invert_spd
   use lme4_family, only : inverse_normal_cdf
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
   implicit none
   private
   public :: wald_confint_lmm, wald_confint_glmm
   public :: profile_confint_lmm_beta, profile_confint_glmm_beta
   public :: parametric_bootstrap_lmm, parametric_bootstrap_glmm
   public :: bootstrap_percentile_confint
   public :: influence_lmm_groups, influence_glmm_groups
   public :: likelihood_ratio_test

contains

   subroutine wald_confint_lmm(fit,level,interval)
      type(lmm_result_t), intent(in) :: fit
      real(dp), intent(in) :: level
      type(profile_result_t), intent(out) :: interval
      call wald_from_beta(fit%beta,fit%vcov_beta,level,interval)
   end subroutine wald_confint_lmm

   subroutine wald_confint_glmm(fit,level,interval)
      type(glmm_result_t), intent(in) :: fit
      real(dp), intent(in) :: level
      type(profile_result_t), intent(out) :: interval
      call wald_from_beta(fit%beta,fit%vcov_beta,level,interval)
   end subroutine wald_confint_glmm

   subroutine wald_from_beta(beta,vcov,level,interval)
      real(dp), intent(in) :: beta(:), vcov(:,:), level
      type(profile_result_t), intent(out) :: interval
      real(dp) :: z
      integer :: i
      interval%level = level
      allocate(interval%estimate(size(beta)),interval%standard_error(size(beta)), &
         interval%lower(size(beta)),interval%upper(size(beta)))
      interval%estimate = beta
      do i = 1, size(beta)
         interval%standard_error(i) = sqrt(max(0.0_dp,vcov(i,i)))
      end do
      z = inverse_normal_cdf(0.5_dp+0.5_dp*level)
      interval%lower = beta-z*interval%standard_error
      interval%upper = beta+z*interval%standard_error
      interval%status = 0
      interval%message = 'Wald confidence interval'
   end subroutine wald_from_beta

   subroutine profile_confint_lmm_beta(y,x,terms,fit,level,interval,weights,control)
      real(dp), intent(in) :: y(:), x(:,:)
      type(random_term_t), intent(in) :: terms(:)
      type(lmm_result_t), intent(in) :: fit
      real(dp), intent(in) :: level
      type(profile_result_t), intent(out) :: interval
      real(dp), intent(in), optional :: weights(:)
      type(lmm_control_t), intent(in), optional :: control

      type(lmm_result_t) :: base_fit
      type(random_term_t), allocatable :: work_terms(:)
      real(dp) :: target, z
      integer :: j, p

      p = size(x,2)
      call wald_confint_lmm(fit,level,interval)
      if (p < 2) then
         interval%status = 1
         interval%message = 'profile intervals require at least two fixed effects; returned Wald intervals'
         return
      end if
      work_terms = terms
      if (present(weights) .and. present(control)) then
         call fit_lmm(y,x,work_terms,base_fit,reml=.false.,weights=weights,control=control)
      else if (present(weights)) then
         call fit_lmm(y,x,work_terms,base_fit,reml=.false.,weights=weights)
      else if (present(control)) then
         call fit_lmm(y,x,work_terms,base_fit,reml=.false.,control=control)
      else
         call fit_lmm(y,x,work_terms,base_fit,reml=.false.)
      end if
      if (.not. base_fit%converged) then
         interval%status = 2
         interval%message = 'ML refit failed; returned Wald intervals'
         return
      end if
      interval%estimate = base_fit%beta
      interval%standard_error = [(sqrt(max(0.0_dp,base_fit%vcov_beta(j,j))),j=1,p)]
      z = inverse_normal_cdf(0.5_dp+0.5_dp*level)
      target = z*z
      do j = 1, p
         call profile_endpoint_lmm(-1,j,target,base_fit%beta(j), &
            max(interval%standard_error(j),1.0e-3_dp),base_fit%deviance,y,x,terms, &
            interval%lower(j),weights,control)
         call profile_endpoint_lmm(1,j,target,base_fit%beta(j), &
            max(interval%standard_error(j),1.0e-3_dp),base_fit%deviance,y,x,terms, &
            interval%upper(j),weights,control)
      end do
      interval%status = 0
      interval%message = 'profile-likelihood confidence interval for fixed effects'
   end subroutine profile_confint_lmm_beta

   subroutine profile_endpoint_lmm(direction,j,target,estimate,se,base_deviance,y,x,terms, &
      endpoint,weights,control)
      integer, intent(in) :: direction,j
      real(dp), intent(in) :: target,estimate,se,base_deviance,y(:),x(:,:)
      type(random_term_t), intent(in) :: terms(:)
      real(dp), intent(out) :: endpoint
      real(dp), intent(in), optional :: weights(:)
      type(lmm_control_t), intent(in), optional :: control
      real(dp) :: inner,outer,mid,finner,fouter,fmid,step
      integer :: k

      inner = estimate
      finner = -target
      step = max(se,0.05_dp*(1.0_dp+abs(estimate)))
      outer = estimate+real(direction,dp)*step
      call profiled_lmm_difference(outer,j,base_deviance,y,x,terms,fouter,weights,control)
      do k = 1, 12
         if (fouter >= target .or. .not. ieee_is_finite(fouter)) exit
         inner = outer
         finner = fouter
         step = 1.7_dp*step
         outer = estimate+real(direction,dp)*step
         call profiled_lmm_difference(outer,j,base_deviance,y,x,terms,fouter,weights,control)
      end do
      if (.not. ieee_is_finite(fouter) .or. fouter < target) then
         endpoint = estimate+real(direction,dp)*1.96_dp*se
         return
      end if
      do k = 1, 28
         mid = 0.5_dp*(inner+outer)
         call profiled_lmm_difference(mid,j,base_deviance,y,x,terms,fmid,weights,control)
         if (.not. ieee_is_finite(fmid)) then
            outer = mid
         else if (fmid >= target) then
            outer = mid
            fouter = fmid
         else
            inner = mid
            finner = fmid
         end if
      end do
      endpoint = 0.5_dp*(inner+outer)
   end subroutine profile_endpoint_lmm

   subroutine profiled_lmm_difference(value,j,base_deviance,y,x,terms,difference,weights,control)
      real(dp), intent(in) :: value,base_deviance,y(:),x(:,:)
      integer, intent(in) :: j
      type(random_term_t), intent(in) :: terms(:)
      real(dp), intent(out) :: difference
      real(dp), intent(in), optional :: weights(:)
      type(lmm_control_t), intent(in), optional :: control
      real(dp), allocatable :: xdrop(:,:), yadjust(:)
      type(random_term_t), allocatable :: work_terms(:)
      type(lmm_result_t) :: refit

      call drop_column(x,j,xdrop)
      yadjust = y-x(:,j)*value
      work_terms = terms
      if (present(weights) .and. present(control)) then
         call fit_lmm(yadjust,xdrop,work_terms,refit,reml=.false.,weights=weights,control=control)
      else if (present(weights)) then
         call fit_lmm(yadjust,xdrop,work_terms,refit,reml=.false.,weights=weights)
      else if (present(control)) then
         call fit_lmm(yadjust,xdrop,work_terms,refit,reml=.false.,control=control)
      else
         call fit_lmm(yadjust,xdrop,work_terms,refit,reml=.false.)
      end if
      if (refit%converged) then
         difference = refit%deviance-base_deviance
      else
         difference = ieee_value(1.0_dp,ieee_quiet_nan)
      end if
   end subroutine profiled_lmm_difference

   subroutine profile_confint_glmm_beta(y,x,terms,family,fit,level,interval,weights, &
      offset,dispersion,control)
      real(dp), intent(in) :: y(:), x(:,:)
      type(random_term_t), intent(in) :: terms(:)
      integer, intent(in) :: family
      type(glmm_result_t), intent(in) :: fit
      real(dp), intent(in) :: level
      type(profile_result_t), intent(out) :: interval
      real(dp), intent(in), optional :: weights(:), offset(:), dispersion
      type(glmm_control_t), intent(in), optional :: control
      real(dp) :: target,z,disp
      integer :: j,p

      call wald_confint_glmm(fit,level,interval)
      p=size(x,2)
      if (p<2 .or. .not. fit%converged) then
         interval%status=1
         interval%message='GLMM profile intervals require a converged fit with at least two fixed effects; returned Wald intervals'
         return
      end if
      disp=fit%dispersion
      if (present(dispersion)) disp=dispersion
      z=inverse_normal_cdf(0.5_dp+0.5_dp*level)
      target=z*z
      do j=1,p
         call profile_endpoint_glmm(-1,j,target,fit%beta(j), &
            max(interval%standard_error(j),1.0e-3_dp),fit%deviance,y,x,terms,family, &
            interval%lower(j),weights,offset,disp,control)
         call profile_endpoint_glmm(1,j,target,fit%beta(j), &
            max(interval%standard_error(j),1.0e-3_dp),fit%deviance,y,x,terms,family, &
            interval%upper(j),weights,offset,disp,control)
      end do
      interval%status=0
      interval%message='profile-likelihood confidence interval for GLMM fixed effects'
   end subroutine profile_confint_glmm_beta

   subroutine profile_endpoint_glmm(direction,j,target,estimate,se,base_deviance,y,x, &
      terms,family,endpoint,weights,offset,dispersion,control)
      integer, intent(in) :: direction,j,family
      real(dp), intent(in) :: target,estimate,se,base_deviance,y(:),x(:,:),dispersion
      type(random_term_t), intent(in) :: terms(:)
      real(dp), intent(out) :: endpoint
      real(dp), intent(in), optional :: weights(:),offset(:)
      type(glmm_control_t), intent(in), optional :: control
      real(dp) :: inner,outer,mid,fouter,fmid,step
      integer :: k

      inner=estimate
      step=max(se,0.05_dp*(1.0_dp+abs(estimate)))
      outer=estimate+real(direction,dp)*step
      call profiled_glmm_difference(outer,j,base_deviance,y,x,terms,family,fouter, &
         weights,offset,dispersion,control)
      do k=1,12
         if (fouter>=target .or. .not. ieee_is_finite(fouter)) exit
         inner=outer
         step=1.7_dp*step
         outer=estimate+real(direction,dp)*step
         call profiled_glmm_difference(outer,j,base_deviance,y,x,terms,family,fouter, &
            weights,offset,dispersion,control)
      end do
      if (.not. ieee_is_finite(fouter) .or. fouter<target) then
         endpoint=estimate+real(direction,dp)*1.96_dp*se
         return
      end if
      do k=1,28
         mid=0.5_dp*(inner+outer)
         call profiled_glmm_difference(mid,j,base_deviance,y,x,terms,family,fmid, &
            weights,offset,dispersion,control)
         if (.not. ieee_is_finite(fmid) .or. fmid>=target) then
            outer=mid
         else
            inner=mid
         end if
      end do
      endpoint=0.5_dp*(inner+outer)
   end subroutine profile_endpoint_glmm

   subroutine profiled_glmm_difference(value,j,base_deviance,y,x,terms,family,difference, &
      weights,offset,dispersion,control)
      real(dp), intent(in) :: value,base_deviance,y(:),x(:,:),dispersion
      integer, intent(in) :: j,family
      type(random_term_t), intent(in) :: terms(:)
      real(dp), intent(out) :: difference
      real(dp), intent(in), optional :: weights(:),offset(:)
      type(glmm_control_t), intent(in), optional :: control
      real(dp), allocatable :: xdrop(:,:), prof_offset(:)
      type(random_term_t), allocatable :: work_terms(:)
      type(glmm_result_t) :: refit

      call drop_column(x,j,xdrop)
      allocate(prof_offset(size(y)),source=x(:,j)*value)
      if (present(offset)) prof_offset=prof_offset+offset
      work_terms=terms
      if (present(weights) .and. present(control)) then
         call fit_glmm(y,xdrop,work_terms,family,refit,weights=weights,offset=prof_offset, &
            control=control,dispersion=dispersion)
      else if (present(weights)) then
         call fit_glmm(y,xdrop,work_terms,family,refit,weights=weights,offset=prof_offset, &
            dispersion=dispersion)
      else if (present(control)) then
         call fit_glmm(y,xdrop,work_terms,family,refit,offset=prof_offset,control=control, &
            dispersion=dispersion)
      else
         call fit_glmm(y,xdrop,work_terms,family,refit,offset=prof_offset,dispersion=dispersion)
      end if
      if (refit%converged) then
         difference=refit%deviance-base_deviance
      else
         difference=ieee_value(1.0_dp,ieee_quiet_nan)
      end if
   end subroutine profiled_glmm_difference

   subroutine parametric_bootstrap_lmm(x,terms,fit,nsim,bootstrap,seed,control)
      real(dp), intent(in) :: x(:,:)
      type(random_term_t), intent(in) :: terms(:)
      type(lmm_result_t), intent(in) :: fit
      integer, intent(in) :: nsim
      type(bootstrap_result_t), intent(out) :: bootstrap
      integer, intent(in), optional :: seed
      type(lmm_control_t), intent(in), optional :: control
      type(random_term_t), allocatable :: work_terms(:)
      type(lmm_result_t) :: refit
      real(dp), allocatable :: ysim(:), usim(:)
      real(dp) :: nan
      integer :: b,base_seed

      nan=ieee_value(1.0_dp,ieee_quiet_nan)
      allocate(bootstrap%beta(size(fit%beta),nsim),bootstrap%theta(size(fit%theta),nsim), &
         bootstrap%scale(nsim),bootstrap%log_likelihood(nsim),bootstrap%converged(nsim))
      bootstrap%beta=nan; bootstrap%theta=nan; bootstrap%scale=nan
      bootstrap%log_likelihood=nan; bootstrap%converged=.false.
      base_seed=13579
      if (present(seed)) base_seed=seed
      do b=1,nsim
         work_terms=terms
         call simulate_lmm(x,work_terms,fit%beta,fit%varcorr,fit%sigma,ysim,usim, &
            seed=base_seed+104729*b)
         if (present(control)) then
            call fit_lmm(ysim,x,work_terms,refit,reml=fit%reml,control=control)
         else
            call fit_lmm(ysim,x,work_terms,refit,reml=fit%reml)
         end if
         if (refit%converged) then
            bootstrap%beta(:,b)=refit%beta
            bootstrap%theta(:,b)=refit%theta
            bootstrap%scale(b)=refit%sigma
            bootstrap%log_likelihood(b)=refit%log_likelihood
            bootstrap%converged(b)=.true.
         end if
      end do
      bootstrap%successful=count(bootstrap%converged)
   end subroutine parametric_bootstrap_lmm

   subroutine parametric_bootstrap_glmm(x,terms,fit,nsim,bootstrap,seed,control)
      real(dp), intent(in) :: x(:,:)
      type(random_term_t), intent(in) :: terms(:)
      type(glmm_result_t), intent(in) :: fit
      integer, intent(in) :: nsim
      type(bootstrap_result_t), intent(out) :: bootstrap
      integer, intent(in), optional :: seed
      type(glmm_control_t), intent(in), optional :: control
      type(random_term_t), allocatable :: work_terms(:)
      type(glmm_result_t) :: refit
      real(dp), allocatable :: ysim(:), usim(:)
      real(dp) :: nan
      integer :: b,base_seed

      nan=ieee_value(1.0_dp,ieee_quiet_nan)
      allocate(bootstrap%beta(size(fit%beta),nsim),bootstrap%theta(size(fit%theta),nsim), &
         bootstrap%scale(nsim),bootstrap%log_likelihood(nsim),bootstrap%converged(nsim))
      bootstrap%beta=nan; bootstrap%theta=nan; bootstrap%scale=nan
      bootstrap%log_likelihood=nan; bootstrap%converged=.false.
      base_seed=24680
      if (present(seed)) base_seed=seed
      do b=1,nsim
         work_terms=terms
         call simulate_glmm(x,work_terms,fit%beta,fit%varcorr,fit%family,ysim,usim, &
            dispersion=fit%dispersion,seed=base_seed+104729*b)
         if (present(control)) then
            call fit_glmm(ysim,x,work_terms,fit%family,refit,control=control, &
               dispersion=fit%dispersion)
         else
            call fit_glmm(ysim,x,work_terms,fit%family,refit,dispersion=fit%dispersion)
         end if
         if (refit%converged) then
            bootstrap%beta(:,b)=refit%beta
            bootstrap%theta(:,b)=refit%theta
            bootstrap%scale(b)=refit%dispersion
            bootstrap%log_likelihood(b)=refit%log_likelihood
            bootstrap%converged(b)=.true.
         end if
      end do
      bootstrap%successful=count(bootstrap%converged)
   end subroutine parametric_bootstrap_glmm

   subroutine bootstrap_percentile_confint(samples,level,lower,upper)
      real(dp), intent(in) :: samples(:,:),level
      real(dp), allocatable, intent(out) :: lower(:),upper(:)
      real(dp), allocatable :: values(:)
      integer :: i,n
      allocate(lower(size(samples,1)),upper(size(samples,1)))
      do i=1,size(samples,1)
         values=pack(samples(i,:),ieee_is_finite(samples(i,:)))
         n=size(values)
         if (n<1) then
            lower(i)=ieee_value(1.0_dp,ieee_quiet_nan)
            upper(i)=lower(i)
         else
            call sort_real(values)
            lower(i)=quantile_sorted(values,0.5_dp*(1.0_dp-level))
            upper(i)=quantile_sorted(values,0.5_dp*(1.0_dp+level))
         end if
      end do
   end subroutine bootstrap_percentile_confint

   subroutine influence_lmm_groups(y,x,terms,fit,term_index,influence,control)
      real(dp), intent(in) :: y(:),x(:,:)
      type(random_term_t), intent(in) :: terms(:)
      type(lmm_result_t), intent(in) :: fit
      integer, intent(in) :: term_index
      type(influence_result_t), intent(out) :: influence
      type(lmm_control_t), intent(in), optional :: control
      type(random_term_t), allocatable :: reduced_terms(:)
      type(lmm_result_t) :: reduced_fit
      real(dp), allocatable :: yr(:),xr(:,:),precision(:,:),delta(:)
      logical, allocatable :: keep(:)
      integer :: g,info,p,nlev

      p=size(fit%beta)
      nlev=terms(term_index)%n_levels
      allocate(influence%dfbeta(p,nlev),influence%cooks_distance(nlev), &
         influence%deleted_log_likelihood(nlev),influence%converged(nlev))
      influence%dfbeta=0.0_dp; influence%cooks_distance=0.0_dp
      influence%deleted_log_likelihood=-huge(1.0_dp); influence%converged=.false.
      call invert_spd(fit%vcov_beta,precision,info)
      if (info/=0) precision=0.0_dp
      allocate(keep(size(y)))
      do g=1,nlev
         keep=terms(term_index)%group/=g
         call subset_model(y,x,terms,keep,term_index,g,yr,xr,reduced_terms)
         if (present(control)) then
            call fit_lmm(yr,xr,reduced_terms,reduced_fit,reml=fit%reml,control=control)
         else
            call fit_lmm(yr,xr,reduced_terms,reduced_fit,reml=fit%reml)
         end if
         if (reduced_fit%converged) then
            delta=reduced_fit%beta-fit%beta
            influence%dfbeta(:,g)=delta
            influence%cooks_distance(g)=dot_product(delta,matmul(precision,delta))/real(p,dp)
            influence%deleted_log_likelihood(g)=reduced_fit%log_likelihood
            influence%converged(g)=.true.
         end if
      end do
   end subroutine influence_lmm_groups

   subroutine influence_glmm_groups(y,x,terms,fit,term_index,influence,control)
      real(dp), intent(in) :: y(:),x(:,:)
      type(random_term_t), intent(in) :: terms(:)
      type(glmm_result_t), intent(in) :: fit
      integer, intent(in) :: term_index
      type(influence_result_t), intent(out) :: influence
      type(glmm_control_t), intent(in), optional :: control
      type(random_term_t), allocatable :: reduced_terms(:)
      type(glmm_result_t) :: reduced_fit
      real(dp), allocatable :: yr(:),xr(:,:),precision(:,:),delta(:)
      logical, allocatable :: keep(:)
      integer :: g,info,p,nlev

      p=size(fit%beta)
      nlev=terms(term_index)%n_levels
      allocate(influence%dfbeta(p,nlev),influence%cooks_distance(nlev), &
         influence%deleted_log_likelihood(nlev),influence%converged(nlev))
      influence%dfbeta=0.0_dp; influence%cooks_distance=0.0_dp
      influence%deleted_log_likelihood=-huge(1.0_dp); influence%converged=.false.
      call invert_spd(fit%vcov_beta,precision,info)
      if (info/=0) precision=0.0_dp
      allocate(keep(size(y)))
      do g=1,nlev
         keep=terms(term_index)%group/=g
         call subset_model(y,x,terms,keep,term_index,g,yr,xr,reduced_terms)
         if (present(control)) then
            call fit_glmm(yr,xr,reduced_terms,fit%family,reduced_fit,control=control, &
               dispersion=fit%dispersion)
         else
            call fit_glmm(yr,xr,reduced_terms,fit%family,reduced_fit,dispersion=fit%dispersion)
         end if
         if (reduced_fit%converged) then
            delta=reduced_fit%beta-fit%beta
            influence%dfbeta(:,g)=delta
            influence%cooks_distance(g)=dot_product(delta,matmul(precision,delta))/real(p,dp)
            influence%deleted_log_likelihood(g)=reduced_fit%log_likelihood
            influence%converged(g)=.true.
         end if
      end do
   end subroutine influence_glmm_groups

   subroutine likelihood_ratio_test(reduced_deviance,full_deviance,degrees_freedom, &
      statistic,p_value)
      real(dp), intent(in) :: reduced_deviance,full_deviance
      integer, intent(in) :: degrees_freedom
      real(dp), intent(out) :: statistic,p_value
      statistic=max(0.0_dp,reduced_deviance-full_deviance)
      if (degrees_freedom<1) then
         p_value=1.0_dp
      else
         p_value=regularized_gamma_q(0.5_dp*real(degrees_freedom,dp),0.5_dp*statistic)
      end if
   end subroutine likelihood_ratio_test

   subroutine subset_model(y,x,terms,keep,deleted_term,deleted_level,yr,xr,reduced_terms)
      real(dp), intent(in) :: y(:),x(:,:)
      type(random_term_t), intent(in) :: terms(:)
      logical, intent(in) :: keep(:)
      integer, intent(in) :: deleted_term,deleted_level
      real(dp), allocatable, intent(out) :: yr(:),xr(:,:)
      type(random_term_t), allocatable, intent(out) :: reduced_terms(:)
      integer :: k
      yr=pack(y,keep)
      allocate(xr(count(keep),size(x,2)))
      do k=1,size(x,2)
         xr(:,k)=pack(x(:,k),keep)
      end do
      reduced_terms=terms
      do k=1,size(terms)
         if (allocated(reduced_terms(k)%z)) deallocate(reduced_terms(k)%z)
         if (allocated(reduced_terms(k)%group)) deallocate(reduced_terms(k)%group)
         allocate(reduced_terms(k)%z(count(keep),size(terms(k)%z,2)),source=0.0_dp)
         reduced_terms(k)%z=reshape(pack(terms(k)%z,spread(keep,2,size(terms(k)%z,2))), &
            [count(keep),size(terms(k)%z,2)])
         reduced_terms(k)%group=pack(terms(k)%group,keep)
         if (k==deleted_term) then
            where (reduced_terms(k)%group>deleted_level)
               reduced_terms(k)%group=reduced_terms(k)%group-1
            end where
            reduced_terms(k)%n_levels=terms(k)%n_levels-1
         end if
      end do
   end subroutine subset_model

   subroutine drop_column(x,j,xdrop)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: j
      real(dp), allocatable, intent(out) :: xdrop(:,:)
      integer :: p
      p=size(x,2)
      allocate(xdrop(size(x,1),p-1))
      if (j>1) xdrop(:,1:j-1)=x(:,1:j-1)
      if (j<p) xdrop(:,j:p-1)=x(:,j+1:p)
   end subroutine drop_column

   subroutine sort_real(x)
      real(dp), intent(inout) :: x(:)
      integer :: i,j
      real(dp) :: key
      do i=2,size(x)
         key=x(i); j=i-1
         do while (j>=1)
            if (x(j)<=key) exit
            x(j+1)=x(j); j=j-1
         end do
         x(j+1)=key
      end do
   end subroutine sort_real

   real(dp) function quantile_sorted(x,p) result(value)
      real(dp), intent(in) :: x(:),p
      real(dp) :: h,fraction
      integer :: lower_index
      if (size(x)==1) then
         value=x(1)
         return
      end if
      h=1.0_dp+real(size(x)-1,dp)*min(1.0_dp,max(0.0_dp,p))
      lower_index=floor(h)
      fraction=h-real(lower_index,dp)
      if (lower_index>=size(x)) then
         value=x(size(x))
      else
         value=(1.0_dp-fraction)*x(lower_index)+fraction*x(lower_index+1)
      end if
   end function quantile_sorted

   real(dp) function regularized_gamma_q(a,x) result(value)
      real(dp), intent(in) :: a,x
      real(dp) :: ap,del,sumv,b,c,d,h,an
      integer :: n
      if (x<=0.0_dp) then
         value=1.0_dp
      else if (x<a+1.0_dp) then
         ap=a; sumv=1.0_dp/a; del=sumv
         do n=1,200
            ap=ap+1.0_dp
            del=del*x/ap
            sumv=sumv+del
            if (abs(del)<abs(sumv)*1.0e-14_dp) exit
         end do
         value=1.0_dp-sumv*exp(-x+a*log(x)-log_gamma(a))
      else
         b=x+1.0_dp-a
         c=1.0_dp/tiny(1.0_dp)
         d=1.0_dp/b
         h=d
         do n=1,200
            an=-real(n,dp)*(real(n,dp)-a)
            b=b+2.0_dp
            d=an*d+b
            if (abs(d)<tiny(1.0_dp)) d=tiny(1.0_dp)
            c=b+an/c
            if (abs(c)<tiny(1.0_dp)) c=tiny(1.0_dp)
            d=1.0_dp/d
            del=d*c
            h=h*del
            if (abs(del-1.0_dp)<1.0e-14_dp) exit
         end do
         value=exp(-x+a*log(x)-log_gamma(a))*h
      end if
      value=min(1.0_dp,max(0.0_dp,value))
   end function regularized_gamma_q

end module lme4_inference
