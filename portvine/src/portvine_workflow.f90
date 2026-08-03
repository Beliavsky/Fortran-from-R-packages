! SPDX-License-Identifier: GPL-3.0-only
module portvine_workflow
   use portvine_kinds, only : dp
   use portvine_types
   use portvine_stats, only : estimate_risk_measure
   use portvine_marginals, only : fit_rolling_marginals, get_marginal_point
   use portvine_dependence, only : fit_vine_windows
   use portvine_conditional, only : conditional_dvine_sample
   use rugarch, only : garch_spec, distribution_quantile
   implicit none
   private
   public :: estimate_risk_roll

contains

   subroutine estimate_risk_roll(returns, marginal_settings, vine_settings, alpha, &
      measures, n_samples, result, weights, cond_indices, cond_u, n_mc_samples, &
      prior_residual_strategy, status)
      real(dp), intent(in) :: returns(:,:)
      type(marginal_settings_type), intent(in) :: marginal_settings
      type(vine_settings_type), intent(in) :: vine_settings
      real(dp), intent(in) :: alpha(:)
      integer, intent(in) :: measures(:), n_samples
      type(portvine_roll_result), intent(out) :: result
      real(dp), intent(in), optional :: weights(:,:)
      integer, intent(in), optional :: cond_indices(:), n_mc_samples
      real(dp), intent(in), optional :: cond_u(:)
      logical, intent(in), optional :: prior_residual_strategy
      integer, intent(out), optional :: status
      integer :: d, nobs, nf, nv, nm, na, nc, ncase, t, f, v, mw, m, c, istat, mc
      integer :: cond_source_index, a, offset
      logical :: conditional, use_prior
      real(dp), allocatable :: w(:,:), u_sample(:,:), portfolio(:), all_portfolio(:)
      real(dp), allocatable :: cond_values(:), risk_value(:), local_cond_u(:)

      result = portvine_roll_result()
      if (present(status)) status = portvine_success
      d = size(returns,1)
      nobs = size(returns,2)
      nm = size(measures)
      na = size(alpha)
      conditional = present(cond_indices)
      nc = 0
      if (conditional) nc = size(cond_indices)
      use_prior = .false.
      if (present(prior_residual_strategy)) use_prior = prior_residual_strategy
      mc = 1000
      if (present(n_mc_samples)) mc = max(1,n_mc_samples)

      if (.not.valid_inputs()) then
         result%status = portvine_invalid_input
         result%message = 'invalid dimensions, window settings, probabilities, or conditioning indices'
         if (present(status)) status = result%status
         return
      end if
      nf = nobs-marginal_settings%train_size
      nv = (nf+vine_settings%refit_size-1)/vine_settings%refit_size
      allocate(w(d,nv))
      call prepare_weights(w,istat)
      if (istat /= portvine_success) then
         result%status = istat
         result%message = 'weights must be nonnegative and conditioning assets must have zero weight'
         if (present(status)) status = result%status
         return
      end if

      call fit_rolling_marginals(returns,marginal_settings,vine_settings%train_size, &
         result%marginal,istat)
      if (istat /= portvine_success) then
         result%status = istat
         result%message = 'one or more rolling marginal fits failed'
         if (present(status)) status = result%status
         return
      end if
      if (conditional) then
         call fit_vine_windows(result%marginal,marginal_settings,vine_settings,nobs, &
            cond_indices,result%dvine,result%cvine,istat)
      else
         call fit_vine_windows(result%marginal,marginal_settings,vine_settings,nobs, &
            dvines=result%dvine,cvines=result%cvine,status=istat)
      end if
      if (istat /= portvine_success) then
         result%status = istat
         result%message = 'vine fitting failed'
         if (present(status)) status = result%status
         return
      end if

      result%vine_type = vine_settings%vine_type
      allocate(result%overall(nm,na,nf),result%realized(nf),result%alpha(na), &
         result%row_index(nf),result%vine_window(nf),result%measure(nm), &
         result%weights(d,nv))
      result%overall = 0.0_dp
      result%alpha = alpha
      result%measure = measures
      result%weights = w
      if (conditional) then
         ncase = size(cond_u)+1
         allocate(result%conditional(nm,na,ncase,nf),result%condition_level(ncase), &
            result%conditional_value(nc,ncase,nf))
         result%conditional = 0.0_dp
         result%condition_level(1:size(cond_u)) = cond_u
         result%condition_level(ncase) = -1.0_dp
         result%conditional_value = 0.0_dp
         allocate(all_portfolio(n_samples*ncase),local_cond_u(nc),cond_values(nc))
      else
         ncase = 0
         allocate(result%conditional(nm,na,0,nf),result%condition_level(0), &
            result%conditional_value(0,0,nf))
         allocate(all_portfolio(n_samples))
      end if
      allocate(u_sample(d,n_samples),portfolio(n_samples),risk_value(na))

      do f = 1, nf
         t = marginal_settings%train_size+f
         v = min(nv,(f-1)/vine_settings%refit_size+1)
         mw = (v*vine_settings%refit_size+marginal_settings%refit_size-1)/ &
            marginal_settings%refit_size
         mw = max(1,min(size(result%marginal(1)%window),mw))
         result%row_index(f) = t
         result%vine_window(f) = v
         result%realized(f) = sum(w(:,v)*returns(:,t))

         if (.not.conditional) then
            if (vine_settings%vine_type == vine_dvine) then
               call result%dvine(v)%simulate(n_samples,u_sample)
            else
               call result%cvine(v)%simulate(n_samples,u_sample)
            end if
            call portfolio_from_uniforms(u_sample,t,mw,w(:,v),portfolio)
            all_portfolio = portfolio
         else
            offset = 0
            do c = 1, size(cond_u)
               local_cond_u = cond_u(c)
               call conditional_dvine_sample(result%dvine(v),n_samples,local_cond_u,.true., &
                  u_sample,istat)
               call portfolio_from_uniforms(u_sample,t,mw,w(:,v),portfolio,cond_indices,cond_values)
               result%conditional_value(:,c,f) = cond_values
               all_portfolio(offset+1:offset+n_samples) = portfolio
               do m = 1, nm
                  call estimate_risk_measure(portfolio,alpha,measures(m),risk_value,mc)
                  result%conditional(m,:,c,f) = risk_value
               end do
               offset = offset+n_samples
            end do
            cond_source_index = t
            if (use_prior) cond_source_index = max(1,t-1)
            do a = 1, nc
               call get_actual_uniform(cond_indices(a),mw,cond_source_index,local_cond_u(a),istat)
               if (istat /= portvine_success) then
                  result%status = istat
                  result%message = 'conditioning residual could not be located'
                  if (present(status)) status = istat
                  return
               end if
            end do
            call conditional_dvine_sample(result%dvine(v),n_samples,local_cond_u,.false., &
               u_sample,istat)
            call portfolio_from_uniforms(u_sample,t,mw,w(:,v),portfolio,cond_indices,cond_values)
            result%conditional_value(:,ncase,f) = cond_values
            all_portfolio(offset+1:offset+n_samples) = portfolio
            do m = 1, nm
               call estimate_risk_measure(portfolio,alpha,measures(m),risk_value,mc)
               result%conditional(m,:,ncase,f) = risk_value
            end do
         end if
         do m = 1, nm
            call estimate_risk_measure(all_portfolio,alpha,measures(m),risk_value,mc)
            result%overall(m,:,f) = risk_value
         end do
      end do
      result%status = portvine_success
      result%message = 'ok'
      if (present(status)) status = result%status

   contains

      logical function valid_inputs() result(ok)
         integer :: i
         ok = d >= 2 .and. nobs > marginal_settings%train_size .and. &
            marginal_settings%train_size >= 30 .and. marginal_settings%refit_size >= 1 .and. &
            vine_settings%train_size <= marginal_settings%train_size .and. &
            vine_settings%refit_size <= marginal_settings%refit_size .and. &
            mod(marginal_settings%refit_size,vine_settings%refit_size) == 0 .and. &
            n_samples >= 1 .and. nm >= 1 .and. na >= 1 .and. &
            all(alpha > 0.0_dp) .and. all(alpha < 1.0_dp)
         if (.not.ok) return
         do i = 1, nm
            if (measures(i) < risk_var .or. measures(i) > risk_es_mc) ok = .false.
         end do
         if (conditional) then
            ok = ok .and. nc >= 1 .and. nc <= 2 .and. &
               vine_settings%vine_type == vine_dvine .and. present(cond_u)
            if (.not.ok) return
            ok = size(cond_u) >= 1 .and. all(cond_u > 0.0_dp) .and. all(cond_u < 1.0_dp)
            do i = 1, nc
               if (cond_indices(i) < 1 .or. cond_indices(i) > d) ok = .false.
               if (count(cond_indices == cond_indices(i)) > 1) ok = .false.
            end do
         end if
      end function valid_inputs

      subroutine prepare_weights(output,istat)
         real(dp), intent(out) :: output(:,:)
         integer, intent(out) :: istat
         integer :: i
         istat = portvine_success
         if (present(weights)) then
            if (size(weights,1) == d .and. size(weights,2) == nv) then
               output = weights
            else if (size(weights,1) == nv .and. size(weights,2) == d) then
               output = transpose(weights)
            else
               istat = portvine_invalid_input
               return
            end if
         else
            output = 1.0_dp
            if (conditional) then
               do i = 1, nc
                  output(cond_indices(i),:) = 0.0_dp
               end do
            end if
         end if
         if (any(output < 0.0_dp)) istat = portvine_invalid_input
         if (conditional) then
            do i = 1, nc
               if (any(abs(output(cond_indices(i),:)) > 1.0e-14_dp)) istat = portvine_invalid_input
            end do
         end if
      end subroutine prepare_weights

      subroutine portfolio_from_uniforms(u,t,mw,weight,port,condition_indices,condition_values)
         real(dp), intent(in) :: u(:,:), weight(:)
         integer, intent(in) :: t, mw
         real(dp), intent(out) :: port(size(u,2))
         integer, intent(in), optional :: condition_indices(:)
         real(dp), intent(out), optional :: condition_values(:)
         type(garch_spec) :: spec
         real(dp) :: mu, sigma, actual_u, value
         integer :: aa, ss, ii, st
         port = 0.0_dp
         if (present(condition_values)) condition_values = 0.0_dp
         do aa = 1, d
            call get_marginal_point(result%marginal(aa),mw,t,mu,sigma,actual_u,spec,st)
            do ss = 1, size(u,2)
               value = mu+sigma*distribution_quantile(u(aa,ss),spec%cond_dist, &
                  spec%shape,spec%skew,spec%lambda)
               port(ss) = port(ss)+weight(aa)*value
               if (present(condition_indices) .and. present(condition_values)) then
                  do ii = 1, size(condition_indices)
                     if (aa == condition_indices(ii) .and. ss == 1) condition_values(ii) = value
                  end do
               end if
            end do
         end do
      end subroutine portfolio_from_uniforms

      subroutine get_actual_uniform(asset_index,mw,index,u,istat)
         integer, intent(in) :: asset_index,mw,index
         real(dp), intent(out) :: u
         integer, intent(out) :: istat
         real(dp) :: mu,sigma
         call get_marginal_point(result%marginal(asset_index),mw,index,mu,sigma,u,status=istat)
      end subroutine get_actual_uniform

   end subroutine estimate_risk_roll

end module portvine_workflow
