! SPDX-License-Identifier: GPL-3.0-only
!
! Experimental modern Fortran translation of computational methods from the
! R package rmgarch, copyright (C) 2008-2025 Alexios Galanos.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 3 only.
module rmgarch_dcc
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use rmgarch_kinds, only : dp
   use rmgarch_math, only : covariance_matrix, outer_product, normalize_covariance, &
      make_positive_definite
   use rmgarch_distributions, only : multivariate_logpdf, random_multivariate
   use rmgarch_optimizer, only : optimizer_result, nelder_mead
   use rmgarch_types, only : dcc_spec, dcc_fit_result, dist_gaussian, dist_student, dist_laplace
   implicit none
   private

   type :: dcc_context
      real(dp), allocatable :: z(:,:)
      integer :: p = 1
      integer :: q = 1
      integer :: g = 0
      integer :: distribution = dist_gaussian
      real(dp) :: shape = 8.0_dp
      logical :: estimate_shape = .false.
   end type dcc_context

   public :: make_dcc_spec, dcc_filter, dcc_log_likelihood
   public :: fit_dcc, fit_dcc11, dcc_forecast, dcc_forecast_history
   public :: simulate_dcc, asymmetric_residuals, dcc_persistence

contains

   function make_dcc_spec(alpha, beta, gamma, distribution, shape) result(spec)
      real(dp), intent(in) :: alpha(:), beta(:)
      real(dp), intent(in), optional :: gamma(:), shape
      integer, intent(in), optional :: distribution
      type(dcc_spec) :: spec

      spec%p = size(alpha)
      spec%q = size(beta)
      allocate(spec%alpha(spec%p),spec%beta(spec%q))
      spec%alpha = alpha
      spec%beta = beta
      if (present(gamma)) then
         spec%g = size(gamma)
         allocate(spec%gamma(spec%g))
         spec%gamma = gamma
      else
         spec%g = 0
         allocate(spec%gamma(0))
      end if
      if (present(distribution)) spec%distribution = distribution
      if (present(shape)) spec%shape = shape
   end function make_dcc_spec

   pure function asymmetric_residuals(z) result(nres)
      real(dp), intent(in) :: z(:,:)
      real(dp) :: nres(size(z,1),size(z,2))
      where (z < 0.0_dp)
         nres = z
      elsewhere
         nres = 0.0_dp
      end where
   end function asymmetric_residuals

   pure function dcc_persistence(spec) result(value)
      type(dcc_spec), intent(in) :: spec
      real(dp) :: value
      value = sum(spec%alpha)+sum(spec%beta)+0.5_dp*sum(spec%gamma)
   end function dcc_persistence

   subroutine dcc_filter(z, spec, q, r, loglikelihoods, qbar_out, nbar_out, valid)
      real(dp), intent(in) :: z(:,:)
      type(dcc_spec), intent(in) :: spec
      real(dp), intent(out) :: q(size(z,2),size(z,2),size(z,1))
      real(dp), intent(out) :: r(size(z,2),size(z,2),size(z,1))
      real(dp), intent(out) :: loglikelihoods(size(z,1))
      real(dp), intent(out), optional :: qbar_out(size(z,2),size(z,2))
      real(dp), intent(out), optional :: nbar_out(size(z,2),size(z,2))
      logical, intent(out), optional :: valid
      real(dp) :: qbar(size(z,2),size(z,2)), nbar(size(z,2),size(z,2))
      real(dp) :: intercept(size(z,2),size(z,2)), nres(size(z,1),size(z,2))
      real(dp) :: suma, sumb, sumg
      integer :: t, i, maxlag, nobs, m
      logical :: ok, step_ok

      nobs = size(z,1)
      m = size(z,2)
      suma = sum(spec%alpha)
      sumb = sum(spec%beta)
      sumg = sum(spec%gamma)
      ok = nobs > 1 .and. m > 0 .and. spec%p == size(spec%alpha) .and. &
         spec%q == size(spec%beta) .and. spec%g == size(spec%gamma)
      ok = ok .and. all_nonnegative(spec%alpha) .and. all_nonnegative(spec%beta) .and. &
         all_nonnegative(spec%gamma) .and. dcc_persistence(spec) < 1.0_dp
      ok = ok .and. any(spec%distribution == [dist_gaussian,dist_student,dist_laplace])
      if (spec%distribution == dist_student) ok = ok .and. spec%shape > 2.0_dp
      q = 0.0_dp
      r = 0.0_dp
      loglikelihoods = -huge(1.0_dp)
      if (.not. ok) then
         if (present(valid)) valid = .false.
         return
      end if

      qbar = make_positive_definite(covariance_matrix(z),1.0e-10_dp)
      nres = asymmetric_residuals(z)
      nbar = make_positive_definite(covariance_matrix(nres),1.0e-12_dp)
      intercept = (1.0_dp-suma-sumb)*qbar-sumg*nbar
      intercept = make_positive_definite(intercept,1.0e-10_dp)
      maxlag = max(1,max(spec%p,max(spec%q,spec%g)))
      do t = 1, min(maxlag,nobs)
         q(:,:,t) = qbar
         r(:,:,t) = normalize_covariance(qbar)
         loglikelihoods(t) = dcc_observation_loglik(z(t,:),r(:,:,t),spec%distribution,spec%shape,step_ok)
         ok = ok .and. step_ok
      end do
      do t = maxlag+1, nobs
         q(:,:,t) = intercept
         do i = 1, spec%p
            q(:,:,t) = q(:,:,t)+spec%alpha(i)*outer_product(z(t-i,:),z(t-i,:))
         end do
         do i = 1, spec%g
            q(:,:,t) = q(:,:,t)+spec%gamma(i)*outer_product(nres(t-i,:),nres(t-i,:))
         end do
         do i = 1, spec%q
            q(:,:,t) = q(:,:,t)+spec%beta(i)*q(:,:,t-i)
         end do
         q(:,:,t) = make_positive_definite(q(:,:,t),1.0e-10_dp)
         r(:,:,t) = normalize_covariance(q(:,:,t))
         loglikelihoods(t) = dcc_observation_loglik(z(t,:),r(:,:,t),spec%distribution,spec%shape,step_ok)
         ok = ok .and. step_ok
      end do
      if (present(qbar_out)) qbar_out = qbar
      if (present(nbar_out)) nbar_out = nbar
      if (present(valid)) valid = ok
   end subroutine dcc_filter

   function dcc_log_likelihood(z, spec, valid) result(value)
      real(dp), intent(in) :: z(:,:)
      type(dcc_spec), intent(in) :: spec
      logical, intent(out), optional :: valid
      real(dp) :: value
      real(dp), allocatable :: q(:,:,:), r(:,:,:), ll(:)
      logical :: ok

      allocate(q(size(z,2),size(z,2),size(z,1)),r(size(z,2),size(z,2),size(z,1)),ll(size(z,1)))
      call dcc_filter(z,spec,q,r,ll,valid=ok)
      if (ok) then
         value = sum(ll)
      else
         value = -huge(1.0_dp)
      end if
      if (present(valid)) valid = ok
   end function dcc_log_likelihood

   function fit_dcc(z, p, q, g, distribution, shape, estimate_shape, max_iterations) result(fit)
      real(dp), intent(in) :: z(:,:)
      integer, intent(in), optional :: p, q, g, distribution, max_iterations
      real(dp), intent(in), optional :: shape
      logical, intent(in), optional :: estimate_shape
      type(dcc_fit_result) :: fit
      type(dcc_context) :: context
      type(optimizer_result) :: opt
      real(dp), allocatable :: x0(:), weights(:), alpha(:), beta(:), gamma(:)
      real(dp) :: remaining, initial_shape
      integer :: maxit, ncoef, k
      logical :: ok

      if (present(p)) context%p = max(0,p)
      if (present(q)) context%q = max(0,q)
      if (present(g)) context%g = max(0,g)
      if (present(distribution)) context%distribution = distribution
      if (present(shape)) context%shape = shape
      if (present(estimate_shape)) context%estimate_shape = estimate_shape
      if (context%distribution /= dist_student) context%estimate_shape = .false.
      allocate(context%z(size(z,1),size(z,2)))
      context%z = z
      ncoef = context%p+context%q+context%g
      if (ncoef == 0 .or. size(z,1) <= max(1,max(context%p,max(context%q,context%g)))) then
         fit%status = 3
         fit%message = 'invalid DCC order or insufficient observations'
         return
      end if
      allocate(x0(ncoef+merge(1,0,context%estimate_shape)),weights(ncoef))
      weights = 0.0_dp
      if (context%p > 0) weights(1:context%p) = 0.05_dp/real(context%p,dp)
      if (context%q > 0) weights(context%p+1:context%p+context%q) = 0.90_dp/real(context%q,dp)
      if (context%g > 0) weights(context%p+context%q+1:ncoef) = 0.02_dp/real(context%g,dp)
      if (sum(weights) >= 0.98_dp) weights = 0.95_dp*weights/sum(weights)
      remaining = max(1.0e-4_dp,1.0_dp-sum(weights))
      x0(1:ncoef) = log(max(weights,1.0e-6_dp)/remaining)
      if (context%estimate_shape) then
         initial_shape = min(49.5_dp,max(2.1_dp,context%shape))
         x0(ncoef+1) = log((initial_shape-2.01_dp)/(50.0_dp-initial_shape))
      end if
      maxit = 800
      if (present(max_iterations)) maxit = max_iterations
      opt = nelder_mead(dcc_objective_general,x0,context,step=0.18_dp,tolerance=1.0e-8_dp, &
         max_iterations=maxit)
      call decode_general_parameters(opt%x,context,alpha,beta,gamma,initial_shape)
      fit%spec = make_dcc_spec(alpha,beta,gamma,context%distribution,initial_shape)
      allocate(fit%qbar(size(z,2),size(z,2)),fit%nbar(size(z,2),size(z,2)))
      allocate(fit%q(size(z,2),size(z,2),size(z,1)),fit%r(size(z,2),size(z,2),size(z,1)))
      allocate(fit%loglikelihoods(size(z,1)),fit%standardized_residuals(size(z,1),size(z,2)))
      fit%standardized_residuals = z
      call dcc_filter(z,fit%spec,fit%q,fit%r,fit%loglikelihoods,fit%qbar,fit%nbar,ok)
      if (ok) then
         fit%log_likelihood = sum(fit%loglikelihoods)
      else
         fit%log_likelihood = -huge(1.0_dp)
      end if
      fit%iterations = opt%iterations
      fit%status = merge(0,2,ok)
      if (opt%status /= 0 .and. fit%status == 0) fit%status = opt%status
      if (fit%status == 0) then
         fit%message = 'converged'
      else if (ok) then
         fit%message = 'valid fit; optimizer stopped before convergence'
      else
         fit%message = 'invalid fitted recursion'
      end if
      k = size(opt%x)
      if (ok) then
         fit%aic = -2.0_dp*fit%log_likelihood+2.0_dp*real(k,dp)
         fit%bic = -2.0_dp*fit%log_likelihood+log(real(size(z,1),dp))*real(k,dp)
      end if
   end function fit_dcc

   function fit_dcc11(z, asymmetric, distribution, shape, max_iterations) result(fit)
      real(dp), intent(in) :: z(:,:)
      logical, intent(in), optional :: asymmetric
      integer, intent(in), optional :: distribution, max_iterations
      real(dp), intent(in), optional :: shape
      type(dcc_fit_result) :: fit
      logical :: use_asym
      integer :: dist, maxit
      real(dp) :: nu

      use_asym = .false.
      if (present(asymmetric)) use_asym = asymmetric
      dist = dist_gaussian
      if (present(distribution)) dist = distribution
      nu = 8.0_dp
      if (present(shape)) nu = shape
      maxit = 600
      if (present(max_iterations)) maxit = max_iterations
      fit = fit_dcc(z,p=1,q=1,g=merge(1,0,use_asym),distribution=dist,shape=nu, &
         max_iterations=maxit)
   end function fit_dcc11

   subroutine dcc_forecast(spec, qbar, nbar, last_q, last_z, horizons, q_forecast, r_forecast)
      type(dcc_spec), intent(in) :: spec
      real(dp), intent(in) :: qbar(:,:), nbar(:,:), last_q(:,:), last_z(:)
      integer, intent(in) :: horizons
      real(dp), intent(out) :: q_forecast(size(qbar,1),size(qbar,2),horizons)
      real(dp), intent(out) :: r_forecast(size(qbar,1),size(qbar,2),horizons)
      real(dp) :: intercept(size(qbar,1),size(qbar,2)), last_n(size(last_z))
      real(dp) :: suma, sumb, sumg
      integer :: h

      suma = sum(spec%alpha)
      sumb = sum(spec%beta)
      sumg = sum(spec%gamma)
      intercept = (1.0_dp-suma-sumb)*qbar-sumg*nbar
      where (last_z < 0.0_dp)
         last_n = last_z
      elsewhere
         last_n = 0.0_dp
      end where
      if (horizons <= 0) return
      q_forecast(:,:,1) = intercept+suma*outer_product(last_z,last_z)+ &
         sumg*outer_product(last_n,last_n)+sumb*last_q
      q_forecast(:,:,1) = make_positive_definite(q_forecast(:,:,1))
      r_forecast(:,:,1) = normalize_covariance(q_forecast(:,:,1))
      do h = 2, horizons
         q_forecast(:,:,h) = intercept+suma*qbar+sumg*nbar+sumb*q_forecast(:,:,h-1)
         q_forecast(:,:,h) = make_positive_definite(q_forecast(:,:,h))
         r_forecast(:,:,h) = normalize_covariance(q_forecast(:,:,h))
      end do
   end subroutine dcc_forecast

   subroutine dcc_forecast_history(spec, qbar, nbar, q_history, z_history, horizons, q_forecast, r_forecast)
      type(dcc_spec), intent(in) :: spec
      real(dp), intent(in) :: qbar(:,:), nbar(:,:), q_history(:,:,:), z_history(:,:)
      integer, intent(in) :: horizons
      real(dp), intent(out) :: q_forecast(size(qbar,1),size(qbar,2),horizons)
      real(dp), intent(out) :: r_forecast(size(qbar,1),size(qbar,2),horizons)
      real(dp) :: intercept(size(qbar,1),size(qbar,2)), nvec(size(qbar,1))
      integer :: h, i, source, nq, nz

      nq = size(q_history,3)
      nz = size(z_history,1)
      intercept = (1.0_dp-sum(spec%alpha)-sum(spec%beta))*qbar-sum(spec%gamma)*nbar
      do h = 1, horizons
         q_forecast(:,:,h) = intercept
         do i = 1, spec%p
            source = nz+h-i
            if (source <= nz .and. source >= 1) then
               q_forecast(:,:,h) = q_forecast(:,:,h)+spec%alpha(i)* &
                  outer_product(z_history(source,:),z_history(source,:))
            else
               q_forecast(:,:,h) = q_forecast(:,:,h)+spec%alpha(i)*qbar
            end if
         end do
         do i = 1, spec%g
            source = nz+h-i
            if (source <= nz .and. source >= 1) then
               where (z_history(source,:) < 0.0_dp)
                  nvec = z_history(source,:)
               elsewhere
                  nvec = 0.0_dp
               end where
               q_forecast(:,:,h) = q_forecast(:,:,h)+spec%gamma(i)*outer_product(nvec,nvec)
            else
               q_forecast(:,:,h) = q_forecast(:,:,h)+spec%gamma(i)*nbar
            end if
         end do
         do i = 1, spec%q
            source = nq+h-i
            if (source <= nq .and. source >= 1) then
               q_forecast(:,:,h) = q_forecast(:,:,h)+spec%beta(i)*q_history(:,:,source)
            else if (source > nq) then
               q_forecast(:,:,h) = q_forecast(:,:,h)+spec%beta(i)*q_forecast(:,:,source-nq)
            else
               q_forecast(:,:,h) = q_forecast(:,:,h)+spec%beta(i)*qbar
            end if
         end do
         q_forecast(:,:,h) = make_positive_definite(q_forecast(:,:,h),1.0e-10_dp)
         r_forecast(:,:,h) = normalize_covariance(q_forecast(:,:,h))
      end do
   end subroutine dcc_forecast_history

   subroutine simulate_dcc(nobs, spec, qbar, z, q, r, burn)
      integer, intent(in) :: nobs
      type(dcc_spec), intent(in) :: spec
      real(dp), intent(in) :: qbar(:,:)
      real(dp), intent(out) :: z(nobs,size(qbar,1))
      real(dp), intent(out) :: q(size(qbar,1),size(qbar,2),nobs)
      real(dp), intent(out) :: r(size(qbar,1),size(qbar,2),nobs)
      integer, intent(in), optional :: burn
      integer :: nburn, ntotal, t, m, i
      real(dp), allocatable :: zall(:,:), qall(:,:,:), rall(:,:,:), nres(:,:), nbar(:,:)
      real(dp) :: intercept(size(qbar,1),size(qbar,2)), zero_mean(size(qbar,1))
      real(dp) :: draw(size(qbar,1))
      real(dp) :: suma, sumb, sumg
      logical :: ok

      nburn = 300
      if (present(burn)) nburn = max(0,burn)
      ntotal = nobs+nburn
      m = size(qbar,1)
      allocate(zall(ntotal,m),qall(m,m,ntotal),rall(m,m,ntotal),nres(ntotal,m),nbar(m,m))
      zall = 0.0_dp
      nres = 0.0_dp
      zero_mean = 0.0_dp
      nbar = 0.5_dp*qbar
      suma = sum(spec%alpha)
      sumb = sum(spec%beta)
      sumg = sum(spec%gamma)
      intercept = (1.0_dp-suma-sumb)*qbar-sumg*nbar
      qall(:,:,1) = make_positive_definite(qbar)
      rall(:,:,1) = normalize_covariance(qall(:,:,1))
      do t = 1, ntotal
         if (t > 1) then
            qall(:,:,t) = intercept
            do i = 1, spec%p
               if (t-i >= 1) then
                  qall(:,:,t) = qall(:,:,t)+spec%alpha(i)*outer_product(zall(t-i,:),zall(t-i,:))
               else
                  qall(:,:,t) = qall(:,:,t)+spec%alpha(i)*qbar
               end if
            end do
            do i = 1, spec%g
               if (t-i >= 1) then
                  qall(:,:,t) = qall(:,:,t)+spec%gamma(i)*outer_product(nres(t-i,:),nres(t-i,:))
               else
                  qall(:,:,t) = qall(:,:,t)+spec%gamma(i)*nbar
               end if
            end do
            do i = 1, spec%q
               if (t-i >= 1) then
                  qall(:,:,t) = qall(:,:,t)+spec%beta(i)*qall(:,:,t-i)
               else
                  qall(:,:,t) = qall(:,:,t)+spec%beta(i)*qbar
               end if
            end do
            qall(:,:,t) = make_positive_definite(qall(:,:,t),1.0e-10_dp)
            rall(:,:,t) = normalize_covariance(qall(:,:,t))
         end if
         call random_multivariate(spec%distribution,zero_mean,rall(:,:,t),draw, &
            shape=spec%shape,valid=ok)
         if (ok) then
            zall(t,:) = draw
         else
            zall(t,:) = 0.0_dp
         end if
         where (zall(t,:) < 0.0_dp)
            nres(t,:) = zall(t,:)
         elsewhere
            nres(t,:) = 0.0_dp
         end where
      end do
      z = zall(nburn+1:ntotal,:)
      q = qall(:,:,nburn+1:ntotal)
      r = rall(:,:,nburn+1:ntotal)
   end subroutine simulate_dcc

   function dcc_observation_loglik(z, r, distribution, shape, ok) result(ll)
      real(dp), intent(in) :: z(:), r(:,:), shape
      integer, intent(in) :: distribution
      logical, intent(out) :: ok
      real(dp) :: ll, zero_mean(size(z))
      zero_mean = 0.0_dp
      ll = multivariate_logpdf(distribution,z,zero_mean,r,shape,ok)
      ok = ok .and. ieee_is_finite(ll)
   end function dcc_observation_loglik

   function dcc_objective_general(x, generic_context) result(value)
      real(dp), intent(in) :: x(:)
      class(*), intent(in) :: generic_context
      real(dp) :: value, shape, ll
      real(dp), allocatable :: alpha(:), beta(:), gamma(:)
      type(dcc_spec) :: spec
      logical :: ok

      select type (context => generic_context)
      type is (dcc_context)
         call decode_general_parameters(x,context,alpha,beta,gamma,shape)
         spec = make_dcc_spec(alpha,beta,gamma,context%distribution,shape)
         ll = dcc_log_likelihood(context%z,spec,ok)
         if (ok) then
            value = -ll
         else
            value = huge(1.0_dp)/100.0_dp
         end if
      class default
         value = huge(1.0_dp)/100.0_dp
      end select
   end function dcc_objective_general

   subroutine decode_general_parameters(x, context, alpha, beta, gamma, shape)
      real(dp), intent(in) :: x(:)
      type(dcc_context), intent(in) :: context
      real(dp), allocatable, intent(out) :: alpha(:), beta(:), gamma(:)
      real(dp), intent(out) :: shape
      real(dp), allocatable :: ex(:), weights(:)
      real(dp) :: denom, logistic
      integer :: ncoef, offset

      ncoef = context%p+context%q+context%g
      allocate(ex(ncoef),weights(ncoef),alpha(context%p),beta(context%q),gamma(context%g))
      ex = exp(max(-30.0_dp,min(30.0_dp,x(1:ncoef))))
      denom = 1.0_dp+sum(ex)
      weights = 0.999_dp*ex/denom
      offset = 0
      if (context%p > 0) alpha = weights(1:context%p)
      offset = context%p
      if (context%q > 0) beta = weights(offset+1:offset+context%q)
      offset = offset+context%q
      if (context%g > 0) gamma = weights(offset+1:offset+context%g)
      shape = context%shape
      if (context%estimate_shape) then
         logistic = 1.0_dp/(1.0_dp+exp(-max(-30.0_dp,min(30.0_dp,x(ncoef+1)))))
         shape = 2.01_dp+47.99_dp*logistic
      end if
   end subroutine decode_general_parameters

   pure function all_nonnegative(x) result(ok)
      real(dp), intent(in) :: x(:)
      logical :: ok
      if (size(x) == 0) then
         ok = .true.
      else
         ok = minval(x) >= 0.0_dp
      end if
   end function all_nonnegative

end module rmgarch_dcc
