! SPDX-License-Identifier: GPL-3.0-only
module fitheavytail_mvt
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
   use fitheavytail_kinds, only: dp
   use fitheavytail_status, only: ht_success, ht_invalid_argument, &
      ht_too_few_observations, ht_singular_matrix, ht_no_convergence
   use fitheavytail_types, only: heavy_tail_fit, clear_fit
   use fitheavytail_linalg, only: clean_complete_rows, column_mean, &
      sample_covariance, quadratic_forms, weighted_covariance, &
      symmetric_eigen, factor_decomposition, identity_matrix, &
      inverse_matrix, logdet_spd, outer_product, symmetrize, trace_matrix
   use fitheavytail_special, only: digamma_dp
   use fitheavytail_tail, only: default_nu_min, default_nu_max, cap_nu, &
      nu_from_average_marginal_kurtosis, nu_from_cross_cumulants, &
      nu_from_all_cumulants, nu_hill_estimator, nu_opp_estimator, &
      nu_pop_estimator
   use fitheavytail_rng, only: rng_state, seed_rng, uniform_random
   implicit none
   private

   public :: fit_mvt, mvt_log_likelihood, nu_mle

contains

   subroutine fit_mvt(x, result, fixed_nu, nu_method, &
      nu_iterative_method, initial_nu, initial_mu, initial_scatter, &
      initial_covariance, na_rm, optimize_mu, observation_weights, &
      scale_covmat, px_em_acceleration, nu_update_start_at_iter, &
      nu_update_every_num_iter, factors, max_iter, ptol, ftol, &
      nu_min, nu_max)
      real(dp), intent(in) :: x(:,:)
      type(heavy_tail_fit), intent(out) :: result
      real(dp), intent(in), optional :: fixed_nu, initial_nu
      character(len=*), intent(in), optional :: nu_method
      character(len=*), intent(in), optional :: nu_iterative_method
      real(dp), intent(in), optional :: initial_mu(:)
      real(dp), intent(in), optional :: initial_scatter(:,:)
      real(dp), intent(in), optional :: initial_covariance(:,:)
      logical, intent(in), optional :: na_rm, optimize_mu
      real(dp), intent(in), optional :: observation_weights(:)
      logical, intent(in), optional :: scale_covmat
      logical, intent(in), optional :: px_em_acceleration
      integer, intent(in), optional :: nu_update_start_at_iter
      integer, intent(in), optional :: nu_update_every_num_iter
      integer, intent(in), optional :: factors, max_iter
      real(dp), intent(in), optional :: ptol, ftol, nu_min, nu_max

      real(dp), allocatable :: data(:,:), complete(:,:), work(:,:)
      real(dp), allocatable :: mu(:), old_mu(:), sigma(:,:), old_sigma(:,:)
      real(dp), allocatable :: covariance(:,:), xc(:,:), q(:), etau(:)
      real(dp), allocatable :: obs_weights(:), b(:,:), psi(:), eigval(:)
      real(dp), allocatable :: eigvec(:,:), s(:,:), ave_etau_x(:)
      real(dp), allocatable :: ave_etau_xx(:,:)
      real(dp) :: nu, old_nu, lo, hi, p_tol, f_tol, alpha
      real(dp) :: ave_etau, ave_elogtau, ll, old_ll, start_time, end_time
      real(dp) :: kappa, gamma_shape, nmse
      integer :: n, t, f, niter, iter, status, nu_start, nu_every
      logical :: drop_na, do_mu, do_scale, do_px, optimize_nu
      logical :: has_missing, factor_structure, params_ok, fun_ok
      logical :: converged_flag
      character(len=32) :: nmethod, iterative_method

      call clear_fit(result)
      n = size(x,2)
      if (n < 1) then
         call set_error(result, ht_invalid_argument, 'X has no columns')
         return
      end if
      drop_na = .true.
      if (present(na_rm)) drop_na = na_rm
      if (n == 1) drop_na = .true.
      if (drop_na) then
         call clean_complete_rows(x,data)
      else
         allocate(data(size(x,1),n))
         data = x
      end if
      t = size(data,1)
      if (t <= n) then
         call set_error(result, ht_too_few_observations, &
            'fit_mvt requires T > N after handling NaNs')
         return
      end if
      has_missing = any_nan(data)

      lo = default_nu_min
      hi = default_nu_max
      if (present(nu_min)) lo = nu_min
      if (present(nu_max)) hi = nu_max
      if (lo <= 2.0_dp .or. hi <= lo) then
         call set_error(result, ht_invalid_argument, 'invalid nu bounds')
         return
      end if
      niter = 100
      if (present(max_iter)) niter = max_iter
      p_tol = 1.0e-3_dp
      if (present(ptol)) p_tol = ptol
      f_tol = huge(1.0_dp)
      if (present(ftol)) f_tol = ftol
      do_mu = .true.
      if (present(optimize_mu)) do_mu = optimize_mu
      do_scale = .false.
      if (present(scale_covmat)) do_scale = scale_covmat
      do_px = .true.
      if (present(px_em_acceleration)) do_px = px_em_acceleration
      nu_start = 1
      if (present(nu_update_start_at_iter)) then
         nu_start = max(1,nu_update_start_at_iter)
      end if
      nu_every = 1
      if (present(nu_update_every_num_iter)) then
         nu_every = max(1,nu_update_every_num_iter)
      end if
      f = n
      if (present(factors)) f = factors
      if (f < 1 .or. f > n .or. niter < 1) then
         call set_error(result, ht_invalid_argument, &
            'factors and max_iter are outside their valid ranges')
         return
      end if
      factor_structure = f /= n

      nmethod = 'iterative'
      if (present(nu_method)) nmethod = adjustl(nu_method)
      optimize_nu = .not.present(fixed_nu) .and. &
         trim(nmethod) == 'iterative'
      if (present(fixed_nu)) then
         if (fixed_nu <= 2.0_dp) then
            call set_error(result, ht_invalid_argument, &
               'fixed_nu must exceed two')
            return
         end if
         nu = fixed_nu
         if (nu >= huge(1.0_dp)/100.0_dp) nu = 1.0e15_dp
      else if (optimize_nu) then
         nu = 4.0_dp
         if (present(initial_nu)) nu = initial_nu
         nu = cap_nu(nu,lo,hi)
      else
         call clean_complete_rows(data,complete)
         select case (trim(nmethod))
         case ('kurtosis')
            nu = nu_from_average_marginal_kurtosis(complete, &
               nu_min=lo,nu_max=hi)
         case ('cross-cumulants')
            nu = nu_from_cross_cumulants(complete,lo,hi)
         case ('all-cumulants')
            nu = nu_from_all_cumulants(complete,lo,hi)
         case ('Hill')
            nu = nu_hill_estimator(complete,lo,hi)
         case ('MLE-diag')
            nu = nu_mle(complete,'MLE-mv-diagcov', &
               nu_min=lo,nu_max=hi)
         case ('MLE-diag-resampled')
            nu = nu_mle(complete,'MLE-mv-diagcov-resampled', &
               nu_min=lo,nu_max=hi)
         case default
            call set_error(result,ht_invalid_argument, &
               'unknown nu_method')
            return
         end select
      end if
      iterative_method = 'POP'
      if (present(nu_iterative_method)) then
         iterative_method = adjustl(nu_iterative_method)
      end if

      allocate(work(t,n))
      work = data
      call initialize_missing(work)
      allocate(mu(n),old_mu(n),sigma(n,n),old_sigma(n,n), &
         covariance(n,n),xc(t,n),q(t),etau(t),obs_weights(t), &
         ave_etau_x(n),ave_etau_xx(n,n),s(n,n))
      if (present(initial_mu)) then
         if (size(initial_mu) /= n) then
            call set_error(result,ht_invalid_argument, &
               'initial_mu has the wrong size')
            return
         end if
         mu = initial_mu
      else if (do_mu) then
         mu = column_mean_nan(data)
      else
         mu = 0.0_dp
      end if
      if (present(initial_scatter)) then
         if (any(shape(initial_scatter) /= [n,n])) then
            call set_error(result,ht_invalid_argument, &
               'initial_scatter has the wrong shape')
            return
         end if
         sigma = initial_scatter
      else if (present(initial_covariance)) then
         if (any(shape(initial_covariance) /= [n,n])) then
            call set_error(result,ht_invalid_argument, &
               'initial_covariance has the wrong shape')
            return
         end if
         sigma = (max(nu,2.1_dp)-2.0_dp)/max(nu,2.1_dp) * &
            initial_covariance
      else
         sigma = (max(nu,2.1_dp)-2.0_dp)/max(nu,2.1_dp) * &
            sample_covariance(work)
      end if
      if (.not.do_mu) mu = 0.0_dp
      if (present(observation_weights)) then
         if (size(observation_weights) /= t .or. &
            any(observation_weights < 0.0_dp) .or. &
            sum(observation_weights) <= 0.0_dp) then
            call set_error(result,ht_invalid_argument, &
               'observation_weights are invalid')
            return
         end if
         obs_weights = observation_weights * real(t,dp) / &
            sum(observation_weights)
      else
         obs_weights = 1.0_dp
      end if

      if (factor_structure) then
         allocate(b(n,f),psi(n),eigval(n),eigvec(n,n))
         call symmetric_eigen(sigma,eigval,eigvec,status)
         if (status /= ht_success) then
            call set_error(result,status,'factor initialization failed')
            return
         end if
         b = 0.0_dp
         do iter=1,f
            b(:,iter) = eigvec(:,iter)*sqrt(max(eigval(iter),0.0_dp))
         end do
         psi = 0.0_dp
         s = matmul(b,transpose(b))
         do iter=1,n
            psi(iter) = max(1.0e-10_dp,sigma(iter,iter)-s(iter,iter))
         end do
         sigma = s
         do iter=1,n
            sigma(iter,iter) = sigma(iter,iter)+psi(iter)
         end do
      end if

      alpha = 1.0_dp
      ll = 0.0_dp
      if (f_tol < huge(1.0_dp)/2.0_dp) then
         ll = mvt_log_likelihood(data,mu,sigma,nu,status)
      end if
      converged_flag = .false.
      call cpu_time(start_time)
      do iter=1,niter
         old_mu = mu
         old_sigma = sigma
         old_nu = nu
         old_ll = ll

         if (has_missing .or. factor_structure) then
            call estep_general(data,mu,sigma,nu,ave_etau, &
               ave_elogtau,ave_etau_x,ave_etau_xx,status)
            if (status /= ht_success) exit
            if (do_mu) mu = ave_etau_x/ave_etau
            s = ave_etau_xx - outer_product(mu,ave_etau_x) - &
               outer_product(ave_etau_x,mu) + &
               ave_etau*outer_product(mu,mu)
            if (do_px) alpha = ave_etau
            s = s/alpha
            if (factor_structure) then
               call factor_decomposition(s,f,psi,b,sigma,status)
            else
               sigma = symmetrize(s)
            end if
            etau = ave_etau
            if (optimize_nu) then
               nu = optimize_ecm_nu(ave_elogtau-ave_etau,lo,hi)
            end if
         else
            call center_rows(data,mu,xc)
            call quadratic_forms(xc,sigma,q,status)
            if (status /= ht_success) exit
            etau = (real(n,dp)+nu)/(nu+q)*obs_weights
            ave_etau = sum(etau)/real(t,dp)

            if (optimize_nu .and. iter >= nu_start .and. &
               modulo(iter,nu_every) == 0) then
               nu = update_nu(iterative_method,xc,sigma,q,nu, &
                  ave_etau,lo,hi,status)
               if (status /= ht_success) exit
            end if
            etau = (real(n,dp)+nu)/(nu+q)*obs_weights
            ave_etau = sum(etau)/real(t,dp)
            ave_etau_x = matmul(etau,data)/real(t,dp)
            if (do_mu) mu = ave_etau_x/ave_etau
            call center_rows(data,mu,xc)
            sigma = weighted_covariance(xc,etau,real(t,dp))
            if (do_px) alpha = ave_etau
            sigma = sigma/alpha
         end if

         if (f_tol < huge(1.0_dp)/2.0_dp) then
            ll = mvt_log_likelihood(data,mu,sigma,nu,status)
            if (status /= ht_success) exit
            fun_ok = relative_scalar_converged(ll,old_ll,f_tol)
         else
            fun_ok = .true.
         end if
         params_ok = relative_vector_converged(mu,old_mu,p_tol) .and. &
            relative_scalar_converged(fnu(nu),fnu(old_nu),p_tol) .and. &
            relative_matrix_converged(sigma,old_sigma,p_tol)
         if (params_ok .and. fun_ok .and. iter >= nu_start) then
            converged_flag = .true.
            exit
         end if
      end do
      call cpu_time(end_time)
      if (status /= ht_success) then
         call set_error(result,status,'numerical failure in fit_mvt')
         return
      end if

      if (do_scale) then
         call center_rows(work,mu,xc)
         kappa = compute_kappa_from_marginals(xc)
         gamma_shape = gamma_s_psi1(sigma,t,1.0_dp+kappa)
         nmse = (kappa*(2.0_dp*gamma_shape+real(n,dp)) + &
            gamma_shape+real(n,dp))/(gamma_shape*real(t,dp))
         sigma = sigma/(1.0_dp+nmse)
      end if
      covariance = nu/(nu-2.0_dp)*sigma
      result%mu = mu
      result%mean = mu
      result%scatter = sigma
      result%covariance = covariance
      result%nu = nu
      result%num_iterations = min(iter,niter)
      result%converged = converged_flag
      result%cpu_time = end_time-start_time
      result%log_likelihood = ll
      result%latent_weights = etau
      if (factor_structure) then
         result%loadings = sqrt(nu/(nu-2.0_dp))*b
         result%psi = nu/(nu-2.0_dp)*psi
      end if
      if (result%converged) then
         result%status = ht_success
         result%message = 'success'
      else
         result%status = ht_no_convergence
         result%message = 'maximum iterations reached'
      end if
   end subroutine fit_mvt

   function update_nu(method,xc,sigma,q,nu,ave_etau,lo,hi,status) &
      result(updated)
      character(len=*),intent(in)::method
      real(dp),intent(in)::xc(:,:),sigma(:,:),q(:),nu,ave_etau,lo,hi
      integer,intent(out)::status
      real(dp)::updated
      real(dp),allocatable::regularized(:,:),q2(:),varx(:)
      real(dp)::c
      integer::n,t,i
      n=size(xc,2)
      t=size(xc,1)
      status=ht_success
      select case(trim(method))
      case('ECM','ECM-diag')
         allocate(q2(t),regularized(n,n))
         if(trim(method)=='ECM-diag') then
            regularized=0.0_dp
            do i=1,n
            regularized(i,i)=max(sigma(i,i),1.0e-12_dp)
            end do
         else
            regularized=0.9_dp*sigma
            do i=1,n
               regularized(i,i)=regularized(i,i)+0.1_dp*sigma(i,i)
            end do
         end if
         call quadratic_forms(xc,regularized,q2,status)
         if(status/=ht_success) then
         updated=nu
         return
         end if
         q2=(real(n,dp)+nu)/(nu+q2)
         c=digamma_dp(0.5_dp*(real(n,dp)+nu)) - &
            log(0.5_dp*(real(n,dp)+nu)) + &
            sum(log(q2)-q2)/real(t,dp)
         updated=optimize_ecm_nu(c,lo,hi)
      case('ECME','ECME-diag')
         if(trim(method)=='ECME-diag') then
            updated=nu_mle(xc,'MLE-mv-diagscat',sigma,lo,hi)
         else
            allocate(regularized(n,n))
            regularized=0.9_dp*sigma
            do i=1,n
               regularized(i,i)=regularized(i,i)+0.1_dp*sigma(i,i)
            end do
            updated=nu_mle(xc,'MLE-mv-scat',regularized,lo,hi)
         end if
      case('OPP')
         allocate(varx(n))
         varx=sum(xc*xc,dim=1)/real(t-1,dp)
         updated=nu_opp_estimator(varx,trace_matrix(sigma), &
            nu_min=lo,nu_max=hi,status=status)
      case('OPP-harmonic')
         allocate(varx(n))
         varx=sum(xc*xc,dim=1)/real(t-1,dp)
         updated=nu_opp_estimator(varx,trace_matrix(sigma),q, &
            'OPP-harmonic',lo,hi,status)
      case('POP','POP-approx-1','POP-approx-2', &
           'POP-approx-4','POP-sigma-corrected', &
           'POP-sigma-corrected-true')
         updated=nu_pop_estimator(q,nu,n,trim(method), &
            alpha=ave_etau,nu_min=lo,nu_max=hi,status=status)
      case('POP-approx-3','POP-exact')
         updated=nu_pop_estimator(q,nu,n,trim(method),xc=xc, &
            sigma=sigma,alpha=ave_etau,nu_min=lo,nu_max=hi, &
            status=status)
      case default
         updated=nu
         status=ht_invalid_argument
      end select
   end function update_nu

   function nu_mle(xc,method,sigma_scatter,nu_min,nu_max) result(nu)
      real(dp),intent(in)::xc(:,:)
      character(len=*),intent(in)::method
      real(dp),intent(in),optional::sigma_scatter(:,:)
      real(dp),intent(in),optional::nu_min,nu_max
      real(dp)::nu,lo,hi
      real(dp),allocatable::q(:),varx(:),delta(:,:),resampled(:)
      real(dp),allocatable::stacked(:)
      real(dp)::nu_sum
      type(rng_state)::rng
      integer::n,t,i,j,k,index,nres,trep,status
      lo=default_nu_min
      hi=default_nu_max
      if(present(nu_min)) lo=nu_min
      if(present(nu_max)) hi=nu_max
      n=size(xc,2)
      t=size(xc,1)
      select case(trim(method))
      case('MLE-mv-cov')
         allocate(q(t))
         if (present(sigma_scatter)) then
            call quadratic_forms(xc,sigma_scatter,q,status)
         else
            block
               real(dp) :: covariance(size(xc,2),size(xc,2))
               covariance = sample_covariance(xc)
               call quadratic_forms(xc,covariance,q,status)
            end block
         end if
         if (status /= ht_success) then
            nu = hi
            return
         end if
         nu = golden_mle_covariance(q,n,t,lo,hi)
      case('MLE-mv-scat')
         if(.not.present(sigma_scatter)) then
         nu=hi
         return
         end if
         allocate(q(t))
         call quadratic_forms(xc,sigma_scatter,q,status)
         if(status/=ht_success) then
         nu=hi
         return
         end if
         nu=golden_mle_scatter(q,n,t,lo,hi)
      case('MLE-mv-diagscat','MLE-mv-diagscat-resampled')
         if(.not.present(sigma_scatter)) then
         nu=hi
         return
         end if
         allocate(q(t))
         q=0.0_dp
         do j=1,n
            q=q+xc(:,j)**2/max(sigma_scatter(j,j),tiny(1.0_dp))
         end do
         nu=golden_mle_scatter(q,n,t,lo,hi)
      case('MLE-mv-diagcov')
         allocate(varx(n),q(t))
         varx=sum(xc*xc,dim=1)/real(t-1,dp)
         q=0.0_dp
         do j=1,n
         q=q+xc(:,j)**2/max(varx(j),tiny(1.0_dp))
         end do
         nu=golden_mle_covariance(q,n,t,lo,hi)
      case('MLE-mv-diagcov-resampled')
         nres=max(1,nint(1.2_dp*real(n,dp)))
         trep=4*t
         allocate(varx(n),delta(t,n),resampled(trep))
         varx=sum(xc*xc,dim=1)/real(t-1,dp)
         do j=1,n
         delta(:,j)=xc(:,j)**2/max(varx(j),tiny(1.0_dp))
         end do
         call seed_rng(rng,13579)
         do i=1,trep
            k=modulo(i-1,t)+1
            resampled(i)=0.0_dp
            do j=1,nres
               index=1+int(uniform_random(rng)*real(n,dp))
               index=min(n,max(1,index))
               resampled(i)=resampled(i)+delta(k,index)
            end do
         end do
         nu=golden_mle_covariance(resampled,n,trep,lo,hi)
      case('MLE-uv-var-ave')
         allocate(varx(n),q(t))
         varx = sum(xc*xc,dim=1)/real(t-1,dp)
         nu_sum = 0.0_dp
         do j = 1, n
            q = xc(:,j)**2/max(varx(j),tiny(1.0_dp))
            nu_sum = nu_sum + golden_mle_covariance(q,1,t,lo,hi)
         end do
         nu = nu_sum/real(n,dp)
      case('MLE-uv-scat-ave')
         if (.not.present(sigma_scatter)) then
            nu = hi
            return
         end if
         allocate(q(t))
         nu_sum = 0.0_dp
         do j = 1, n
            q = xc(:,j)**2/max(sigma_scatter(j,j),tiny(1.0_dp))
            nu_sum = nu_sum + golden_mle_scatter(q,1,t,lo,hi)
         end do
         nu = nu_sum/real(n,dp)
      case('MLE-uv-var-stacked')
         allocate(varx(n),stacked(t*n))
         varx = sum(xc*xc,dim=1)/real(t-1,dp)
         do j = 1, n
            stacked((j-1)*t+1:j*t) = &
               xc(:,j)**2/max(varx(j),tiny(1.0_dp))
         end do
         nu = golden_mle_covariance(stacked,1,t*n,lo,hi)
      case default
         nu=hi
      end select
   end function nu_mle

   function golden_mle_scatter(q,n,t,lo,hi) result(root)
      real(dp),intent(in)::q(:),lo,hi
      integer,intent(in)::n,t
      real(dp)::root,a,b,c,d,fc,fd
      real(dp),parameter::gr=0.6180339887498948482_dp
      integer::iter
      a=lo
      b=hi
      c=b-gr*(b-a)
      d=a+gr*(b-a)
      fc=mle_scatter_objective(c,q,n,t)
      fd=mle_scatter_objective(d,q,n,t)
      do iter=1,120
         if(abs(b-a)<=1.0e-8_dp*(1.0_dp+abs(a)+abs(b))) exit
         if(fc<fd) then
            b=d
            d=c
            fd=fc
            c=b-gr*(b-a)
            fc=mle_scatter_objective(c,q,n,t)
         else
            a=c
            c=d
            fc=fd
            d=a+gr*(b-a)
            fd=mle_scatter_objective(d,q,n,t)
         end if
      end do
      root=0.5_dp*(a+b)
   end function golden_mle_scatter

   function golden_mle_covariance(q,n,t,lo,hi) result(root)
      real(dp),intent(in)::q(:),lo,hi
      integer,intent(in)::n,t
      real(dp)::root,a,b,c,d,fc,fd
      real(dp),parameter::gr=0.6180339887498948482_dp
      integer::iter
      a=lo
      b=hi
      c=b-gr*(b-a)
      d=a+gr*(b-a)
      fc=mle_cov_objective(c,q,n,t)
      fd=mle_cov_objective(d,q,n,t)
      do iter=1,120
         if(abs(b-a)<=1.0e-8_dp*(1.0_dp+abs(a)+abs(b))) exit
         if(fc<fd) then
            b=d
            d=c
            fd=fc
            c=b-gr*(b-a)
            fc=mle_cov_objective(c,q,n,t)
         else
            a=c
            c=d
            fc=fd
            d=a+gr*(b-a)
            fd=mle_cov_objective(d,q,n,t)
         end if
      end do
      root=0.5_dp*(a+b)
   end function golden_mle_covariance

   pure function mle_scatter_objective(nu,q,n,t) result(value)
      real(dp),intent(in)::nu,q(:)
      integer,intent(in)::n,t
      real(dp)::value
      value=0.5_dp*(real(n,dp)+nu)*sum(log(nu+q)) - &
         real(t,dp)*log_gamma(0.5_dp*(real(n,dp)+nu)) + &
         real(t,dp)*log_gamma(0.5_dp*nu) - &
         0.5_dp*nu*real(t,dp)*log(nu)
   end function mle_scatter_objective

   pure function mle_cov_objective(nu,q,n,t) result(value)
      real(dp),intent(in)::nu,q(:)
      integer,intent(in)::n,t
      real(dp)::value,ratio
      ratio=nu/(nu-2.0_dp)
      value=0.5_dp*real(n*t,dp)*log((nu-2.0_dp)/nu) + &
         0.5_dp*(real(n,dp)+nu)*sum(log(nu+ratio*q)) - &
         real(t,dp)*log_gamma(0.5_dp*(real(n,dp)+nu)) + &
         real(t,dp)*log_gamma(0.5_dp*nu) - &
         0.5_dp*nu*real(t,dp)*log(nu)
   end function mle_cov_objective

   function optimize_ecm_nu(cvalue,lo,hi) result(root)
      real(dp),intent(in)::cvalue,lo,hi
      real(dp)::root,a,b,c,d,fc,fd
      real(dp),parameter::gr=0.6180339887498948482_dp
      integer::iter
      a=lo
      b=hi
      c=b-gr*(b-a)
      d=a+gr*(b-a)
      fc=ecm_objective(c,cvalue)
      fd=ecm_objective(d,cvalue)
      do iter=1,120
         if(abs(b-a)<=1.0e-8_dp*(1.0_dp+abs(a)+abs(b))) exit
         if(fc<fd) then
            b=d
            d=c
            fd=fc
            c=b-gr*(b-a)
            fc=ecm_objective(c,cvalue)
         else
            a=c
            c=d
            fc=fd
            d=a+gr*(b-a)
            fd=ecm_objective(d,cvalue)
         end if
      end do
      root=0.5_dp*(a+b)
   end function optimize_ecm_nu

   pure function ecm_objective(nu,cvalue) result(value)
      real(dp),intent(in)::nu,cvalue
      real(dp)::value
      value=-0.5_dp*nu*log(0.5_dp*nu)+log_gamma(0.5_dp*nu) - &
         0.5_dp*nu*cvalue
   end function ecm_objective

   function mvt_log_likelihood(x,mu,sigma,nu,status) result(value)
      real(dp),intent(in)::x(:,:),mu(:),sigma(:,:),nu
      integer,intent(out),optional::status
      real(dp)::value
      real(dp),allocatable::subsigma(:,:),delta(:),solution(:)
      integer,allocatable::index(:)
      integer::i,j,p,istat,n
      real(dp)::q,logdet,constant
      n=size(x,2)
      value=0.0_dp
      istat=ht_success
      do i=1,size(x,1)
         p=count([( .not.ieee_is_nan(x(i,j)),j=1,n )])
         if(p==0) cycle
         allocate(index(p),subsigma(p,p),delta(p),solution(p))
         p=0
         do j=1,n
            if(.not.ieee_is_nan(x(i,j))) then
               p=p+1
               index(p)=j
               delta(p)=x(i,j)-mu(j)
            end if
         end do
         do j=1,p
            subsigma(j,:)=sigma(index(j),index)
         end do
         call solve_vector_local(subsigma,delta,solution,istat)
         if(istat/=ht_success) then
            value=-huge(1.0_dp)
            exit
         end if
         q=dot_product(delta,solution)
         logdet=logdet_spd(subsigma,istat)
         if(istat/=ht_success) then
            value=-huge(1.0_dp)
            exit
         end if
         constant=log_gamma(0.5_dp*(nu+real(p,dp))) - &
            log_gamma(0.5_dp*nu) - 0.5_dp*(real(p,dp)* &
            log(nu*acos(-1.0_dp))+logdet)
         value=value+constant-0.5_dp*(nu+real(p,dp))* &
            log(1.0_dp+q/nu)
         deallocate(index,subsigma,delta,solution)
      end do
      if(present(status)) status=istat
   end function mvt_log_likelihood

   subroutine estep_general(x,mu,sigma,nu,ave_etau,ave_elogtau, &
      ave_etau_x,ave_etau_xx,status)
      real(dp),intent(in)::x(:,:),mu(:),sigma(:,:),nu
      real(dp),intent(out)::ave_etau,ave_elogtau,ave_etau_x(:)
      real(dp),intent(out)::ave_etau_xx(:,:)
      integer,intent(out)::status
      integer::t,n,i,j,k,p,m,istat
      integer,allocatable::obs(:),mis(:)
      real(dp),allocatable::soo(:,:),invsoo(:,:),delta(:),completed(:)
      real(dp),allocatable::condcov(:,:),smo(:,:),tmpmat(:,:)
      real(dp)::tmp,etau_i,elog_i
      t=size(x,1)
      n=size(x,2)
      status=ht_success
      ave_etau=0.0_dp
      ave_elogtau=0.0_dp
      ave_etau_x=0.0_dp
      ave_etau_xx=0.0_dp
      do i=1,t
         p=count([( .not.ieee_is_nan(x(i,j)),j=1,n )])
         m=n-p
         if(p==0) cycle
         allocate(obs(p),mis(m),soo(p,p),invsoo(p,p),delta(p),completed(n))
         p=0
         m=0
         do j=1,n
            if(ieee_is_nan(x(i,j))) then
               m=m+1
               mis(m)=j
            else
               p=p+1
               obs(p)=j
            end if
         end do
         do j=1,p
         soo(j,:)=sigma(obs(j),obs)
         end do
         call inverse_matrix(soo,invsoo,istat)
         if(istat/=ht_success) then
         status=istat
         return
         end if
         delta=x(i,obs)-mu(obs)
         tmp=nu+dot_product(delta,matmul(invsoo,delta))
         etau_i=(nu+real(p,dp))/tmp
         elog_i=digamma_dp(0.5_dp*(nu+real(p,dp)))-log(0.5_dp*tmp)
         completed=mu
         completed(obs)=x(i,obs)
         if(m>0) then
            allocate(smo(m,p))
            do j=1,m
            smo(j,:)=sigma(mis(j),obs)
            end do
            completed(mis)=mu(mis)+matmul(smo,matmul(invsoo,delta))
         end if
         ave_etau=ave_etau+etau_i/real(t,dp)
         ave_elogtau=ave_elogtau+elog_i/real(t,dp)
         ave_etau_x=ave_etau_x+etau_i*completed/real(t,dp)
         ave_etau_xx=ave_etau_xx+etau_i* &
            outer_product(completed,completed)
         if(m>0) then
            allocate(condcov(m,m),tmpmat(m,p))
            tmpmat=matmul(smo,invsoo)
            do j=1,m
               do k=1,m
                  condcov(j,k)=sigma(mis(j),mis(k)) - &
                     dot_product(tmpmat(j,:),sigma(obs,mis(k)))
                  ave_etau_xx(mis(j),mis(k)) = &
                     ave_etau_xx(mis(j),mis(k))+condcov(j,k)
               end do
            end do
         end if
         deallocate(obs,mis,soo,invsoo,delta,completed)
         if(allocated(smo)) deallocate(smo)
         if(allocated(condcov)) deallocate(condcov)
         if(allocated(tmpmat)) deallocate(tmpmat)
      end do
      ave_etau_xx=symmetrize(ave_etau_xx/real(t,dp))
   end subroutine estep_general

   subroutine initialize_missing(x)
      real(dp),intent(inout)::x(:,:)
      real(dp)::mu(size(x,2))
      integer::i,j,count_valid
      do j=1,size(x,2)
         mu(j)=0.0_dp
         count_valid=0
         do i=1,size(x,1)
            if(.not.ieee_is_nan(x(i,j))) then
               mu(j)=mu(j)+x(i,j)
               count_valid=count_valid+1
            end if
         end do
         if(count_valid>0) mu(j)=mu(j)/real(count_valid,dp)
         do i=1,size(x,1)
            if(ieee_is_nan(x(i,j))) x(i,j)=mu(j)
         end do
      end do
   end subroutine initialize_missing

   function column_mean_nan(x) result(mu)
      real(dp),intent(in)::x(:,:)
      real(dp)::mu(size(x,2))
      integer::i,j,count_valid
      do j=1,size(x,2)
         mu(j)=0.0_dp
         count_valid=0
         do i=1,size(x,1)
            if(.not.ieee_is_nan(x(i,j))) then
               mu(j)=mu(j)+x(i,j)
               count_valid=count_valid+1
            end if
         end do
         if(count_valid>0) mu(j)=mu(j)/real(count_valid,dp)
      end do
   end function column_mean_nan

   function any_nan(x) result(found)
      real(dp),intent(in)::x(:,:)
      logical::found
      integer::i,j
      found=.false.
      do j=1,size(x,2)
      do i=1,size(x,1)
         if(ieee_is_nan(x(i,j))) then
         found=.true.
         return
         end if
      end do
      end do
   end function any_nan

   subroutine center_rows(x,mu,xc)
      real(dp),intent(in)::x(:,:),mu(:)
      real(dp),intent(out)::xc(:,:)
      integer::i
      do i=1,size(x,1)
      xc(i,:)=x(i,:)-mu
      end do
   end subroutine center_rows

   function compute_kappa_from_marginals(xc) result(kappa)
      real(dp),intent(in)::xc(:,:)
      real(dp)::kappa,ki(size(xc,2)),m2,m4,k
      integer::j,t
      t=size(xc,1)
      if(t<4) then
      kappa=0.0_dp
      return
      end if
      do j=1,size(xc,2)
         m2=sum(xc(:,j)**2)/real(t,dp)
         m4=sum(xc(:,j)**4)/real(t,dp)
         if(m2<=tiny(1.0_dp)) then
            ki(j)=0.0_dp
         else
            ki(j)=real((t+1)*(t-1),dp)/real((t-2)*(t-3),dp) * &
               (m4/(m2*m2)-3.0_dp*real(t-1,dp)/real(t+1,dp))
         end if
      end do
      k=sum(ki)/real(size(ki),dp)
      kappa=real(t-1,dp)/real((t-2)*(t-3),dp) * &
         (real(t+1,dp)*k+6.0_dp)/3.0_dp
   end function compute_kappa_from_marginals

   function gamma_s_psi1(s,t,psi1) result(gamma)
      real(dp),intent(in)::s(:,:),psi1
      integer,intent(in)::t
      real(dp)::gamma,snorm(size(s,1),size(s,2)),naive,a,b,tr
      integer::n
      n=size(s,1)
      tr=trace_matrix(s)
      snorm=s/tr
      naive=real(n,dp)*sum(snorm*snorm)
      a=real(t,dp)/real(t,dp)+0.0_dp
      a=real(t,dp)/(real(t,dp)+psi1-1.0_dp)*psi1
      b=real(t,dp)/real(t-1,dp) * &
         (real(t,dp)+psi1-1.0_dp)/(real(t,dp)+3.0_dp*psi1-1.0_dp)
      gamma=min(real(n,dp),max(1.0_dp,b*(naive-a*real(n,dp)/real(t,dp))))
   end function gamma_s_psi1

   pure function fnu(nu) result(value)
      real(dp),intent(in)::nu
      real(dp)::value
      value=nu/(nu-2.0_dp)
   end function fnu

   pure function relative_vector_converged(a,b,tol) result(ok)
      real(dp),intent(in)::a(:),b(:),tol
      logical::ok
      ok=all(abs(a-b)<=0.5_dp*tol*(abs(a)+abs(b)+tiny(1.0_dp)))
   end function relative_vector_converged

   pure function relative_matrix_converged(a,b,tol) result(ok)
      real(dp),intent(in)::a(:,:),b(:,:),tol
      logical::ok
      ok=all(abs(a-b)<=0.5_dp*tol*(abs(a)+abs(b)+tiny(1.0_dp)))
   end function relative_matrix_converged

   pure function relative_scalar_converged(a,b,tol) result(ok)
      real(dp),intent(in)::a,b,tol
      logical::ok
      ok=abs(a-b)<=0.5_dp*tol*(abs(a)+abs(b)+tiny(1.0_dp))
   end function relative_scalar_converged

   subroutine solve_vector_local(a,b,x,status)
      use fitheavytail_linalg, only: solve_linear
      real(dp),intent(in)::a(:,:),b(:)
      real(dp),intent(out)::x(:)
      integer,intent(out)::status
      call solve_linear(a,b,x,status)
   end subroutine solve_vector_local

   subroutine set_error(result,status,message)
      type(heavy_tail_fit),intent(inout)::result
      integer,intent(in)::status
      character(len=*),intent(in)::message
      result%status=status
      result%message=message
      result%converged=.false.
   end subroutine set_error

end module fitheavytail_mvt
