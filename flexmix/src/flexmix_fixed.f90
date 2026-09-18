! SPDX-License-Identifier: GPL-2.0-or-later
module flexmix_fixed
   use flexmix_kinds, only : dp
   use flexmix_types
   use flexmix_components, only : fit_gaussian_regression_component, gaussian_regression_log_density
   use flexmix_components, only : fit_poisson_regression_component, poisson_regression_log_density
   use flexmix_components, only : fit_binomial_regression_component, binomial_regression_log_density
   use flexmix_components, only : fit_gamma_regression_component, gamma_regression_log_density
   use flexmix_em_utils, only : initialize_posteriors, classify_posteriors
   use flexmix_numeric, only : normalize_log_probabilities
   implicit none
   private

   public :: flexmix_glmfix_gaussian
   public :: flexmix_glmfix_poisson
   public :: flexmix_glmfix_binomial
   public :: flexmix_glmfix_gamma

contains

   subroutine flexmix_glmfix_gaussian(x, y, design, result, control, case_weights, initial_cluster, &
                                      initial_posterior, offset, variance_group)
      real(dp), intent(in) :: x(:,:,:) !! Component-specific global design arrays, shape `(n,q,k)`; inactive columns are masked.
      real(dp), intent(in) :: y(:) !! Gaussian response vector, size `n`, shared by all mixture components.
      logical, intent(in) :: design(:,:) !! Coefficient-incidence matrix, shape `(k,q)`; true entries share one global.
      type(flexmix_result), intent(out) :: result !! Fitted fixed/shared-coefficient Gaussian mixture result.
      type(flexmix_control), intent(in), optional :: control !! Optional EM iteration, tolerance, and classification controls.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`; defaults to one.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based initial component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional initial membership matrix, shape `(n,k)`.
      real(dp), intent(in), optional :: offset(:,:) !! Optional component-specific additive offsets, shape `(n,k)`.
      integer, intent(in), optional :: variance_group(:) !! Optional positive variance-group labels, size `k`; default gives.
      call flexmix_fit_glmfix(x, y, design, flexmix_model_gaussian, result, control, case_weights, &
                              initial_cluster, initial_posterior, offset, variance_group=variance_group)
   end subroutine flexmix_glmfix_gaussian

   subroutine flexmix_glmfix_poisson(x, y, design, result, control, case_weights, initial_cluster, initial_posterior, offset)
      real(dp), intent(in) :: x(:,:,:) !! Component-specific global design arrays, shape `(n,q,k)`; inactive columns are masked.
      real(dp), intent(in) :: y(:) !! Nonnegative Poisson count response, size `n`.
      logical, intent(in) :: design(:,:) !! Coefficient-incidence matrix, shape `(k,q)`, identifying globally shared coefficient.
      type(flexmix_result), intent(out) :: result !! Fitted fixed/shared-coefficient Poisson mixture result.
      type(flexmix_control), intent(in), optional :: control !! Optional EM iteration, tolerance, and classification controls.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`; defaults to one.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based initial component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional initial membership matrix, shape `(n,k)`.
      real(dp), intent(in), optional :: offset(:,:) !! Optional component-specific additive log-link offsets, shape `(n,k)`.
      call flexmix_fit_glmfix(x, y, design, flexmix_model_poisson, result, control, case_weights, &
                              initial_cluster, initial_posterior, offset)
   end subroutine flexmix_glmfix_poisson

   subroutine flexmix_glmfix_binomial(x, success, failure, design, result, control, case_weights, initial_cluster, &
                                      initial_posterior, offset)
      real(dp), intent(in) :: x(:,:,:) !! Component-specific global design arrays, shape `(n,q,k)`; inactive columns are masked.
      real(dp), intent(in) :: success(:) !! Nonnegative binomial success counts, size `n`.
      real(dp), intent(in) :: failure(:) !! Nonnegative binomial failure counts, size `n`.
      logical, intent(in) :: design(:,:) !! Coefficient-incidence matrix, shape `(k,q)`, identifying globally shared coefficient.
      type(flexmix_result), intent(out) :: result !! Fitted fixed/shared-coefficient binomial mixture result.
      type(flexmix_control), intent(in), optional :: control !! Optional EM iteration, tolerance, and classification controls.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`; defaults to one.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based initial component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional initial membership matrix, shape `(n,k)`.
      real(dp), intent(in), optional :: offset(:,:) !! Optional component-specific additive logit offsets, shape `(n,k)`.
      real(dp), allocatable :: total(:)

      if (size(success) /= size(failure) .or. any(success < 0.0_dp) .or. any(failure < 0.0_dp)) then
         result%status = -1
         return
      end if
      allocate(total(size(success)))
      total = success + failure
      call flexmix_fit_glmfix(x, success, design, flexmix_model_binomial, result, control, case_weights, &
                              initial_cluster, initial_posterior, offset, trials=total)
   end subroutine flexmix_glmfix_binomial

   subroutine flexmix_glmfix_gamma(x, y, design, result, control, case_weights, initial_cluster, initial_posterior, offset)
      real(dp), intent(in) :: x(:,:,:) !! Component-specific global design arrays, shape `(n,q,k)`; inactive columns are masked.
      real(dp), intent(in) :: y(:) !! Strictly positive Gamma response vector, size `n`.
      logical, intent(in) :: design(:,:) !! Coefficient-incidence matrix, shape `(k,q)`, identifying globally shared coefficient.
      type(flexmix_result), intent(out) :: result !! Fitted fixed/shared-coefficient Gamma inverse-link mixture result.
      type(flexmix_control), intent(in), optional :: control !! Optional EM iteration, tolerance, and classification controls.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`; defaults to one.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based initial component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional initial membership matrix, shape `(n,k)`.
      real(dp), intent(in), optional :: offset(:,:) !! Optional component-specific additive inverse-link offsets, shape `(n,k)`.
      call flexmix_fit_glmfix(x, y, design, flexmix_model_gamma_regression, result, control, case_weights, &
                              initial_cluster, initial_posterior, offset)
   end subroutine flexmix_glmfix_gamma

   subroutine flexmix_fit_glmfix(x, y, design, model_kind, result, control, case_weights, initial_cluster, &
                                 initial_posterior, offset, trials, variance_group)
      real(dp), intent(in) :: x(:,:,:) !! Component-specific global design arrays, shape `(n,q,k)`.
      real(dp), intent(in) :: y(:) !! Response vector, size `n`; for binomial models this contains success counts.
      logical, intent(in) :: design(:,:) !! Coefficient-incidence matrix `(k,q)` for shared/component-specific terms.
      integer, intent(in) :: model_kind !! Regression family identifier for Gaussian, Poisson, binomial, or Gamma inverse-link.
      type(flexmix_result), intent(out) :: result !! Fitted fixed/shared-coefficient mixture result.
      type(flexmix_control), intent(in), optional :: control !! Optional EM controls; `minprior` deletion is intentionally not.
      real(dp), intent(in), optional :: case_weights(:) !! Optional nonnegative case/frequency weights, size `n`.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based initial component labels, size `n`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional initial membership matrix, shape `(n,k)`.
      real(dp), intent(in), optional :: offset(:,:) !! Optional component-specific family-scale offsets, shape `(n,k)`.
      real(dp), intent(in), optional :: trials(:) !! Optional binomial trial counts, size `n`, required only for binomial.
      integer, intent(in), optional :: variance_group(:) !! Optional positive Gaussian variance-group labels, size `k`.
      type(flexmix_control) :: ctl
      real(dp), allocatable :: weights(:), posterior(:,:), classified(:,:), prior(:), log_density(:,:), logp(:,:), logsum(:)
      real(dp), allocatable :: xstack(:,:), ystack(:), wstack(:), offset_stack(:), trials_stack(:), beta_global(:), old_beta(:)
      real(dp), allocatable :: sigma(:), shape_value(:), beta_component(:)
      integer, allocatable :: var_group(:), cluster(:), component_size(:)
      real(dp) :: llh, old_llh, dummy_sigma, dummy_shape
      integer :: n, q, k, i, j, c, pos, iter, fit_info, df_component, info, ngroup

      ctl = flexmix_control()
      if (present(control)) ctl = control
      n = size(y)
      q = size(x,2)
      k = size(x,3)
      result%status = 0
      result%model_kind = model_kind
      result%n = n
      result%p = q
      result%d = 1
      result%k = k
      result%k0 = k
      if (size(x,1) /= n .or. size(design,1) /= k .or. size(design,2) /= q .or. k < 1 .or. q < 1) then
         result%status = -1
         return
      end if
      if (any(.not. any(design,dim=1))) then
         result%status = -2
         return
      end if
      if (model_kind == flexmix_model_poisson .and. any(y < 0.0_dp)) then
         result%status = -3
         return
      end if
      if (model_kind == flexmix_model_gamma_regression .and. any(y <= 0.0_dp)) then
         result%status = -4
         return
      end if
      if (model_kind == flexmix_model_binomial) then
         if (.not. present(trials)) then
            result%status = -5
            return
         end if
         if (size(trials) /= n .or. any(trials < y) .or. any(y < 0.0_dp)) then
            result%status = -6
            return
         end if
      end if
      if (present(offset)) then
         if (size(offset,1) /= n .or. size(offset,2) /= k) then
            result%status = -7
            return
         end if
      end if
      allocate(weights(n), posterior(n,k), classified(n,k), prior(k), log_density(n,k), logp(n,k), logsum(n))
      allocate(xstack(n*k,q), ystack(n*k), wstack(n*k), offset_stack(n*k), beta_global(q), old_beta(q))
      allocate(sigma(k), shape_value(k), beta_component(q))
      allocate(var_group(k), cluster(n), component_size(k))
      if (present(trials)) allocate(trials_stack(n*k))
      weights = 1.0_dp
      if (present(case_weights)) then
         if (size(case_weights) /= n .or. any(case_weights < 0.0_dp)) then
            result%status = -8
            return
         end if
         weights = case_weights
      end if
      var_group = [(j,j=1,k)]
      if (present(variance_group)) then
         if (size(variance_group) /= k .or. any(variance_group < 1)) then
            result%status = -9
            return
         end if
         var_group = variance_group
      end if
      ngroup = maxval(var_group)
      do c = 1, ngroup
         if (.not. any(var_group == c)) then
            result%status = -10
            return
         end if
      end do
      call initialize_posteriors(y, k, posterior, initial_cluster, initial_posterior, info)
      if (info /= 0) then
         result%status = info
         return
      end if
      beta_global = 0.0_dp
      sigma = 1.0_dp
      shape_value = 1.0_dp
      prior = 1.0_dp / real(k,dp)
      old_llh = -huge(1.0_dp)
      llh = old_llh
      do iter = 1, max(1,ctl%iter_max)
         call classify_posteriors(posterior, ctl%classify, classified)
         do j = 1, k
            prior(j) = sum(weights * classified(:,j))
         end do
         if (sum(prior) > 0.0_dp) then
            prior = prior / sum(prior)
         else
            prior = 1.0_dp / real(k,dp)
         end if
         pos = 0
         do j = 1, k
            do i = 1, n
               pos = pos + 1
               xstack(pos,:) = 0.0_dp
               do c = 1, q
                  if (design(j,c)) xstack(pos,c) = x(i,c,j)
               end do
               ystack(pos) = y(i)
               wstack(pos) = weights(i) * classified(i,j)
               if (present(offset)) then
                  offset_stack(pos) = offset(i,j)
               else
                  offset_stack(pos) = 0.0_dp
               end if
               if (present(trials)) trials_stack(pos) = trials(i)
            end do
         end do
         old_beta = beta_global
         select case (model_kind)
         case (flexmix_model_gaussian)
            call fit_gaussian_regression_component(xstack, ystack, wstack, beta_global, dummy_sigma, df_component, fit_info, &
                                                   offset_stack)
         case (flexmix_model_poisson)
            call fit_poisson_regression_component(xstack, ystack, wstack, beta_global, df_component, fit_info, &
                                                  old_beta, offset_stack)
         case (flexmix_model_binomial)
            call fit_binomial_regression_component(xstack, ystack, trials_stack, wstack, beta_global, df_component, fit_info, &
                                                   old_beta, offset_stack)
         case (flexmix_model_gamma_regression)
            call fit_gamma_regression_component(xstack, ystack, wstack, beta_global, dummy_shape, df_component, fit_info, &
                                                 old_beta, offset_stack)
         case default
            fit_info = 1
         end select
         if (fit_info /= 0) then
            result%status = 20 + fit_info
            exit
         end if
         if (model_kind == flexmix_model_gaussian) then
            call update_grouped_sigma(x, y, design, beta_global, classified, weights, offset, var_group, sigma)
         else if (model_kind == flexmix_model_gamma_regression) then
            shape_value = dummy_shape
         end if
         do j = 1, k
            beta_component = 0.0_dp
            do c = 1, q
               if (design(j,c)) beta_component(c) = beta_global(c)
            end do
            select case (model_kind)
            case (flexmix_model_gaussian)
               if (present(offset)) then
                  call gaussian_regression_log_density(x(:,:,j), y, beta_component, sigma(j), log_density(:,j), offset(:,j))
               else
                  call gaussian_regression_log_density(x(:,:,j), y, beta_component, sigma(j), log_density(:,j))
               end if
            case (flexmix_model_poisson)
               if (present(offset)) then
                  call poisson_regression_log_density(x(:,:,j), y, beta_component, log_density(:,j), offset(:,j))
               else
                  call poisson_regression_log_density(x(:,:,j), y, beta_component, log_density(:,j))
               end if
            case (flexmix_model_binomial)
               if (present(offset)) then
                  call binomial_regression_log_density(x(:,:,j), y, trials, beta_component, log_density(:,j), offset(:,j))
               else
                  call binomial_regression_log_density(x(:,:,j), y, trials, beta_component, log_density(:,j))
               end if
            case (flexmix_model_gamma_regression)
               if (present(offset)) then
                  call gamma_regression_log_density(x(:,:,j), y, beta_component, shape_value(j), log_density(:,j), offset(:,j))
               else
                  call gamma_regression_log_density(x(:,:,j), y, beta_component, shape_value(j), log_density(:,j))
               end if
            end select
            logp(:,j) = log_density(:,j) + log(max(prior(j),tiny(1.0_dp)))
         end do
         call normalize_log_probabilities(logp, posterior, logsum)
         llh = sum(weights * logsum)
         if (iter > 1) then
            if (abs(llh-old_llh)/(abs(llh)+0.1_dp) < ctl%tolerance) then
               result%converged = .true.
               exit
            end if
         end if
         old_llh = llh
      end do
      result%iterations = min(iter,max(1,ctl%iter_max))
      result%loglik = llh
      result%prior = prior
      result%posterior = posterior
      result%posterior_unscaled = exp(logp)
      result%log_posterior_unscaled = logp
      result%case_weights = weights
      result%parameter_design = design
      result%variance_group = var_group
      allocate(result%component_role(k))
      result%component_role = flexmix_role_regular
      allocate(result%beta(q,k))
      result%beta = 0.0_dp
      do j = 1, k
         do c = 1, q
            if (design(j,c)) result%beta(c,j) = beta_global(c)
         end do
      end do
      if (model_kind == flexmix_model_gaussian) result%sigma = sigma
      if (model_kind == flexmix_model_gamma_regression) result%shape = shape_value
      component_size = 0
      do i = 1, n
         cluster(i) = maxloc(posterior(i,:),dim=1)
         component_size(cluster(i)) = component_size(cluster(i)) + nint(weights(i))
      end do
      result%cluster = cluster
      result%size = component_size
      result%df = count(any(design,dim=1)) + k - 1
      if (model_kind == flexmix_model_gaussian) result%df = result%df + ngroup
      if (model_kind == flexmix_model_gamma_regression) result%df = result%df + 1
      if (result%status == 0 .and. .not. result%converged .and. ctl%iter_max > 0) result%status = 1
   end subroutine flexmix_fit_glmfix

   subroutine update_grouped_sigma(x, y, design, beta_global, classified, case_weights, offset, variance_group, sigma)
      real(dp), intent(in) :: x(:,:,:) !! Component-specific global design arrays, shape `(n,q,k)`.
      real(dp), intent(in) :: y(:) !! Gaussian response vector, size `n`.
      logical, intent(in) :: design(:,:) !! Coefficient-incidence matrix, shape `(k,q)`.
      real(dp), intent(in) :: beta_global(:) !! Current unique global coefficient vector, size `q`.
      real(dp), intent(in) :: classified(:,:) !! M-step component weights, shape `(n,k)`.
      real(dp), intent(in) :: case_weights(:) !! Nonnegative case/frequency weights, size `n`.
      real(dp), intent(in), optional :: offset(:,:) !! Optional component-specific Gaussian offsets, shape `(n,k)`.
      integer, intent(in) :: variance_group(:) !! Positive variance-group labels, size `k`.
      real(dp), intent(out) :: sigma(:) !! Updated component standard deviations, size `k`, shared within equal variance groups.
      real(dp), allocatable :: beta_component(:), residual(:), w(:)
      real(dp) :: numerator, denominator, mean_weight, sigma_group
      integer :: g, j, active_parameters

      allocate(beta_component(size(beta_global)), residual(size(y)), w(size(y)))
      do g = 1, maxval(variance_group)
         numerator = 0.0_dp
         denominator = 0.0_dp
         active_parameters = 0
         do j = 1, size(variance_group)
            if (variance_group(j) /= g) cycle
            beta_component = 0.0_dp
            where (design(j,:)) beta_component = beta_global
            residual = y - matmul(x(:,:,j),beta_component)
            if (present(offset)) residual = residual - offset(:,j)
            w = case_weights * classified(:,j)
            mean_weight = sum(w) / real(size(w),dp)
            if (mean_weight > tiny(1.0_dp)) numerator = numerator + sum(w*residual**2) / mean_weight
            denominator = denominator + real(size(y),dp)
            active_parameters = max(active_parameters,count(design(j,:)))
         end do
         denominator = denominator - real(active_parameters,dp)
         if (denominator > 0.0_dp) then
            sigma_group = sqrt(max(tiny(1.0_dp),numerator/denominator))
         else
            sigma_group = 1.0_dp
         end if
         do j = 1, size(variance_group)
            if (variance_group(j) == g) sigma(j) = sigma_group
         end do
      end do
   end subroutine update_grouped_sigma

end module flexmix_fixed
