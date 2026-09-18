! SPDX-License-Identifier: GPL-2.0-or-later
module flexmix_mixed
   use flexmix_kinds, only : dp
   use flexmix_types, only : flexmix_control, flexmix_result
   use flexmix_types, only : flexmix_class_hard, flexmix_model_lmm, flexmix_model_lmer, flexmix_model_lmmc, flexmix_model_lmc
   use flexmix_em_utils, only : initialize_posteriors
   use flexmix_numeric, only : normalize_log_probabilities, weighted_least_squares, inverse_logdet_spd
   implicit none
   private
   real(dp), parameter :: pi = acos(-1.0_dp)
   public :: flexmix_lmm
   public :: flexmix_lmer
   public :: flexmix_lmmc
   public :: flexmix_lmc

contains

   subroutine flexmix_lmm(x, y, z, group, k, result, control, group_weights, initial_cluster, initial_posterior, &
                          varfix_random, varfix_residual)
      real(dp), intent(in) :: x(:,:) !! Fixed-effect design matrix `(n,p)`, including any desired intercept column.
      real(dp), intent(in) :: y(:) !! Gaussian response vector, size `n`.
      real(dp), intent(in) :: z(:,:) !! Random-effect design matrix `(n,q)` aligned row-for-row with `x` and `y`.
      integer, intent(in) :: group(:) !! Integer subject/group labels, size `n`; each group is one mixture-classification unit.
      integer, intent(in) :: k !! Requested initial number of mixture components.
      type(flexmix_result), intent(out) :: result !! Fitted Gaussian linear mixed-model mixture using upstream-style.
      type(flexmix_control), intent(in), optional :: control !! Optional EM iteration, minprior, tolerance, and classification.
      real(dp), intent(in), optional :: group_weights(:) !! Optional nonnegative weights for unique groups in first-occurrence.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based starting labels for the unique groups.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional starting group posterior probabilities `(ngroup,k)`.
      logical, intent(in), optional :: varfix_random !! Share the random-effect covariance across mixture components when true.
      logical, intent(in), optional :: varfix_residual !! Share residual variance across mixture components when true.
      call flexmix_mixed_common(x, y, z, group, k, flexmix_model_lmm, result, control, group_weights, &
                                initial_cluster, initial_posterior, varfix_random, varfix_residual)
   end subroutine flexmix_lmm

   subroutine flexmix_lmer(x, y, z, group, k, result, control, group_weights, initial_cluster, initial_posterior)
      real(dp), intent(in) :: x(:,:) !! Fixed-effect design matrix `(n,p)`, including any desired intercept column.
      real(dp), intent(in) :: y(:) !! Gaussian response vector, size `n`.
      real(dp), intent(in) :: z(:,:) !! Random-effect design matrix `(n,q)` aligned row-for-row with `x` and `y`.
      integer, intent(in) :: group(:) !! Integer subject/group labels, size `n`; each group is one mixture-classification unit.
      integer, intent(in) :: k !! Requested initial number of mixture components.
      type(flexmix_result), intent(out) :: result !! Fitted marginal Gaussian mixed-model mixture using translated EM updates.
      type(flexmix_control), intent(in), optional :: control !! Optional EM iteration, minprior, tolerance, and classification.
      real(dp), intent(in), optional :: group_weights(:) !! Optional nonnegative weights for unique groups in first-occurrence.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based starting labels for the unique groups.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional starting group posterior probabilities `(ngroup,k)`.
      call flexmix_mixed_common(x, y, z, group, k, flexmix_model_lmer, result, control, group_weights, &
                                initial_cluster, initial_posterior)
   end subroutine flexmix_lmer

   subroutine flexmix_lmmc(x, y, z, group, censored, k, result, control, group_weights, initial_cluster, &
                            initial_posterior, varfix_random, varfix_residual)
      real(dp), intent(in) :: x(:,:) !! Fixed-effect design matrix `(n,p)`, including any desired intercept column.
      real(dp), intent(in) :: y(:) !! Observed Gaussian responses; censored rows contain their left-censoring upper limits.
      real(dp), intent(in) :: z(:,:) !! Random-effect design matrix `(n,q)` aligned row-for-row with `x` and `y`.
      integer, intent(in) :: group(:) !! Integer subject/group labels, size `n`; each group is one mixture unit.
      logical, intent(in) :: censored(:) !! Left-censoring mask, size `n`; true rows store their upper censoring limits in `y`.
      integer, intent(in) :: k !! Requested initial number of mixture components.
      type(flexmix_result), intent(out) :: result !! Fitted censored Gaussian linear mixed-model mixture.
      type(flexmix_control), intent(in), optional :: control !! Optional EM iteration, minprior, tolerance, and classification.
      real(dp), intent(in), optional :: group_weights(:) !! Optional nonnegative weights for unique groups in group order.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based starting labels for the unique groups.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional starting group posterior probabilities `(ngroup,k)`.
      logical, intent(in), optional :: varfix_random !! Share the random-effect covariance across mixture components when true.
      logical, intent(in), optional :: varfix_residual !! Share residual variance across mixture components when true.
      type(flexmix_control) :: ctl
      integer, allocatable :: group_index(:), group_label(:), group_size(:), component_df(:)
      logical, allocatable :: active(:), first(:)
      real(dp), allocatable :: gw(:), group_signal(:), posterior(:,:), classified(:,:), prior(:), logp(:,:), logsum(:)
      real(dp), allocatable :: log_density(:,:), beta(:,:), residual_var(:), random_cov(:,:,:), old_beta(:)
      real(dp), allocatable :: latent_y(:), latent_var(:), cross_yb(:,:), random_mean(:,:), random_post_cov(:,:,:)
      real(dp) :: llh, old_llh
      logical :: fix_random, fix_residual
      integer :: fit_info, g, info, iter, j, n, ng, p, q, rank_value, wls_info

      ctl = flexmix_control()
      if (present(control)) ctl = control
      n = size(y)
      p = size(x,2)
      q = size(z,2)
      if (size(x,1) /= n .or. size(z,1) /= n .or. size(group) /= n .or. size(censored) /= n .or. &
          k < 1 .or. p < 1 .or. q < 1) then
         call mixed_invalid_result(result, flexmix_model_lmmc, -1)
         return
      end if
      call index_groups(group, group_index, group_label, group_size, first)
      ng = size(group_label)
      if (ng < 1) then
         call mixed_invalid_result(result, flexmix_model_lmmc, -2)
         return
      end if
      allocate(gw(ng), group_signal(ng), posterior(ng,k), classified(ng,k), prior(k), active(k))
      allocate(log_density(ng,k), logp(ng,k), logsum(ng), beta(p,k), residual_var(k), random_cov(q,q,k))
      allocate(latent_y(n), latent_var(n), cross_yb(n,q), random_mean(q,ng), random_post_cov(q,q,ng))
      allocate(component_df(k), old_beta(p))
      gw = 1.0_dp
      if (present(group_weights)) then
         if (size(group_weights) /= ng .or. any(group_weights < 0.0_dp)) then
            call mixed_invalid_result(result, flexmix_model_lmmc, -3)
            return
         end if
         gw = group_weights
      end if
      if (sum(gw) <= tiny(1.0_dp)) then
         call mixed_invalid_result(result, flexmix_model_lmmc, -4)
         return
      end if
      do g = 1, ng
         group_signal(g) = sum(pack(y,group_index == g)) / real(group_size(g),dp)
      end do
      call initialize_posteriors(group_signal, k, posterior, initial_cluster, initial_posterior, info)
      if (info /= 0) then
         call mixed_invalid_result(result, flexmix_model_lmmc, info)
         return
      end if
      prior = 1.0_dp / real(k,dp)
      active = .true.
      beta = 0.0_dp
      residual_var = max(1.0e-6_dp, sum((y - sum(y)/real(n,dp))**2) / real(max(1,n-1),dp))
      random_cov = 0.0_dp
      do j = 1, k
         do g = 1, q
            random_cov(g,g,j) = max(1.0e-4_dp, 0.5_dp * residual_var(j))
         end do
         do info = 1, n
            latent_var(info) = gw(group_index(info)) * posterior(group_index(info),j)
         end do
         call weighted_least_squares(x, y, latent_var, beta(:,j), rank_value, latent_y, wls_info)
         if (wls_info /= 0 .or. rank_value < 1) beta(:,j) = 0.0_dp
      end do
      fix_random = .false.
      fix_residual = .false.
      if (present(varfix_random)) fix_random = varfix_random
      if (present(varfix_residual)) fix_residual = varfix_residual
      old_llh = -huge(1.0_dp)
      llh = old_llh
      do iter = 1, max(1,ctl%iter_max)
         call classify_group_posterior(posterior, ctl%classify, classified)
         call update_group_prior(classified, gw, active, ctl%minprior, prior)
         call renormalize_active_groups(classified,active)
         do j = 1, k
            if (.not. active(j)) cycle
            old_beta = beta(:,j)
            call censored_random_moments(x, y, z, group_index, censored, beta(:,j), random_cov(:,:,j), &
                                         residual_var(j), latent_y, latent_var, random_mean, random_post_cov, &
                                         cross_yb, fit_info)
            if (fit_info /= 0) then
               active(j) = .false.
               cycle
            end if
            call fit_lmmc_component(x, z, group_index, gw * classified(:,j), latent_y, latent_var, random_mean, &
                                    random_post_cov, cross_yb, beta(:,j), random_cov(:,:,j), residual_var(j), &
                                    component_df(j), fit_info)
            if (fit_info /= 0) then
               active(j) = .false.
               beta(:,j) = old_beta
            end if
         end do
         call ensure_active_component(active,prior)
         if (fix_random) call share_random_covariance(random_cov,prior,active)
         if (fix_residual) call share_residual_variance(residual_var,prior,active)
         call evaluate_censored_group_likelihood(x, y, z, group_index, censored, active, beta, random_cov, &
                                                 residual_var, log_density)
         do j = 1, k
            if (active(j)) then
               logp(:,j) = log_density(:,j) + log(max(prior(j),tiny(1.0_dp)))
            else
               logp(:,j) = -huge(1.0_dp)
            end if
         end do
         call normalize_log_probabilities(logp,posterior,logsum)
         llh = dot_product(gw,logsum)
         if (iter > 1) then
            if (abs(llh-old_llh)/(abs(llh)+0.1_dp) < ctl%tolerance) exit
         end if
         old_llh = llh
      end do
      call pack_mixed_result(result, flexmix_model_lmmc, x, group_index, first, gw, active, prior, posterior, &
                             logp, beta, random_cov, residual_var, component_df, llh, iter, ctl)
   end subroutine flexmix_lmmc

   subroutine flexmix_lmc(x, y, group, censored, k, result, control, group_weights, initial_cluster, &
                           initial_posterior, varfix_residual)
      real(dp), intent(in) :: x(:,:) !! Fixed-effect design matrix `(n,p)`, including any desired intercept column.
      real(dp), intent(in) :: y(:) !! Observed Gaussian responses; censored rows contain their left-censoring upper limits.
      integer, intent(in) :: group(:) !! Integer subject/group labels, size `n`; each group is one mixture unit.
      logical, intent(in) :: censored(:) !! Left-censoring mask, size `n`; true rows store their upper censoring limits in `y`.
      integer, intent(in) :: k !! Requested initial number of mixture components.
      type(flexmix_result), intent(out) :: result !! Fitted censored Gaussian regression mixture without random effects.
      type(flexmix_control), intent(in), optional :: control !! Optional EM iteration, minprior, tolerance, and classification.
      real(dp), intent(in), optional :: group_weights(:) !! Optional nonnegative weights for unique groups in group order.
      integer, intent(in), optional :: initial_cluster(:) !! Optional one-based starting labels for the unique groups.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional starting group posterior probabilities `(ngroup,k)`.
      logical, intent(in), optional :: varfix_residual !! Share residual variance across mixture components when true.
      type(flexmix_control) :: ctl
      integer, allocatable :: group_index(:), group_label(:), group_size(:), component_df(:)
      logical, allocatable :: active(:), first(:)
      real(dp), allocatable :: gw(:), group_signal(:), posterior(:,:), classified(:,:), prior(:), logp(:,:), logsum(:)
      real(dp), allocatable :: log_density(:,:), beta(:,:), residual_var(:), latent_y(:), latent_var(:), row_weight(:)
      real(dp), allocatable :: fitted(:), old_beta(:)
      real(dp) :: llh, old_llh, denom
      logical :: fix_residual
      integer :: g, info, iter, j, n, ng, p, rank_value, wls_info

      ctl = flexmix_control()
      if (present(control)) ctl = control
      n = size(y)
      p = size(x,2)
      if (size(x,1) /= n .or. size(group) /= n .or. size(censored) /= n .or. k < 1 .or. p < 1) then
         call mixed_invalid_result(result, flexmix_model_lmc, -1)
         return
      end if
      call index_groups(group, group_index, group_label, group_size, first)
      ng = size(group_label)
      if (ng < 1) then
         call mixed_invalid_result(result, flexmix_model_lmc, -2)
         return
      end if
      allocate(gw(ng), group_signal(ng), posterior(ng,k), classified(ng,k), prior(k), active(k))
      allocate(log_density(ng,k), logp(ng,k), logsum(ng), beta(p,k), residual_var(k), component_df(k))
      allocate(latent_y(n), latent_var(n), row_weight(n), fitted(n), old_beta(p))
      gw = 1.0_dp
      if (present(group_weights)) then
         if (size(group_weights) /= ng .or. any(group_weights < 0.0_dp)) then
            call mixed_invalid_result(result, flexmix_model_lmc, -3)
            return
         end if
         gw = group_weights
      end if
      if (sum(gw) <= tiny(1.0_dp)) then
         call mixed_invalid_result(result, flexmix_model_lmc, -4)
         return
      end if
      do g = 1, ng
         group_signal(g) = sum(pack(y,group_index == g)) / real(group_size(g),dp)
      end do
      call initialize_posteriors(group_signal, k, posterior, initial_cluster, initial_posterior, info)
      if (info /= 0) then
         call mixed_invalid_result(result, flexmix_model_lmc, info)
         return
      end if
      prior = 1.0_dp / real(k,dp)
      active = .true.
      beta = 0.0_dp
      residual_var = max(1.0e-6_dp, sum((y - sum(y)/real(n,dp))**2) / real(max(1,n),dp))
      do j = 1, k
         do info = 1, n
            row_weight(info) = gw(group_index(info)) * posterior(group_index(info),j)
         end do
         call weighted_least_squares(x, y, row_weight, beta(:,j), rank_value, fitted, wls_info)
         if (wls_info /= 0 .or. rank_value < 1) then
            beta(:,j) = 0.0_dp
         else
            denom = sum(row_weight)
            if (denom > tiny(1.0_dp)) then
               residual_var(j) = max(1.0e-8_dp, dot_product(row_weight,fitted**2) / denom)
            end if
         end if
         component_df(j) = p + 1
      end do
      fix_residual = .false.
      if (present(varfix_residual)) fix_residual = varfix_residual
      old_llh = -huge(1.0_dp)
      llh = old_llh
      do iter = 1, max(1,ctl%iter_max)
         call classify_group_posterior(posterior, ctl%classify, classified)
         call update_group_prior(classified, gw, active, ctl%minprior, prior)
         call renormalize_active_groups(classified,active)
         do j = 1, k
            if (.not. active(j)) cycle
            old_beta = beta(:,j)
            call censored_fixed_moments(x, y, censored, beta(:,j), residual_var(j), latent_y, latent_var)
            do info = 1, n
               row_weight(info) = gw(group_index(info)) * classified(group_index(info),j)
            end do
            call weighted_least_squares(x, latent_y, row_weight, beta(:,j), rank_value, fitted, wls_info)
            if (wls_info /= 0 .or. rank_value < 1) then
               active(j) = .false.
               beta(:,j) = old_beta
               cycle
            end if
            denom = sum(row_weight)
            if (denom <= tiny(1.0_dp)) then
               active(j) = .false.
               beta(:,j) = old_beta
               cycle
            end if
            residual_var(j) = max(1.0e-10_dp, &
                                  dot_product(row_weight,fitted**2 + latent_var) / denom)
            component_df(j) = rank_value + 1
         end do
         call ensure_active_component(active,prior)
         if (fix_residual) call share_residual_variance(residual_var,prior,active)
         call evaluate_censored_fixed_likelihood(x, y, group_index, censored, active, beta, residual_var, log_density)
         do j = 1, k
            if (active(j)) then
               logp(:,j) = log_density(:,j) + log(max(prior(j),tiny(1.0_dp)))
            else
               logp(:,j) = -huge(1.0_dp)
            end if
         end do
         call normalize_log_probabilities(logp,posterior,logsum)
         llh = dot_product(gw,logsum)
         if (iter > 1) then
            if (abs(llh-old_llh)/(abs(llh)+0.1_dp) < ctl%tolerance) exit
         end if
         old_llh = llh
      end do
      call pack_lmc_result(result, x, group_index, first, gw, active, prior, posterior, logp, beta, residual_var, &
                           component_df, llh, iter, ctl)
   end subroutine flexmix_lmc

   subroutine flexmix_mixed_common(x, y, z, group, k, model_kind, result, control, group_weights, initial_cluster, &
                                   initial_posterior, varfix_random, varfix_residual)
      real(dp), intent(in) :: x(:,:) !! Fixed-effect design matrix `(n,p)`.
      real(dp), intent(in) :: y(:) !! Gaussian response vector, size `n`.
      real(dp), intent(in) :: z(:,:) !! Random-effect design matrix `(n,q)`.
      integer, intent(in) :: group(:) !! Integer group labels, size `n`.
      integer, intent(in) :: k !! Requested initial number of mixture components.
      integer, intent(in) :: model_kind !! Result family identifier distinguishing translated `lmm` and `lmer` entry points.
      type(flexmix_result), intent(out) :: result !! Fitted grouped mixed-model mixture result.
      type(flexmix_control), intent(in), optional :: control !! Optional EM controls.
      real(dp), intent(in), optional :: group_weights(:) !! Optional unique-group weights, size `ngroup`.
      integer, intent(in), optional :: initial_cluster(:) !! Optional starting labels for unique groups, size `ngroup`.
      real(dp), intent(in), optional :: initial_posterior(:,:) !! Optional starting posterior matrix `(ngroup,k)`.
      logical, intent(in), optional :: varfix_random !! Share random covariance across components when true.
      logical, intent(in), optional :: varfix_residual !! Share residual variance across components when true.
      type(flexmix_control) :: ctl
      integer, allocatable :: group_index(:), group_label(:), group_size(:), component_df(:)
      logical, allocatable :: active(:), first(:)
      real(dp), allocatable :: gw(:), group_signal(:), posterior(:,:), classified(:,:), prior(:), logp(:,:), logsum(:)
      real(dp), allocatable :: log_density(:,:), beta(:,:), residual_var(:), random_cov(:,:,:)
      real(dp), allocatable :: random_mean(:,:,:), random_post_cov(:,:,:,:), old_beta(:)
      real(dp) :: llh, old_llh
      logical :: fix_random, fix_residual
      integer :: fit_info, g, info, iter, j, n, ng, p, q

      ctl = flexmix_control()
      if (present(control)) ctl = control
      n = size(y)
      p = size(x,2)
      q = size(z,2)
      if (size(x,1) /= n .or. size(z,1) /= n .or. size(group) /= n .or. k < 1 .or. p < 1 .or. q < 1) then
         call mixed_invalid_result(result, model_kind, -1)
         return
      end if
      call index_groups(group, group_index, group_label, group_size, first)
      ng = size(group_label)
      if (ng < 1) then
         call mixed_invalid_result(result, model_kind, -2)
         return
      end if
      allocate(gw(ng), group_signal(ng), posterior(ng,k), classified(ng,k), prior(k), active(k))
      allocate(log_density(ng,k), logp(ng,k), logsum(ng), beta(p,k), residual_var(k), random_cov(q,q,k))
      allocate(random_mean(q,ng,k), random_post_cov(q,q,ng,k), component_df(k), old_beta(p))
      gw = 1.0_dp
      if (present(group_weights)) then
         if (size(group_weights) /= ng .or. any(group_weights < 0.0_dp)) then
            call mixed_invalid_result(result, model_kind, -3)
            return
         end if
         gw = group_weights
      end if
      if (sum(gw) <= tiny(1.0_dp)) then
         call mixed_invalid_result(result, model_kind, -4)
         return
      end if
      do g = 1, ng
         group_signal(g) = sum(pack(y,group_index == g)) / real(group_size(g),dp)
      end do
      call initialize_posteriors(group_signal, k, posterior, initial_cluster, initial_posterior, info)
      if (info /= 0) then
         call mixed_invalid_result(result, model_kind, info)
         return
      end if
      prior = 1.0_dp / real(k,dp)
      active = .true.
      beta = 0.0_dp
      residual_var = max(1.0e-6_dp, sum((y - sum(y)/real(n,dp))**2) / real(max(1,n-1),dp))
      random_cov = 0.0_dp
      random_mean = 0.0_dp
      random_post_cov = 0.0_dp
      do j = 1, k
         do g = 1, q
            random_cov(g,g,j) = max(1.0e-4_dp, 0.5_dp * residual_var(j))
         end do
         do g = 1, ng
            random_post_cov(:,:,g,j) = random_cov(:,:,j)
         end do
      end do
      fix_random = .false.
      fix_residual = .false.
      if (present(varfix_random)) fix_random = varfix_random
      if (present(varfix_residual)) fix_residual = varfix_residual
      old_llh = -huge(1.0_dp)
      llh = old_llh
      do iter = 1, max(1,ctl%iter_max)
         call classify_group_posterior(posterior, ctl%classify, classified)
         call update_group_prior(classified, gw, active, ctl%minprior, prior)
         call renormalize_active_groups(classified,active)
         do j = 1, k
            if (.not. active(j)) cycle
            old_beta = beta(:,j)
            call fit_lmm_component(x, y, z, group_index, gw * classified(:,j), random_mean(:,:,j), &
                                   random_post_cov(:,:,:,j), beta(:,j), random_cov(:,:,j), residual_var(j), &
                                   component_df(j), fit_info)
            if (fit_info /= 0) then
               active(j) = .false.
               beta(:,j) = old_beta
            end if
         end do
         call ensure_active_component(active,prior)
         if (fix_random) call share_random_covariance(random_cov,prior,active)
         if (fix_residual) call share_residual_variance(residual_var,prior,active)
         do j = 1, k
            if (.not. active(j)) cycle
            call update_random_posteriors(x, y, z, group_index, beta(:,j), random_cov(:,:,j), residual_var(j), &
                                          random_mean(:,:,j), random_post_cov(:,:,:,j), fit_info)
            if (fit_info /= 0) active(j) = .false.
         end do
         call ensure_active_component(active,prior)
         call evaluate_group_likelihood(x, y, z, group_index, active, beta, random_cov, residual_var, log_density)
         do j = 1, k
            if (active(j)) then
               logp(:,j) = log_density(:,j) + log(max(prior(j),tiny(1.0_dp)))
            else
               logp(:,j) = -huge(1.0_dp)
            end if
         end do
         call normalize_log_probabilities(logp,posterior,logsum)
         llh = dot_product(gw,logsum)
         if (iter > 1) then
            if (abs(llh-old_llh)/(abs(llh)+0.1_dp) < ctl%tolerance) then
               exit
            end if
         end if
         old_llh = llh
      end do
      call pack_mixed_result(result, model_kind, x, group_index, first, gw, active, prior, posterior, logp, beta, &
                             random_cov, residual_var, component_df, llh, iter, ctl)
   end subroutine flexmix_mixed_common

   subroutine index_groups(group, group_index, group_label, group_size, first)
      integer, intent(in) :: group(:) !! Row-level arbitrary integer group labels, size `n`.
      integer, allocatable, intent(out) :: group_index(:) !! Compact one-based group index for each row, size `n`.
      integer, allocatable, intent(out) :: group_label(:) !! Unique original labels in first-occurrence order.
      integer, allocatable, intent(out) :: group_size(:) !! Number of rows belonging to each compact group.
      logical, allocatable, intent(out) :: first(:) !! Row mask marking the first observation in every compact group.
      integer, allocatable :: labels_work(:), size_work(:)
      integer :: g, i, ng, n
      n = size(group)
      allocate(group_index(n), labels_work(n), size_work(n), first(n))
      labels_work = 0
      size_work = 0
      first = .false.
      ng = 0
      do i = 1, n
         g = 0
         if (ng > 0) then
            do g = 1, ng
               if (labels_work(g) == group(i)) exit
            end do
            if (g > ng) g = 0
         end if
         if (g == 0) then
            ng = ng + 1
            labels_work(ng) = group(i)
            g = ng
            first(i) = .true.
         end if
         group_index(i) = g
         size_work(g) = size_work(g) + 1
      end do
      allocate(group_label(ng), group_size(ng))
      group_label = labels_work(1:ng)
      group_size = size_work(1:ng)
   end subroutine index_groups

   subroutine classify_group_posterior(posterior, classify, classified)
      real(dp), intent(in) :: posterior(:,:) !! Current group posterior probabilities, shape `(ngroup,k)`.
      integer, intent(in) :: classify !! Classification mode; weighted unless equal to `flexmix_class_hard`.
      real(dp), intent(out) :: classified(:,:) !! Group M-step weights, shape `(ngroup,k)`.
      integer :: g, j
      classified = posterior
      if (classify /= flexmix_class_hard) return
      classified = 0.0_dp
      do g = 1, size(posterior,1)
         j = maxloc(posterior(g,:),dim=1)
         classified(g,j) = 1.0_dp
      end do
   end subroutine classify_group_posterior

   subroutine update_group_prior(classified, weights, active, minprior, prior)
      real(dp), intent(in) :: classified(:,:) !! Group classification weights, shape `(ngroup,k)`.
      real(dp), intent(in) :: weights(:) !! Nonnegative unique-group case weights, size `ngroup`.
      logical, intent(inout) :: active(:) !! Active-component mask updated for `minprior` deletion.
      real(dp), intent(in) :: minprior !! Minimum fitted component prior below which a component is deleted.
      real(dp), intent(out) :: prior(:) !! Updated mixture prior vector, size `k`.
      real(dp) :: total
      integer :: j
      total = sum(weights)
      if (total <= tiny(1.0_dp)) then
         prior = 0.0_dp
         prior(1) = 1.0_dp
         active = .false.
         active(1) = .true.
         return
      end if
      do j = 1, size(prior)
         if (active(j)) then
            prior(j) = dot_product(weights,classified(:,j)) / total
            if (prior(j) < minprior) active(j) = .false.
         else
            prior(j) = 0.0_dp
         end if
      end do
      call ensure_active_component(active,prior)
      if (sum(prior,mask=active) > tiny(1.0_dp)) then
         where (active) prior = prior / sum(prior,mask=active)
      end if
      where (.not. active) prior = 0.0_dp
   end subroutine update_group_prior

   subroutine renormalize_active_groups(classified, active)
      real(dp), intent(inout) :: classified(:,:) !! Group M-step weights renormalized over active components.
      logical, intent(in) :: active(:) !! Active-component mask, size `k`.
      real(dp) :: s
      integer :: g, j
      do g = 1, size(classified,1)
         do j = 1, size(active)
            if (.not. active(j)) classified(g,j) = 0.0_dp
         end do
         s = sum(classified(g,:))
         if (s > tiny(1.0_dp)) classified(g,:) = classified(g,:) / s
      end do
   end subroutine renormalize_active_groups

   subroutine ensure_active_component(active, prior)
      logical, intent(inout) :: active(:) !! Active-component mask guaranteed to retain at least one component.
      real(dp), intent(inout) :: prior(:) !! Component prior vector repaired consistently with `active`.
      integer :: j
      if (any(active)) return
      j = maxloc(prior,dim=1)
      if (j < 1) j = 1
      active(j) = .true.
      prior = 0.0_dp
      prior(j) = 1.0_dp
   end subroutine ensure_active_component

   subroutine fit_lmm_component(x, y, z, group_index, group_weight, random_mean, random_post_cov, beta, random_cov, &
                                residual_var, component_df, info)
      real(dp), intent(in) :: x(:,:) !! Fixed-effect design matrix `(n,p)`.
      real(dp), intent(in) :: y(:) !! Gaussian response vector, size `n`.
      real(dp), intent(in) :: z(:,:) !! Random-effect design matrix `(n,q)`.
      integer, intent(in) :: group_index(:) !! Compact one-based group index for each row, size `n`.
      real(dp), intent(in) :: group_weight(:) !! Current component weight for each unique group, size `ngroup`.
      real(dp), intent(in) :: random_mean(:,:) !! Previous posterior random-effect means, shape `(q,ngroup)`.
      real(dp), intent(in) :: random_post_cov(:,:,:) !! Previous posterior random-effect covariances, shape `(q,q,ngroup)`.
      real(dp), intent(out) :: beta(:) !! Updated fixed-effect coefficient vector, size `p`.
      real(dp), intent(out) :: random_cov(:,:) !! Updated random-effect covariance matrix, shape `(q,q)`.
      real(dp), intent(out) :: residual_var !! Updated positive residual variance.
      integer, intent(out) :: component_df !! Parameter count `p + q(q+1)/2 + 1`.
      integer, intent(out) :: info !! Zero on success; nonzero for degenerate weighted fits.
      real(dp), allocatable :: adjusted(:), row_weight(:), residual(:), covariance_sum(:,:)
      real(dp) :: denom, extra, trace_term
      integer :: g, i, q, rank_value, wls_info
      q = size(z,2)
      allocate(adjusted(size(y)), row_weight(size(y)), residual(size(y)), covariance_sum(q,q))
      do i = 1, size(y)
         g = group_index(i)
         adjusted(i) = y(i) - dot_product(z(i,:),random_mean(:,g))
         row_weight(i) = group_weight(g)
      end do
      call weighted_least_squares(x, adjusted, row_weight, beta, rank_value, residual, wls_info)
      if (wls_info /= 0 .or. rank_value < 1) then
         info = 1
         residual_var = 1.0_dp
         random_cov = 0.0_dp
         component_df = 0
         return
      end if
      residual = y - matmul(x,beta)
      do i = 1, size(y)
         residual(i) = residual(i) - dot_product(z(i,:),random_mean(:,group_index(i)))
      end do
      denom = sum(row_weight)
      extra = 0.0_dp
      do g = 1, size(group_weight)
         if (group_weight(g) <= 0.0_dp) cycle
         trace_term = trace_ztz_sigma(z, group_index, g, random_post_cov(:,:,g))
         extra = extra + group_weight(g) * trace_term
      end do
      if (denom <= tiny(1.0_dp)) then
         info = 2
         residual_var = 1.0_dp
         random_cov = 0.0_dp
         component_df = 0
         return
      end if
      residual_var = max(1.0e-10_dp, (sum(row_weight * residual**2) + extra) / denom)
      covariance_sum = 0.0_dp
      denom = sum(group_weight)
      if (denom <= tiny(1.0_dp)) then
         info = 3
         random_cov = 0.0_dp
         component_df = 0
         return
      end if
      do g = 1, size(group_weight)
         covariance_sum = covariance_sum + group_weight(g) * &
                          (random_post_cov(:,:,g) + outer_product(random_mean(:,g),random_mean(:,g)))
      end do
      random_cov = covariance_sum / denom
      call stabilize_covariance(random_cov)
      component_df = size(x,2) + q * (q + 1) / 2 + 1
      info = 0
   end subroutine fit_lmm_component

   function trace_ztz_sigma(z, group_index, group_value, sigma) result(value)
      real(dp), intent(in) :: z(:,:) !! Random-effect design matrix `(n,q)`.
      integer, intent(in) :: group_index(:) !! Compact row-level group indices, size `n`.
      integer, intent(in) :: group_value !! Unique-group index whose design contribution is requested.
      real(dp), intent(in) :: sigma(:,:) !! Posterior random-effect covariance matrix `(q,q)` for the group.
      real(dp) :: value
      real(dp), allocatable :: ztz(:,:)
      integer :: i
      allocate(ztz(size(z,2),size(z,2)))
      ztz = 0.0_dp
      do i = 1, size(group_index)
         if (group_index(i) == group_value) ztz = ztz + outer_product(z(i,:),z(i,:))
      end do
      value = sum(ztz * transpose(sigma))
   end function trace_ztz_sigma

   pure function outer_product(a, b) result(c)
      real(dp), intent(in) :: a(:) !! Left vector of the outer product.
      real(dp), intent(in) :: b(:) !! Right vector of the outer product.
      real(dp) :: c(size(a),size(b))
      integer :: i
      do i = 1, size(a)
         c(i,:) = a(i) * b
      end do
   end function outer_product

   subroutine stabilize_covariance(a)
      real(dp), intent(inout) :: a(:,:) !! Symmetric covariance matrix regularized in place to remain positive definite.
      real(dp), allocatable :: inverse(:,:)
      real(dp) :: jitter, logdet
      integer :: attempt, i, info
      a = 0.5_dp * (a + transpose(a))
      allocate(inverse(size(a,1),size(a,2)))
      jitter = 1.0e-10_dp * max(1.0_dp,maxval(abs(a)))
      do attempt = 1, 10
         call inverse_logdet_spd(a,inverse,logdet,info)
         if (info == 0) return
         do i = 1, size(a,1)
            a(i,i) = a(i,i) + jitter
         end do
         jitter = 10.0_dp * jitter
      end do
   end subroutine stabilize_covariance

   subroutine update_random_posteriors(x, y, z, group_index, beta, random_cov, residual_var, random_mean, &
                                       random_post_cov, info)
      real(dp), intent(in) :: x(:,:) !! Fixed-effect design matrix `(n,p)`.
      real(dp), intent(in) :: y(:) !! Gaussian response vector, size `n`.
      real(dp), intent(in) :: z(:,:) !! Random-effect design matrix `(n,q)`.
      integer, intent(in) :: group_index(:) !! Compact row-level group indices, size `n`.
      real(dp), intent(in) :: beta(:) !! Current fixed-effect coefficient vector, size `p`.
      real(dp), intent(in) :: random_cov(:,:) !! Current random-effect covariance matrix `(q,q)`.
      real(dp), intent(in) :: residual_var !! Current positive residual variance.
      real(dp), intent(out) :: random_mean(:,:) !! Updated posterior random-effect means, shape `(q,ngroup)`.
      real(dp), intent(out) :: random_post_cov(:,:,:) !! Updated posterior random-effect covariances, shape `(q,q,ngroup)`.
      integer, intent(out) :: info !! Zero on success; nonzero if covariance inversion fails.
      real(dp), allocatable :: psi_inv(:,:), precision(:,:), sigma(:,:), rhs(:), resid(:)
      real(dp) :: logdet
      integer :: g, i, inv_info, ng, q
      q = size(z,2)
      ng = size(random_mean,2)
      allocate(psi_inv(q,q), precision(q,q), sigma(q,q), rhs(q), resid(size(y)))
      call inverse_logdet_spd(random_cov,psi_inv,logdet,inv_info)
      if (inv_info /= 0) then
         info = 1
         return
      end if
      resid = y - matmul(x,beta)
      do g = 1, ng
         precision = psi_inv
         rhs = 0.0_dp
         do i = 1, size(y)
            if (group_index(i) /= g) cycle
            precision = precision + outer_product(z(i,:),z(i,:)) / residual_var
            rhs = rhs + z(i,:) * resid(i) / residual_var
         end do
         call inverse_logdet_spd(precision,sigma,logdet,inv_info)
         if (inv_info /= 0) then
            info = 2
            return
         end if
         random_post_cov(:,:,g) = sigma
         random_mean(:,g) = matmul(sigma,rhs)
      end do
      info = 0
   end subroutine update_random_posteriors

   subroutine share_random_covariance(random_cov, prior, active)
      real(dp), intent(inout) :: random_cov(:,:,:) !! Component random-effect covariance matrices, shape `(q,q,k)`.
      real(dp), intent(in) :: prior(:) !! Current component prior probabilities, size `k`.
      logical, intent(in) :: active(:) !! Active-component mask, size `k`.
      real(dp), allocatable :: shared(:,:)
      real(dp) :: denom
      integer :: j
      allocate(shared(size(random_cov,1),size(random_cov,2)))
      shared = 0.0_dp
      denom = 0.0_dp
      do j = 1, size(active)
         if (.not. active(j)) cycle
         shared = shared + prior(j) * random_cov(:,:,j)
         denom = denom + prior(j)
      end do
      if (denom > tiny(1.0_dp)) shared = shared / denom
      call stabilize_covariance(shared)
      do j = 1, size(active)
         if (active(j)) random_cov(:,:,j) = shared
      end do
   end subroutine share_random_covariance

   subroutine share_residual_variance(residual_var, prior, active)
      real(dp), intent(inout) :: residual_var(:) !! Component residual variances, size `k`.
      real(dp), intent(in) :: prior(:) !! Current component prior probabilities, size `k`.
      logical, intent(in) :: active(:) !! Active-component mask, size `k`.
      real(dp) :: denom, shared
      integer :: j
      shared = 0.0_dp
      denom = 0.0_dp
      do j = 1, size(active)
         if (.not. active(j)) cycle
         shared = shared + prior(j) * residual_var(j)
         denom = denom + prior(j)
      end do
      if (denom > tiny(1.0_dp)) shared = shared / denom
      shared = max(shared,1.0e-10_dp)
      do j = 1, size(active)
         if (active(j)) residual_var(j) = shared
      end do
   end subroutine share_residual_variance

   subroutine evaluate_group_likelihood(x, y, z, group_index, active, beta, random_cov, residual_var, log_density)
      real(dp), intent(in) :: x(:,:) !! Fixed-effect design matrix `(n,p)`.
      real(dp), intent(in) :: y(:) !! Gaussian response vector, size `n`.
      real(dp), intent(in) :: z(:,:) !! Random-effect design matrix `(n,q)`.
      integer, intent(in) :: group_index(:) !! Compact row-level group indices, size `n`.
      logical, intent(in) :: active(:) !! Active-component mask, size `k`.
      real(dp), intent(in) :: beta(:,:) !! Component fixed-effect coefficients, shape `(p,k)`.
      real(dp), intent(in) :: random_cov(:,:,:) !! Component random-effect covariance matrices, shape `(q,q,k)`.
      real(dp), intent(in) :: residual_var(:) !! Component residual variances, size `k`.
      real(dp), intent(out) :: log_density(:,:) !! Marginal group log densities, shape `(ngroup,k)`.
      real(dp), allocatable :: vg(:,:), inv(:,:), rg(:), zg(:,:), xg(:,:), yg(:)
      real(dp) :: logdet, quad
      integer :: g, i, info, j, m, ng
      ng = size(log_density,1)
      log_density = -huge(1.0_dp)
      do g = 1, ng
         m = count(group_index == g)
         allocate(vg(m,m),inv(m,m),rg(m),zg(m,size(z,2)),xg(m,size(x,2)),yg(m))
         i = 0
         do j = 1, size(y)
            if (group_index(j) /= g) cycle
            i = i + 1
            zg(i,:) = z(j,:)
            xg(i,:) = x(j,:)
            yg(i) = y(j)
         end do
         do j = 1, size(active)
            if (.not. active(j)) cycle
            vg = matmul(matmul(zg,random_cov(:,:,j)),transpose(zg))
            do i = 1, m
               vg(i,i) = vg(i,i) + residual_var(j)
            end do
            call inverse_logdet_spd(vg,inv,logdet,info)
            if (info /= 0) cycle
            rg = yg - matmul(xg,beta(:,j))
            quad = dot_product(rg,matmul(inv,rg))
            log_density(g,j) = -0.5_dp * (real(m,dp)*log(2.0_dp*pi) + logdet + quad)
         end do
         deallocate(vg,inv,rg,zg,xg,yg)
      end do
   end subroutine evaluate_group_likelihood

   subroutine censored_random_moments(x, y, z, group_index, censored, beta, random_cov, residual_var, &
                                      latent_y, latent_var, random_mean, random_post_cov, cross_yb, info)
      real(dp), intent(in) :: x(:,:) !! Fixed-effect design matrix `(n,p)`.
      real(dp), intent(in) :: y(:) !! Observed responses with censoring thresholds in censored rows.
      real(dp), intent(in) :: z(:,:) !! Random-effect design matrix `(n,q)`.
      integer, intent(in) :: group_index(:) !! Compact one-based group index for each row, size `n`.
      logical, intent(in) :: censored(:) !! Left-censoring mask, size `n`; multiple censored rows per group are allowed.
      real(dp), intent(in) :: beta(:) !! Current fixed-effect coefficient vector, size `p`.
      real(dp), intent(in) :: random_cov(:,:) !! Current random-effect covariance matrix `(q,q)`.
      real(dp), intent(in) :: residual_var !! Current positive residual variance.
      real(dp), intent(out) :: latent_y(:) !! Conditional response means, equal to observed values on uncensored rows.
      real(dp), intent(out) :: latent_var(:) !! Conditional response marginal variances; zero on uncensored rows.
      real(dp), intent(out) :: random_mean(:,:) !! Conditional random-effect means, shape `(q,ngroup)`.
      real(dp), intent(out) :: random_post_cov(:,:,:) !! Conditional random-effect covariance matrices `(q,q,ngroup)`.
      real(dp), intent(out) :: cross_yb(:,:) !! Conditional covariance of each response row with random effects `(n,q)`.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid covariance or numerical truncation failure.
      real(dp), allocatable :: vg(:,:), invv(:,:), zg(:,:), xg(:,:), yg(:), mug(:), gain(:,:), base_cov(:,:)
      real(dp), allocatable :: centered(:), voo(:,:), invoo(:,:), vco(:,:), ro(:), cond_mu(:), cond_cov(:,:)
      real(dp), allocatable :: trunc_mean(:), trunc_cov(:,:), upper(:), local_cov(:,:), cross_local(:,:)
      integer, allocatable :: idx(:), obs(:), cidx(:)
      real(dp) :: logcdf, logdet
      integer :: g, i, inv_info, j, m, mo, nc, ng, q, trunc_info

      ng = size(random_mean,2)
      q = size(z,2)
      latent_y = y
      latent_var = 0.0_dp
      random_mean = 0.0_dp
      random_post_cov = 0.0_dp
      cross_yb = 0.0_dp
      info = 0
      do g = 1, ng
         m = count(group_index == g)
         allocate(idx(m), vg(m,m), invv(m,m), zg(m,q), xg(m,size(x,2)), yg(m), mug(m), gain(q,m), &
                  base_cov(q,q), centered(m), local_cov(m,m), cross_local(m,q))
         j = 0
         do i = 1, size(y)
            if (group_index(i) /= g) cycle
            j = j + 1
            idx(j) = i
            zg(j,:) = z(i,:)
            xg(j,:) = x(i,:)
            yg(j) = y(i)
         end do
         mug = matmul(xg,beta)
         vg = matmul(matmul(zg,random_cov),transpose(zg))
         do i = 1, m
            vg(i,i) = vg(i,i) + residual_var
         end do
         call inverse_logdet_spd(vg,invv,logdet,inv_info)
         if (inv_info /= 0) then
            info = 1
            return
         end if
         nc = count(censored(idx))
         mo = m - nc
         allocate(cidx(nc),obs(mo))
         j = 0
         do i = 1, m
            if (.not. censored(idx(i))) cycle
            j = j + 1
            cidx(j) = i
         end do
         j = 0
         do i = 1, m
            if (censored(idx(i))) cycle
            j = j + 1
            obs(j) = i
         end do
         local_cov = 0.0_dp
         if (nc > 0) then
            allocate(cond_mu(nc), cond_cov(nc,nc), trunc_mean(nc), trunc_cov(nc,nc), upper(nc))
            upper = yg(cidx)
            if (mo > 0) then
               allocate(voo(mo,mo), invoo(mo,mo), vco(nc,mo), ro(mo))
               voo = vg(obs,obs)
               call inverse_logdet_spd(voo,invoo,logdet,inv_info)
               if (inv_info /= 0) then
                  info = 2
                  return
               end if
               vco = vg(cidx,obs)
               ro = yg(obs) - mug(obs)
               cond_mu = mug(cidx) + matmul(vco,matmul(invoo,ro))
               cond_cov = vg(cidx,cidx) - matmul(matmul(vco,invoo),transpose(vco))
            else
               cond_mu = mug(cidx)
               cond_cov = vg(cidx,cidx)
            end if
            call truncated_normal_upper_moments(cond_mu,cond_cov,upper,trunc_mean,trunc_cov,logcdf,trunc_info)
            if (trunc_info /= 0) then
               info = 3
               return
            end if
            do i = 1, nc
               latent_y(idx(cidx(i))) = trunc_mean(i)
               latent_var(idx(cidx(i))) = max(0.0_dp,trunc_cov(i,i))
            end do
            local_cov(cidx,cidx) = trunc_cov
         end if
         gain = matmul(matmul(random_cov,transpose(zg)),invv)
         base_cov = random_cov - matmul(matmul(gain,zg),random_cov)
         base_cov = 0.5_dp * (base_cov + transpose(base_cov))
         centered = latent_y(idx) - mug
         random_mean(:,g) = matmul(gain,centered)
         random_post_cov(:,:,g) = base_cov + matmul(matmul(gain,local_cov),transpose(gain))
         cross_local = matmul(local_cov,transpose(gain))
         do i = 1, m
            cross_yb(idx(i),:) = cross_local(i,:)
         end do
         call stabilize_covariance(random_post_cov(:,:,g))
         deallocate(idx,vg,invv,zg,xg,yg,mug,gain,base_cov,centered,local_cov,cross_local)
         deallocate(obs,cidx)
         if (allocated(voo)) deallocate(voo)
         if (allocated(invoo)) deallocate(invoo)
         if (allocated(vco)) deallocate(vco)
         if (allocated(ro)) deallocate(ro)
         if (allocated(cond_mu)) deallocate(cond_mu)
         if (allocated(cond_cov)) deallocate(cond_cov)
         if (allocated(trunc_mean)) deallocate(trunc_mean)
         if (allocated(trunc_cov)) deallocate(trunc_cov)
         if (allocated(upper)) deallocate(upper)
      end do
   end subroutine censored_random_moments

   subroutine fit_lmmc_component(x, z, group_index, group_weight, latent_y, latent_var, random_mean, &
                                 random_post_cov, cross_yb, beta, random_cov, residual_var, component_df, info)
      real(dp), intent(in) :: x(:,:) !! Fixed-effect design matrix `(n,p)`.
      real(dp), intent(in) :: z(:,:) !! Random-effect design matrix `(n,q)`.
      integer, intent(in) :: group_index(:) !! Compact one-based group index for each row, size `n`.
      real(dp), intent(in) :: group_weight(:) !! Current component weight for each unique group, size `ngroup`.
      real(dp), intent(in) :: latent_y(:) !! Conditional response means from the censoring E-step, size `n`.
      real(dp), intent(in) :: latent_var(:) !! Conditional response marginal variances, size `n`.
      real(dp), intent(in) :: random_mean(:,:) !! Conditional random-effect means, shape `(q,ngroup)`.
      real(dp), intent(in) :: random_post_cov(:,:,:) !! Conditional random-effect covariances `(q,q,ngroup)`.
      real(dp), intent(in) :: cross_yb(:,:) !! Conditional response/random-effect covariances `(n,q)`.
      real(dp), intent(out) :: beta(:) !! Updated fixed-effect coefficient vector, size `p`.
      real(dp), intent(out) :: random_cov(:,:) !! Updated random-effect covariance matrix `(q,q)`.
      real(dp), intent(out) :: residual_var !! Updated positive residual variance.
      integer, intent(out) :: component_df !! Parameter count `p + q(q+1)/2 + 1`.
      integer, intent(out) :: info !! Zero on success; nonzero for a degenerate weighted fit.
      real(dp), allocatable :: adjusted(:), row_weight(:), residual(:), covariance_sum(:,:)
      real(dp) :: denom_group, denom_row, correction, variance_term
      integer :: g, i, q, rank_value, wls_info

      q = size(z,2)
      allocate(adjusted(size(latent_y)), row_weight(size(latent_y)), residual(size(latent_y)), covariance_sum(q,q))
      do i = 1, size(latent_y)
         g = group_index(i)
         adjusted(i) = latent_y(i) - dot_product(z(i,:),random_mean(:,g))
         row_weight(i) = group_weight(g)
      end do
      call weighted_least_squares(x, adjusted, row_weight, beta, rank_value, residual, wls_info)
      if (wls_info /= 0 .or. rank_value < 1) then
         info = 1
         random_cov = 0.0_dp
         residual_var = 1.0_dp
         component_df = 0
         return
      end if
      denom_group = sum(group_weight)
      denom_row = sum(row_weight)
      if (denom_group <= tiny(1.0_dp) .or. denom_row <= tiny(1.0_dp)) then
         info = 2
         random_cov = 0.0_dp
         residual_var = 1.0_dp
         component_df = 0
         return
      end if
      covariance_sum = 0.0_dp
      do g = 1, size(group_weight)
         covariance_sum = covariance_sum + group_weight(g) * &
                          (random_post_cov(:,:,g) + outer_product(random_mean(:,g),random_mean(:,g)))
      end do
      random_cov = covariance_sum / denom_group
      call stabilize_covariance(random_cov)
      residual = latent_y - matmul(x,beta)
      variance_term = 0.0_dp
      correction = 0.0_dp
      do i = 1, size(latent_y)
         g = group_index(i)
         residual(i) = residual(i) - dot_product(z(i,:),random_mean(:,g))
         variance_term = variance_term + group_weight(g) * &
                         (latent_var(i) + dot_product(z(i,:),matmul(random_post_cov(:,:,g),z(i,:))))
         correction = correction + group_weight(g) * dot_product(cross_yb(i,:),z(i,:))
      end do
      residual_var = max(1.0e-10_dp, &
                         (sum(row_weight * residual**2) + variance_term - 2.0_dp * correction) / denom_row)
      component_df = size(x,2) + q * (q + 1) / 2 + 1
      info = 0
   end subroutine fit_lmmc_component

   subroutine evaluate_censored_group_likelihood(x, y, z, group_index, censored, active, beta, random_cov, &
                                                 residual_var, log_density)
      real(dp), intent(in) :: x(:,:) !! Fixed-effect design matrix `(n,p)`.
      real(dp), intent(in) :: y(:) !! Observed responses with censoring thresholds in censored rows.
      real(dp), intent(in) :: z(:,:) !! Random-effect design matrix `(n,q)`.
      integer, intent(in) :: group_index(:) !! Compact one-based group index for each row, size `n`.
      logical, intent(in) :: censored(:) !! Left-censoring mask, size `n`; multiple censored rows per group are allowed.
      logical, intent(in) :: active(:) !! Active-component mask, size `k`.
      real(dp), intent(in) :: beta(:,:) !! Component fixed-effect coefficients, shape `(p,k)`.
      real(dp), intent(in) :: random_cov(:,:,:) !! Component random-effect covariance matrices `(q,q,k)`.
      real(dp), intent(in) :: residual_var(:) !! Component residual variances, size `k`.
      real(dp), intent(out) :: log_density(:,:) !! Group observed-data log likelihoods under censoring.
      real(dp), allocatable :: vg(:,:), xg(:,:), zg(:,:), yg(:), mug(:), voo(:,:), invoo(:,:), ro(:), vco(:,:)
      real(dp), allocatable :: cond_mu(:), cond_cov(:,:), trunc_mean(:), trunc_cov(:,:), upper(:)
      integer, allocatable :: idx(:), obs(:), cidx(:)
      real(dp) :: logcdf, logdet, quad
      integer :: g, i, info, j, m, mo, nc, ng, trunc_info

      ng = size(log_density,1)
      log_density = -huge(1.0_dp)
      do g = 1, ng
         m = count(group_index == g)
         allocate(idx(m), xg(m,size(x,2)), zg(m,size(z,2)), yg(m), mug(m), vg(m,m))
         i = 0
         do j = 1, size(y)
            if (group_index(j) /= g) cycle
            i = i + 1
            idx(i) = j
            xg(i,:) = x(j,:)
            zg(i,:) = z(j,:)
            yg(i) = y(j)
         end do
         nc = count(censored(idx))
         mo = m - nc
         allocate(cidx(nc),obs(mo))
         i = 0
         do j = 1, m
            if (.not. censored(idx(j))) cycle
            i = i + 1
            cidx(i) = j
         end do
         i = 0
         do j = 1, m
            if (censored(idx(j))) cycle
            i = i + 1
            obs(i) = j
         end do
         do j = 1, size(active)
            if (.not. active(j)) cycle
            mug = matmul(xg,beta(:,j))
            vg = matmul(matmul(zg,random_cov(:,:,j)),transpose(zg))
            do i = 1, m
               vg(i,i) = vg(i,i) + residual_var(j)
            end do
            if (nc == 0) then
               allocate(invoo(m,m), ro(m))
               call inverse_logdet_spd(vg,invoo,logdet,info)
               if (info == 0) then
                  ro = yg - mug
                  quad = dot_product(ro,matmul(invoo,ro))
                  log_density(g,j) = -0.5_dp * (real(m,dp)*log(2.0_dp*pi) + logdet + quad)
               end if
               deallocate(invoo,ro)
            else
               allocate(cond_mu(nc), cond_cov(nc,nc), trunc_mean(nc), trunc_cov(nc,nc), upper(nc))
               upper = yg(cidx)
               if (mo > 0) then
                  allocate(voo(mo,mo), invoo(mo,mo), ro(mo), vco(nc,mo))
                  voo = vg(obs,obs)
                  call inverse_logdet_spd(voo,invoo,logdet,info)
                  if (info == 0) then
                     ro = yg(obs) - mug(obs)
                     quad = dot_product(ro,matmul(invoo,ro))
                     vco = vg(cidx,obs)
                     cond_mu = mug(cidx) + matmul(vco,matmul(invoo,ro))
                     cond_cov = vg(cidx,cidx) - matmul(matmul(vco,invoo),transpose(vco))
                     call truncated_normal_upper_moments(cond_mu,cond_cov,upper,trunc_mean,trunc_cov, &
                                                         logcdf,trunc_info)
                     if (trunc_info == 0) then
                        log_density(g,j) = -0.5_dp * (real(mo,dp)*log(2.0_dp*pi) + logdet + quad) + logcdf
                     end if
                  end if
                  deallocate(voo,invoo,ro,vco)
               else
                  cond_mu = mug(cidx)
                  cond_cov = vg(cidx,cidx)
                  call truncated_normal_upper_moments(cond_mu,cond_cov,upper,trunc_mean,trunc_cov,logcdf,trunc_info)
                  if (trunc_info == 0) log_density(g,j) = logcdf
               end if
               deallocate(cond_mu,cond_cov,trunc_mean,trunc_cov,upper)
            end if
         end do
         deallocate(idx,xg,zg,yg,mug,vg)
         deallocate(obs,cidx)
      end do
   end subroutine evaluate_censored_group_likelihood

   subroutine truncated_normal_upper_moments(mu, covariance, upper, mean_value, covariance_value, log_probability, info)
      real(dp), intent(in) :: mu(:) !! Mean vector of the Gaussian distribution being upper-truncated.
      real(dp), intent(in) :: covariance(:,:) !! Positive-definite covariance matrix of the Gaussian distribution.
      real(dp), intent(in) :: upper(:) !! Componentwise upper truncation limits, size equal to `mu`.
      real(dp), intent(out) :: mean_value(:) !! Conditional mean vector given all coordinates are below `upper`.
      real(dp), intent(out) :: covariance_value(:,:) !! Conditional covariance matrix under the upper truncation event.
      real(dp), intent(out) :: log_probability !! Log probability of the multivariate upper-truncation event.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid shape, covariance, or unsupported dimension.
      integer, parameter :: nsample = 1024
      integer, parameter :: primes(32) = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53, &
                                          59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131]
      real(dp), allocatable :: lower(:,:), zscore(:), sample(:), sum1(:), sum2(:,:)
      real(dp) :: alpha, cdf, logcdf, mills, p, total_weight, u, variance_factor, weight
      integer :: chol_info, d, h, i

      d = size(mu)
      if (d < 1 .or. size(upper) /= d .or. size(covariance,1) /= d .or. size(covariance,2) /= d .or. &
          size(mean_value) /= d .or. size(covariance_value,1) /= d .or. size(covariance_value,2) /= d) then
         info = 1
         return
      end if
      if (d == 1) then
         if (covariance(1,1) <= tiny(1.0_dp)) then
            info = 2
            return
         end if
         alpha = (upper(1)-mu(1))/sqrt(covariance(1,1))
         call left_truncated_standard_moments(alpha,logcdf,mills,variance_factor)
         mean_value(1) = mu(1) - sqrt(covariance(1,1))*mills
         covariance_value(1,1) = covariance(1,1)*variance_factor
         log_probability = logcdf
         info = 0
         return
      end if
      if (d > size(primes)) then
         info = 3
         return
      end if
      allocate(lower(d,d), zscore(d), sample(d), sum1(d), sum2(d,d))
      call cholesky_lower(covariance,lower,chol_info)
      if (chol_info /= 0) then
         info = 4
         return
      end if
      sum1 = 0.0_dp
      sum2 = 0.0_dp
      total_weight = 0.0_dp
      do h = 1, nsample
         zscore = 0.0_dp
         weight = 1.0_dp
         do i = 1, d
            alpha = upper(i) - mu(i)
            if (i > 1) alpha = alpha - dot_product(lower(i,1:i-1),zscore(1:i-1))
            alpha = alpha / lower(i,i)
            call left_truncated_standard_moments(alpha,logcdf,mills,variance_factor)
            if (logcdf < log(tiny(1.0_dp))) then
               cdf = tiny(1.0_dp)
            else
               cdf = exp(logcdf)
            end if
            weight = weight * cdf
            u = halton_value(h+31,primes(i))
            p = max(tiny(1.0_dp),min(1.0_dp-epsilon(1.0_dp),u*cdf))
            zscore(i) = inverse_standard_normal(p)
         end do
         sample = mu + matmul(lower,zscore)
         total_weight = total_weight + weight
         sum1 = sum1 + weight*sample
         sum2 = sum2 + weight*outer_product(sample,sample)
      end do
      if (total_weight <= tiny(1.0_dp)) then
         info = 5
         return
      end if
      mean_value = sum1 / total_weight
      covariance_value = sum2 / total_weight - outer_product(mean_value,mean_value)
      covariance_value = 0.5_dp * (covariance_value + transpose(covariance_value))
      do i = 1, d
         covariance_value(i,i) = max(0.0_dp,covariance_value(i,i))
      end do
      log_probability = log(total_weight) - log(real(nsample,dp))
      info = 0
   end subroutine truncated_normal_upper_moments

   subroutine cholesky_lower(a, lower, info)
      real(dp), intent(in) :: a(:,:) !! Symmetric positive-definite matrix to factor.
      real(dp), intent(out) :: lower(:,:) !! Lower-triangular Cholesky factor satisfying `a = lower*transpose(lower)`.
      integer, intent(out) :: info !! Zero on success; positive pivot index when the matrix is not positive definite.
      real(dp) :: value
      integer :: i, j, k, n

      n = size(a,1)
      lower = 0.0_dp
      info = 0
      if (size(a,2) /= n .or. size(lower,1) /= n .or. size(lower,2) /= n) then
         info = -1
         return
      end if
      do i = 1, n
         do j = 1, i
            value = a(i,j)
            do k = 1, j-1
               value = value - lower(i,k)*lower(j,k)
            end do
            if (i == j) then
               if (value <= tiny(1.0_dp)) then
                  info = i
                  return
               end if
               lower(i,j) = sqrt(value)
            else
               lower(i,j) = value / lower(j,j)
            end if
         end do
      end do
   end subroutine cholesky_lower

   pure function halton_value(index_value, base) result(value)
      integer, intent(in) :: index_value !! Positive sequence index for the radical-inverse Halton coordinate.
      integer, intent(in) :: base !! Prime radix used for this Halton coordinate; must exceed one.
      real(dp) :: value
      real(dp) :: factor
      integer :: index_work

      value = 0.0_dp
      if (index_value < 1 .or. base < 2) return
      factor = 1.0_dp / real(base,dp)
      index_work = index_value
      do while (index_work > 0)
         value = value + factor*real(mod(index_work,base),dp)
         index_work = index_work / base
         factor = factor / real(base,dp)
      end do
   end function halton_value

   pure function inverse_standard_normal(probability) result(value)
      real(dp), intent(in) :: probability !! Probability in `(0,1)` whose standard-normal quantile is requested.
      real(dp) :: value
      real(dp), parameter :: a1 = -3.969683028665376e1_dp
      real(dp), parameter :: a2 = 2.209460984245205e2_dp
      real(dp), parameter :: a3 = -2.759285104469687e2_dp
      real(dp), parameter :: a4 = 1.383577518672690e2_dp
      real(dp), parameter :: a5 = -3.066479806614716e1_dp
      real(dp), parameter :: a6 = 2.506628277459239_dp
      real(dp), parameter :: b1 = -5.447609879822406e1_dp
      real(dp), parameter :: b2 = 1.615858368580409e2_dp
      real(dp), parameter :: b3 = -1.556989798598866e2_dp
      real(dp), parameter :: b4 = 6.680131188771972e1_dp
      real(dp), parameter :: b5 = -1.328068155288572e1_dp
      real(dp), parameter :: c1 = -7.784894002430293e-3_dp
      real(dp), parameter :: c2 = -3.223964580411365e-1_dp
      real(dp), parameter :: c3 = -2.400758277161838_dp
      real(dp), parameter :: c4 = -2.549732539343734_dp
      real(dp), parameter :: c5 = 4.374664141464968_dp
      real(dp), parameter :: c6 = 2.938163982698783_dp
      real(dp), parameter :: d1 = 7.784695709041462e-3_dp
      real(dp), parameter :: d2 = 3.224671290700398e-1_dp
      real(dp), parameter :: d3 = 2.445134137142996_dp
      real(dp), parameter :: d4 = 3.754408661907416_dp
      real(dp), parameter :: plow = 2.425e-2_dp
      real(dp), parameter :: phigh = 1.0_dp - plow
      real(dp) :: p, q, r

      p = max(tiny(1.0_dp),min(1.0_dp-epsilon(1.0_dp),probability))
      if (p < plow) then
         q = sqrt(-2.0_dp*log(p))
         value = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
      else if (p > phigh) then
         q = sqrt(-2.0_dp*log(1.0_dp-p))
         value = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
      else
         q = p - 0.5_dp
         r = q*q
         value = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q / &
                 (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0_dp)
      end if
   end function inverse_standard_normal

   subroutine left_truncated_standard_moments(alpha, logcdf, mills, variance_factor)
      real(dp), intent(in) :: alpha !! Standardized upper truncation point `(T-mu)/sigma`.
      real(dp), intent(out) :: logcdf !! Natural logarithm of the standard-normal CDF at `alpha`.
      real(dp), intent(out) :: mills !! Inverse Mills ratio `phi(alpha)/Phi(alpha)`.
      real(dp), intent(out) :: variance_factor !! Variance of a standard normal conditional on `Z <= alpha`.
      real(dp) :: cdf, invx, x
      real(dp), parameter :: log_sqrt_2pi = 0.5_dp * log(2.0_dp*pi)
      real(dp), parameter :: sqrt2 = sqrt(2.0_dp)

      if (alpha > -8.0_dp) then
         cdf = max(0.5_dp * erfc(-alpha/sqrt2), tiny(1.0_dp))
         logcdf = log(cdf)
         mills = exp(-0.5_dp*alpha*alpha - log_sqrt_2pi - logcdf)
      else
         x = -alpha
         invx = 1.0_dp / x
         mills = x + invx - 2.0_dp*invx**3 + 10.0_dp*invx**5 - 74.0_dp*invx**7
         logcdf = -0.5_dp*alpha*alpha - log_sqrt_2pi - log(mills)
      end if
      variance_factor = max(0.0_dp,1.0_dp - alpha*mills - mills*mills)
   end subroutine left_truncated_standard_moments

   subroutine censored_fixed_moments(x, y, censored, beta, residual_var, latent_y, latent_var)
      real(dp), intent(in) :: x(:,:) !! Fixed-effect design matrix `(n,p)`.
      real(dp), intent(in) :: y(:) !! Observed responses or left-censoring limits, size `n`.
      logical, intent(in) :: censored(:) !! Left-censoring mask, size `n`.
      real(dp), intent(in) :: beta(:) !! Current fixed-effect coefficient vector, size `p`.
      real(dp), intent(in) :: residual_var !! Current positive Gaussian residual variance.
      real(dp), intent(out) :: latent_y(:) !! Conditional first moments of the latent responses, size `n`.
      real(dp), intent(out) :: latent_var(:) !! Conditional variances of the latent responses, size `n`.
      real(dp) :: alpha, logcdf, mills, mu, sigma, variance_factor
      integer :: i

      sigma = sqrt(max(residual_var,tiny(1.0_dp)))
      do i = 1, size(y)
         if (censored(i)) then
            mu = dot_product(x(i,:),beta)
            alpha = (y(i)-mu) / sigma
            call left_truncated_standard_moments(alpha,logcdf,mills,variance_factor)
            latent_y(i) = mu - sigma*mills
            latent_var(i) = residual_var*variance_factor
         else
            latent_y(i) = y(i)
            latent_var(i) = 0.0_dp
         end if
      end do
   end subroutine censored_fixed_moments

   subroutine evaluate_censored_fixed_likelihood(x, y, group_index, censored, active, beta, residual_var, log_density)
      real(dp), intent(in) :: x(:,:) !! Fixed-effect design matrix `(n,p)`.
      real(dp), intent(in) :: y(:) !! Observed responses or left-censoring limits, size `n`.
      integer, intent(in) :: group_index(:) !! Compact one-based group index for each row, size `n`.
      logical, intent(in) :: censored(:) !! Left-censoring mask, size `n`.
      logical, intent(in) :: active(:) !! Active-component mask, size `k`.
      real(dp), intent(in) :: beta(:,:) !! Fixed-effect coefficients `(p,k)`.
      real(dp), intent(in) :: residual_var(:) !! Positive residual variances, size `k`.
      real(dp), intent(out) :: log_density(:,:) !! Group observed-data log likelihoods, shape `(ngroup,k)`.
      real(dp) :: alpha, logcdf, mills, mu, sigma2, variance_factor
      integer :: g, i, j

      log_density = 0.0_dp
      do j = 1, size(active)
         if (.not. active(j)) then
            log_density(:,j) = -huge(1.0_dp)
            cycle
         end if
         sigma2 = max(residual_var(j),tiny(1.0_dp))
         do i = 1, size(y)
            g = group_index(i)
            mu = dot_product(x(i,:),beta(:,j))
            if (censored(i)) then
               alpha = (y(i)-mu) / sqrt(sigma2)
               call left_truncated_standard_moments(alpha,logcdf,mills,variance_factor)
               log_density(g,j) = log_density(g,j) + logcdf
            else
               log_density(g,j) = log_density(g,j) - 0.5_dp * &
                                  (log(2.0_dp*pi*sigma2) + (y(i)-mu)**2/sigma2)
            end if
         end do
      end do
   end subroutine evaluate_censored_fixed_likelihood

   subroutine pack_lmc_result(result, x, group_index, first, group_weights, active, prior, posterior, logp, beta, &
                              residual_var, component_df, loglik, iterations, control)
      type(flexmix_result), intent(out) :: result !! Destination fitted censored fixed-effect mixture result.
      real(dp), intent(in) :: x(:,:) !! Fixed-effect design matrix `(n,p)` used by the fit.
      integer, intent(in) :: group_index(:) !! Compact row-level group index, size `n`.
      logical, intent(in) :: first(:) !! Row mask marking each group's first observation.
      real(dp), intent(in) :: group_weights(:) !! Unique-group observation weights, size `ngroup`.
      logical, intent(in) :: active(:) !! Active-component mask over the original `k` components.
      real(dp), intent(in) :: prior(:) !! Original component prior vector, size `k`.
      real(dp), intent(in) :: posterior(:,:) !! Final group posterior matrix `(ngroup,k)`.
      real(dp), intent(in) :: logp(:,:) !! Final unnormalized group log posterior weights `(ngroup,k)`.
      real(dp), intent(in) :: beta(:,:) !! Original fixed-effect coefficients `(p,k)`.
      real(dp), intent(in) :: residual_var(:) !! Original residual variances, size `k`.
      integer, intent(in) :: component_df(:) !! Per-component numerical parameter counts before packing.
      real(dp), intent(in) :: loglik !! Maximized grouped observed-data log likelihood.
      integer, intent(in) :: iterations !! Number of EM iterations executed.
      type(flexmix_control), intent(in) :: control !! Controls used by the fit, including convergence tolerance.
      integer :: a, g, i, j, kk, n

      n = size(x,1)
      kk = count(active)
      result%model_kind = flexmix_model_lmc
      result%n = n
      result%p = size(x,2)
      result%d = 1
      result%k0 = size(active)
      result%k = kk
      result%iterations = min(iterations,max(1,control%iter_max))
      result%loglik = loglik
      result%status = 0
      result%converged = iterations <= max(1,control%iter_max)
      allocate(result%prior(kk),result%posterior(n,kk),result%posterior_unscaled(n,kk))
      allocate(result%log_posterior_unscaled(n,kk),result%cluster(n),result%size(kk),result%component_role(kk))
      allocate(result%beta(size(beta,1),kk),result%sigma(kk),result%residual_variance(kk))
      allocate(result%group_first(n),result%case_weights(n))
      result%component_role = 0
      result%group_first = first
      do i = 1, n
         result%case_weights(i) = group_weights(group_index(i))
      end do
      a = 0
      do j = 1, size(active)
         if (.not. active(j)) cycle
         a = a + 1
         result%prior(a) = prior(j)
         result%beta(:,a) = beta(:,j)
         result%residual_variance(a) = residual_var(j)
         result%sigma(a) = sqrt(residual_var(j))
         do i = 1, n
            g = group_index(i)
            result%posterior(i,a) = posterior(g,j)
            result%log_posterior_unscaled(i,a) = logp(g,j)
            if (logp(g,j) > log(tiny(1.0_dp))) then
               result%posterior_unscaled(i,a) = exp(logp(g,j))
            else
               result%posterior_unscaled(i,a) = 0.0_dp
            end if
         end do
      end do
      do i = 1, n
         result%cluster(i) = maxloc(result%posterior(i,:),dim=1)
      end do
      do j = 1, kk
         result%size(j) = count(result%cluster == j)
      end do
      result%df = max(0,kk-1)
      do j = 1, size(active)
         if (active(j)) result%df = result%df + component_df(j)
      end do
   end subroutine pack_lmc_result

   subroutine pack_mixed_result(result, model_kind, x, group_index, first, group_weights, active, prior, posterior, logp, beta, &
                                random_cov, residual_var, component_df, loglik, iterations, control)
      type(flexmix_result), intent(out) :: result !! Destination fitted mixed-model result.
      integer, intent(in) :: model_kind !! Mixed-model family identifier stored in the result.
      real(dp), intent(in) :: x(:,:) !! Fixed-effect design matrix `(n,p)` used by the fit.
      integer, intent(in) :: group_index(:) !! Compact row-level group index, size `n`.
      logical, intent(in) :: first(:) !! Row mask marking each group's first observation.
      real(dp), intent(in) :: group_weights(:) !! Unique-group observation weights, size `ngroup`.
      logical, intent(in) :: active(:) !! Active-component mask over the original `k` components.
      real(dp), intent(in) :: prior(:) !! Original component prior vector, size `k`.
      real(dp), intent(in) :: posterior(:,:) !! Final group posterior matrix `(ngroup,k)`.
      real(dp), intent(in) :: logp(:,:) !! Final unnormalized group log posterior weights `(ngroup,k)`.
      real(dp), intent(in) :: beta(:,:) !! Original fixed-effect coefficients `(p,k)`.
      real(dp), intent(in) :: random_cov(:,:,:) !! Original random-effect covariance matrices `(q,q,k)`.
      real(dp), intent(in) :: residual_var(:) !! Original residual variances, size `k`.
      integer, intent(in) :: component_df(:) !! Per-component numerical parameter counts before packing.
      real(dp), intent(in) :: loglik !! Maximized grouped observed-data log likelihood.
      integer, intent(in) :: iterations !! Number of EM iterations executed.
      type(flexmix_control), intent(in) :: control !! Controls used by the fit, including convergence tolerance.
      integer :: a, g, i, j, kk, n
      n = size(x,1)
      kk = count(active)
      result%model_kind = model_kind
      result%n = n
      result%p = size(x,2)
      result%d = 1
      result%k0 = size(active)
      result%k = kk
      result%iterations = min(iterations,max(1,control%iter_max))
      result%loglik = loglik
      result%status = 0
      result%converged = iterations <= max(1,control%iter_max)
      allocate(result%prior(kk),result%posterior(n,kk),result%posterior_unscaled(n,kk))
      allocate(result%log_posterior_unscaled(n,kk),result%cluster(n),result%size(kk),result%component_role(kk))
      allocate(result%beta(size(beta,1),kk),result%sigma(kk),result%random_covariance(size(random_cov,1),size(random_cov,2),kk))
      allocate(result%residual_variance(kk),result%group_first(n),result%case_weights(n))
      result%component_role = 0
      result%group_first = first
      do i = 1, n
         result%case_weights(i) = group_weights(group_index(i))
      end do
      a = 0
      do j = 1, size(active)
         if (.not. active(j)) cycle
         a = a + 1
         result%prior(a) = prior(j)
         result%beta(:,a) = beta(:,j)
         result%random_covariance(:,:,a) = random_cov(:,:,j)
         result%residual_variance(a) = residual_var(j)
         result%sigma(a) = sqrt(residual_var(j))
         do i = 1, n
            g = group_index(i)
            result%posterior(i,a) = posterior(g,j)
            result%log_posterior_unscaled(i,a) = logp(g,j)
            if (logp(g,j) > log(tiny(1.0_dp))) then
               result%posterior_unscaled(i,a) = exp(logp(g,j))
            else
               result%posterior_unscaled(i,a) = 0.0_dp
            end if
         end do
      end do
      do i = 1, n
         result%cluster(i) = maxloc(result%posterior(i,:),dim=1)
      end do
      do j = 1, kk
         result%size(j) = count(result%cluster == j)
      end do
      result%df = 0
      a = 0
      do j = 1, size(active)
         if (.not. active(j)) cycle
         a = a + 1
         result%df = result%df + component_df(j)
      end do
      result%df = result%df + max(0,kk-1)
   end subroutine pack_mixed_result

   subroutine mixed_invalid_result(result, model_kind, status)
      type(flexmix_result), intent(out) :: result !! Result object initialized to an invalid-fit state.
      integer, intent(in) :: model_kind !! Mixed-model family identifier to preserve in the failed result.
      integer, intent(in) :: status !! Negative validation or numerical status code.
      result%model_kind = model_kind
      result%status = status
      result%converged = .false.
   end subroutine mixed_invalid_result

end module flexmix_mixed
