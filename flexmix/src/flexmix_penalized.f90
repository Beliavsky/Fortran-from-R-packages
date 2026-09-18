! SPDX-License-Identifier: GPL-2.0-or-later
module flexmix_penalized
   use flexmix_kinds, only : dp
   use flexmix_types, only : flexmix_model_gaussian, flexmix_model_poisson, flexmix_model_binomial
   use flexmix_components, only : fit_gaussian_regression_component, fit_poisson_regression_component
   use flexmix_components, only : fit_binomial_regression_component
   implicit none
   private
   public :: fit_glmnet_component

contains

   subroutine fit_glmnet_component(x, y, trials, weights, model_kind, beta, sigma, component_df, selected_lambda, info, &
                                   adaptive, select, alpha, nfolds, fold_id, offset, lambda_grid)
      real(dp), intent(in) :: x(:,:) !! Regression design matrix `(n,p)` whose first column must be an intercept of ones.
      real(dp), intent(in) :: y(:) !! Gaussian/Poisson response, or binomial successes when `model_kind` is binomial.
      real(dp), intent(in), optional :: trials(:) !! Binomial trial totals, size `n`; required only for binomial regression.
      real(dp), intent(in) :: weights(:) !! Nonnegative observation weights, size `n`, including current EM responsibilities.
      integer, intent(in) :: model_kind !! Regression family identifier: Gaussian, Poisson, or binomial.
      real(dp), intent(out) :: beta(:) !! Penalized fitted coefficient vector, size `p`, including the unpenalized intercept.
      real(dp), intent(out) :: sigma !! Gaussian residual scale; set to one for non-Gaussian families.
      integer, intent(out) :: component_df !! Number of nonzero coefficients plus one Gaussian scale parameter when applicable.
      real(dp), intent(out) :: selected_lambda !! Cross-validation penalty minimizing weighted held-out deviance.
      integer, intent(out) :: info !! Zero on success; negative for invalid inputs and positive for numerical failure.
      logical, intent(in), optional :: adaptive !! Use adaptive-lasso penalty factors from the unpenalized weighted fit; default.
      logical, intent(in), optional :: select(:) !! Selection mask for non-intercept coefficients, size `p-1`; false means.
      real(dp), intent(in), optional :: alpha !! Elastic-net mixing value in `[0,1]`; one is lasso and zero is ridge.
      integer, intent(in), optional :: nfolds !! Number of deterministic CV folds when `fold_id` is absent; default five.
      integer, intent(in), optional :: fold_id(:) !! Optional fixed one-based CV fold labels, size `n`, reused at every EM.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive linear-predictor offset, size `n`.
      real(dp), intent(in), optional :: lambda_grid(:) !! Optional nonnegative penalty grid; generated when absent.
      real(dp), allocatable :: penalty(:), lambdas(:), cv(:), beta0(:), residuals(:), work_trials(:)
      real(dp), allocatable :: train_weights(:), candidate(:), best_beta(:)
      real(dp) :: alpha_value, best, sw, meanw, eps_coef
      logical :: adaptive_value
      logical, allocatable :: select_value(:)
      integer, allocatable :: folds(:)
      integer :: f, fit_info, fold_count, j, l, n, p, unpen_df

      n = size(y)
      p = size(x,2)
      beta = 0.0_dp
      sigma = 1.0_dp
      component_df = 0
      selected_lambda = 0.0_dp
      info = 0
      if (size(x,1) /= n .or. size(weights) /= n .or. size(beta) /= p .or. p < 1 .or. n < 2) then
         info = -1
         return
      end if
      if (any(weights < 0.0_dp) .or. sum(weights) <= tiny(1.0_dp)) then
         info = -2
         return
      end if
      if (maxval(abs(x(:,1) - 1.0_dp)) > 1.0e-10_dp) then
         info = -3
         return
      end if
      if (model_kind /= flexmix_model_gaussian .and. model_kind /= flexmix_model_poisson .and. &
          model_kind /= flexmix_model_binomial) then
         info = -4
         return
      end if
      if (model_kind == flexmix_model_binomial) then
         if (.not. present(trials)) then
            info = -5
            return
         end if
         if (size(trials) /= n .or. any(trials < y) .or. any(y < 0.0_dp)) then
            info = -6
            return
         end if
      end if
      if (present(offset)) then
         if (size(offset) /= n) then
            info = -7
            return
         end if
      end if

      adaptive_value = .true.
      if (present(adaptive)) adaptive_value = adaptive
      alpha_value = 1.0_dp
      if (present(alpha)) alpha_value = min(1.0_dp, max(0.0_dp, alpha))
      allocate(select_value(max(0,p-1)))
      select_value = .true.
      if (present(select)) then
         if (size(select) /= p - 1) then
            info = -8
            return
         end if
         select_value = select
      end if
      allocate(penalty(p), beta0(p), residuals(n), work_trials(n))
      penalty = 0.0_dp
      beta0 = 0.0_dp
      work_trials = 1.0_dp
      if (present(trials)) work_trials = trials

      if (adaptive_value .and. p > 1 .and. any(select_value)) then
         select case (model_kind)
         case (flexmix_model_gaussian)
            call fit_gaussian_regression_component(x, y, weights, beta0, sigma, unpen_df, fit_info, offset=offset)
         case (flexmix_model_poisson)
            call fit_poisson_regression_component(x, y, weights, beta0, unpen_df, fit_info, offset=offset)
         case (flexmix_model_binomial)
            call fit_binomial_regression_component(x, y, work_trials, weights, beta0, unpen_df, fit_info, offset=offset)
         end select
         if (fit_info /= 0) beta0 = 0.0_dp
      end if
      meanw = sum(weights) / real(n,dp)
      eps_coef = sqrt(epsilon(1.0_dp))
      do j = 2, p
         if (.not. select_value(j-1)) then
            penalty(j) = 0.0_dp
         else if (adaptive_value) then
            if (abs(beta0(j)) <= eps_coef) then
               penalty(j) = 1.0e12_dp
            else
               penalty(j) = meanw / abs(beta0(j))
            end if
         else
            penalty(j) = 1.0_dp
         end if
      end do

      if (.not. any(select_value)) then
         call fit_unpenalized(x, y, work_trials, weights, model_kind, beta, sigma, component_df, info, offset)
         return
      end if

      if (present(lambda_grid)) then
         if (size(lambda_grid) < 1 .or. any(lambda_grid < 0.0_dp)) then
            info = -9
            return
         end if
         allocate(lambdas(size(lambda_grid)))
         lambdas = lambda_grid
      else
         call make_lambda_grid(x, y, work_trials, weights, model_kind, penalty, alpha_value, lambdas, info, offset)
         if (info /= 0) return
      end if
      fold_count = 5
      if (present(nfolds)) fold_count = max(2, nfolds)
      allocate(folds(n))
      if (present(fold_id)) then
         if (size(fold_id) /= n .or. minval(fold_id) < 1) then
            info = -10
            return
         end if
         folds = fold_id
         fold_count = maxval(folds)
      else
         do j = 1, n
            folds(j) = 1 + modulo(j - 1, min(fold_count,n))
         end do
         fold_count = min(fold_count,n)
      end if
      if (fold_count < 2) then
         info = -11
         return
      end if

      allocate(cv(size(lambdas)), train_weights(n), candidate(p), best_beta(p))
      cv = 0.0_dp
      do f = 1, fold_count
         train_weights = weights
         where (folds == f) train_weights = 0.0_dp
         if (sum(train_weights) <= tiny(1.0_dp)) cycle
         candidate = 0.0_dp
         do l = 1, size(lambdas)
            call fit_penalized_at_lambda(x, y, work_trials, train_weights, model_kind, penalty, alpha_value, &
                                         lambdas(l), candidate, fit_info, offset)
            if (fit_info /= 0) then
               cv(l) = huge(1.0_dp)
            else if (cv(l) < 0.5_dp * huge(1.0_dp)) then
               cv(l) = cv(l) + validation_loss(x, y, work_trials, weights, folds, f, model_kind, candidate, offset)
            end if
         end do
      end do
      best = huge(1.0_dp)
      l = 1
      do j = 1, size(lambdas)
         if (cv(j) < best) then
            best = cv(j)
            l = j
         end if
      end do
      selected_lambda = lambdas(l)
      best_beta = 0.0_dp
      do j = 1, l
         call fit_penalized_at_lambda(x, y, work_trials, weights, model_kind, penalty, alpha_value, &
                                      lambdas(j), best_beta, fit_info, offset)
         if (fit_info /= 0) then
            info = 20 + fit_info
            return
         end if
      end do
      beta = best_beta
      component_df = count(abs(beta) > 1.0e-9_dp)
      if (model_kind == flexmix_model_gaussian) then
         residuals = y - matmul(x,beta)
         if (present(offset)) residuals = residuals - offset
         sw = sum(weights)
         meanw = sw / real(n,dp)
         if (n > component_df .and. meanw > tiny(1.0_dp)) then
            sigma = sqrt(max(tiny(1.0_dp), sum(weights * residuals**2 / meanw) / real(n - component_df,dp)))
         else
            sigma = sqrt(max(tiny(1.0_dp), sum(weights * residuals**2) / max(sw,tiny(1.0_dp))))
         end if
         component_df = component_df + 1
      else
         sigma = 1.0_dp
      end if
   end subroutine fit_glmnet_component

   subroutine fit_unpenalized(x, y, trials, weights, model_kind, beta, sigma, component_df, info, offset)
      real(dp), intent(in) :: x(:,:) !! Regression design matrix `(n,p)` including the intercept column.
      real(dp), intent(in) :: y(:) !! Family response vector, or binomial successes.
      real(dp), intent(in) :: trials(:) !! Binomial trial totals; ignored for Gaussian and Poisson families.
      real(dp), intent(in) :: weights(:) !! Nonnegative fitting weights, size `n`.
      integer, intent(in) :: model_kind !! Regression family identifier.
      real(dp), intent(out) :: beta(:) !! Unpenalized coefficient vector, size `p`.
      real(dp), intent(out) :: sigma !! Gaussian residual scale or one for non-Gaussian models.
      integer, intent(out) :: component_df !! Upstream-style component parameter count.
      integer, intent(out) :: info !! Zero on success; otherwise component-fit error code.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed linear-predictor offset, size `n`.
      select case (model_kind)
      case (flexmix_model_gaussian)
         call fit_gaussian_regression_component(x, y, weights, beta, sigma, component_df, info, offset=offset)
      case (flexmix_model_poisson)
         call fit_poisson_regression_component(x, y, weights, beta, component_df, info, offset=offset)
         sigma = 1.0_dp
      case (flexmix_model_binomial)
         call fit_binomial_regression_component(x, y, trials, weights, beta, component_df, info, offset=offset)
         sigma = 1.0_dp
      case default
         info = -1
         beta = 0.0_dp
         sigma = 1.0_dp
         component_df = 0
      end select
   end subroutine fit_unpenalized

   subroutine make_lambda_grid(x, y, trials, weights, model_kind, penalty, alpha, lambdas, info, offset)
      real(dp), intent(in) :: x(:,:) !! Regression design matrix `(n,p)` including an intercept in column one.
      real(dp), intent(in) :: y(:) !! Family response vector or binomial successes.
      real(dp), intent(in) :: trials(:) !! Binomial trial totals; ignored for Gaussian and Poisson families.
      real(dp), intent(in) :: weights(:) !! Nonnegative observation weights used for the component fit.
      integer, intent(in) :: model_kind !! Regression family identifier.
      real(dp), intent(in) :: penalty(:) !! Nonnegative coefficient-specific penalty factors, size `p`.
      real(dp), intent(in) :: alpha !! Elastic-net mixing fraction in `[0,1]`.
      real(dp), allocatable, intent(out) :: lambdas(:) !! Generated decreasing penalty grid.
      integer, intent(out) :: info !! Zero on success; nonzero if the null fit cannot be formed.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed linear-predictor offset, size `n`.
      real(dp), allocatable :: beta(:), grad(:), eta(:), mu(:), residual(:)
      real(dp) :: lam_max, ratio, sw, pbar
      integer :: j, m, n, p

      n = size(y)
      p = size(x,2)
      allocate(beta(p), grad(p), eta(n), mu(n), residual(n))
      beta = 0.0_dp
      sw = sum(weights)
      if (sw <= tiny(1.0_dp)) then
         info = 1
         allocate(lambdas(1))
         lambdas = 0.0_dp
         return
      end if
      select case (model_kind)
      case (flexmix_model_gaussian)
         beta(1) = dot_product(weights, y) / sw
         if (present(offset)) beta(1) = beta(1) - dot_product(weights,offset) / sw
      case (flexmix_model_poisson)
         if (present(offset)) then
            beta(1) = log(max(tiny(1.0_dp), dot_product(weights,y) / &
                      max(tiny(1.0_dp), dot_product(weights,exp(min(30.0_dp,max(-30.0_dp,offset)))))))
         else
            beta(1) = log(max(tiny(1.0_dp), dot_product(weights,y) / sw))
         end if
      case (flexmix_model_binomial)
         pbar = min(1.0_dp - 1.0e-10_dp, max(1.0e-10_dp, dot_product(weights,y) / &
                max(tiny(1.0_dp),dot_product(weights,trials))))
         beta(1) = log(pbar / (1.0_dp - pbar))
         call optimize_intercept_offset(y, trials, weights, beta(1), offset)
      case default
         info = 2
         allocate(lambdas(1))
         lambdas = 0.0_dp
         return
      end select
      call smooth_gradient(x, y, trials, weights, model_kind, beta, 0.0_dp, penalty, grad, eta, mu, offset)
      lam_max = 0.0_dp
      do j = 2, p
         if (penalty(j) > 0.0_dp .and. penalty(j) < 0.5_dp * huge(1.0_dp)) then
            lam_max = max(lam_max, abs(grad(j)) / max(1.0e-6_dp, alpha * penalty(j)))
         end if
      end do
      if (lam_max <= 100.0_dp * epsilon(1.0_dp)) lam_max = 1.0_dp
      m = 32
      ratio = 1.0e-4_dp
      if (n < p) ratio = 1.0e-2_dp
      allocate(lambdas(m))
      do j = 1, m
         lambdas(j) = lam_max * exp(log(ratio) * real(j - 1,dp) / real(m - 1,dp))
      end do
      info = 0
   end subroutine make_lambda_grid

   subroutine optimize_intercept_offset(success, trials, weights, intercept, offset)
      real(dp), intent(in) :: success(:) !! Binomial success counts, size `n`.
      real(dp), intent(in) :: trials(:) !! Binomial trial totals, size `n`.
      real(dp), intent(in) :: weights(:) !! Nonnegative observation weights, size `n`.
      real(dp), intent(inout) :: intercept !! Intercept updated by scalar Newton iterations.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed logit offset, size `n`.
      real(dp) :: eta, p, score, hess, delta
      integer :: i, iter
      if (.not. present(offset)) return
      do iter = 1, 50
         score = 0.0_dp
         hess = 0.0_dp
         do i = 1, size(success)
            eta = min(30.0_dp, max(-30.0_dp, intercept + offset(i)))
            p = 1.0_dp / (1.0_dp + exp(-eta))
            score = score + weights(i) * (success(i) - trials(i) * p)
            hess = hess - weights(i) * trials(i) * p * (1.0_dp - p)
         end do
         if (abs(hess) <= tiny(1.0_dp)) exit
         delta = score / hess
         intercept = intercept - delta
         if (abs(delta) <= 1.0e-10_dp * (1.0_dp + abs(intercept))) exit
      end do
   end subroutine optimize_intercept_offset

   subroutine fit_penalized_at_lambda(x, y, trials, weights, model_kind, penalty, alpha, lambda_value, beta, info, offset)
      real(dp), intent(in) :: x(:,:) !! Regression design matrix `(n,p)` including intercept column one.
      real(dp), intent(in) :: y(:) !! Family response vector or binomial successes.
      real(dp), intent(in) :: trials(:) !! Binomial trial totals; ignored for Gaussian and Poisson families.
      real(dp), intent(in) :: weights(:) !! Nonnegative fitting weights, size `n`.
      integer, intent(in) :: model_kind !! Regression family identifier.
      real(dp), intent(in) :: penalty(:) !! Coefficient penalty factors, size `p`, with zero for the intercept/unselected terms.
      real(dp), intent(in) :: alpha !! Elastic-net mixing fraction in `[0,1]`.
      real(dp), intent(in) :: lambda_value !! Nonnegative penalty strength.
      real(dp), intent(inout) :: beta(:) !! Warm-start coefficient vector on entry and fitted coefficients on return.
      integer, intent(out) :: info !! Zero on convergence; positive if iterations or line search fail.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive linear-predictor offset, size `n`.
      real(dp), allocatable :: grad(:), eta(:), mu(:), candidate(:), eta2(:), mu2(:)
      real(dp) :: old_obj, new_obj, prev_obj, step, threshold, delta_max
      integer :: iter, j, line_iter, p

      p = size(beta)
      allocate(grad(p), eta(size(y)), mu(size(y)), candidate(p), eta2(size(y)), mu2(size(y)))
      call smooth_gradient(x, y, trials, weights, model_kind, beta, lambda_value*(1.0_dp-alpha), penalty, &
                           grad, eta, mu, offset)
      old_obj = total_objective(y, trials, weights, model_kind, eta, beta, penalty, alpha, lambda_value)
      info = 1
      do iter = 1, 1000
         step = 1.0_dp
         do line_iter = 1, 40
            candidate = beta - step * grad
            do j = 2, p
               threshold = step * lambda_value * alpha * penalty(j)
               if (candidate(j) > threshold) then
                  candidate(j) = candidate(j) - threshold
               else if (candidate(j) < -threshold) then
                  candidate(j) = candidate(j) + threshold
               else
                  candidate(j) = 0.0_dp
               end if
            end do
            call linear_predictor(x, candidate, eta2, offset)
            call response_mean(eta2, model_kind, mu2)
            new_obj = total_objective(y, trials, weights, model_kind, eta2, candidate, penalty, alpha, lambda_value)
            if (new_obj <= old_obj + 1.0e-12_dp * (1.0_dp + abs(old_obj))) exit
            step = 0.5_dp * step
         end do
         if (line_iter > 40) then
            info = 2
            return
         end if
         delta_max = maxval(abs(candidate - beta))
         prev_obj = old_obj
         beta = candidate
         old_obj = new_obj
         call smooth_gradient(x, y, trials, weights, model_kind, beta, lambda_value*(1.0_dp-alpha), penalty, &
                              grad, eta, mu, offset)
         if (delta_max <= 2.0e-8_dp * (1.0_dp + maxval(abs(beta))) .or. &
             abs(new_obj - prev_obj) <= 1.0e-11_dp * (1.0_dp + abs(prev_obj))) then
            info = 0
            return
         end if
      end do
   end subroutine fit_penalized_at_lambda

   subroutine smooth_gradient(x, y, trials, weights, model_kind, beta, ridge_lambda, penalty, grad, eta, mu, offset)
      real(dp), intent(in) :: x(:,:) !! Regression design matrix `(n,p)`.
      real(dp), intent(in) :: y(:) !! Family response vector or binomial successes.
      real(dp), intent(in) :: trials(:) !! Binomial trial totals; ignored for Gaussian and Poisson families.
      real(dp), intent(in) :: weights(:) !! Nonnegative fitting weights, size `n`.
      integer, intent(in) :: model_kind !! Regression family identifier.
      real(dp), intent(in) :: beta(:) !! Current coefficient vector, size `p`.
      real(dp), intent(in) :: ridge_lambda !! Ridge multiplier already including `(1-alpha)`.
      real(dp), intent(in) :: penalty(:) !! Coefficient penalty factors, size `p`.
      real(dp), intent(out) :: grad(:) !! Gradient of the smooth normalized objective, size `p`.
      real(dp), intent(out) :: eta(:) !! Current linear predictor, size `n`.
      real(dp), intent(out) :: mu(:) !! Current response mean or Gaussian predictor, size `n`.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive linear-predictor offset, size `n`.
      real(dp), allocatable :: score(:)
      real(dp) :: sw
      integer :: j
      call linear_predictor(x, beta, eta, offset)
      call response_mean(eta, model_kind, mu)
      allocate(score(size(y)))
      select case (model_kind)
      case (flexmix_model_gaussian)
         score = weights * (mu - y)
      case (flexmix_model_poisson)
         score = weights * (mu - y)
      case (flexmix_model_binomial)
         score = weights * (trials * mu - y)
      end select
      sw = max(sum(weights), tiny(1.0_dp))
      grad = matmul(transpose(x), score) / sw
      do j = 2, size(beta)
         grad(j) = grad(j) + ridge_lambda * penalty(j) * beta(j)
      end do
   end subroutine smooth_gradient

   pure subroutine linear_predictor(x, beta, eta, offset)
      real(dp), intent(in) :: x(:,:) !! Regression design matrix `(n,p)`.
      real(dp), intent(in) :: beta(:) !! Coefficient vector, size `p`.
      real(dp), intent(out) :: eta(:) !! Linear predictor, size `n`.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive predictor offset, size `n`.
      eta = matmul(x,beta)
      if (present(offset)) eta = eta + offset
   end subroutine linear_predictor

   pure subroutine response_mean(eta, model_kind, mu)
      real(dp), intent(in) :: eta(:) !! Family linear predictor, size `n`.
      integer, intent(in) :: model_kind !! Regression family identifier.
      real(dp), intent(out) :: mu(:) !! Family response mean, size `n`.
      integer :: i
      select case (model_kind)
      case (flexmix_model_gaussian)
         mu = eta
      case (flexmix_model_poisson)
         do i = 1, size(eta)
            mu(i) = exp(min(30.0_dp, max(-30.0_dp, eta(i))))
         end do
      case (flexmix_model_binomial)
         do i = 1, size(eta)
            if (eta(i) >= 0.0_dp) then
               mu(i) = 1.0_dp / (1.0_dp + exp(-min(30.0_dp,eta(i))))
            else
               mu(i) = exp(max(-30.0_dp,eta(i))) / (1.0_dp + exp(max(-30.0_dp,eta(i))))
            end if
         end do
      end select
   end subroutine response_mean

   pure function total_objective(y, trials, weights, model_kind, eta, beta, penalty, alpha, lambda_value) result(value)
      real(dp), intent(in) :: y(:) !! Family response vector or binomial successes.
      real(dp), intent(in) :: trials(:) !! Binomial trial totals; ignored for Gaussian and Poisson families.
      real(dp), intent(in) :: weights(:) !! Nonnegative fitting weights, size `n`.
      integer, intent(in) :: model_kind !! Regression family identifier.
      real(dp), intent(in) :: eta(:) !! Current linear predictor, size `n`.
      real(dp), intent(in) :: beta(:) !! Current coefficient vector, size `p`.
      real(dp), intent(in) :: penalty(:) !! Coefficient penalty factors, size `p`.
      real(dp), intent(in) :: alpha !! Elastic-net mixing fraction.
      real(dp), intent(in) :: lambda_value !! Nonnegative penalty strength.
      real(dp) :: value
      real(dp) :: sw, smooth, l1, ridge, e
      integer :: i, j
      sw = max(sum(weights), tiny(1.0_dp))
      smooth = 0.0_dp
      select case (model_kind)
      case (flexmix_model_gaussian)
         smooth = 0.5_dp * sum(weights * (y - eta)**2) / sw
      case (flexmix_model_poisson)
         do i = 1, size(y)
            e = min(30.0_dp, max(-30.0_dp,eta(i)))
            smooth = smooth + weights(i) * (exp(e) - y(i) * e)
         end do
         smooth = smooth / sw
      case (flexmix_model_binomial)
         do i = 1, size(y)
            e = eta(i)
            if (e > 0.0_dp) then
               smooth = smooth + weights(i) * (trials(i) * (e + log(1.0_dp + exp(-e))) - y(i) * e)
            else
               smooth = smooth + weights(i) * (trials(i) * log(1.0_dp + exp(e)) - y(i) * e)
            end if
         end do
         smooth = smooth / sw
      end select
      l1 = 0.0_dp
      ridge = 0.0_dp
      do j = 2, size(beta)
         l1 = l1 + penalty(j) * abs(beta(j))
         ridge = ridge + penalty(j) * beta(j)**2
      end do
      value = smooth + lambda_value * (alpha * l1 + 0.5_dp * (1.0_dp - alpha) * ridge)
   end function total_objective

   function validation_loss(x, y, trials, weights, folds, fold_value, model_kind, beta, offset) result(value)
      real(dp), intent(in) :: x(:,:) !! Regression design matrix `(n,p)`.
      real(dp), intent(in) :: y(:) !! Family response vector or binomial successes.
      real(dp), intent(in) :: trials(:) !! Binomial trial totals; ignored for Gaussian and Poisson families.
      real(dp), intent(in) :: weights(:) !! Base observation weights, size `n`.
      integer, intent(in) :: folds(:) !! One-based fold labels, size `n`.
      integer, intent(in) :: fold_value !! Held-out fold label to evaluate.
      integer, intent(in) :: model_kind !! Regression family identifier.
      real(dp), intent(in) :: beta(:) !! Candidate coefficient vector, size `p`.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive linear-predictor offset, size `n`.
      real(dp) :: value
      real(dp), allocatable :: eta(:), mu(:)
      real(dp) :: sw, e
      integer :: i
      allocate(eta(size(y)), mu(size(y)))
      call linear_predictor(x,beta,eta,offset)
      call response_mean(eta,model_kind,mu)
      value = 0.0_dp
      sw = 0.0_dp
      do i = 1, size(y)
         if (folds(i) /= fold_value .or. weights(i) <= 0.0_dp) cycle
         sw = sw + weights(i)
         select case (model_kind)
         case (flexmix_model_gaussian)
            value = value + weights(i) * (y(i) - eta(i))**2
         case (flexmix_model_poisson)
            e = min(30.0_dp,max(-30.0_dp,eta(i)))
            value = value + 2.0_dp * weights(i) * (exp(e) - y(i) * e)
         case (flexmix_model_binomial)
            e = eta(i)
            if (e > 0.0_dp) then
               value = value + 2.0_dp * weights(i) * (trials(i) * (e + log(1.0_dp + exp(-e))) - y(i) * e)
            else
               value = value + 2.0_dp * weights(i) * (trials(i) * log(1.0_dp + exp(e)) - y(i) * e)
            end if
         end select
      end do
      if (sw > tiny(1.0_dp)) then
         value = value / sw
      else
         value = huge(1.0_dp)
      end if
   end function validation_loss

end module flexmix_penalized
