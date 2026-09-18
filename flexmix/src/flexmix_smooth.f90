! SPDX-License-Identifier: GPL-2.0-or-later
module flexmix_smooth
   use flexmix_kinds, only : dp
   use flexmix_types, only : flexmix_model_gaussian, flexmix_model_poisson, flexmix_model_binomial
   use flexmix_numeric, only : solve_linear_system
   implicit none
   private
   public :: fit_smooth_component

contains

   subroutine fit_smooth_component(x, y, trials, weights, model_kind, penalty_matrix, beta, sigma, effective_df, &
                                   selected_lambda, info, lambda_grid, offset)
      real(dp), intent(in) :: x(:,:) !! Prefit GAM/design matrix `(n,p)` containing parametric and smooth-basis columns.
      real(dp), intent(in) :: y(:) !! Gaussian/Poisson response, or binomial successes for a binomial family.
      real(dp), intent(in), optional :: trials(:) !! Binomial trial totals, size `n`; required for binomial models.
      real(dp), intent(in) :: weights(:) !! Nonnegative fitting weights, size `n`, including current EM responsibilities.
      integer, intent(in) :: model_kind !! Family identifier: Gaussian, Poisson, or binomial.
      real(dp), intent(in) :: penalty_matrix(:,:) !! Symmetric positive-semidefinite smoothing penalty matrix `(p,p)`.
      real(dp), intent(out) :: beta(:) !! Penalized fitted coefficients, size `p`.
      real(dp), intent(out) :: sigma !! Gaussian residual scale; set to one for Poisson/binomial families.
      real(dp), intent(out) :: effective_df !! Component effective degrees of freedom, including Gaussian scale when applicable.
      real(dp), intent(out) :: selected_lambda !! Smoothing multiplier selected from the supplied/generated grid.
      integer, intent(out) :: info !! Zero on success; negative for invalid inputs and positive for numerical failure.
      real(dp), intent(in), optional :: lambda_grid(:) !! Optional nonnegative smoothing grid; logarithmic defaults if absent.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive linear-predictor offset, size `n`.
      real(dp), allocatable :: lambdas(:), candidate(:), best_beta(:), work_trials(:)
      real(dp) :: criterion, best_criterion, edf_candidate, sigma_candidate
      integer :: fit_info, j, n, p

      n = size(y)
      p = size(x,2)
      beta = 0.0_dp
      sigma = 1.0_dp
      effective_df = 0.0_dp
      selected_lambda = 0.0_dp
      info = 0
      if (size(x,1) /= n .or. size(weights) /= n .or. size(beta) /= p .or. &
          any(shape(penalty_matrix) /= [p,p]) .or. n < 2 .or. p < 1) then
         info = -1
         return
      end if
      if (any(weights < 0.0_dp) .or. sum(weights) <= tiny(1.0_dp)) then
         info = -2
         return
      end if
      if (maxval(abs(penalty_matrix - transpose(penalty_matrix))) > 1.0e-9_dp * &
          max(1.0_dp,maxval(abs(penalty_matrix)))) then
         info = -3
         return
      end if
      if (model_kind /= flexmix_model_gaussian .and. model_kind /= flexmix_model_poisson .and. &
          model_kind /= flexmix_model_binomial) then
         info = -4
         return
      end if
      allocate(work_trials(n))
      work_trials = 1.0_dp
      if (model_kind == flexmix_model_binomial) then
         if (.not. present(trials)) then
            info = -5
            return
         end if
         if (size(trials) /= n .or. any(trials < y) .or. any(y < 0.0_dp)) then
            info = -6
            return
         end if
         work_trials = trials
      end if
      if (present(offset)) then
         if (size(offset) /= n) then
            info = -7
            return
         end if
      end if
      if (present(lambda_grid)) then
         if (size(lambda_grid) < 1 .or. any(lambda_grid < 0.0_dp)) then
            info = -8
            return
         end if
         allocate(lambdas(size(lambda_grid)))
         lambdas = lambda_grid
      else
         allocate(lambdas(25))
         do j = 1, size(lambdas)
            lambdas(j) = exp(log(1.0e-6_dp) + real(j-1,dp) * &
                         (log(1.0e6_dp)-log(1.0e-6_dp)) / real(size(lambdas)-1,dp))
         end do
      end if
      allocate(candidate(p),best_beta(p))
      candidate = 0.0_dp
      best_beta = 0.0_dp
      best_criterion = huge(1.0_dp)
      do j = 1, size(lambdas)
         call fit_at_lambda(x, y, work_trials, weights, model_kind, penalty_matrix, lambdas(j), candidate, &
                            sigma_candidate, edf_candidate, criterion, fit_info, offset)
         if (fit_info /= 0) cycle
         if (criterion < best_criterion) then
            best_criterion = criterion
            best_beta = candidate
            sigma = sigma_candidate
            effective_df = edf_candidate
            selected_lambda = lambdas(j)
         end if
      end do
      if (best_criterion >= 0.5_dp * huge(1.0_dp)) then
         info = 1
         return
      end if
      beta = best_beta
   end subroutine fit_smooth_component

   subroutine fit_at_lambda(x, y, trials, weights, model_kind, penalty_matrix, lambda_value, beta, sigma, &
                            effective_df, criterion, info, offset)
      real(dp), intent(in) :: x(:,:) !! Prefit GAM/design matrix `(n,p)`.
      real(dp), intent(in) :: y(:) !! Family response or binomial successes, size `n`.
      real(dp), intent(in) :: trials(:) !! Binomial trials; ignored for Gaussian/Poisson families.
      real(dp), intent(in) :: weights(:) !! Nonnegative fitting weights, size `n`.
      integer, intent(in) :: model_kind !! Family identifier.
      real(dp), intent(in) :: penalty_matrix(:,:) !! Symmetric smoothing penalty matrix `(p,p)`.
      real(dp), intent(in) :: lambda_value !! Nonnegative smoothing multiplier.
      real(dp), intent(inout) :: beta(:) !! Warm-start coefficient vector on entry and fitted coefficients on return.
      real(dp), intent(out) :: sigma !! Gaussian residual scale or one for non-Gaussian families.
      real(dp), intent(out) :: effective_df !! Effective component parameter count, including Gaussian scale when applicable.
      real(dp), intent(out) :: criterion !! GCV-like Gaussian or AIC-like GLM smoothing-selection criterion.
      integer, intent(out) :: info !! Zero on success; positive when the penalized solve/IRLS fails.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive predictor offset, size `n`.
      real(dp), allocatable :: eta(:), mu(:), work_weight(:), working_response(:), adjusted_response(:), old_beta(:)
      real(dp), allocatable :: xtwx(:,:), system(:,:), rhs(:), residual(:)
      real(dp) :: sw, meanw, denom, deviance, delta, e, prob, var
      integer :: i, iter, p, solve_info

      p = size(x,2)
      allocate(eta(size(y)),mu(size(y)),work_weight(size(y)),working_response(size(y)),adjusted_response(size(y)))
      allocate(old_beta(p),xtwx(p,p),system(p,p),rhs(p),residual(size(y)))
      sigma = 1.0_dp
      info = 0
      sw = sum(weights)
      meanw = sw / real(size(y),dp)
      if (model_kind == flexmix_model_gaussian) then
         adjusted_response = y
         if (present(offset)) adjusted_response = adjusted_response - offset
         call weighted_crossproducts(x, adjusted_response, weights, xtwx, rhs)
         system = xtwx + lambda_value * penalty_matrix
         call solve_linear_system(system,rhs,beta,solve_info)
         if (solve_info /= 0) then
            info = 1
            return
         end if
         call effective_df_from_system(system,xtwx,effective_df,solve_info)
         if (solve_info /= 0) then
            info = 2
            return
         end if
         eta = matmul(x,beta)
         if (present(offset)) eta = eta + offset
         residual = y - eta
         denom = real(size(y),dp) - effective_df
         if (denom <= 1.0e-8_dp .or. meanw <= tiny(1.0_dp)) then
            info = 3
            return
         end if
         sigma = sqrt(max(tiny(1.0_dp),sum(weights * residual**2 / meanw) / denom))
         criterion = real(size(y),dp) * sum(weights*residual**2) / &
                     max(1.0e-12_dp,(sw - meanw*effective_df)**2)
         effective_df = effective_df + 1.0_dp
         return
      end if

      if (maxval(abs(beta)) <= tiny(1.0_dp)) call initialize_glm_beta(y,trials,weights,model_kind,beta,offset)
      do iter = 1, 100
         old_beta = beta
         eta = matmul(x,beta)
         if (present(offset)) eta = eta + offset
         select case (model_kind)
         case (flexmix_model_poisson)
            do i = 1, size(y)
               e = min(30.0_dp,max(-30.0_dp,eta(i)))
               mu(i) = exp(e)
               work_weight(i) = weights(i) * mu(i)
               working_response(i) = e + (y(i)-mu(i))/max(mu(i),1.0e-12_dp)
            end do
         case (flexmix_model_binomial)
            do i = 1, size(y)
               e = min(30.0_dp,max(-30.0_dp,eta(i)))
               if (e >= 0.0_dp) then
                  prob = 1.0_dp / (1.0_dp + exp(-e))
               else
                  prob = exp(e) / (1.0_dp + exp(e))
               end if
               mu(i) = trials(i) * prob
               var = max(1.0e-10_dp,trials(i)*prob*(1.0_dp-prob))
               work_weight(i) = weights(i) * var
               working_response(i) = e + (y(i)-mu(i))/var
            end do
         end select
         adjusted_response = working_response
         if (present(offset)) adjusted_response = adjusted_response - offset
         call weighted_crossproducts(x, adjusted_response, work_weight, xtwx, rhs)
         system = xtwx + lambda_value * penalty_matrix
         call solve_linear_system(system,rhs,beta,solve_info)
         if (solve_info /= 0) then
            info = 4
            return
         end if
         delta = maxval(abs(beta-old_beta))
         if (delta <= 1.0e-9_dp * (1.0_dp + maxval(abs(beta)))) exit
      end do
      if (iter > 100) then
         info = 5
         return
      end if
      eta = matmul(x,beta)
      if (present(offset)) eta = eta + offset
      select case (model_kind)
      case (flexmix_model_poisson)
         do i = 1, size(y)
            e = min(30.0_dp,max(-30.0_dp,eta(i)))
            mu(i) = exp(e)
            work_weight(i) = weights(i) * mu(i)
         end do
      case (flexmix_model_binomial)
         do i = 1, size(y)
            e = min(30.0_dp,max(-30.0_dp,eta(i)))
            if (e >= 0.0_dp) then
               prob = 1.0_dp / (1.0_dp + exp(-e))
            else
               prob = exp(e) / (1.0_dp + exp(e))
            end if
            mu(i) = trials(i) * prob
            work_weight(i) = weights(i) * max(1.0e-10_dp,trials(i)*prob*(1.0_dp-prob))
         end do
      end select
      call weighted_crossproducts(x, y*0.0_dp, work_weight, xtwx, rhs)
      system = xtwx + lambda_value * penalty_matrix
      call effective_df_from_system(system,xtwx,effective_df,solve_info)
      if (solve_info /= 0) then
         info = 6
         return
      end if
      deviance = family_deviance(y,trials,weights,eta,model_kind)
      criterion = deviance + 2.0_dp * effective_df
   end subroutine fit_at_lambda

   subroutine weighted_crossproducts(x, response, weights, xtwx, rhs)
      real(dp), intent(in) :: x(:,:) !! Design matrix `(n,p)`.
      real(dp), intent(in) :: response(:) !! Working response, size `n`.
      real(dp), intent(in) :: weights(:) !! Nonnegative working weights, size `n`.
      real(dp), intent(out) :: xtwx(:,:) !! Weighted normal-equation matrix `(p,p)`.
      real(dp), intent(out) :: rhs(:) !! Weighted design-response crossproduct, size `p`.
      integer :: i
      xtwx = 0.0_dp
      rhs = 0.0_dp
      do i = 1, size(response)
         xtwx = xtwx + weights(i) * outer_product(x(i,:),x(i,:))
         rhs = rhs + weights(i) * x(i,:) * response(i)
      end do
   end subroutine weighted_crossproducts

   subroutine effective_df_from_system(system, xtwx, edf, info)
      real(dp), intent(in) :: system(:,:) !! Penalized normal-equation matrix `(p,p)`.
      real(dp), intent(in) :: xtwx(:,:) !! Unpenalized weighted information matrix `(p,p)`.
      real(dp), intent(out) :: edf !! Trace of `(X'WX + lambda S)^(-1) X'WX`.
      integer, intent(out) :: info !! Zero on success; nonzero if a system column cannot be solved.
      real(dp), allocatable :: column(:), solution(:)
      integer :: j, p, solve_info
      p = size(system,1)
      allocate(column(p),solution(p))
      edf = 0.0_dp
      do j = 1, p
         column = xtwx(:,j)
         call solve_linear_system(system,column,solution,solve_info)
         if (solve_info /= 0) then
            info = solve_info
            return
         end if
         edf = edf + solution(j)
      end do
      info = 0
   end subroutine effective_df_from_system

   subroutine initialize_glm_beta(y, trials, weights, model_kind, beta, offset)
      real(dp), intent(in) :: y(:) !! Poisson response or binomial successes, size `n`.
      real(dp), intent(in) :: trials(:) !! Binomial trial totals; ignored for Poisson.
      real(dp), intent(in) :: weights(:) !! Nonnegative fitting weights, size `n`.
      integer, intent(in) :: model_kind !! Poisson or binomial family identifier.
      real(dp), intent(out) :: beta(:) !! Initial coefficient vector with an intercept-only starting value.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed predictor offset, size `n`.
      real(dp) :: sw, pbar
      beta = 0.0_dp
      sw = max(sum(weights),tiny(1.0_dp))
      select case (model_kind)
      case (flexmix_model_poisson)
         beta(1) = log(max(1.0e-10_dp,dot_product(weights,y)/sw))
         if (present(offset)) beta(1) = beta(1) - dot_product(weights,offset)/sw
      case (flexmix_model_binomial)
         pbar = dot_product(weights,y) / max(1.0e-10_dp,dot_product(weights,trials))
         pbar = min(1.0_dp-1.0e-8_dp,max(1.0e-8_dp,pbar))
         beta(1) = log(pbar/(1.0_dp-pbar))
         if (present(offset)) beta(1) = beta(1) - dot_product(weights,offset)/sw
      end select
   end subroutine initialize_glm_beta

   function family_deviance(y, trials, weights, eta, model_kind) result(value)
      real(dp), intent(in) :: y(:) !! Poisson response or binomial successes, size `n`.
      real(dp), intent(in) :: trials(:) !! Binomial trial totals; ignored for Poisson.
      real(dp), intent(in) :: weights(:) !! Nonnegative observation weights, size `n`.
      real(dp), intent(in) :: eta(:) !! Fitted linear predictor, size `n`.
      integer, intent(in) :: model_kind !! Poisson or binomial family identifier.
      real(dp) :: value
      real(dp) :: mu, p, e, term
      integer :: i
      value = 0.0_dp
      select case (model_kind)
      case (flexmix_model_poisson)
         do i = 1, size(y)
            e = min(30.0_dp,max(-30.0_dp,eta(i)))
            mu = exp(e)
            term = mu - y(i)
            if (y(i) > 0.0_dp) term = term + y(i)*log(y(i)/mu)
            value = value + 2.0_dp*weights(i)*term
         end do
      case (flexmix_model_binomial)
         do i = 1, size(y)
            e = min(30.0_dp,max(-30.0_dp,eta(i)))
            if (e >= 0.0_dp) then
               p = 1.0_dp/(1.0_dp+exp(-e))
            else
               p = exp(e)/(1.0_dp+exp(e))
            end if
            p = min(1.0_dp-1.0e-12_dp,max(1.0e-12_dp,p))
            term = 0.0_dp
            if (y(i) > 0.0_dp) term = term + y(i)*log(y(i)/(trials(i)*p))
            if (trials(i)-y(i) > 0.0_dp) then
               term = term + (trials(i)-y(i))*log((trials(i)-y(i))/(trials(i)*(1.0_dp-p)))
            end if
            value = value + 2.0_dp*weights(i)*term
         end do
      end select
   end function family_deviance

   pure function outer_product(a, b) result(c)
      real(dp), intent(in) :: a(:) !! Left vector in an outer product.
      real(dp), intent(in) :: b(:) !! Right vector in an outer product.
      real(dp) :: c(size(a),size(b))
      integer :: i
      do i = 1, size(a)
         c(i,:) = a(i)*b
      end do
   end function outer_product

end module flexmix_smooth
