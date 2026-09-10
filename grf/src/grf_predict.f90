module grf_predict
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_nan
   use grf_kinds, only : dp
   use grf_types, only : grf_forest, grf_boosted_forest
   use grf_core, only : compute_forest_weights, find_leaf_node
   use grf_utils, only : weighted_quantile, weighted_covariance, weighted_variance
   use grf_utils, only : solve_linear_system, objective_bayes_debias
   implicit none
   private

   public :: predict_regression_forest
   public :: predict_multi_regression_forest
   public :: predict_probability_forest
   public :: predict_quantile_forest
   public :: predict_causal_forest
   public :: predict_instrumental_forest
   public :: predict_survival_forest
   public :: predict_causal_survival_forest
   public :: predict_lm_forest
   public :: predict_multi_arm_causal_forest
   public :: predict_ll_regression_forest
   public :: predict_boosted_regression_forest

contains

   pure subroutine predict_regression_forest(forest, x_new, prediction, oob, variance)
      type(grf_forest), intent(in) :: forest !! Trained scalar regression forest supplying trees and training responses.
      real(dp), intent(in) :: x_new(:,:) !! Predictor rows at which to estimate the conditional mean.
      real(dp), allocatable, intent(out) :: prediction(:) !! Conditional-mean estimate for each query row.
      logical, intent(in), optional :: oob !! If true and x_new represents training rows, exclude trees containing query row i in-bag.
      real(dp), allocatable, intent(out), optional :: variance(:) !! Objective-Bayes-debiased little-bag variance estimates; NaN when CI groups are unavailable.
      real(dp), allocatable :: leaf_outcome(:)
      real(dp), allocatable :: leaf_weight(:)
      logical, allocatable :: valid(:)
      real(dp), allocatable :: rho(:)
      real(dp) :: average_outcome
      real(dp) :: average_weight
      real(dp) :: nan_value
      logical :: use_oob
      integer :: i
      integer :: n_valid

      nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
      use_oob = .false.
      if (present(oob)) use_oob = oob
      allocate(prediction(size(x_new,1)))
      if (present(variance)) allocate(variance(size(x_new,1)))
      allocate(leaf_outcome(size(forest%trees)), leaf_weight(size(forest%trees)), valid(size(forest%trees)))
      allocate(rho(size(forest%trees)))
      do i = 1, size(x_new,1)
         call collect_regression_leaf_values(forest, x_new(i,:), i, use_oob, leaf_outcome, leaf_weight, valid)
         n_valid = count(valid)
         if (n_valid == 0) then
            prediction(i) = nan_value
            if (present(variance)) variance(i) = nan_value
            cycle
         end if
         average_outcome = sum(leaf_outcome, mask=valid) / real(n_valid, dp)
         average_weight = sum(leaf_weight, mask=valid) / real(n_valid, dp)
         if (abs(average_weight) <= 1.0e-16_dp) then
            prediction(i) = nan_value
            if (present(variance)) variance(i) = nan_value
            cycle
         end if
         prediction(i) = average_outcome / average_weight
         if (present(variance)) then
            rho = 0.0_dp
            where (valid)
               rho = (leaf_outcome - prediction(i) * leaf_weight) / average_weight
            end where
            call grouped_rho_variance(rho, valid, forest%options%ci_group_size, variance(i))
         end if
      end do
   end subroutine predict_regression_forest

   pure subroutine predict_multi_regression_forest(forest, x_new, prediction, oob)
      type(grf_forest), intent(in) :: forest !! Trained multi-response regression forest supplying trees and training outcomes.
      real(dp), intent(in) :: x_new(:,:) !! Predictor rows at which to estimate all conditional response means.
      real(dp), allocatable, intent(out) :: prediction(:,:) !! Matrix of predictions with rows matching x_new and columns matching outcomes.
      logical, intent(in), optional :: oob !! If true for training rows, use only out-of-bag trees for each observation.
      real(dp), allocatable :: kernel(:,:)

      call compute_forest_weights(forest, x_new, kernel, oob)
      allocate(prediction(size(x_new,1), forest%n_outputs))
      prediction = matmul(kernel, forest%y(:,1:forest%n_outputs))
   end subroutine predict_multi_regression_forest

   pure subroutine predict_probability_forest(forest, x_new, prediction, oob)
      type(grf_forest), intent(in) :: forest !! Trained probability forest with one-based training class labels.
      real(dp), intent(in) :: x_new(:,:) !! Predictor rows at which class probabilities are requested.
      real(dp), allocatable, intent(out) :: prediction(:,:) !! Class-probability matrix; columns correspond to labels 1 through n_classes.
      logical, intent(in), optional :: oob !! If true for training rows, use only out-of-bag trees for each observation.
      real(dp), allocatable :: kernel(:,:)
      integer :: i
      integer :: c

      call compute_forest_weights(forest, x_new, kernel, oob)
      allocate(prediction(size(x_new,1), forest%n_classes))
      prediction = 0.0_dp
      do i = 1, size(forest%classes)
         c = forest%classes(i)
         if (c >= 1 .and. c <= forest%n_classes) prediction(:,c) = prediction(:,c) + kernel(:,i)
      end do
   end subroutine predict_probability_forest

   pure subroutine predict_quantile_forest(forest, x_new, prediction, quantiles, oob)
      type(grf_forest), intent(in) :: forest !! Trained quantile forest retaining its training outcomes and default quantile grid.
      real(dp), intent(in) :: x_new(:,:) !! Predictor rows at which conditional quantiles are requested.
      real(dp), allocatable, intent(out) :: prediction(:,:) !! Conditional quantiles with one row per query and one column per requested probability.
      real(dp), intent(in), optional :: quantiles(:) !! Optional probabilities in (0,1); defaults to the grid used when training the forest.
      logical, intent(in), optional :: oob !! If true for training rows, use only out-of-bag trees for each observation.
      real(dp), allocatable :: kernel(:,:)
      real(dp), allocatable :: probs(:)
      integer :: i
      integer :: j

      if (present(quantiles)) then
         allocate(probs(size(quantiles)))
         probs = quantiles
      else
         allocate(probs(size(forest%quantiles)))
         probs = forest%quantiles
      end if
      call compute_forest_weights(forest, x_new, kernel, oob)
      allocate(prediction(size(x_new,1), size(probs)))
      do i = 1, size(x_new,1)
         do j = 1, size(probs)
            call weighted_quantile(forest%y(:,1), kernel(i,:), probs(j), prediction(i,j))
         end do
      end do
   end subroutine predict_quantile_forest

   pure subroutine predict_causal_forest(forest, x_new, prediction, oob, variance)
      type(grf_forest), intent(in) :: forest !! Trained scalar causal forest with outcome and treatment nuisance estimates.
      real(dp), intent(in) :: x_new(:,:) !! Predictor rows at which conditional treatment effects are requested.
      real(dp), allocatable, intent(out) :: prediction(:) !! Local residual-on-residual treatment-effect estimates for each query row.
      logical, intent(in), optional :: oob !! If true for training rows, use only trees for which query row i is out-of-bag.
      real(dp), allocatable, intent(out), optional :: variance(:) !! Objective-Bayes-debiased little-bag treatment-effect variance estimates.
      real(dp), allocatable :: yr(:)
      real(dp), allocatable :: wr(:)

      yr = forest%y(:,1) - forest%y_hat(:,1)
      wr = forest%w(:,1) - forest%w_hat(:,1)
      call predict_instrumental_moments(forest, x_new, yr, wr, wr, prediction, oob, variance)
   end subroutine predict_causal_forest

   pure subroutine predict_instrumental_forest(forest, x_new, prediction, oob, variance)
      type(grf_forest), intent(in) :: forest !! Trained instrumental forest retaining residualization values for Y, W, and Z.
      real(dp), intent(in) :: x_new(:,:) !! Predictor rows at which local IV treatment effects are requested.
      real(dp), allocatable, intent(out) :: prediction(:) !! Local IV treatment-effect estimate for each query row.
      logical, intent(in), optional :: oob !! If true for training rows, use only trees for which query row i is out-of-bag.
      real(dp), allocatable, intent(out), optional :: variance(:) !! Objective-Bayes-debiased little-bag IV variance estimates.
      real(dp), allocatable :: yr(:)
      real(dp), allocatable :: wr(:)
      real(dp), allocatable :: zr(:)

      yr = forest%y(:,1) - forest%y_hat(:,1)
      wr = forest%w(:,1) - forest%w_hat(:,1)
      zr = forest%z - forest%z_hat
      call predict_instrumental_moments(forest, x_new, yr, wr, zr, prediction, oob, variance)
   end subroutine predict_instrumental_forest

   pure subroutine predict_survival_forest(forest, x_new, prediction, cumulative_hazard, oob, failure_times, &
                                           nelson_aalen, individual_times)
      type(grf_forest), intent(in) :: forest !! Trained survival forest with event indicators and sorted failure-time grid.
      real(dp), intent(in) :: x_new(:,:) !! Predictor rows at which conditional survival or cumulative hazard is requested.
      real(dp), allocatable, intent(out) :: prediction(:,:) !! Query-by-time survival/hazard matrix, or query-by-one for individual_times mode.
      logical, intent(in), optional :: cumulative_hazard !! If true return Nelson-Aalen cumulative hazard instead of survival probability.
      logical, intent(in), optional :: oob !! If true for training rows, use only out-of-bag trees for localization.
      real(dp), intent(in), optional :: failure_times(:) !! Evaluation grid, or one time per query row when individual_times is true.
      logical, intent(in), optional :: nelson_aalen !! If true return exp(-Nelson-Aalen cumulative hazard); otherwise Kaplan-Meier survival.
      logical, intent(in), optional :: individual_times !! If true, failure_times must match query rows and one value is returned per row.
      real(dp), allocatable :: kernel(:,:)
      real(dp), allocatable :: times(:)
      real(dp) :: risk
      real(dp) :: deaths
      real(dp) :: hazard
      real(dp) :: survival
      real(dp) :: event_time
      real(dp) :: target_time
      logical :: want_hazard
      logical :: want_nelson
      logical :: per_row
      integer :: i
      integer :: j
      integer :: k
      integer :: cursor

      want_hazard = .false.
      if (present(cumulative_hazard)) want_hazard = cumulative_hazard
      want_nelson = .false.
      if (present(nelson_aalen)) want_nelson = nelson_aalen
      per_row = .false.
      if (present(individual_times)) per_row = individual_times
      if (per_row) then
         if (.not. present(failure_times)) then
            allocate(prediction(size(x_new,1),1))
            prediction = ieee_value(0.0_dp, ieee_quiet_nan)
            return
         end if
         if (size(failure_times) /= size(x_new,1)) then
            allocate(prediction(size(x_new,1),1))
            prediction = ieee_value(0.0_dp, ieee_quiet_nan)
            return
         end if
      else if (present(failure_times)) then
         allocate(times(size(failure_times)))
         times = failure_times
      else
         allocate(times(size(forest%failure_times)))
         times = forest%failure_times
      end if
      call compute_forest_weights(forest, x_new, kernel, oob)
      if (per_row) then
         allocate(prediction(size(x_new,1),1))
      else
         allocate(prediction(size(x_new,1),size(times)))
      end if
      do i = 1, size(x_new,1)
         if (per_row) then
            target_time = failure_times(i)
            hazard = 0.0_dp
            survival = 1.0_dp
            do cursor = 1, forest%n_times
               event_time = forest%failure_times(cursor)
               if (event_time > target_time) exit
               risk = 0.0_dp
               deaths = 0.0_dp
               do k = 1, forest%n_train
                  if (forest%y(k,1) >= event_time) risk = risk + kernel(i,k)
                  if (forest%event(k) == 1 .and. same_time(forest%y(k,1), event_time)) deaths = deaths + kernel(i,k)
               end do
               if (risk > 100.0_dp * epsilon(1.0_dp)) then
                  hazard = hazard + deaths / risk
                  survival = survival * max(0.0_dp, 1.0_dp - deaths / risk)
               end if
            end do
            if (want_hazard) then
               prediction(i,1) = hazard
            else if (want_nelson) then
               prediction(i,1) = exp(-hazard)
            else
               prediction(i,1) = survival
            end if
         else
            hazard = 0.0_dp
            survival = 1.0_dp
            cursor = 1
            do j = 1, size(times)
               do while (cursor <= forest%n_times)
                  event_time = forest%failure_times(cursor)
                  if (event_time > times(j)) exit
                  risk = 0.0_dp
                  deaths = 0.0_dp
                  do k = 1, forest%n_train
                     if (forest%y(k,1) >= event_time) risk = risk + kernel(i,k)
                     if (forest%event(k) == 1 .and. same_time(forest%y(k,1), event_time)) deaths = deaths + kernel(i,k)
                  end do
                  if (risk > 100.0_dp * epsilon(1.0_dp)) then
                     hazard = hazard + deaths / risk
                     survival = survival * max(0.0_dp, 1.0_dp - deaths / risk)
                  end if
                  cursor = cursor + 1
               end do
               if (want_hazard) then
                  prediction(i,j) = hazard
               else if (want_nelson) then
                  prediction(i,j) = exp(-hazard)
               else
                  prediction(i,j) = survival
               end if
            end do
         end if
      end do
   end subroutine predict_survival_forest

   pure subroutine predict_causal_survival_forest(forest, x_new, prediction, oob, variance)
      type(grf_forest), intent(in) :: forest !! Trained causal-survival forest containing doubly robust numerator and denominator scores.
      real(dp), intent(in) :: x_new(:,:) !! Predictor rows at which horizon-specific treatment effects are requested.
      real(dp), allocatable, intent(out) :: prediction(:) !! Local ratio of causal-survival numerator to denominator moments.
      logical, intent(in), optional :: oob !! If true for training rows, use only trees for which query row i is out-of-bag.
      real(dp), allocatable, intent(out), optional :: variance(:) !! Objective-Bayes-debiased little-bag variance for the horizon-specific effect.

      call predict_ratio_moments(forest, x_new, forest%causal_survival_numerator, &
         forest%causal_survival_denominator, prediction, oob, variance)
   end subroutine predict_causal_survival_forest

   pure subroutine predict_lm_forest(forest, x_new, prediction, ridge, oob)
      type(grf_forest), intent(in) :: forest !! Trained local-coefficient forest with residualized outcome and regressor matrices.
      real(dp), intent(in) :: x_new(:,:) !! Predictor rows at which local coefficient matrices are requested.
      real(dp), allocatable, intent(out) :: prediction(:,:,:) !! Coefficients shaped query by treatment/regressor by outcome.
      real(dp), intent(in), optional :: ridge !! Nonnegative diagonal regularization added to local treatment moment matrices; default 1e-8.
      logical, intent(in), optional :: oob !! If true for training rows, use only out-of-bag trees for localization.
      real(dp), allocatable :: kernel(:,:)
      real(dp), allocatable :: wr(:,:)
      real(dp), allocatable :: yr(:,:)
      real(dp), allocatable :: gram(:,:)
      real(dp), allocatable :: rhs(:)
      real(dp), allocatable :: beta(:)
      real(dp) :: lambda
      integer :: i
      integer :: j
      integer :: k
      integer :: m
      integer :: info

      lambda = 1.0e-8_dp
      if (present(ridge)) lambda = max(0.0_dp, ridge)
      call compute_forest_weights(forest, x_new, kernel, oob)
      wr = forest%w(:,1:forest%n_treatments) - forest%w_hat(:,1:forest%n_treatments)
      yr = forest%y(:,1:forest%n_outputs) - forest%y_hat(:,1:forest%n_outputs)
      allocate(prediction(size(x_new,1), forest%n_treatments, forest%n_outputs))
      allocate(gram(forest%n_treatments,forest%n_treatments), rhs(forest%n_treatments), beta(forest%n_treatments))
      do i = 1, size(x_new,1)
         gram = 0.0_dp
         do j = 1, forest%n_treatments
            do k = 1, forest%n_treatments
               gram(j,k) = sum(kernel(i,:) * wr(:,j) * wr(:,k))
            end do
            gram(j,j) = gram(j,j) + lambda
         end do
         do m = 1, forest%n_outputs
            do j = 1, forest%n_treatments
               rhs(j) = sum(kernel(i,:) * wr(:,j) * yr(:,m))
            end do
            call solve_linear_system(gram, rhs, beta, info)
            if (info == 0) then
               prediction(i,:,m) = beta
            else
               prediction(i,:,m) = ieee_value(0.0_dp, ieee_quiet_nan)
            end if
         end do
      end do
   end subroutine predict_lm_forest

   pure subroutine predict_multi_arm_causal_forest(forest, x_new, prediction, ridge, oob)
      type(grf_forest), intent(in) :: forest !! Trained multi-treatment forest represented by the shared GRF local-linear family.
      real(dp), intent(in) :: x_new(:,:) !! Predictor rows at which treatment contrasts relative to baseline are requested.
      real(dp), allocatable, intent(out) :: prediction(:,:) !! Treatment-effect matrix with query rows and one column per treatment contrast.
      real(dp), intent(in), optional :: ridge !! Nonnegative regularization for local treatment moment inversion; default 1e-8.
      logical, intent(in), optional :: oob !! If true for training rows, use only out-of-bag trees for localization.
      real(dp), allocatable :: all_prediction(:,:,:)

      call predict_lm_forest(forest, x_new, all_prediction, ridge, oob)
      allocate(prediction(size(x_new,1), forest%n_treatments))
      prediction = all_prediction(:,:,1)
   end subroutine predict_multi_arm_causal_forest

   pure subroutine predict_ll_regression_forest(forest, x_new, prediction, linear_variables, ridge, oob, variance, &
                                                 weight_penalty, lambda_path, selected_ridge)
      type(grf_forest), intent(in) :: forest !! Regression forest whose leaf co-membership weights define local neighborhoods.
      real(dp), intent(in) :: x_new(:,:) !! Predictor rows at which local-linear corrected regression predictions are requested.
      real(dp), allocatable, intent(out) :: prediction(:) !! Local ridge-regression intercept at each query row.
      integer, intent(in), optional :: linear_variables(:) !! Optional one-based predictor columns entering the local linear correction; defaults to all columns.
      real(dp), intent(in), optional :: ridge !! Nonnegative ridge multiplier; when absent it is selected by OOB cross-validation on the upstream default path.
      logical, intent(in), optional :: oob !! If true for training rows, use only out-of-bag trees for localization.
      real(dp), allocatable, intent(out), optional :: variance(:) !! Objective-Bayes-debiased little-bag variance of the local intercept estimate.
      logical, intent(in), optional :: weight_penalty !! If true, scale each slope penalty by its local second moment; default false uses trace normalization.
      real(dp), intent(in), optional :: lambda_path(:) !! Optional nonnegative ridge multipliers searched when ridge is absent; defaults to the upstream ten-value path.
      real(dp), intent(out), optional :: selected_ridge !! Effective ridge multiplier used for prediction, including an automatically selected value when ridge is absent.
      real(dp), allocatable :: kernel(:,:)
      real(dp), allocatable :: design(:,:)
      real(dp), allocatable :: gram(:,:)
      real(dp), allocatable :: rhs(:)
      real(dp), allocatable :: beta(:)
      real(dp), allocatable :: zeta(:)
      real(dp), allocatable :: pseudo_residual(:)
      integer, allocatable :: variables(:)
      real(dp) :: lambda
      real(dp) :: nan_value
      logical :: use_oob
      logical :: use_weight_penalty
      integer :: i
      integer :: j
      integer :: k
      integer :: info

      if (present(linear_variables)) then
         allocate(variables(size(linear_variables)))
         variables = linear_variables
      else
         allocate(variables(forest%p))
         variables = [(j, j = 1, forest%p)]
      end if
      use_weight_penalty = .false.
      if (present(weight_penalty)) use_weight_penalty = weight_penalty
      if (present(ridge)) then
         lambda = max(0.0_dp, ridge)
      else
         if (present(lambda_path)) then
            call tune_ll_ridge(forest, variables, use_weight_penalty, lambda, lambda_path)
         else
            call tune_ll_ridge(forest, variables, use_weight_penalty, lambda)
         end if
      end if
      if (present(selected_ridge)) selected_ridge = lambda
      use_oob = .false.
      if (present(oob)) use_oob = oob
      nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
      call compute_forest_weights(forest, x_new, kernel, oob)
      allocate(prediction(size(x_new,1)))
      if (present(variance)) allocate(variance(size(x_new,1)))
      allocate(design(forest%n_train,size(variables)+1))
      allocate(gram(size(variables)+1,size(variables)+1), rhs(size(variables)+1), beta(size(variables)+1))
      allocate(zeta(size(variables)+1), pseudo_residual(forest%n_train))
      do i = 1, size(x_new,1)
         call build_local_linear_system(forest, x_new(i,:), kernel(i,:), variables, design, gram, rhs)
         call apply_ll_penalty(gram, lambda, use_weight_penalty)
         call solve_linear_system(gram, rhs, beta, info)
         if (info == 0) then
            prediction(i) = beta(1)
         else
            prediction(i) = nan_value
            if (present(variance)) variance(i) = nan_value
            cycle
         end if
         if (present(variance)) then
            rhs = 0.0_dp
            rhs(1) = 1.0_dp
            call solve_linear_system(gram, rhs, zeta, info)
            if (info /= 0) then
               variance(i) = nan_value
               cycle
            end if
            do k = 1, forest%n_train
               pseudo_residual(k) = dot_product(design(k,:), zeta) * &
                  (forest%y(k,1) - dot_product(design(k,:), beta))
            end do
            call grouped_leaf_score_variance(forest, x_new(i,:), i, use_oob, pseudo_residual, variance(i))
         end if
      end do
   end subroutine predict_ll_regression_forest

   pure subroutine tune_ll_ridge(forest, variables, weight_penalty, ridge, lambda_path)
      type(grf_forest), intent(in) :: forest !! Regression forest whose OOB kernel weights are used to cross-validate local-linear ridge multipliers.
      integer, intent(in) :: variables(:) !! One-based predictor columns entering each local linear correction.
      logical, intent(in) :: weight_penalty !! True for covariance-scaled slope penalties; false for the upstream trace-normalized standard ridge penalty.
      real(dp), intent(out) :: ridge !! Ridge multiplier attaining the smallest finite OOB mean squared error on the candidate path.
      real(dp), intent(in), optional :: lambda_path(:) !! Optional nonnegative candidate multipliers; defaults to the upstream GRF path.
      real(dp), parameter :: default_path(10) = [0.0_dp, 0.001_dp, 0.01_dp, 0.05_dp, 0.1_dp, &
                                                0.3_dp, 0.5_dp, 0.7_dp, 1.0_dp, 10.0_dp]
      real(dp), allocatable :: path(:)
      real(dp), allocatable :: kernel(:,:)
      real(dp), allocatable :: design(:,:)
      real(dp), allocatable :: gram(:,:)
      real(dp), allocatable :: penalized(:,:)
      real(dp), allocatable :: rhs(:)
      real(dp), allocatable :: beta(:)
      real(dp), allocatable :: error_sum(:)
      integer, allocatable :: error_count(:)
      real(dp) :: prediction
      real(dp) :: error
      real(dp) :: best_error
      integer :: i
      integer :: j
      integer :: info
      integer :: best

      if (present(lambda_path)) then
         if (size(lambda_path) > 0) then
            allocate(path(size(lambda_path)))
            path = max(lambda_path, 0.0_dp)
         else
            allocate(path(size(default_path)))
            path = default_path
         end if
      else
         allocate(path(size(default_path)))
         path = default_path
      end if
      call compute_forest_weights(forest, forest%x, kernel, oob=.true.)
      allocate(design(forest%n_train,size(variables)+1))
      allocate(gram(size(variables)+1,size(variables)+1), penalized(size(variables)+1,size(variables)+1))
      allocate(rhs(size(variables)+1), beta(size(variables)+1))
      allocate(error_sum(size(path)), error_count(size(path)))
      error_sum = 0.0_dp
      error_count = 0
      do i = 1, forest%n_train
         if (sum(kernel(i,:)) <= 0.0_dp) cycle
         call build_local_linear_system(forest, forest%x(i,:), kernel(i,:), variables, design, gram, rhs)
         do j = 1, size(path)
            penalized = gram
            call apply_ll_penalty(penalized, path(j), weight_penalty)
            call solve_linear_system(penalized, rhs, beta, info)
            if (info /= 0) cycle
            prediction = beta(1)
            if (ieee_is_nan(prediction) .or. ieee_is_nan(forest%y(i,1))) cycle
            error = prediction - forest%y(i,1)
            error_sum(j) = error_sum(j) + error * error
            error_count(j) = error_count(j) + 1
         end do
      end do
      ridge = path(1)
      best = 0
      best_error = huge(1.0_dp)
      do j = 1, size(path)
         if (error_count(j) <= 0) cycle
         error = error_sum(j) / real(error_count(j), dp)
         if (error < best_error) then
            best_error = error
            best = j
         end if
      end do
      if (best > 0) ridge = path(best)
   end subroutine tune_ll_ridge

   pure subroutine build_local_linear_system(forest, query, weights, variables, design, gram, rhs)
      type(grf_forest), intent(in) :: forest !! Regression forest supplying training predictors and outcomes for the local weighted regression.
      real(dp), intent(in) :: query(:) !! Predictor row defining the centering point of the local linear design matrix.
      real(dp), intent(in) :: weights(:) !! Normalized forest kernel weights over training rows for this query.
      integer, intent(in) :: variables(:) !! One-based predictor columns included as centered local slopes.
      real(dp), intent(out) :: design(:,:) !! Local design matrix with an intercept followed by centered selected predictors.
      real(dp), intent(out) :: gram(:,:) !! Unpenalized weighted cross-product matrix X'WX.
      real(dp), intent(out) :: rhs(:) !! Weighted outcome cross-product vector X'WY.
      integer :: j
      integer :: k
      integer :: m

      design(:,1) = 1.0_dp
      do j = 1, size(variables)
         do k = 1, forest%n_train
            if (ieee_is_nan(forest%x(k,variables(j))) .or. ieee_is_nan(query(variables(j)))) then
               design(k,j+1) = 0.0_dp
            else
               design(k,j+1) = forest%x(k,variables(j)) - query(variables(j))
            end if
         end do
      end do
      gram = 0.0_dp
      rhs = 0.0_dp
      do j = 1, size(variables) + 1
         rhs(j) = sum(weights * design(:,j) * forest%y(:,1))
         do m = 1, size(variables) + 1
            gram(j,m) = sum(weights * design(:,j) * design(:,m))
         end do
      end do
   end subroutine build_local_linear_system

   pure subroutine apply_ll_penalty(gram, ridge, weight_penalty)
      real(dp), intent(inout) :: gram(:,:) !! Local weighted cross-product matrix modified in place by the requested slope penalty.
      real(dp), intent(in) :: ridge !! Nonnegative local-linear ridge multiplier.
      logical, intent(in) :: weight_penalty !! True to scale each slope by its own local second moment; false to use average diagonal normalization.
      real(dp) :: normalization
      real(dp) :: diagonal
      integer :: j

      if (ridge <= 0.0_dp .or. size(gram,1) <= 1) return
      if (weight_penalty) then
         do j = 2, size(gram,1)
            diagonal = gram(j,j)
            gram(j,j) = gram(j,j) + ridge * diagonal
         end do
      else
         normalization = 0.0_dp
         do j = 1, size(gram,1)
            normalization = normalization + gram(j,j)
         end do
         normalization = normalization / real(size(gram,1), dp)
         do j = 2, size(gram,1)
            gram(j,j) = gram(j,j) + ridge * normalization
         end do
      end if
   end subroutine apply_ll_penalty

   pure subroutine predict_boosted_regression_forest(forest, x_new, prediction)
      type(grf_boosted_forest), intent(in) :: forest !! Additive boosted forest containing residual stages, intercept, and learning rate.
      real(dp), intent(in) :: x_new(:,:) !! Predictor rows at which the additive boosted prediction is requested.
      real(dp), allocatable, intent(out) :: prediction(:) !! Intercept plus the learning-rate-scaled sum of stage predictions.
      real(dp), allocatable :: stage_prediction(:)
      integer :: s

      allocate(prediction(size(x_new,1)))
      prediction = forest%intercept
      do s = 1, forest%n_stages
         call predict_regression_forest(forest%stages(s), x_new, stage_prediction)
         prediction = prediction + forest%learning_rate * stage_prediction
      end do
   end subroutine predict_boosted_regression_forest

   pure subroutine predict_ratio_moments(forest, x_new, numerator, denominator, prediction, oob, variance)
      type(grf_forest), intent(in) :: forest !! Trained forest whose honest leaves contain ratio-estimator sufficient statistics.
      real(dp), intent(in) :: x_new(:,:) !! Predictor rows at which local numerator/denominator ratios are requested.
      real(dp), intent(in) :: numerator(:) !! Per-training-row numerator pseudo-outcomes.
      real(dp), intent(in) :: denominator(:) !! Per-training-row denominator pseudo-outcomes.
      real(dp), allocatable, intent(out) :: prediction(:) !! Local numerator divided by local denominator for each query row.
      logical, intent(in), optional :: oob !! If true, exclude trees containing query row i in-bag for training-row prediction.
      real(dp), allocatable, intent(out), optional :: variance(:) !! Objective-Bayes-debiased little-bag delta-method variance estimates.
      real(dp), allocatable :: leaf_numerator(:)
      real(dp), allocatable :: leaf_denominator(:)
      real(dp), allocatable :: rho(:)
      logical, allocatable :: valid(:)
      real(dp) :: average_numerator
      real(dp) :: average_denominator
      real(dp) :: nan_value
      logical :: use_oob
      integer :: i
      integer :: n_valid

      nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
      use_oob = .false.
      if (present(oob)) use_oob = oob
      allocate(prediction(size(x_new,1)))
      if (present(variance)) allocate(variance(size(x_new,1)))
      allocate(leaf_numerator(size(forest%trees)), leaf_denominator(size(forest%trees)))
      allocate(rho(size(forest%trees)), valid(size(forest%trees)))
      do i = 1, size(x_new,1)
         call collect_ratio_leaf_values(forest, x_new(i,:), i, use_oob, numerator, denominator, &
            leaf_numerator, leaf_denominator, valid)
         n_valid = count(valid)
         if (n_valid == 0) then
            prediction(i) = nan_value
            if (present(variance)) variance(i) = nan_value
            cycle
         end if
         average_numerator = sum(leaf_numerator, mask=valid) / real(n_valid, dp)
         average_denominator = sum(leaf_denominator, mask=valid) / real(n_valid, dp)
         if (abs(average_denominator) <= 1.0e-16_dp) then
            prediction(i) = nan_value
            if (present(variance)) variance(i) = nan_value
            cycle
         end if
         prediction(i) = average_numerator / average_denominator
         if (present(variance)) then
            rho = 0.0_dp
            where (valid)
               rho = leaf_numerator - leaf_denominator * prediction(i)
            end where
            call grouped_rho_variance(rho, valid, forest%options%ci_group_size, variance(i))
            if (.not. ieee_is_nan(variance(i))) variance(i) = variance(i) / (average_denominator * average_denominator)
         end if
      end do
   end subroutine predict_ratio_moments

   pure subroutine collect_ratio_leaf_values(forest, row, query_index, use_oob, numerator, denominator, &
      leaf_numerator, leaf_denominator, valid)
      type(grf_forest), intent(in) :: forest !! Trained forest supplying honest leaf memberships and observation weights.
      real(dp), intent(in) :: row(:) !! Single predictor row routed through every eligible tree.
      integer, intent(in) :: query_index !! One-based query-row index used to determine out-of-bag tree eligibility.
      logical, intent(in) :: use_oob !! If true, trees containing query_index in-bag are excluded.
      real(dp), intent(in) :: numerator(:) !! Per-training-row numerator pseudo-outcomes.
      real(dp), intent(in) :: denominator(:) !! Per-training-row denominator pseudo-outcomes.
      real(dp), intent(out) :: leaf_numerator(:) !! Weighted numerator sum divided by leaf sample count for each tree.
      real(dp), intent(out) :: leaf_denominator(:) !! Weighted denominator sum divided by leaf sample count for each tree.
      logical, intent(out) :: valid(:) !! True for eligible nonempty tree leaves having positive total sample weight.
      integer :: t
      integer :: node
      integer :: start
      integer :: count_leaf
      integer :: j
      integer :: sample
      real(dp) :: sum_weight

      leaf_numerator = 0.0_dp
      leaf_denominator = 0.0_dp
      valid = .false.
      do t = 1, size(forest%trees)
         if (use_oob .and. query_index <= forest%n_train) then
            if (forest%trees(t)%inbag(query_index)) cycle
         end if
         node = find_leaf_node(forest%trees(t), row)
         if (node < 1 .or. node > forest%trees(t)%n_nodes) cycle
         count_leaf = forest%trees(t)%leaf_count(node)
         if (count_leaf <= 0) cycle
         start = forest%trees(t)%leaf_start(node)
         sum_weight = 0.0_dp
         do j = start, start + count_leaf - 1
            sample = forest%trees(t)%leaf_samples(j)
            leaf_numerator(t) = leaf_numerator(t) + forest%sample_weights(sample) * numerator(sample)
            leaf_denominator(t) = leaf_denominator(t) + forest%sample_weights(sample) * denominator(sample)
            sum_weight = sum_weight + forest%sample_weights(sample)
         end do
         if (abs(sum_weight) <= 1.0e-16_dp) cycle
         leaf_numerator(t) = leaf_numerator(t) / real(count_leaf, dp)
         leaf_denominator(t) = leaf_denominator(t) / real(count_leaf, dp)
         valid(t) = .true.
      end do
   end subroutine collect_ratio_leaf_values

   pure subroutine collect_regression_leaf_values(forest, row, query_index, use_oob, leaf_outcome, leaf_weight, valid)
      type(grf_forest), intent(in) :: forest !! Trained regression forest supplying honest leaf memberships and sample weights.
      real(dp), intent(in) :: row(:) !! Single predictor row routed through every eligible tree.
      integer, intent(in) :: query_index !! One-based query-row index used to identify the corresponding training row for OOB prediction.
      logical, intent(in) :: use_oob !! If true, trees containing query_index in-bag are excluded when that index is a training row.
      real(dp), intent(out) :: leaf_outcome(:) !! Per-tree weighted outcome sum divided by leaf sample count.
      real(dp), intent(out) :: leaf_weight(:) !! Per-tree sample-weight sum divided by leaf sample count.
      logical, intent(out) :: valid(:) !! True for eligible trees with nonempty leaves and positive total sample weight.
      integer :: t
      integer :: node
      integer :: start
      integer :: count_leaf
      integer :: j
      integer :: sample
      real(dp) :: sum_outcome
      real(dp) :: sum_weight

      leaf_outcome = 0.0_dp
      leaf_weight = 0.0_dp
      valid = .false.
      do t = 1, size(forest%trees)
         if (use_oob .and. query_index <= forest%n_train) then
            if (forest%trees(t)%inbag(query_index)) cycle
         end if
         node = find_leaf_node(forest%trees(t), row)
         if (node < 1 .or. node > forest%trees(t)%n_nodes) cycle
         count_leaf = forest%trees(t)%leaf_count(node)
         if (count_leaf <= 0) cycle
         start = forest%trees(t)%leaf_start(node)
         sum_outcome = 0.0_dp
         sum_weight = 0.0_dp
         do j = start, start + count_leaf - 1
            sample = forest%trees(t)%leaf_samples(j)
            sum_outcome = sum_outcome + forest%sample_weights(sample) * forest%y(sample,1)
            sum_weight = sum_weight + forest%sample_weights(sample)
         end do
         if (abs(sum_weight) <= 1.0e-16_dp) cycle
         leaf_outcome(t) = sum_outcome / real(count_leaf, dp)
         leaf_weight(t) = sum_weight / real(count_leaf, dp)
         valid(t) = .true.
      end do
   end subroutine collect_regression_leaf_values

   pure subroutine predict_instrumental_moments(forest, x_new, outcome, treatment, instrument, prediction, oob, variance)
      type(grf_forest), intent(in) :: forest !! Trained causal, IV, or causal-survival forest whose honest leaves define local moments.
      real(dp), intent(in) :: x_new(:,:) !! Predictor rows at which local moment estimates are requested.
      real(dp), intent(in) :: outcome(:) !! Residualized or transformed outcome used in the estimating equation.
      real(dp), intent(in) :: treatment(:) !! Residualized treatment aligned with the training rows.
      real(dp), intent(in) :: instrument(:) !! Residualized instrument; equal to treatment for causal forests.
      real(dp), allocatable, intent(out) :: prediction(:) !! Local treatment-effect estimate for each query row.
      logical, intent(in), optional :: oob !! If true, exclude trees containing query row i in-bag for training-row prediction.
      real(dp), allocatable, intent(out), optional :: variance(:) !! Objective-Bayes-debiased little-bag variance estimate for each effect.
      integer, parameter :: outcome_index = 1
      integer, parameter :: treatment_index = 2
      integer, parameter :: instrument_index = 3
      integer, parameter :: outcome_instrument_index = 4
      integer, parameter :: treatment_instrument_index = 5
      integer, parameter :: weight_index = 7
      real(dp), allocatable :: leaf_value(:,:)
      real(dp), allocatable :: average(:)
      real(dp), allocatable :: rho(:)
      logical, allocatable :: valid(:)
      real(dp) :: numerator
      real(dp) :: denominator
      real(dp) :: main_effect
      real(dp) :: psi_1
      real(dp) :: psi_2
      real(dp) :: nan_value
      logical :: use_oob
      integer :: i
      integer :: t
      integer :: n_valid

      nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
      use_oob = .false.
      if (present(oob)) use_oob = oob
      allocate(prediction(size(x_new,1)))
      if (present(variance)) allocate(variance(size(x_new,1)))
      allocate(leaf_value(7,size(forest%trees)), average(7), rho(size(forest%trees)), valid(size(forest%trees)))
      do i = 1, size(x_new,1)
         call collect_instrumental_leaf_values(forest, x_new(i,:), i, use_oob, outcome, treatment, instrument, leaf_value, valid)
         n_valid = count(valid)
         if (n_valid == 0) then
            prediction(i) = nan_value
            if (present(variance)) variance(i) = nan_value
            cycle
         end if
         do t = 1, 7
            average(t) = sum(leaf_value(t,:), mask=valid) / real(n_valid, dp)
         end do
         numerator = average(outcome_instrument_index) * average(weight_index) - &
            average(outcome_index) * average(instrument_index)
         denominator = average(treatment_instrument_index) * average(weight_index) - &
            average(treatment_index) * average(instrument_index)
         if (abs(denominator) <= 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(numerator))) then
            prediction(i) = nan_value
            if (present(variance)) variance(i) = nan_value
            cycle
         end if
         prediction(i) = numerator / denominator
         if (present(variance)) then
            if (abs(average(weight_index)) <= 1.0e-16_dp) then
               variance(i) = nan_value
               cycle
            end if
            main_effect = (average(outcome_index) - average(treatment_index) * prediction(i)) / average(weight_index)
            rho = 0.0_dp
            do t = 1, size(forest%trees)
               if (.not. valid(t)) cycle
               psi_1 = leaf_value(outcome_instrument_index,t) - &
                  leaf_value(treatment_instrument_index,t) * prediction(i) - &
                  leaf_value(instrument_index,t) * main_effect
               psi_2 = leaf_value(outcome_index,t) - leaf_value(treatment_index,t) * prediction(i) - &
                  leaf_value(weight_index,t) * main_effect
               rho(t) = (average(weight_index) * psi_1 - average(instrument_index) * psi_2) / denominator
            end do
            call grouped_rho_variance(rho, valid, forest%options%ci_group_size, variance(i))
         end if
      end do
   end subroutine predict_instrumental_moments

   pure subroutine collect_instrumental_leaf_values(forest, row, query_index, use_oob, outcome, treatment, &
      instrument, leaf_value, valid)
      type(grf_forest), intent(in) :: forest !! Trained forest supplying honest leaf memberships and observation weights.
      real(dp), intent(in) :: row(:) !! Single predictor row routed through every eligible tree.
      integer, intent(in) :: query_index !! One-based query-row index used to select out-of-bag trees for training-row prediction.
      logical, intent(in) :: use_oob !! If true, trees containing query_index in-bag are excluded.
      real(dp), intent(in) :: outcome(:) !! Residualized or transformed outcome values for training observations.
      real(dp), intent(in) :: treatment(:) !! Residualized treatment values for training observations.
      real(dp), intent(in) :: instrument(:) !! Residualized instrument values for training observations.
      real(dp), intent(out) :: leaf_value(:,:) !! Seven sufficient statistics by tree: Y, W, Z, YZ, WZ, ZZ, and sample weight.
      logical, intent(out) :: valid(:) !! True for eligible trees with a nonempty honest leaf and positive total sample weight.
      integer :: t
      integer :: node
      integer :: start
      integer :: count_leaf
      integer :: j
      integer :: sample
      real(dp) :: g
      real(dp) :: y
      real(dp) :: w
      real(dp) :: z
      real(dp) :: sum_weight

      leaf_value = 0.0_dp
      valid = .false.
      do t = 1, size(forest%trees)
         if (use_oob .and. query_index <= forest%n_train) then
            if (forest%trees(t)%inbag(query_index)) cycle
         end if
         node = find_leaf_node(forest%trees(t), row)
         if (node < 1 .or. node > forest%trees(t)%n_nodes) cycle
         count_leaf = forest%trees(t)%leaf_count(node)
         if (count_leaf <= 0) cycle
         start = forest%trees(t)%leaf_start(node)
         sum_weight = 0.0_dp
         do j = start, start + count_leaf - 1
            sample = forest%trees(t)%leaf_samples(j)
            g = forest%sample_weights(sample)
            y = outcome(sample)
            w = treatment(sample)
            z = instrument(sample)
            leaf_value(1,t) = leaf_value(1,t) + g * y
            leaf_value(2,t) = leaf_value(2,t) + g * w
            leaf_value(3,t) = leaf_value(3,t) + g * z
            leaf_value(4,t) = leaf_value(4,t) + g * y * z
            leaf_value(5,t) = leaf_value(5,t) + g * w * z
            leaf_value(6,t) = leaf_value(6,t) + g * z * z
            sum_weight = sum_weight + g
         end do
         if (abs(sum_weight) <= 1.0e-16_dp) cycle
         leaf_value(1:6,t) = leaf_value(1:6,t) / real(count_leaf, dp)
         leaf_value(7,t) = sum_weight / real(count_leaf, dp)
         valid(t) = .true.
      end do
   end subroutine collect_instrumental_leaf_values

   pure subroutine grouped_rho_variance(rho, valid, ci_group_size, variance)
      real(dp), intent(in) :: rho(:) !! Per-tree delta-method pseudo-outcomes ordered in the same CI groups as forest trees.
      logical, intent(in) :: valid(:) !! Tree eligibility flags; a group contributes only when all its trees are valid.
      integer, intent(in) :: ci_group_size !! Number of trees sharing each half-sample; must exceed one for variance estimation.
      real(dp), intent(out) :: variance !! Objective-Bayes-debiased between-group variance, or NaN when no complete group contributes.
      real(dp) :: rho_squared
      real(dp) :: rho_grouped_squared
      real(dp) :: group_rho
      real(dp) :: var_between
      real(dp) :: var_total
      real(dp) :: group_noise
      real(dp) :: num_good_groups
      integer :: group
      integer :: first
      integer :: last

      if (ci_group_size <= 1) then
         variance = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      num_good_groups = 0.0_dp
      rho_squared = 0.0_dp
      rho_grouped_squared = 0.0_dp
      do group = 1, size(rho) / ci_group_size
         first = (group - 1) * ci_group_size + 1
         last = first + ci_group_size - 1
         if (.not. all(valid(first:last))) cycle
         num_good_groups = num_good_groups + 1.0_dp
         group_rho = sum(rho(first:last)) / real(ci_group_size, dp)
         rho_squared = rho_squared + sum(rho(first:last) * rho(first:last))
         rho_grouped_squared = rho_grouped_squared + group_rho * group_rho
      end do
      if (num_good_groups <= 0.0_dp) then
         variance = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      var_between = rho_grouped_squared / num_good_groups
      var_total = rho_squared / (num_good_groups * real(ci_group_size, dp))
      group_noise = (var_total - var_between) / real(ci_group_size - 1, dp)
      variance = objective_bayes_debias(var_between, group_noise, num_good_groups)
   end subroutine grouped_rho_variance

   pure subroutine grouped_leaf_score_variance(forest, row, query_index, use_oob, pseudo_residual, variance)
      type(grf_forest), intent(in) :: forest !! Trained forest whose grouped honest leaves define the bootstrap-of-little-bags replicates.
      real(dp), intent(in) :: row(:) !! Single predictor row used to identify one honest leaf in each tree.
      integer, intent(in) :: query_index !! One-based query-row index used for OOB eligibility when predicting training rows.
      logical, intent(in) :: use_oob !! If true, exclude trees containing query_index in-bag.
      real(dp), intent(in) :: pseudo_residual(:) !! Per-training-observation influence residual from the local estimating equation.
      real(dp), intent(out) :: variance !! Objective-Bayes-debiased variance of grouped leaf-average influence scores.
      real(dp), allocatable :: psi(:)
      logical, allocatable :: valid(:)
      real(dp) :: avg_score
      real(dp) :: psi_squared
      real(dp) :: psi_grouped_squared
      real(dp) :: group_psi
      real(dp) :: var_between
      real(dp) :: var_total
      real(dp) :: group_noise
      real(dp) :: num_good_groups
      integer :: t
      integer :: node
      integer :: start
      integer :: count_leaf
      integer :: j
      integer :: group
      integer :: first
      integer :: last

      if (forest%options%ci_group_size <= 1) then
         variance = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      allocate(psi(size(forest%trees)), valid(size(forest%trees)))
      psi = 0.0_dp
      valid = .false.
      do t = 1, size(forest%trees)
         if (use_oob .and. query_index <= forest%n_train) then
            if (forest%trees(t)%inbag(query_index)) cycle
         end if
         node = find_leaf_node(forest%trees(t), row)
         if (node < 1 .or. node > forest%trees(t)%n_nodes) cycle
         count_leaf = forest%trees(t)%leaf_count(node)
         if (count_leaf <= 0) cycle
         start = forest%trees(t)%leaf_start(node)
         do j = start, start + count_leaf - 1
            psi(t) = psi(t) + pseudo_residual(forest%trees(t)%leaf_samples(j))
         end do
         psi(t) = psi(t) / real(count_leaf, dp)
         valid(t) = .true.
      end do
      num_good_groups = 0.0_dp
      avg_score = 0.0_dp
      psi_squared = 0.0_dp
      psi_grouped_squared = 0.0_dp
      do group = 1, size(psi) / forest%options%ci_group_size
         first = (group - 1) * forest%options%ci_group_size + 1
         last = first + forest%options%ci_group_size - 1
         if (.not. all(valid(first:last))) cycle
         num_good_groups = num_good_groups + 1.0_dp
         group_psi = sum(psi(first:last)) / real(forest%options%ci_group_size, dp)
         psi_squared = psi_squared + sum(psi(first:last) * psi(first:last))
         psi_grouped_squared = psi_grouped_squared + group_psi * group_psi
         avg_score = avg_score + group_psi
      end do
      if (num_good_groups <= 0.0_dp) then
         variance = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      avg_score = avg_score / num_good_groups
      var_between = psi_grouped_squared / num_good_groups - avg_score * avg_score
      var_total = psi_squared / (num_good_groups * real(forest%options%ci_group_size, dp)) - avg_score * avg_score
      group_noise = (var_total - var_between) / real(forest%options%ci_group_size - 1, dp)
      variance = objective_bayes_debias(var_between, group_noise, num_good_groups)
   end subroutine grouped_leaf_score_variance

   pure elemental logical function same_time(a, b) result(equal)
      real(dp), intent(in) :: a !! First nonnegative survival time being compared.
      real(dp), intent(in) :: b !! Second nonnegative survival time being compared.

      equal = abs(a - b) <= 10.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(a), abs(b))
   end function same_time

end module grf_predict
