module grf_train
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_nan
   use grf_kinds, only : dp
   use grf_rng, only : grf_rng_state
   use grf_types, only : grf_options, grf_forest, grf_boosted_forest
   use grf_types, only : grf_regression, grf_multi_regression, grf_probability, grf_quantile
   use grf_types, only : grf_causal, grf_instrumental, grf_survival, grf_causal_survival, grf_lm
   use grf_core, only : train_forest_core, compute_forest_weights, find_leaf_node
   use grf_utils, only : weighted_mean, unique_sorted, solve_linear_system
   implicit none
   private

   public :: regression_forest
   public :: multi_regression_forest
   public :: probability_forest
   public :: quantile_forest
   public :: causal_forest
   public :: instrumental_forest
   public :: survival_forest
   public :: causal_survival_forest
   public :: lm_forest
   public :: multi_arm_causal_forest
   public :: ll_regression_forest
   public :: boosted_regression_forest

contains

   subroutine regression_forest(x, y, forest, info, options, sample_weights)
      real(dp), intent(in) :: x(:,:) !! Predictor matrix with observations in rows and variables in columns; NaNs are routed by learned splits.
      real(dp), intent(in) :: y(:) !! Numeric response vector with one value for each row of x.
      type(grf_forest), intent(out) :: forest !! Trained honest regression forest, including copied training data and tree structure.
      integer, intent(out) :: info !! Zero on success; negative values indicate invalid dimensions or training options.
      type(grf_options), intent(in), optional :: options !! Optional forest controls; default values match the package-level Fortran defaults.
      real(dp), intent(in), optional :: sample_weights(:) !! Optional nonnegative observation weights; defaults to one for every training row.
      type(grf_options) :: local_options
      type(grf_options) :: tuned_options

      call initialize_forest(x, forest, info, options, sample_weights)
      if (info /= 0) return
      if (size(y) /= size(x,1)) then
         info = -101
         return
      end if
      allocate(forest%y(size(y),1))
      forest%y(:,1) = y
      forest%family = grf_regression
      forest%n_outputs = 1
      call choose_options(options, local_options)
      if (local_options%tune_parameters) then
         call tune_regression_options(x, y, forest%sample_weights, local_options, tuned_options, info)
         if (info /= 0) return
         local_options = tuned_options
         local_options%tune_parameters = .false.
      end if
      call train_forest_core(forest, local_options, info)
   end subroutine regression_forest

   subroutine multi_regression_forest(x, y, forest, info, options, sample_weights)
      real(dp), intent(in) :: x(:,:) !! Predictor matrix with observations in rows and variables in columns.
      real(dp), intent(in) :: y(:,:) !! Multi-response outcome matrix with observations in rows and responses in columns.
      type(grf_forest), intent(out) :: forest !! Trained multi-response regression forest and copied training data.
      integer, intent(out) :: info !! Zero on success; negative values indicate invalid dimensions or options.
      type(grf_options), intent(in), optional :: options !! Optional forest controls; defaults are used when absent.
      real(dp), intent(in), optional :: sample_weights(:) !! Optional nonnegative observation weights matching rows of x.
      type(grf_options) :: local_options

      call initialize_forest(x, forest, info, options, sample_weights)
      if (info /= 0) return
      if (size(y,1) /= size(x,1) .or. size(y,2) < 1) then
         info = -102
         return
      end if
      forest%y = y
      forest%family = grf_multi_regression
      forest%n_outputs = size(y,2)
      call choose_options(options, local_options)
      call train_forest_core(forest, local_options, info)
   end subroutine multi_regression_forest

   subroutine probability_forest(x, classes, forest, info, options, sample_weights)
      real(dp), intent(in) :: x(:,:) !! Predictor matrix with observations in rows and variables in columns.
      integer, intent(in) :: classes(:) !! One-based class labels; labels must lie in 1 through maxval(classes).
      type(grf_forest), intent(out) :: forest !! Trained probability forest and copied class labels.
      integer, intent(out) :: info !! Zero on success; negative values indicate invalid class labels, dimensions, or options.
      type(grf_options), intent(in), optional :: options !! Optional forest controls; defaults are used when absent.
      real(dp), intent(in), optional :: sample_weights(:) !! Optional nonnegative observation weights matching rows of x.
      type(grf_options) :: local_options

      call initialize_forest(x, forest, info, options, sample_weights)
      if (info /= 0) return
      if (size(classes) /= size(x,1) .or. size(classes) == 0) then
         info = -103
         return
      end if
      if (minval(classes) < 1) then
         info = -104
         return
      end if
      forest%classes = classes
      forest%n_classes = maxval(classes)
      forest%family = grf_probability
      forest%n_outputs = forest%n_classes
      call choose_options(options, local_options)
      call train_forest_core(forest, local_options, info)
   end subroutine probability_forest

   subroutine quantile_forest(x, y, quantiles, forest, info, options)
      real(dp), intent(in) :: x(:,:) !! Predictor matrix with observations in rows and variables in columns.
      real(dp), intent(in) :: y(:) !! Numeric outcome vector whose conditional quantiles are estimated.
      real(dp), intent(in) :: quantiles(:) !! Requested probabilities strictly between zero and one, such as [0.1,0.5,0.9].
      type(grf_forest), intent(out) :: forest !! Trained quantile forest with the requested probability grid retained.
      integer, intent(out) :: info !! Zero on success; negative values indicate invalid dimensions, quantiles, or options.
      type(grf_options), intent(in), optional :: options !! Optional forest controls; defaults are used when absent.
      type(grf_options) :: local_options

      call initialize_forest(x, forest, info, options)
      if (info /= 0) return
      if (size(y) /= size(x,1) .or. size(quantiles) < 1) then
         info = -105
         return
      end if
      if (any(quantiles <= 0.0_dp) .or. any(quantiles >= 1.0_dp)) then
         info = -106
         return
      end if
      allocate(forest%y(size(y),1))
      forest%y(:,1) = y
      forest%quantiles = quantiles
      forest%family = grf_quantile
      forest%n_outputs = 1
      call choose_options(options, local_options)
      call train_forest_core(forest, local_options, info)
   end subroutine quantile_forest

   subroutine causal_forest(x, y, w, forest, info, options, y_hat, w_hat, sample_weights)
      real(dp), intent(in) :: x(:,:) !! Predictor matrix used to discover treatment-effect heterogeneity.
      real(dp), intent(in) :: y(:) !! Scalar observed outcome for each training observation.
      real(dp), intent(in) :: w(:) !! Scalar treatment or exposure for each training observation.
      type(grf_forest), intent(out) :: forest !! Trained causal forest with residualization values retained for prediction.
      integer, intent(out) :: info !! Zero on success; negative values indicate invalid dimensions or options.
      type(grf_options), intent(in), optional :: options !! Optional forest controls; defaults are used when absent.
      real(dp), intent(in), optional :: y_hat(:) !! Optional nuisance estimate E[Y|X]; defaults to out-of-bag regression-forest predictions.
      real(dp), intent(in), optional :: w_hat(:) !! Optional nuisance estimate E[W|X]; defaults to out-of-bag regression-forest predictions.
      real(dp), intent(in), optional :: sample_weights(:) !! Optional nonnegative observation weights matching rows of x.
      type(grf_options) :: local_options

      call initialize_forest(x, forest, info, options, sample_weights)
      if (info /= 0) return
      if (size(y) /= size(x,1) .or. size(w) /= size(x,1)) then
         info = -107
         return
      end if
      allocate(forest%y(size(y),1), forest%w(size(w),1))
      forest%y(:,1) = y
      forest%w(:,1) = w
      allocate(forest%y_hat(size(y),1), forest%w_hat(size(w),1))
      call choose_options(options, local_options)
      call resolve_nuisance_vector(x, y, forest%sample_weights, y_hat, local_options, 1009, &
                                   forest%y_hat(:,1), info)
      if (info /= 0) return
      call resolve_nuisance_vector(x, w, forest%sample_weights, w_hat, local_options, 2017, &
                                   forest%w_hat(:,1), info)
      if (info /= 0) return
      allocate(forest%treatment_variance_hat(size(w)))
      call estimate_treatment_variance(x, w, forest%w_hat(:,1), forest%sample_weights, local_options, &
         2503, forest%treatment_variance_hat, info)
      if (info /= 0) return
      forest%family = grf_causal
      forest%n_outputs = 1
      forest%n_treatments = 1
      call train_forest_core(forest, local_options, info)
   end subroutine causal_forest

   subroutine instrumental_forest(x, y, w, z, forest, info, options, y_hat, w_hat, z_hat, &
                                  reduced_form_weight, sample_weights)
      real(dp), intent(in) :: x(:,:) !! Predictor matrix used to discover heterogeneous local IV effects.
      real(dp), intent(in) :: y(:) !! Scalar observed outcome for each training observation.
      real(dp), intent(in) :: w(:) !! Scalar treatment for each training observation.
      real(dp), intent(in) :: z(:) !! Scalar instrument for each training observation.
      type(grf_forest), intent(out) :: forest !! Trained instrumental forest with nuisance values retained.
      integer, intent(out) :: info !! Zero on success; negative values indicate invalid dimensions, nuisance vectors, or options.
      type(grf_options), intent(in), optional :: options !! Optional forest controls; defaults are used when absent.
      real(dp), intent(in), optional :: y_hat(:) !! Optional E[Y|X] nuisance estimates; defaults to out-of-bag regression-forest predictions.
      real(dp), intent(in), optional :: w_hat(:) !! Optional E[W|X] nuisance estimates; defaults to out-of-bag regression-forest predictions.
      real(dp), intent(in), optional :: z_hat(:) !! Optional E[Z|X] nuisance estimates; defaults to out-of-bag regression-forest predictions.
      real(dp), intent(in), optional :: reduced_form_weight !! Optional mixing weight in [0,1] for regularizing Z toward W; default zero.
      real(dp), intent(in), optional :: sample_weights(:) !! Optional nonnegative observation weights matching rows of x.
      type(grf_options) :: local_options
      type(grf_options) :: compliance_options
      type(grf_forest) :: compliance_forest
      real(dp), allocatable :: compliance_kernel(:,:)
      real(dp), allocatable :: compliance_y_residual(:)
      real(dp), allocatable :: compliance_z_residual(:)
      real(dp) :: local_y
      real(dp) :: local_z
      real(dp) :: local_yz
      real(dp) :: local_z2
      real(dp) :: local_denominator
      logical :: binary_z
      integer :: i

      call initialize_forest(x, forest, info, options, sample_weights)
      if (info /= 0) return
      if (size(y) /= size(x,1) .or. size(w) /= size(x,1) .or. size(z) /= size(x,1)) then
         info = -108
         return
      end if
      allocate(forest%y(size(y),1), forest%w(size(w),1), forest%z(size(z)))
      forest%y(:,1) = y
      forest%w(:,1) = w
      forest%z = z
      allocate(forest%y_hat(size(y),1), forest%w_hat(size(w),1), forest%z_hat(size(z)))
      call choose_options(options, local_options)
      call resolve_nuisance_vector(x, y, forest%sample_weights, y_hat, local_options, 1009, &
                                   forest%y_hat(:,1), info)
      if (info /= 0) return
      call resolve_nuisance_vector(x, w, forest%sample_weights, w_hat, local_options, 2017, &
                                   forest%w_hat(:,1), info)
      if (info /= 0) return
      call resolve_nuisance_vector(x, z, forest%sample_weights, z_hat, local_options, 3011, forest%z_hat, info)
      if (info /= 0) return
      binary_z = all((abs(z) <= 1.0e-12_dp) .or. (abs(z - 1.0_dp) <= 1.0e-12_dp))
      if (binary_z) then
         compliance_options = local_options
         compliance_options%num_trees = max(50, min(max(1, local_options%num_trees), 500))
         compliance_options%ci_group_size = 1
         compliance_options%seed = local_options%seed + 4013
         if (cluster_equalized(compliance_options)) then
            call causal_forest(x, w, z, compliance_forest, info, compliance_options, &
               y_hat=forest%w_hat(:,1), w_hat=forest%z_hat)
         else
            call causal_forest(x, w, z, compliance_forest, info, compliance_options, &
               y_hat=forest%w_hat(:,1), w_hat=forest%z_hat, sample_weights=forest%sample_weights)
         end if
         if (info /= 0) return
         call compute_forest_weights(compliance_forest, x, compliance_kernel, oob=.true.)
         compliance_y_residual = w - forest%w_hat(:,1)
         compliance_z_residual = z - forest%z_hat
         allocate(forest%compliance_hat(size(z)))
         do i = 1, size(z)
            local_y = sum(compliance_kernel(i,:) * compliance_y_residual)
            local_z = sum(compliance_kernel(i,:) * compliance_z_residual)
            local_yz = sum(compliance_kernel(i,:) * compliance_y_residual * compliance_z_residual)
            local_z2 = sum(compliance_kernel(i,:) * compliance_z_residual * compliance_z_residual)
            local_denominator = local_z2 - local_z * local_z
            if (abs(local_denominator) <= 100.0_dp * epsilon(1.0_dp)) then
               forest%compliance_hat(i) = 0.0_dp
            else
               forest%compliance_hat(i) = (local_yz - local_y * local_z) / local_denominator
            end if
         end do
      end if
      forest%reduced_form_weight = 0.0_dp
      if (present(reduced_form_weight)) forest%reduced_form_weight = reduced_form_weight
      if (forest%reduced_form_weight < 0.0_dp .or. forest%reduced_form_weight > 1.0_dp) then
         info = -109
         return
      end if
      forest%family = grf_instrumental
      forest%n_outputs = 1
      forest%n_treatments = 1
      call train_forest_core(forest, local_options, info)
   end subroutine instrumental_forest

   subroutine survival_forest(x, y, event, forest, info, options, sample_weights, failure_times)
      real(dp), intent(in) :: x(:,:) !! Predictor matrix for conditional survival estimation.
      real(dp), intent(in) :: y(:) !! Nonnegative observed event or censoring times.
      integer, intent(in) :: event(:) !! Event indicators equal to one for failures and zero for right-censoring.
      type(grf_forest), intent(out) :: forest !! Trained survival forest retaining the sorted failure-time grid and relabeled time indices.
      integer, intent(out) :: info !! Zero on success; negative values indicate invalid dimensions, event indicators, grid, or options.
      type(grf_options), intent(in), optional :: options !! Optional forest controls, including the fast_logrank approximation flag.
      real(dp), intent(in), optional :: sample_weights(:) !! Optional nonnegative observation weights matching rows of x.
      real(dp), intent(in), optional :: failure_times(:) !! Optional strictly increasing training event grid; defaults to all distinct observed failure times.
      type(grf_options) :: local_options
      real(dp), allocatable :: event_times(:)
      integer :: i

      call initialize_forest(x, forest, info, options, sample_weights)
      if (info /= 0) return
      if (size(y) /= size(x,1) .or. size(event) /= size(x,1)) then
         info = -110
         return
      end if
      if (any(event < 0) .or. any(event > 1) .or. any(y < 0.0_dp)) then
         info = -111
         return
      end if
      if (count(event == 1) == 0) then
         info = -112
         return
      end if
      allocate(forest%y(size(y),1))
      forest%y(:,1) = y
      forest%event = event
      if (present(failure_times)) then
         if (size(failure_times) < 1 .or. .not. strictly_increasing(failure_times)) then
            info = -120
            return
         end if
         forest%failure_times = failure_times
      else
         event_times = pack(y, event == 1)
         call unique_sorted(event_times, forest%failure_times)
      end if
      forest%n_times = size(forest%failure_times)
      allocate(forest%survival_time_index(size(y)))
      do i = 1, size(y)
         forest%survival_time_index(i) = failure_interval_index(y(i), forest%failure_times)
      end do
      forest%family = grf_survival
      forest%n_outputs = forest%n_times
      call choose_options(options, local_options)
      call train_forest_core(forest, local_options, info)
   end subroutine survival_forest

   subroutine causal_survival_forest(x, y, w, event, horizon, forest, info, options, &
                                     survival_probability, y_hat, w_hat, sample_weights, failure_times)
      real(dp), intent(in) :: x(:,:) !! Predictor matrix used to discover heterogeneous treatment effects on survival outcomes.
      real(dp), intent(in) :: y(:) !! Nonnegative observed event or censoring times.
      real(dp), intent(in) :: w(:) !! Scalar treatment for each observation.
      integer, intent(in) :: event(:) !! Event indicators equal to one for failures and zero for censoring.
      real(dp), intent(in) :: horizon !! Positive time horizon defining the restricted-mean or survival-probability target.
      type(grf_forest), intent(out) :: forest !! Trained causal-survival forest with doubly robust pseudo-outcomes and nuisance estimates retained.
      integer, intent(out) :: info !! Zero on success; negative values indicate invalid dimensions, indicators, horizon, grid, or options.
      type(grf_options), intent(in), optional :: options !! Optional main-forest controls; nuisance forests disable CI grouping and use honest fitting.
      logical, intent(in), optional :: survival_probability !! If true target survival probability at horizon; otherwise target restricted mean survival time.
      real(dp), intent(in), optional :: y_hat(:) !! Optional E[f(T)|X] nuisance vector; when absent it is estimated from auxiliary survival forests.
      real(dp), intent(in), optional :: w_hat(:) !! Optional E[W|X] nuisance vector; when absent it is estimated by an OOB regression forest.
      real(dp), intent(in), optional :: sample_weights(:) !! Optional nonnegative observation weights matching rows of x.
      real(dp), intent(in), optional :: failure_times(:) !! Optional increasing nuisance-survival grid beginning no later than min(y); defaults to unique observed times.
      type(grf_options) :: local_options
      type(grf_options) :: nuisance_options
      type(grf_forest) :: survival_nuisance
      type(grf_forest) :: marginal_survival_nuisance
      type(grf_forest) :: censor_nuisance
      real(dp), allocatable :: xw(:,:)
      real(dp), allocatable :: xw_treated(:,:)
      real(dp), allocatable :: xw_control(:,:)
      real(dp), allocatable :: y_work(:)
      real(dp), allocatable :: f_y(:)
      real(dp), allocatable :: grid(:)
      real(dp), allocatable :: s_hat(:,:)
      real(dp), allocatable :: s1_hat(:,:)
      real(dp), allocatable :: s0_hat(:,:)
      real(dp), allocatable :: sy_hat(:,:)
      real(dp), allocatable :: c_hat(:,:)
      real(dp), allocatable :: nuisance_y(:)
      real(dp), allocatable :: nuisance_control(:)
      real(dp), allocatable :: nuisance_w(:)
      real(dp), allocatable :: numerator(:)
      real(dp), allocatable :: denominator(:)
      real(dp), allocatable :: c_y_hat(:)
      integer, allocatable :: event_work(:)
      integer, allocatable :: censor_event(:)
      integer :: i
      integer :: p
      logical :: target_survival_probability
      logical :: binary_w

      call choose_options(options, local_options)
      call initialize_forest(x, forest, info, options, sample_weights)
      if (info /= 0) return
      if (size(y) /= size(x,1) .or. size(w) /= size(x,1) .or. size(event) /= size(x,1)) then
         info = -113
         return
      end if
      if (horizon <= 0.0_dp .or. any(y < 0.0_dp) .or. any(event < 0) .or. any(event > 1)) then
         info = -114
         return
      end if
      if (count(event == 1) == 0) then
         info = -121
         return
      end if
      target_survival_probability = .false.
      if (present(survival_probability)) target_survival_probability = survival_probability
      allocate(y_work(size(y)), event_work(size(event)), f_y(size(y)))
      y_work = y
      event_work = event
      if (target_survival_probability) then
         do i = 1, size(y)
            if (y(i) > horizon) then
               f_y(i) = 1.0_dp
            else
               f_y(i) = 0.0_dp
            end if
         end do
      else
         do i = 1, size(y)
            if (y_work(i) >= horizon) then
               y_work(i) = horizon
               event_work(i) = 1
            end if
            f_y(i) = y_work(i)
         end do
      end if
      if (present(failure_times)) then
         if (size(failure_times) <= 2 .or. any(failure_times(2:) <= failure_times(:size(failure_times)-1))) then
            info = -122
            return
         end if
         allocate(grid(size(failure_times)))
         grid = failure_times
      else
         call unique_sorted(y_work, grid)
      end if
      if (size(grid) <= 2 .or. minval(y_work) < grid(1) .or. horizon < grid(1)) then
         info = -123
         return
      end if

      nuisance_options = local_options
      nuisance_options%num_trees = max(50, min(max(1, local_options%num_trees / 4), 500))
      nuisance_options%ci_group_size = 1
      nuisance_options%min_node_size = 5
      nuisance_options%honesty = .true.
      nuisance_options%honesty_fraction = 0.5_dp
      nuisance_options%honesty_prune_leaves = .true.
      nuisance_options%seed = local_options%seed + 7001

      allocate(nuisance_w(size(w)))
      call resolve_nuisance_vector(x, w, forest%sample_weights, w_hat, nuisance_options, 101, nuisance_w, info)
      if (info /= 0) return
      allocate(forest%treatment_variance_hat(size(w)))
      call estimate_treatment_variance(x, w, nuisance_w, forest%sample_weights, nuisance_options, &
         151, forest%treatment_variance_hat, info)
      if (info /= 0) return
      nuisance_options%min_node_size = 15

      p = size(x,2)
      allocate(xw(size(x,1),p+1))
      xw(:,1:p) = x
      xw(:,p+1) = w
      if (cluster_equalized(nuisance_options)) then
         call survival_forest(xw, y_work, event_work, survival_nuisance, info, nuisance_options)
      else
         call survival_forest(xw, y_work, event_work, survival_nuisance, info, nuisance_options, forest%sample_weights)
      end if
      if (info /= 0) return
      call predict_survival_nuisance(survival_nuisance, xw, grid, s_hat, oob=.true.)

      allocate(nuisance_y(size(y)))
      if (present(y_hat)) then
         if (size(y_hat) /= size(y)) then
            info = -124
            return
         end if
         nuisance_y = y_hat
      else
         binary_w = all((abs(w) <= 1.0e-12_dp) .or. (abs(w - 1.0_dp) <= 1.0e-12_dp))
         if (binary_w) then
            allocate(xw_treated(size(x,1),p+1), xw_control(size(x,1),p+1))
            xw_treated = xw
            xw_control = xw
            xw_treated(:,p+1) = 1.0_dp
            xw_control(:,p+1) = 0.0_dp
            call predict_survival_nuisance(survival_nuisance, xw_treated, grid, s1_hat, oob=.true.)
            call predict_survival_nuisance(survival_nuisance, xw_control, grid, s0_hat, oob=.true.)
            allocate(nuisance_control(size(y)))
            if (target_survival_probability) then
               call survival_at_horizon(s1_hat, grid, horizon, nuisance_y)
               call survival_at_horizon(s0_hat, grid, horizon, nuisance_control)
               nuisance_y = nuisance_w * nuisance_y + (1.0_dp - nuisance_w) * nuisance_control
            else
               call expected_survival_curve(s1_hat, grid, nuisance_y)
               call expected_survival_curve(s0_hat, grid, nuisance_control)
               nuisance_y = nuisance_w * nuisance_y + (1.0_dp - nuisance_w) * nuisance_control
            end if
         else
            nuisance_options%seed = local_options%seed + 8009
            if (cluster_equalized(nuisance_options)) then
               call survival_forest(x, y_work, event_work, marginal_survival_nuisance, info, nuisance_options)
            else
               call survival_forest(x, y_work, event_work, marginal_survival_nuisance, info, nuisance_options, &
                                    forest%sample_weights)
            end if
            if (info /= 0) return
            call predict_survival_nuisance(marginal_survival_nuisance, x, grid, sy_hat, oob=.true.)
            if (target_survival_probability) then
               call survival_at_horizon(sy_hat, grid, horizon, nuisance_y)
            else
               call expected_survival_curve(sy_hat, grid, nuisance_y)
            end if
         end if
      end if

      allocate(censor_event(size(event_work)))
      censor_event = 1 - event_work
      if (count(censor_event == 1) > 0) then
         nuisance_options%seed = local_options%seed + 9001
         if (cluster_equalized(nuisance_options)) then
            call survival_forest(xw, y_work, censor_event, censor_nuisance, info, nuisance_options)
         else
            call survival_forest(xw, y_work, censor_event, censor_nuisance, info, nuisance_options, forest%sample_weights)
         end if
         if (info /= 0) return
         call predict_survival_nuisance(censor_nuisance, xw, grid, c_hat, oob=.true.)
      else
         allocate(c_hat(size(x,1),size(grid)))
         c_hat = 1.0_dp
      end if

      if (target_survival_probability) then
         do i = 1, size(y_work)
            if (y_work(i) > horizon) then
               y_work(i) = horizon
               event_work(i) = 1
            end if
         end do
      end if
      allocate(numerator(size(y)), denominator(size(y)), c_y_hat(size(y)))
      call compute_causal_survival_psi(s_hat, c_hat, nuisance_y, w - nuisance_w, event_work, f_y, &
         y_work, grid, target_survival_probability, horizon, numerator, denominator, c_y_hat, info)
      if (info /= 0) return

      allocate(forest%y(size(y),1), forest%w(size(w),1))
      forest%y(:,1) = y_work
      forest%w(:,1) = w
      forest%event = event_work
      forest%horizon = horizon
      forest%survival_probability_target = target_survival_probability
      allocate(forest%y_hat(size(y),1), forest%w_hat(size(w),1))
      forest%y_hat(:,1) = nuisance_y
      forest%w_hat(:,1) = nuisance_w
      forest%causal_survival_numerator = numerator
      forest%causal_survival_denominator = denominator
      forest%censor_survival_at_y = c_y_hat
      forest%failure_times = grid
      forest%n_times = size(grid)
      forest%family = grf_causal_survival
      forest%n_outputs = 1
      forest%n_treatments = 1
      call train_forest_core(forest, local_options, info)
   end subroutine causal_survival_forest

   subroutine lm_forest(x, y, w, forest, info, options, y_hat, w_hat, sample_weights, stabilize_splits)
      real(dp), intent(in) :: x(:,:) !! Predictor matrix used to discover heterogeneity in local linear coefficients.
      real(dp), intent(in) :: y(:,:) !! Outcome matrix with observations in rows and one or more outcome columns.
      real(dp), intent(in) :: w(:,:) !! Regressor/treatment matrix with observations in rows and coefficient variables in columns.
      type(grf_forest), intent(out) :: forest !! Trained local-linear coefficient forest with nuisance matrices retained.
      integer, intent(out) :: info !! Zero on success; negative values indicate invalid dimensions, nuisance matrices, or options.
      type(grf_options), intent(in), optional :: options !! Optional forest controls; defaults are used when absent.
      real(dp), intent(in), optional :: y_hat(:,:) !! Optional conditional outcome nuisance matrix; defaults to out-of-bag multi-response regression-forest predictions.
      real(dp), intent(in), optional :: w_hat(:,:) !! Optional conditional regressor nuisance matrix; defaults to out-of-bag multi-response regression-forest predictions.
      real(dp), intent(in), optional :: sample_weights(:) !! Optional nonnegative observation weights matching rows of x.
      logical, intent(in), optional :: stabilize_splits !! If true request treatment-balance split stabilization; upstream lm_forest default is false.
      type(grf_options) :: local_options

      call choose_options(options, local_options)
      local_options%stabilize_splits = .false.
      if (present(stabilize_splits)) local_options%stabilize_splits = stabilize_splits
      call initialize_forest(x, forest, info, local_options, sample_weights)
      if (info /= 0) return
      if (size(y,1) /= size(x,1) .or. size(w,1) /= size(x,1) .or. size(y,2) < 1 .or. size(w,2) < 1) then
         info = -115
         return
      end if
      forest%y = y
      forest%w = w
      forest%n_outputs = size(y,2)
      forest%n_treatments = size(w,2)
      allocate(forest%y_hat(size(y,1),size(y,2)), forest%w_hat(size(w,1),size(w,2)))
      call resolve_nuisance_matrix(x, y, forest%sample_weights, y_hat, local_options, 4001, forest%y_hat, info)
      if (info /= 0) return
      call resolve_nuisance_matrix(x, w, forest%sample_weights, w_hat, local_options, 5003, forest%w_hat, info)
      if (info /= 0) return
      forest%family = grf_lm
      call train_forest_core(forest, local_options, info)
   end subroutine lm_forest

   subroutine multi_arm_causal_forest(x, y, w, forest, info, options, y_hat, w_hat, sample_weights, stabilize_splits)
      real(dp), intent(in) :: x(:,:) !! Predictor matrix for multi-arm or multi-treatment heterogeneous-effect estimation.
      real(dp), intent(in) :: y(:) !! Scalar outcome vector for each training observation.
      real(dp), intent(in) :: w(:,:) !! Non-baseline one-hot treatment contrasts for categorical arms, or a general multivariate treatment matrix.
      type(grf_forest), intent(out) :: forest !! Trained multi-treatment forest using GRF multi-causal relabeling mechanics.
      integer, intent(out) :: info !! Zero on success; negative values indicate invalid dimensions or nuisance estimates.
      type(grf_options), intent(in), optional :: options !! Optional forest controls; defaults are used when absent.
      real(dp), intent(in), optional :: y_hat(:) !! Optional E[Y|X] nuisance vector; defaults to an out-of-bag regression-forest prediction.
      real(dp), intent(in), optional :: w_hat(:,:) !! Optional propensity matrix: all categorical arms or non-baseline contrasts; general treatments require the same shape as w.
      real(dp), intent(in), optional :: sample_weights(:) !! Optional nonnegative observation weights matching rows of x.
      logical, intent(in), optional :: stabilize_splits !! If true request GRF treatment-balance split constraints; default true for categorical multi-arm forests.
      type(grf_options) :: local_options
      real(dp), allocatable :: yy(:,:)
      real(dp), allocatable :: yy_hat(:,:)
      real(dp), allocatable :: resolved_w_hat(:,:)
      real(dp), allocatable :: propensity_hat(:,:)
      logical :: categorical
      logical :: use_stabilized_splits
      integer :: nt

      if (size(y) /= size(x,1) .or. size(w,1) /= size(x,1) .or. size(w,2) < 1) then
         info = -116
         return
      end if
      nt = size(w,2)
      categorical = is_categorical_contrast_matrix(w)
      call choose_options(options, local_options)
      use_stabilized_splits = .true.
      if (present(stabilize_splits)) use_stabilized_splits = stabilize_splits
      local_options%stabilize_splits = use_stabilized_splits
      allocate(yy(size(y),1), yy_hat(size(y),1))
      yy(:,1) = y

      if (categorical) then
         allocate(resolved_w_hat(size(w,1),nt), propensity_hat(size(w,1),nt+1))
         if (present(w_hat)) then
            if (size(w_hat,1) /= size(w,1)) then
               info = -117
               return
            end if
            if (size(w_hat,2) == nt + 1) then
               propensity_hat = w_hat
               resolved_w_hat = w_hat(:,2:nt+1)
            else if (size(w_hat,2) == nt) then
               resolved_w_hat = w_hat
               propensity_hat(:,2:nt+1) = w_hat
               propensity_hat(:,1) = 1.0_dp - sum(w_hat, dim=2)
            else
               info = -118
               return
            end if
         else
            call resolve_categorical_propensity(x, w, local_options, sample_weights, propensity_hat, info)
            if (info /= 0) return
            resolved_w_hat = propensity_hat(:,2:nt+1)
         end if
         if (present(y_hat)) then
            if (size(y_hat) /= size(y)) then
               info = -119
               return
            end if
            yy_hat(:,1) = y_hat
            call lm_forest(x, yy, w, forest, info, local_options, yy_hat, resolved_w_hat, sample_weights, use_stabilized_splits)
         else
            call lm_forest(x, yy, w, forest, info, local_options, w_hat=resolved_w_hat, sample_weights=sample_weights, &
                           stabilize_splits=use_stabilized_splits)
         end if
         if (info /= 0) return
         forest%multi_arm_categorical = .true.
         forest%arm_propensity_hat = propensity_hat
         return
      end if

      if (present(w_hat)) then
         if (any(shape(w_hat) /= shape(w))) then
            info = -118
            return
         end if
      end if
      if (present(y_hat)) then
         if (size(y_hat) /= size(y)) then
            info = -117
            return
         end if
         yy_hat(:,1) = y_hat
         call lm_forest(x, yy, w, forest, info, local_options, yy_hat, w_hat, sample_weights, use_stabilized_splits)
      else
         call lm_forest(x, yy, w, forest, info, local_options, w_hat=w_hat, sample_weights=sample_weights, &
                        stabilize_splits=use_stabilized_splits)
      end if
      if (info == 0) forest%multi_arm_categorical = .false.
   end subroutine multi_arm_causal_forest

   subroutine ll_regression_forest(x, y, forest, info, options, sample_weights, enable_ll_split, &
                                  split_weight_penalty, split_ridge, split_variables, split_cutoff)
      real(dp), intent(in) :: x(:,:) !! Predictor matrix used both for forest localization and later local-linear correction.
      real(dp), intent(in) :: y(:) !! Scalar numeric outcome vector matching rows of x.
      type(grf_forest), intent(out) :: forest !! Local-linear forest; optionally trained with node-wise ridge-residual relabeling.
      integer, intent(out) :: info !! Zero on success; negative values indicate invalid dimensions, local-linear controls, or training options.
      type(grf_options), intent(in), optional :: options !! Optional forest controls; automatic forest tuning is unavailable when residual splitting is enabled.
      real(dp), intent(in), optional :: sample_weights(:) !! Optional nonnegative observation weights matching rows of x; retained as a Fortran extension with residual splitting.
      logical, intent(in), optional :: enable_ll_split !! If true, split on local ridge residuals as in upstream enable.ll.split; default false.
      logical, intent(in), optional :: split_weight_penalty !! If true, scale split ridge penalties by predictor second moments; default false.
      real(dp), intent(in), optional :: split_ridge !! Nonnegative ridge multiplier for residual splitting; default 0.1.
      integer, intent(in), optional :: split_variables(:) !! One-based predictor columns used in split-time local ridge regressions; defaults to all columns.
      integer, intent(in), optional :: split_cutoff !! Node size below which the overall ridge fit is reused; default floor(sqrt(n)).
      type(grf_options) :: local_options
      real(dp), allocatable :: design(:,:)
      real(dp), allocatable :: gram(:,:)
      real(dp), allocatable :: rhs(:)
      real(dp) :: normalization
      logical :: use_ll_split
      integer :: j
      integer :: solve_info

      use_ll_split = .false.
      if (present(enable_ll_split)) use_ll_split = enable_ll_split
      if (.not. use_ll_split) then
         call regression_forest(x, y, forest, info, options, sample_weights)
         return
      end if
      call initialize_forest(x, forest, info, options, sample_weights)
      if (info /= 0) return
      if (size(y) /= size(x,1)) then
         info = -130
         return
      end if
      call choose_options(options, local_options)
      if (local_options%tune_parameters) then
         info = -131
         return
      end if
      allocate(forest%y(size(y),1))
      forest%y(:,1) = y
      forest%family = grf_regression
      forest%n_outputs = 1
      forest%ll_split_enabled = .true.
      forest%ll_split_lambda = 0.1_dp
      if (present(split_ridge)) forest%ll_split_lambda = split_ridge
      if (forest%ll_split_lambda < 0.0_dp) then
         info = -132
         return
      end if
      forest%ll_split_weight_penalty = .false.
      if (present(split_weight_penalty)) forest%ll_split_weight_penalty = split_weight_penalty
      forest%ll_split_cutoff = floor(sqrt(real(size(x,1), dp)))
      if (present(split_cutoff)) forest%ll_split_cutoff = split_cutoff
      if (forest%ll_split_cutoff < 0 .or. forest%ll_split_cutoff > size(x,1)) then
         info = -133
         return
      end if
      if (present(split_variables)) then
         if (size(split_variables) < 1 .or. minval(split_variables) < 1 .or. maxval(split_variables) > size(x,2)) then
            info = -134
            return
         end if
         forest%ll_split_variables = split_variables
      else
         allocate(forest%ll_split_variables(size(x,2)))
         forest%ll_split_variables = [(j, j = 1, size(x,2))]
      end if
      if (forest%ll_split_cutoff > 0) then
         allocate(design(size(x,1),size(forest%ll_split_variables)+1))
         call build_ll_split_design(x, forest%ll_split_variables, design)
         allocate(gram(size(design,2),size(design,2)), rhs(size(design,2)))
         gram = matmul(transpose(design), design)
         rhs = matmul(transpose(design), y)
         if (forest%ll_split_weight_penalty) then
            do j = 2, size(gram,1)
               gram(j,j) = gram(j,j) + forest%ll_split_lambda * gram(j,j)
            end do
         else
            normalization = 0.0_dp
            do j = 1, size(gram,1)
               normalization = normalization + gram(j,j)
            end do
            normalization = normalization / real(size(gram,1), dp)
            do j = 2, size(gram,1)
               gram(j,j) = gram(j,j) + forest%ll_split_lambda * normalization
            end do
         end if
         allocate(forest%ll_split_overall_beta(size(design,2)))
         call solve_linear_system(gram, rhs, forest%ll_split_overall_beta, solve_info)
         if (solve_info /= 0) then
            forest%ll_split_overall_beta = 0.0_dp
            forest%ll_split_overall_beta(1) = weighted_mean(y, forest%sample_weights)
         end if
      else
         allocate(forest%ll_split_overall_beta(0))
      end if
      call train_forest_core(forest, local_options, info)
   end subroutine ll_regression_forest

   pure subroutine build_ll_split_design(x, variables, design)
      real(dp), intent(in) :: x(:,:) !! Predictor matrix from which selected split-time ridge covariates are extracted.
      integer, intent(in) :: variables(:) !! One-based predictor columns used in local-linear split relabeling.
      real(dp), intent(out) :: design(:,:) !! Design matrix with an intercept and selected predictors; missing predictor values are replaced by zero.
      integer :: i
      integer :: j

      design(:,1) = 1.0_dp
      do j = 1, size(variables)
         do i = 1, size(x,1)
            if (ieee_is_nan(x(i,variables(j)))) then
               design(i,j+1) = 0.0_dp
            else
               design(i,j+1) = x(i,variables(j))
            end if
         end do
      end do
   end subroutine build_ll_split_design

   subroutine boosted_regression_forest(x, y, forest, info, num_stages, learning_rate, options, sample_weights, &
                                        error_reduction, max_stages, tune_trees)
      real(dp), intent(in) :: x(:,:) !! Predictor matrix used by every residual-fitting forest stage.
      real(dp), intent(in) :: y(:) !! Scalar response vector whose conditional mean is approximated additively.
      type(grf_boosted_forest), intent(out) :: forest !! Sequence of residual regression forests with OOB training predictions and stage errors.
      integer, intent(out) :: info !! Zero on success; negative values indicate invalid dimensions, stage controls, shrinkage, or forest training.
      integer, intent(in), optional :: num_stages !! Fixed number of stages; when absent the stage count is selected by an OOB error-reduction check.
      real(dp), intent(in), optional :: learning_rate !! Positive stage multiplier; default one matches upstream unshrunk additive boosting.
      type(grf_options), intent(in), optional :: options !! Optional controls shared by full residual forests; stage seeds are offset deterministically.
      real(dp), intent(in), optional :: sample_weights(:) !! Optional nonnegative training weights passed to every full and tuning forest.
      real(dp), intent(in), optional :: error_reduction !! Required next-stage OOB error fraction for automatic stopping; default 0.97 in (0,1].
      integer, intent(in), optional :: max_stages !! Maximum automatic stage count when num_stages is absent; default five.
      integer, intent(in), optional :: tune_trees !! Number of trees in each small forest used to test an additional stage; default ten.
      type(grf_options) :: base_options
      type(grf_options) :: stage_options
      type(grf_options) :: small_options
      type(grf_forest), allocatable :: stage_buffer(:)
      type(grf_forest), allocatable :: final_stages(:)
      real(dp), allocatable :: residual(:)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: small_fitted(:)
      real(dp), allocatable :: weights(:)
      real(dp), allocatable :: error_buffer(:)
      real(dp), allocatable :: final_errors(:)
      real(dp) :: rate
      real(dp) :: threshold
      real(dp) :: previous_error
      real(dp) :: candidate_error
      integer :: limit
      integer :: small_tree_count
      integer :: actual_stages
      integer :: s
      logical :: fixed_stages

      info = 0
      if (size(y) /= size(x,1)) then
         info = -118
         return
      end if
      fixed_stages = present(num_stages)
      limit = 5
      if (present(max_stages)) limit = max_stages
      if (fixed_stages) limit = num_stages
      rate = 1.0_dp
      if (present(learning_rate)) rate = learning_rate
      threshold = 0.97_dp
      if (present(error_reduction)) threshold = error_reduction
      small_tree_count = 10
      if (present(tune_trees)) small_tree_count = tune_trees
      if (limit < 1 .or. rate <= 0.0_dp .or. threshold <= 0.0_dp .or. threshold > 1.0_dp .or. &
          small_tree_count < 2) then
         info = -119
         return
      end if
      allocate(weights(size(y)))
      weights = 1.0_dp
      if (present(sample_weights)) then
         if (size(sample_weights) /= size(y) .or. any(sample_weights < 0.0_dp) .or. sum(sample_weights) <= 0.0_dp) then
            info = -120
            return
         end if
         weights = sample_weights
      end if
      allocate(stage_buffer(limit), error_buffer(limit), residual(size(y)), fitted(size(y)))
      allocate(forest%training_prediction(size(y)))
      forest%training_prediction = 0.0_dp
      forest%intercept = 0.0_dp
      forest%learning_rate = rate
      residual = y
      error_buffer = 0.0_dp
      actual_stages = 0
      previous_error = huge(1.0_dp)
      call choose_options(options, base_options)
      do s = 1, limit
         if (s > 1 .and. .not. fixed_stages) then
            small_options = base_options
            small_options%num_trees = small_tree_count
            small_options%seed = base_options%seed + 65537 * s + 17
            if (cluster_equalized(small_options)) then
               call regression_forest(x, residual, stage_buffer(s), info, small_options)
            else
               call regression_forest(x, residual, stage_buffer(s), info, small_options, weights)
            end if
            if (info /= 0) return
            allocate(small_fitted(size(y)))
            call regression_oob_debiased_error(stage_buffer(s), residual, small_fitted, candidate_error)
            deallocate(small_fitted)
            if (.not. ieee_is_nan(candidate_error)) then
               if (candidate_error > threshold * previous_error) exit
            end if
         end if
         stage_options = base_options
         stage_options%seed = base_options%seed + 104729 * (s - 1)
         if (cluster_equalized(stage_options)) then
            call regression_forest(x, residual, stage_buffer(s), info, stage_options)
         else
            call regression_forest(x, residual, stage_buffer(s), info, stage_options, weights)
         end if
         if (info /= 0) return
         call regression_oob_debiased_error(stage_buffer(s), residual, fitted, error_buffer(s))
         if (s == 1 .and. base_options%tune_parameters) then
            stage_options = stage_buffer(s)%options
            stage_options%seed = base_options%seed
            stage_options%tune_parameters = .false.
            base_options = stage_options
         end if
         forest%training_prediction = forest%training_prediction + rate * fitted
         residual = y - forest%training_prediction
         actual_stages = s
         previous_error = error_buffer(s)
      end do
      if (actual_stages < 1) then
         info = -119
         return
      end if
      allocate(final_stages(actual_stages), final_errors(actual_stages))
      final_stages = stage_buffer(1:actual_stages)
      final_errors = error_buffer(1:actual_stages)
      call move_alloc(final_stages, forest%stages)
      call move_alloc(final_errors, forest%stage_error)
      forest%n_stages = actual_stages
   end subroutine boosted_regression_forest

   pure subroutine regression_oob_debiased_error(forest, outcome, prediction, mean_error)
      type(grf_forest), intent(in) :: forest !! Regression forest whose OOB tree leaves are used for prediction and Monte-Carlo error debiasing.
      real(dp), intent(in) :: outcome(:) !! Training response or residual vector aligned with forest rows.
      real(dp), intent(out) :: prediction(:) !! Out-of-bag forest prediction for each training row, with NaN when fewer than two trees contribute.
      real(dp), intent(out) :: mean_error !! Mean finite per-row debiased squared error, or NaN when no row has a usable estimate.
      real(dp), allocatable :: leaf_outcome(:)
      real(dp), allocatable :: leaf_weight(:)
      real(dp) :: sum_outcome
      real(dp) :: sum_weight
      real(dp) :: average_outcome_stat
      real(dp) :: average_weight_stat
      real(dp) :: rho
      real(dp) :: bias
      real(dp) :: row_error
      real(dp) :: error_sum
      integer :: i
      integer :: t
      integer :: node
      integer :: first
      integer :: count_leaf
      integer :: j
      integer :: sample
      integer :: valid_trees
      integer :: valid_rows

      allocate(leaf_outcome(size(forest%trees)), leaf_weight(size(forest%trees)))
      prediction = ieee_value(0.0_dp, ieee_quiet_nan)
      error_sum = 0.0_dp
      valid_rows = 0
      do i = 1, forest%n_train
         leaf_outcome = 0.0_dp
         leaf_weight = 0.0_dp
         valid_trees = 0
         do t = 1, size(forest%trees)
            if (forest%trees(t)%inbag(i)) cycle
            node = find_leaf_node(forest%trees(t), forest%x(i,:))
            if (node < 1 .or. node > forest%trees(t)%n_nodes) cycle
            count_leaf = forest%trees(t)%leaf_count(node)
            if (count_leaf <= 0) cycle
            first = forest%trees(t)%leaf_start(node)
            sum_outcome = 0.0_dp
            sum_weight = 0.0_dp
            do j = first, first + count_leaf - 1
               sample = forest%trees(t)%leaf_samples(j)
               sum_outcome = sum_outcome + forest%sample_weights(sample) * outcome(sample)
               sum_weight = sum_weight + forest%sample_weights(sample)
            end do
            if (sum_weight <= 1.0e-16_dp) cycle
            valid_trees = valid_trees + 1
            leaf_outcome(valid_trees) = sum_outcome / real(count_leaf, dp)
            leaf_weight(valid_trees) = sum_weight / real(count_leaf, dp)
         end do
         if (valid_trees <= 1) cycle
         average_outcome_stat = sum(leaf_outcome(1:valid_trees)) / real(valid_trees, dp)
         average_weight_stat = sum(leaf_weight(1:valid_trees)) / real(valid_trees, dp)
         if (average_weight_stat <= 1.0e-16_dp) cycle
         prediction(i) = average_outcome_stat / average_weight_stat
         bias = 0.0_dp
         do t = 1, valid_trees
            rho = (leaf_outcome(t) - prediction(i) * leaf_weight(t)) / average_weight_stat
            bias = bias + rho * rho
         end do
         bias = bias / real(valid_trees * (valid_trees - 1), dp)
         row_error = (prediction(i) - outcome(i))**2 - bias
         if (.not. ieee_is_nan(row_error)) then
            error_sum = error_sum + row_error
            valid_rows = valid_rows + 1
         end if
      end do
      if (valid_rows > 0) then
         mean_error = error_sum / real(valid_rows, dp)
      else
         mean_error = ieee_value(0.0_dp, ieee_quiet_nan)
      end if
   end subroutine regression_oob_debiased_error

   subroutine tune_regression_options(x, y, weights, input_options, tuned_options, info)
      real(dp), intent(in) :: x(:,:) !! Predictor matrix used to train mini forests across randomized tuning draws.
      real(dp), intent(in) :: y(:) !! Scalar regression target used for out-of-bag tuning error evaluation.
      real(dp), intent(in) :: weights(:) !! Nonnegative observation weights passed to every tuning forest.
      type(grf_options), intent(in) :: input_options !! Baseline options whose seven upstream tunable fields are optimized.
      type(grf_options), intent(out) :: tuned_options !! Selected options, reverting to the baseline when tuning cannot improve OOB error.
      integer, intent(out) :: info !! Zero on success; negative values indicate invalid tuning controls or an auxiliary forest failure.
      integer, parameter :: n_parameters = 7
      type(grf_rng_state) :: rng
      type(grf_options) :: candidate
      type(grf_options) :: trial_options
      type(grf_options) :: default_options
      type(grf_forest) :: mini_forest
      type(grf_forest) :: trial_forest
      type(grf_forest) :: default_forest
      real(dp), allocatable :: draws(:,:)
      real(dp), allocatable :: fit_draws(:,:)
      real(dp), allocatable :: fit_errors(:)
      real(dp), allocatable :: all_errors(:)
      real(dp), allocatable :: prediction(:)
      real(dp), allocatable :: covariance(:,:)
      real(dp), allocatable :: rhs(:)
      real(dp), allocatable :: alpha(:)
      real(dp), allocatable :: draw(:)
      real(dp) :: error
      real(dp) :: response_mean
      real(dp) :: response_variance
      real(dp) :: distance_squared
      real(dp) :: kernel_value
      real(dp) :: surface_value
      real(dp) :: best_surface
      real(dp) :: trial_error
      real(dp) :: default_error
      real(dp) :: length_scale
      real(dp) :: noise_ratio
      integer :: reps
      integer :: draws_count
      integer :: mini_trees
      integer :: keep_count
      integer :: minimum_keep
      integer :: r
      integer :: j
      integer :: k
      integer :: d
      integer :: solve_info
      integer :: best_index

      info = 0
      tuned_options = input_options
      reps = input_options%tune_num_reps
      draws_count = input_options%tune_num_draws
      mini_trees = input_options%tune_num_trees
      if (reps < 2 .or. draws_count < 1 .or. mini_trees < 2) then
         info = -128
         return
      end if
      allocate(draws(reps,n_parameters), all_errors(reps), prediction(size(y)))
      call rng%seed(input_options%seed)
      do r = 1, reps
         do j = 1, n_parameters
            draws(r,j) = rng%uniform()
         end do
         call options_from_tuning_draw(input_options, draws(r,:), size(x,1), size(x,2), candidate)
         candidate%num_trees = mini_trees
         candidate%ci_group_size = 1
         candidate%tune_parameters = .false.
         call train_regression_direct(x, y, weights, candidate, mini_forest, info)
         if (info /= 0) then
            all_errors(r) = ieee_value(0.0_dp, ieee_quiet_nan)
            info = 0
            cycle
         end if
         call regression_oob_debiased_error(mini_forest, y, prediction, error)
         all_errors(r) = error
      end do
      keep_count = count(.not. ieee_is_nan(all_errors))
      minimum_keep = min(10, max(2, reps / 2))
      if (keep_count < minimum_keep) return
      allocate(fit_draws(keep_count,n_parameters), fit_errors(keep_count))
      k = 0
      do r = 1, reps
         if (ieee_is_nan(all_errors(r))) cycle
         k = k + 1
         fit_draws(k,:) = draws(r,:)
         fit_errors(k) = all_errors(r)
      end do
      response_mean = sum(fit_errors) / real(keep_count, dp)
      if (keep_count > 1) then
         response_variance = sum((fit_errors - response_mean)**2) / real(keep_count - 1, dp)
      else
         response_variance = 0.0_dp
      end if
      if (response_variance <= 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, response_mean * response_mean)) return

      length_scale = 0.30_dp
      noise_ratio = 0.50_dp
      allocate(covariance(keep_count,keep_count), rhs(keep_count), alpha(keep_count))
      do r = 1, keep_count
         do k = 1, keep_count
            distance_squared = sum((fit_draws(r,:) - fit_draws(k,:))**2)
            covariance(r,k) = exp(-0.5_dp * distance_squared / (length_scale * length_scale))
         end do
         covariance(r,r) = covariance(r,r) + noise_ratio + 1.0e-10_dp
      end do
      rhs = (fit_errors - response_mean) / sqrt(response_variance)
      call solve_linear_system(covariance, rhs, alpha, solve_info)
      if (solve_info /= 0) return

      allocate(draw(n_parameters))
      call rng%seed(input_options%seed)
      best_surface = huge(1.0_dp)
      best_index = 0
      do d = 1, draws_count
         do j = 1, n_parameters
            draw(j) = rng%uniform()
         end do
         surface_value = response_mean
         do r = 1, keep_count
            distance_squared = sum((draw - fit_draws(r,:))**2)
            kernel_value = exp(-0.5_dp * distance_squared / (length_scale * length_scale))
            surface_value = surface_value + sqrt(response_variance) * kernel_value * alpha(r)
         end do
         if (surface_value < best_surface) then
            best_surface = surface_value
            best_index = d
            draws(1,:) = draw
         end if
      end do
      if (best_index == 0) return
      call options_from_tuning_draw(input_options, draws(1,:), size(x,1), size(x,2), candidate)

      trial_options = candidate
      trial_options%num_trees = 4 * mini_trees
      trial_options%ci_group_size = 1
      trial_options%tune_parameters = .false.
      call train_regression_direct(x, y, weights, trial_options, trial_forest, info)
      if (info /= 0) then
         info = 0
         return
      end if
      call regression_oob_debiased_error(trial_forest, y, prediction, trial_error)

      default_options = input_options
      default_options%num_trees = 4 * mini_trees
      default_options%ci_group_size = 1
      default_options%tune_parameters = .false.
      call train_regression_direct(x, y, weights, default_options, default_forest, info)
      if (info /= 0) then
         info = 0
         return
      end if
      call regression_oob_debiased_error(default_forest, y, prediction, default_error)
      if (.not. ieee_is_nan(trial_error)) then
         if (ieee_is_nan(default_error) .or. trial_error <= default_error) tuned_options = candidate
      end if
   end subroutine tune_regression_options

   pure subroutine options_from_tuning_draw(base, draw, n, p, options)
      type(grf_options), intent(in) :: base !! Baseline forest controls from which untuned settings are retained.
      real(dp), intent(in) :: draw(:) !! Seven independent values in [0,1] mapped to upstream GRF tuning parameter ranges.
      integer, intent(in) :: n !! Number of training observations used to scale the minimum-node-size draw.
      integer, intent(in) :: p !! Number of predictor columns used to scale the mtry draw.
      type(grf_options), intent(out) :: options !! Forest controls obtained by applying the GRF randomized tuning map.
      real(dp) :: log2_n
      real(dp) :: mtry_max

      options = base
      log2_n = log(real(max(n,2), dp)) / log(2.0_dp)
      options%min_node_size = max(1, floor(2.0_dp**(draw(1) * (log2_n - 4.0_dp))))
      options%sample_fraction = 0.05_dp + 0.45_dp * draw(2)
      mtry_max = min(real(p,dp), sqrt(real(p,dp)) + 20.0_dp)
      options%mtry = max(1, ceiling(mtry_max * draw(3)))
      options%alpha = 0.25_dp * draw(4)
      options%imbalance_penalty = -log(max(draw(5), epsilon(1.0_dp)))
      if (base%honesty) then
         options%honesty_fraction = 0.5_dp + 0.3_dp * draw(6)
         options%honesty_prune_leaves = draw(7) < 0.5_dp
      end if
   end subroutine options_from_tuning_draw

   subroutine train_regression_direct(x, y, weights, options, forest, info)
      real(dp), intent(in) :: x(:,:) !! Predictor matrix copied into the tuning forest.
      real(dp), intent(in) :: y(:) !! Scalar response vector copied into the tuning forest.
      real(dp), intent(in) :: weights(:) !! Nonnegative observation weights used by the tuning forest.
      type(grf_options), intent(in) :: options !! Fully resolved mini-forest controls with recursive tuning disabled.
      type(grf_forest), intent(out) :: forest !! Trained regression mini forest used only for tuning-error evaluation.
      integer, intent(out) :: info !! Zero on success or a negative code propagated from forest initialization/training.

      call initialize_forest(x, forest, info, options, weights)
      if (info /= 0) return
      allocate(forest%y(size(y),1))
      forest%y(:,1) = y
      forest%family = grf_regression
      forest%n_outputs = 1
      call train_forest_core(forest, options, info)
   end subroutine train_regression_direct

   pure logical function is_categorical_contrast_matrix(w) result(categorical)
      real(dp), intent(in) :: w(:,:) !! Candidate non-baseline one-hot treatment matrix with a zero row representing the baseline arm.
      real(dp) :: row_sum
      real(dp) :: tolerance
      integer :: i
      integer :: j

      tolerance = 100.0_dp * epsilon(1.0_dp)
      categorical = size(w,2) >= 1
      if (.not. categorical) return
      do i = 1, size(w,1)
         row_sum = 0.0_dp
         do j = 1, size(w,2)
            if (abs(w(i,j)) <= tolerance) then
               cycle
            else if (abs(w(i,j) - 1.0_dp) <= tolerance) then
               row_sum = row_sum + 1.0_dp
            else
               categorical = .false.
               return
            end if
         end do
         if (row_sum > 1.0_dp + tolerance) then
            categorical = .false.
            return
         end if
      end do
   end function is_categorical_contrast_matrix

   subroutine resolve_categorical_propensity(x, w, options, sample_weights, propensity, info)
      real(dp), intent(in) :: x(:,:) !! Predictor matrix used to fit the auxiliary probability forest.
      real(dp), intent(in) :: w(:,:) !! Non-baseline one-hot treatment matrix; all-zero rows identify the baseline treatment.
      type(grf_options), intent(in) :: options !! Main-forest controls inherited by the nuisance probability forest.
      real(dp), intent(in), optional :: sample_weights(:) !! Optional nonnegative observation weights; omitted for equalized-cluster training.
      real(dp), allocatable, intent(out) :: propensity(:,:) !! OOB class probabilities for baseline plus every non-baseline treatment arm.
      integer, intent(out) :: info !! Zero on success or a negative code propagated from the auxiliary probability forest.
      type(grf_options) :: nuisance_options
      type(grf_forest) :: nuisance_forest
      real(dp), allocatable :: kernel(:,:)
      real(dp), allocatable :: fallback_weight(:)
      integer, allocatable :: classes(:)
      real(dp) :: total_weight
      integer :: c
      integer :: i
      integer :: j

      allocate(classes(size(w,1)))
      classes = 1
      do i = 1, size(w,1)
         do j = 1, size(w,2)
            if (w(i,j) > 0.5_dp) then
               classes(i) = j + 1
               exit
            end if
         end do
      end do
      nuisance_options = options
      nuisance_options%ci_group_size = 1
      nuisance_options%num_trees = max(50, max(1, options%num_trees / 4))
      nuisance_options%seed = options%seed + 5501
      nuisance_options%tune_parameters = .false.
      if (cluster_equalized(nuisance_options)) then
         call probability_forest(x, classes, nuisance_forest, info, nuisance_options)
      else if (present(sample_weights)) then
         call probability_forest(x, classes, nuisance_forest, info, nuisance_options, sample_weights)
      else
         call probability_forest(x, classes, nuisance_forest, info, nuisance_options)
      end if
      if (info /= 0) return
      call compute_forest_weights(nuisance_forest, x, kernel, oob=.true.)
      allocate(propensity(size(w,1),size(w,2)+1))
      propensity = 0.0_dp
      do i = 1, size(classes)
         c = classes(i)
         propensity(:,c) = propensity(:,c) + kernel(:,i)
      end do

      allocate(fallback_weight(size(w,1)))
      if (present(sample_weights)) then
         fallback_weight = sample_weights
      else
         fallback_weight = 1.0_dp
      end if
      total_weight = sum(fallback_weight)
      do i = 1, size(w,1)
         if (sum(kernel(i,:)) > 0.0_dp) cycle
         do c = 1, size(propensity,2)
            propensity(i,c) = sum(fallback_weight, mask=classes == c) / total_weight
         end do
      end do
   end subroutine resolve_categorical_propensity

   pure subroutine initialize_forest(x, forest, info, options, sample_weights)
      real(dp), intent(in) :: x(:,:) !! Predictor matrix copied into the forest training object.
      type(grf_forest), intent(out) :: forest !! Fresh forest object initialized with predictors, dimensions, weights, and optional cluster metadata.
      integer, intent(out) :: info !! Zero on success; negative values indicate an empty matrix, invalid weights, or invalid clustering controls.
      type(grf_options), intent(in), optional :: options !! Optional controls, including cluster labels and equalized-cluster sampling settings.
      real(dp), intent(in), optional :: sample_weights(:) !! Optional nonnegative observation weights matching rows of x; disallowed with equalized cluster weights.
      integer, allocatable :: cluster_counts(:)

      info = 0
      if (size(x,1) < 2 .or. size(x,2) < 1) then
         info = -90
         return
      end if
      forest%x = x
      forest%n_train = size(x,1)
      forest%p = size(x,2)
      allocate(forest%sample_weights(size(x,1)))
      forest%sample_weights = 1.0_dp
      if (present(options)) then
         if (allocated(options%clusters)) then
            if (size(options%clusters) /= size(x,1)) then
               info = -92
               return
            end if
            call canonicalize_clusters(options%clusters, forest%clusters, forest%n_clusters)
            if (forest%n_clusters < 1) then
               info = -93
               return
            end if
            allocate(cluster_counts(forest%n_clusters))
            call count_cluster_members(forest%clusters, forest%n_clusters, cluster_counts)
            forest%equalize_cluster_weights = options%equalize_cluster_weights
            if (forest%equalize_cluster_weights) then
               if (present(sample_weights)) then
                  info = -94
                  return
               end if
               forest%samples_per_cluster = minval(cluster_counts)
            else
               forest%samples_per_cluster = maxval(cluster_counts)
            end if
         end if
      end if
      if (present(sample_weights)) then
         if (size(sample_weights) /= size(x,1) .or. any(sample_weights < 0.0_dp) .or. sum(sample_weights) <= 0.0_dp) then
            info = -91
            return
         end if
         forest%sample_weights = sample_weights
      end if
   end subroutine initialize_forest

   pure subroutine canonicalize_clusters(labels, clusters, n_clusters)
      integer, intent(in) :: labels(:) !! Caller-supplied integer cluster labels; equal values identify observations from the same sampling cluster.
      integer, allocatable, intent(out) :: clusters(:) !! Canonical one-based cluster index for every observation, preserving first-appearance order.
      integer, intent(out) :: n_clusters !! Number of distinct cluster labels represented in the training rows.
      integer, allocatable :: seen(:)
      integer :: i
      integer :: j
      integer :: found

      allocate(clusters(size(labels)), seen(size(labels)))
      n_clusters = 0
      do i = 1, size(labels)
         found = 0
         do j = 1, n_clusters
            if (labels(i) == seen(j)) then
               found = j
               exit
            end if
         end do
         if (found == 0) then
            n_clusters = n_clusters + 1
            seen(n_clusters) = labels(i)
            found = n_clusters
         end if
         clusters(i) = found
      end do
   end subroutine canonicalize_clusters

   pure subroutine count_cluster_members(clusters, n_clusters, counts)
      integer, intent(in) :: clusters(:) !! Canonical one-based cluster index for every training observation.
      integer, intent(in) :: n_clusters !! Number of distinct cluster identifiers represented by clusters.
      integer, intent(out) :: counts(:) !! Observation count for each cluster; must have length at least n_clusters.
      integer :: i

      counts = 0
      do i = 1, size(clusters)
         if (clusters(i) >= 1 .and. clusters(i) <= n_clusters) counts(clusters(i)) = counts(clusters(i)) + 1
      end do
   end subroutine count_cluster_members

   pure subroutine choose_options(options, local_options)
      type(grf_options), intent(in), optional :: options !! Optional caller-supplied forest controls.
      type(grf_options), intent(out) :: local_options !! Effective controls, equal to defaults when options is absent.

      local_options = grf_options()
      if (present(options)) local_options = options
   end subroutine choose_options

   pure subroutine predict_survival_nuisance(forest, x_new, times, prediction, oob)
      type(grf_forest), intent(in) :: forest !! Auxiliary survival forest whose local weights define Nelson-Aalen survival estimates.
      real(dp), intent(in) :: x_new(:,:) !! Predictor rows at which nuisance survival curves are evaluated.
      real(dp), intent(in) :: times(:) !! Increasing time grid on which survival probabilities are requested.
      real(dp), allocatable, intent(out) :: prediction(:,:) !! Nelson-Aalen survival probabilities with query rows and time-grid columns.
      logical, intent(in), optional :: oob !! If true for training-row queries, exclude trees containing query row i in-bag.
      real(dp), allocatable :: kernel(:,:)
      real(dp) :: hazard
      real(dp) :: risk
      real(dp) :: deaths
      real(dp) :: event_time
      integer :: i
      integer :: j
      integer :: k
      integer :: cursor

      call compute_forest_weights(forest, x_new, kernel, oob)
      allocate(prediction(size(x_new,1),size(times)))
      do i = 1, size(x_new,1)
         hazard = 0.0_dp
         cursor = 1
         do j = 1, size(times)
            do while (cursor <= forest%n_times)
               event_time = forest%failure_times(cursor)
               if (event_time > times(j)) exit
               risk = 0.0_dp
               deaths = 0.0_dp
               do k = 1, forest%n_train
                  if (forest%y(k,1) >= event_time) risk = risk + kernel(i,k)
                  if (forest%event(k) == 1 .and. nearly_equal_time(forest%y(k,1), event_time)) then
                     deaths = deaths + kernel(i,k)
                  end if
               end do
               if (risk > 100.0_dp * epsilon(1.0_dp)) hazard = hazard + deaths / risk
               cursor = cursor + 1
            end do
            prediction(i,j) = exp(-hazard)
         end do
      end do
   end subroutine predict_survival_nuisance

   pure subroutine expected_survival_curve(survival, times, expectation)
      real(dp), intent(in) :: survival(:,:) !! Conditional survival curves evaluated on an increasing finite time grid.
      real(dp), intent(in) :: times(:) !! Increasing time grid corresponding to survival columns.
      real(dp), intent(out) :: expectation(:) !! Trapezoid-free step-function integral E[min(T,max(times))] for each row.
      integer :: i
      integer :: j

      do i = 1, size(survival,1)
         expectation(i) = times(1)
         do j = 1, size(times) - 1
            expectation(i) = expectation(i) + survival(i,j) * (times(j+1) - times(j))
         end do
      end do
   end subroutine expected_survival_curve

   pure subroutine survival_at_horizon(survival, times, horizon, value)
      real(dp), intent(in) :: survival(:,:) !! Conditional survival curves evaluated on an increasing time grid.
      real(dp), intent(in) :: times(:) !! Increasing time grid corresponding to survival columns.
      real(dp), intent(in) :: horizon !! Target time at which survival probability is requested.
      real(dp), intent(out) :: value(:) !! Survival probability at the greatest grid time not exceeding horizon, or one before the first time.
      integer :: index
      integer :: j

      index = 0
      do j = 1, size(times)
         if (times(j) <= horizon) index = j
      end do
      if (index == 0) then
         value = 1.0_dp
      else
         value = survival(:,index)
      end if
   end subroutine survival_at_horizon

   pure subroutine compute_causal_survival_psi(s_hat, c_hat, y_hat, w_centered, event, f_y, y, times, &
      survival_probability, horizon, numerator, denominator, c_y_hat, info)
      real(dp), intent(in) :: s_hat(:,:) !! Estimated conditional event survival S(t|X,W) on the nuisance time grid.
      real(dp), intent(in) :: c_hat(:,:) !! Estimated conditional censoring survival C(t|X,W) on the nuisance time grid.
      real(dp), intent(in) :: y_hat(:) !! Estimated marginal target mean E[f(T)|X].
      real(dp), intent(in) :: w_centered(:) !! Residualized treatment W-E[W|X].
      integer, intent(in) :: event(:) !! Event indicators after any target-horizon truncation.
      real(dp), intent(in) :: f_y(:) !! Observed target transform: truncated time for RMST or I(T>horizon) for survival probability.
      real(dp), intent(in) :: y(:) !! Event/censoring times after target-horizon truncation where applicable.
      real(dp), intent(in) :: times(:) !! Increasing nuisance time grid used by the survival and censoring curves.
      logical, intent(in) :: survival_probability !! Selects survival-probability rather than RMST Q-function construction.
      real(dp), intent(in) :: horizon !! Target horizon used for the survival-probability conditional Q function.
      real(dp), intent(out) :: numerator(:) !! Doubly robust causal-survival numerator pseudo-outcome for each observation.
      real(dp), intent(out) :: denominator(:) !! Causal-survival denominator pseudo-outcome, equal to residualized treatment squared.
      real(dp), intent(out) :: c_y_hat(:) !! Estimated censoring survival evaluated at each observed/truncated time.
      integer, intent(out) :: info !! Zero on success or a negative value for incompatible dimensions or a time outside the grid.
      real(dp), allocatable :: q_hat(:,:)
      real(dp), allocatable :: dot_products(:,:)
      real(dp) :: q_y
      real(dp) :: numerator_one
      real(dp) :: numerator_two
      real(dp) :: previous_c
      real(dp) :: current_c
      real(dp) :: dlambda
      real(dp) :: integrand
      real(dp) :: s_safe
      integer :: n
      integer :: m
      integer :: i
      integer :: j
      integer :: y_index
      integer :: horizon_index

      info = 0
      n = size(y)
      m = size(times)
      if (size(s_hat,1) /= n .or. size(c_hat,1) /= n .or. size(s_hat,2) /= m .or. size(c_hat,2) /= m) then
         info = -125
         return
      end if
      allocate(q_hat(n,m))
      if (survival_probability) then
         horizon_index = 0
         do j = 1, m
            if (times(j) <= horizon) horizon_index = j
         end do
         if (horizon_index == 0) then
            info = -126
            return
         end if
         do i = 1, n
            do j = 1, m
               s_safe = max(s_hat(i,j), 1.0e-12_dp)
               q_hat(i,j) = s_hat(i,horizon_index) / s_safe
            end do
            q_hat(i,horizon_index:m) = 1.0_dp
         end do
      else
         allocate(dot_products(n,m-1))
         do j = 1, m - 1
            dot_products(:,j) = s_hat(:,j) * (times(j+1) - times(j))
         end do
         q_hat(:,1) = sum(dot_products, dim=2)
         do j = 2, m - 1
            q_hat(:,j) = q_hat(:,j-1) - dot_products(:,j-1)
         end do
         do j = 1, m - 1
            q_hat(:,j) = times(j) + q_hat(:,j) / max(s_hat(:,j), 1.0e-12_dp)
         end do
         q_hat(:,m) = times(m)
      end if

      do i = 1, n
         y_index = 0
         do j = 1, m
            if (times(j) <= y(i)) y_index = j
         end do
         if (y_index == 0) then
            info = -127
            return
         end if
         current_c = max(c_hat(i,y_index), 1.0e-12_dp)
         c_y_hat(i) = current_c
         q_y = q_hat(i,y_index)
         numerator_one = (real(event(i),dp) * (f_y(i) - y_hat(i)) + &
            real(1 - event(i),dp) * (q_y - y_hat(i))) * w_centered(i) / current_c
         numerator_two = 0.0_dp
         previous_c = 1.0_dp
         do j = 1, y_index
            current_c = max(c_hat(i,j), 1.0e-12_dp)
            dlambda = log(max(previous_c, 1.0e-12_dp) / current_c)
            integrand = dlambda / current_c * (q_hat(i,j) - y_hat(i))
            numerator_two = numerator_two + integrand * w_centered(i)
            previous_c = current_c
         end do
         numerator(i) = numerator_one - numerator_two
         denominator(i) = w_centered(i) * w_centered(i)
      end do
   end subroutine compute_causal_survival_psi

   pure logical function strictly_increasing(values) result(ok)
      real(dp), intent(in) :: values(:) !! Candidate numeric grid whose ordering is checked without modifying it.
      integer :: i

      ok = size(values) >= 1
      do i = 2, size(values)
         if (values(i) <= values(i-1)) then
            ok = .false.
            return
         end if
      end do
   end function strictly_increasing

   pure integer function failure_interval_index(value, grid) result(index)
      real(dp), intent(in) :: value !! Observed event or censoring time to map onto the training failure grid.
      real(dp), intent(in) :: grid(:) !! Strictly increasing survival-grid values.
      integer :: j

      index = 0
      do j = 1, size(grid)
         if (grid(j) > value) exit
         index = j
      end do
   end function failure_interval_index

   pure logical function cluster_equalized(options) result(equalized)
      type(grf_options), intent(in) :: options !! Forest controls whose cluster metadata determines whether explicit unit weights must be omitted.

      equalized = allocated(options%clusters) .and. options%equalize_cluster_weights
   end function cluster_equalized

   pure elemental logical function nearly_equal_time(a, b) result(equal)
      real(dp), intent(in) :: a !! First nonnegative survival time being compared.
      real(dp), intent(in) :: b !! Second nonnegative survival time being compared.

      equal = abs(a - b) <= 10.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(a), abs(b))
   end function nearly_equal_time

   subroutine estimate_treatment_variance(x, w, w_hat, weights, options, seed_offset, result, info)
      real(dp), intent(in) :: x(:,:) !! Predictor matrix used when conditional treatment variance requires an auxiliary forest.
      real(dp), intent(in) :: w(:) !! Observed scalar treatment or exposure values.
      real(dp), intent(in) :: w_hat(:) !! Estimated conditional treatment mean E[W|X] aligned with w.
      real(dp), intent(in) :: weights(:) !! Nonnegative observation weights used by any auxiliary variance forest.
      type(grf_options), intent(in) :: options !! Base forest controls; CI grouping is disabled for the auxiliary variance forest.
      integer, intent(in) :: seed_offset !! Deterministic seed offset separating the variance nuisance forest from other auxiliaries.
      real(dp), intent(out) :: result(:) !! Estimated Var[W|X], using e(X)(1-e(X)) for binary W or OOB squared-residual regression otherwise.
      integer, intent(out) :: info !! Zero on success or a negative value propagated from auxiliary regression training.
      type(grf_options) :: variance_options
      type(grf_forest) :: variance_forest
      real(dp), allocatable :: residual_squared(:)
      real(dp), allocatable :: kernel(:,:)
      real(dp) :: fallback
      logical :: binary_w
      integer :: i

      info = 0
      binary_w = all((abs(w) <= 1.0e-12_dp) .or. (abs(w - 1.0_dp) <= 1.0e-12_dp))
      if (binary_w) then
         result = max(w_hat * (1.0_dp - w_hat), 0.0_dp)
         return
      end if
      residual_squared = (w - w_hat)**2
      variance_options = options
      variance_options%ci_group_size = 1
      variance_options%seed = options%seed + seed_offset
      if (cluster_equalized(variance_options)) then
         call regression_forest(x, residual_squared, variance_forest, info, variance_options)
      else
         call regression_forest(x, residual_squared, variance_forest, info, variance_options, weights)
      end if
      if (info /= 0) return
      call compute_forest_weights(variance_forest, x, kernel, oob=.true.)
      result = matmul(kernel, residual_squared)
      fallback = weighted_mean(residual_squared, weights)
      do i = 1, size(result)
         if (sum(kernel(i,:)) <= 0.0_dp) result(i) = fallback
         result(i) = max(result(i), 0.0_dp)
      end do
   end subroutine estimate_treatment_variance

   subroutine resolve_nuisance_vector(x, values, weights, supplied, options, seed_offset, result, info)
      real(dp), intent(in) :: x(:,:) !! Predictor matrix used to fit the auxiliary nuisance regression forest.
      real(dp), intent(in) :: values(:) !! Scalar nuisance target aligned with rows of x.
      real(dp), intent(in) :: weights(:) !! Nonnegative observation weights passed to the auxiliary forest.
      real(dp), intent(in), optional :: supplied(:) !! Optional caller-supplied nuisance estimates; when present no auxiliary forest is fit.
      type(grf_options), intent(in) :: options !! Main-forest options inherited by the auxiliary forest except grouped-CI training.
      integer, intent(in) :: seed_offset !! Deterministic offset separating the auxiliary forest RNG stream from the main forest.
      real(dp), intent(out) :: result(:) !! Row-specific nuisance estimates, preferably out-of-bag auxiliary-forest predictions.
      integer, intent(out) :: info !! Zero on success or a negative code propagated from validation or auxiliary training.
      type(grf_forest) :: nuisance_forest
      type(grf_options) :: nuisance_options
      real(dp), allocatable :: kernel(:,:)
      real(dp) :: fallback
      integer :: i

      info = 0
      if (present(supplied)) then
         if (size(supplied) /= size(values)) then
            info = -92
            return
         end if
         result = supplied
         return
      end if
      nuisance_options = options
      nuisance_options%ci_group_size = 1
      nuisance_options%seed = options%seed + seed_offset
      if (cluster_equalized(nuisance_options)) then
         call regression_forest(x, values, nuisance_forest, info, nuisance_options)
      else
         call regression_forest(x, values, nuisance_forest, info, nuisance_options, weights)
      end if
      if (info /= 0) return
      call compute_forest_weights(nuisance_forest, x, kernel, oob=.true.)
      result = matmul(kernel, values)
      fallback = weighted_mean(values, weights)
      do i = 1, size(values)
         if (sum(kernel(i,:)) <= 0.0_dp) result(i) = fallback
      end do
   end subroutine resolve_nuisance_vector

   subroutine resolve_nuisance_matrix(x, values, weights, supplied, options, seed_offset, result, info)
      real(dp), intent(in) :: x(:,:) !! Predictor matrix used to fit the auxiliary multi-response nuisance forest.
      real(dp), intent(in) :: values(:,:) !! Nuisance target matrix with observations in rows.
      real(dp), intent(in) :: weights(:) !! Nonnegative observation weights passed to the auxiliary forest.
      real(dp), intent(in), optional :: supplied(:,:) !! Optional caller-supplied nuisance matrix; when present no auxiliary forest is fit.
      type(grf_options), intent(in) :: options !! Main-forest options inherited by the auxiliary forest except grouped-CI training.
      integer, intent(in) :: seed_offset !! Deterministic offset separating the auxiliary forest RNG stream from the main forest.
      real(dp), intent(out) :: result(:,:) !! Row-specific nuisance estimates from out-of-bag multi-response forest predictions.
      integer, intent(out) :: info !! Zero on success or a negative code propagated from validation or auxiliary training.
      type(grf_forest) :: nuisance_forest
      type(grf_options) :: nuisance_options
      real(dp), allocatable :: kernel(:,:)
      real(dp) :: fallback
      integer :: i
      integer :: j

      info = 0
      if (present(supplied)) then
         if (any(shape(supplied) /= shape(values))) then
            info = -93
            return
         end if
         result = supplied
         return
      end if
      nuisance_options = options
      nuisance_options%ci_group_size = 1
      nuisance_options%seed = options%seed + seed_offset
      if (cluster_equalized(nuisance_options)) then
         call multi_regression_forest(x, values, nuisance_forest, info, nuisance_options)
      else
         call multi_regression_forest(x, values, nuisance_forest, info, nuisance_options, weights)
      end if
      if (info /= 0) return
      call compute_forest_weights(nuisance_forest, x, kernel, oob=.true.)
      result = matmul(kernel, values)
      do i = 1, size(values,1)
         if (sum(kernel(i,:)) <= 0.0_dp) then
            do j = 1, size(values,2)
               fallback = weighted_mean(values(:,j), weights)
               result(i,j) = fallback
            end do
         end if
      end do
   end subroutine resolve_nuisance_matrix


end module grf_train
