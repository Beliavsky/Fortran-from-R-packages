module grf_analysis
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_nan
   use grf_kinds, only : dp
   use grf_rng, only : grf_rng_state
   use grf_types, only : grf_tree, grf_forest, grf_ate_result, grf_linear_result, grf_rate_result
   use grf_types, only : grf_causal, grf_instrumental, grf_causal_survival, grf_lm
   use grf_core, only : find_leaf_node, compute_forest_weights
   use grf_predict, only : predict_causal_forest, predict_instrumental_forest
   use grf_predict, only : predict_causal_survival_forest, predict_multi_arm_causal_forest
   use grf_utils, only : weighted_mean, standard_error_mean, weighted_linear_regression, invert_matrix
   use grf_utils, only : cluster_standard_error_mean, cluster_robust_linear_regression
   use grf_utils, only : sort_indices_by_values
   implicit none
   private

   public :: get_tree
   public :: split_frequencies
   public :: variable_importance
   public :: get_forest_weights
   public :: get_leaf_node
   public :: merge_forests
   public :: get_scores_causal_forest
   public :: get_scores_instrumental_forest
   public :: get_scores_multi_arm_causal_forest
   public :: get_scores_causal_survival_forest
   public :: average_treatment_effect
   public :: test_calibration
   public :: best_linear_projection
   public :: rank_average_treatment_effect_fit
   public :: rank_average_treatment_effect

contains

   pure subroutine get_tree(forest, index, tree, info)
      type(grf_forest), intent(in) :: forest !! Trained forest containing the requested tree.
      integer, intent(in) :: index !! One-based tree index in the inclusive range 1:size(forest%trees).
      type(grf_tree), intent(out) :: tree !! Deep copy of the selected tree structure and honest leaf membership.
      integer, intent(out) :: info !! Zero on success or -1 when index is outside the forest.

      info = 0
      if (.not. allocated(forest%trees) .or. index < 1 .or. index > size(forest%trees)) then
         info = -1
         return
      end if
      tree = forest%trees(index)
   end subroutine get_tree

   pure subroutine split_frequencies(forest, frequencies, max_depth)
      type(grf_forest), intent(in) :: forest !! Trained forest whose internal split-variable counts are summarized by depth.
      integer, allocatable, intent(out) :: frequencies(:,:) !! Depth-by-variable count matrix; row one corresponds to root splits at depth zero.
      integer, intent(in), optional :: max_depth !! Maximum number of depth levels to include; default four.
      integer :: depth_count
      integer :: t
      integer :: node
      integer :: d
      integer :: v

      depth_count = 4
      if (present(max_depth)) depth_count = max(1, max_depth)
      allocate(frequencies(depth_count,forest%p))
      frequencies = 0
      if (.not. allocated(forest%trees)) return
      do t = 1, size(forest%trees)
         do node = 1, forest%trees(t)%n_nodes
            v = forest%trees(t)%split_var(node)
            if (v <= 0 .or. v > forest%p) cycle
            d = forest%trees(t)%depth(node) + 1
            if (d >= 1 .and. d <= depth_count) frequencies(d,v) = frequencies(d,v) + 1
         end do
      end do
   end subroutine split_frequencies

   pure subroutine variable_importance(forest, importance, decay_exponent, max_depth)
      type(grf_forest), intent(in) :: forest !! Trained forest whose depth-weighted split counts define importance.
      real(dp), allocatable, intent(out) :: importance(:) !! Nonnegative variable-importance scores normalized to sum to one when any split occurs.
      real(dp), intent(in), optional :: decay_exponent !! Positive exponent controlling depth decay; default two.
      integer, intent(in), optional :: max_depth !! Number of depth levels included in the summary; default four.
      integer, allocatable :: frequencies(:,:)
      real(dp) :: decay
      real(dp) :: total
      integer :: depth_count
      integer :: d

      decay = 2.0_dp
      if (present(decay_exponent)) decay = decay_exponent
      depth_count = 4
      if (present(max_depth)) depth_count = max(1, max_depth)
      call split_frequencies(forest, frequencies, depth_count)
      allocate(importance(forest%p))
      importance = 0.0_dp
      do d = 1, depth_count
         importance = importance + real(frequencies(d,:),dp) / real(d,dp)**decay
      end do
      total = sum(importance)
      if (total > 0.0_dp) importance = importance / total
   end subroutine variable_importance

   pure subroutine get_forest_weights(forest, x_new, weights, oob)
      type(grf_forest), intent(in) :: forest !! Trained forest defining the adaptive nearest-neighbor kernel.
      real(dp), intent(in) :: x_new(:,:) !! Predictor rows whose weights over training observations are requested.
      real(dp), allocatable, intent(out) :: weights(:,:) !! Dense query-by-training matrix of normalized forest kernel weights.
      logical, intent(in), optional :: oob !! If true for training rows, exclude in-bag trees observation by observation.

      call compute_forest_weights(forest, x_new, weights, oob)
   end subroutine get_forest_weights

   pure subroutine get_leaf_node(tree, x_new, node_id)
      type(grf_tree), intent(in) :: tree !! Tree structure returned by get_tree or extracted directly from a trained forest.
      real(dp), intent(in) :: x_new(:,:) !! Predictor rows routed through the tree; columns must include every split variable used.
      integer, allocatable, intent(out) :: node_id(:) !! One-based terminal node reached by each row of x_new.
      integer :: i

      allocate(node_id(size(x_new,1)))
      do i = 1, size(x_new,1)
         node_id(i) = find_leaf_node(tree, x_new(i,:))
      end do
   end subroutine get_leaf_node

   pure subroutine merge_forests(forests, merged, info)
      type(grf_forest), intent(in) :: forests(:) !! Forests trained on identical data and with the same family to concatenate.
      type(grf_forest), intent(out) :: merged !! Copy of the first forest with all tree arrays concatenated in input order.
      integer, intent(out) :: info !! Zero on success; negative values indicate no forests or incompatible family/training dimensions.
      type(grf_tree), allocatable :: combined(:)
      integer :: total
      integer :: i
      integer :: offset

      info = 0
      if (size(forests) < 1) then
         info = -1
         return
      end if
      total = 0
      do i = 1, size(forests)
         if (forests(i)%family /= forests(1)%family .or. forests(i)%n_train /= forests(1)%n_train .or. &
             forests(i)%p /= forests(1)%p) then
            info = -2
            return
         end if
         if (.not. allocated(forests(i)%trees)) then
            info = -3
            return
         end if
         total = total + size(forests(i)%trees)
      end do
      merged = forests(1)
      if (allocated(merged%trees)) deallocate(merged%trees)
      allocate(combined(total))
      offset = 0
      do i = 1, size(forests)
         combined(offset+1:offset+size(forests(i)%trees)) = forests(i)%trees
         offset = offset + size(forests(i)%trees)
      end do
      call move_alloc(combined, merged%trees)
      merged%options%num_trees = total
   end subroutine merge_forests

   pure subroutine get_scores_causal_forest(forest, scores)
      type(grf_forest), intent(in) :: forest !! Trained scalar causal forest used to form out-of-bag augmented inverse-weight scores.
      real(dp), allocatable, intent(out) :: scores(:) !! Doubly robust treatment-effect score for each training observation.
      real(dp), allocatable :: tau(:)
      real(dp), allocatable :: residual(:)
      real(dp), allocatable :: wr(:)
      real(dp), allocatable :: variance_hat(:)
      real(dp) :: tolerance
      integer :: i

      call predict_causal_forest(forest, forest%x, tau, oob=.true.)
      wr = forest%w(:,1) - forest%w_hat(:,1)
      residual = forest%y(:,1) - (forest%y_hat(:,1) + tau * wr)
      allocate(variance_hat(forest%n_train), scores(forest%n_train))
      if (allocated(forest%treatment_variance_hat)) then
         variance_hat = forest%treatment_variance_hat
      else
         variance_hat = weighted_mean(wr * wr, forest%sample_weights)
      end if
      tolerance = 100.0_dp * epsilon(1.0_dp)
      do i = 1, forest%n_train
         if (variance_hat(i) <= tolerance .or. ieee_is_nan(tau(i))) then
            scores(i) = ieee_value(0.0_dp, ieee_quiet_nan)
         else
            scores(i) = tau(i) + wr(i) * residual(i) / variance_hat(i)
         end if
      end do
   end subroutine get_scores_causal_forest

   pure subroutine get_scores_instrumental_forest(forest, scores)
      type(grf_forest), intent(in) :: forest !! Trained instrumental forest used to form out-of-bag locally robust IV scores.
      real(dp), allocatable, intent(out) :: scores(:) !! Doubly robust local-IV score for each training observation.
      real(dp), allocatable :: tau(:)
      real(dp), allocatable :: residual(:)
      real(dp), allocatable :: wr(:)
      real(dp), allocatable :: zr(:)
      real(dp) :: denom
      real(dp) :: instrument_variance
      real(dp) :: tolerance
      logical :: binary_z
      integer :: i

      call predict_instrumental_forest(forest, forest%x, tau, oob=.true.)
      wr = forest%w(:,1) - forest%w_hat(:,1)
      zr = forest%z - forest%z_hat
      residual = forest%y(:,1) - (forest%y_hat(:,1) + tau * wr)
      allocate(scores(forest%n_train))
      binary_z = all((abs(forest%z) <= 1.0e-12_dp) .or. (abs(forest%z - 1.0_dp) <= 1.0e-12_dp))
      tolerance = 100.0_dp * epsilon(1.0_dp)
      if (binary_z .and. allocated(forest%compliance_hat)) then
         do i = 1, forest%n_train
            instrument_variance = forest%z_hat(i) * (1.0_dp - forest%z_hat(i))
            denom = instrument_variance * forest%compliance_hat(i)
            if (abs(denom) <= tolerance .or. ieee_is_nan(tau(i))) then
               scores(i) = ieee_value(0.0_dp, ieee_quiet_nan)
            else
               scores(i) = tau(i) + zr(i) * residual(i) / denom
            end if
         end do
      else
         denom = weighted_mean(zr * wr, forest%sample_weights)
         if (abs(denom) <= tolerance) then
            scores = ieee_value(0.0_dp, ieee_quiet_nan)
         else
            scores = tau + zr * residual / denom
         end if
      end if
   end subroutine get_scores_instrumental_forest

   pure subroutine get_scores_multi_arm_causal_forest(forest, scores)
      type(grf_forest), intent(in) :: forest !! Trained multi-arm forest supplying OOB contrasts, nuisances, and categorical propensities when applicable.
      real(dp), allocatable, intent(out) :: scores(:,:) !! Training-row by treatment-contrast matrix of doubly robust or generalized orthogonal scores.
      real(dp), allocatable :: tau(:,:)
      real(dp), allocatable :: wr(:,:)
      real(dp), allocatable :: yr(:)
      real(dp), allocatable :: gram(:,:)
      real(dp), allocatable :: gram_inv(:,:)
      real(dp), allocatable :: influence(:)
      real(dp) :: baseline_mu
      real(dp) :: residual
      real(dp) :: propensity
      real(dp) :: tolerance
      integer :: observed_arm
      integer :: i
      integer :: j
      integer :: k
      integer :: info

      call predict_multi_arm_causal_forest(forest, forest%x, tau, oob=.true.)
      tolerance = 100.0_dp * epsilon(1.0_dp)
      if (forest%multi_arm_categorical .and. allocated(forest%arm_propensity_hat)) then
         allocate(scores(forest%n_train,forest%n_treatments))
         scores = tau
         do i = 1, forest%n_train
            if (any([(ieee_is_nan(tau(i,j)), j=1,forest%n_treatments)])) then
               scores(i,:) = ieee_value(0.0_dp, ieee_quiet_nan)
               cycle
            end if
            observed_arm = 1
            do j = 1, forest%n_treatments
               if (forest%w(i,j) > 0.5_dp) then
                  observed_arm = j + 1
                  exit
               end if
            end do
            baseline_mu = forest%y_hat(i,1) - dot_product(forest%arm_propensity_hat(i,2:), tau(i,:))
            if (observed_arm == 1) then
               residual = forest%y(i,1) - baseline_mu
               propensity = forest%arm_propensity_hat(i,1)
               if (propensity <= tolerance) then
                  scores(i,:) = ieee_value(0.0_dp, ieee_quiet_nan)
               else
                  scores(i,:) = tau(i,:) - residual / propensity
               end if
            else
               residual = forest%y(i,1) - (baseline_mu + tau(i,observed_arm-1))
               propensity = forest%arm_propensity_hat(i,observed_arm)
               if (propensity <= tolerance) then
                  scores(i,:) = ieee_value(0.0_dp, ieee_quiet_nan)
               else
                  scores(i,observed_arm-1) = tau(i,observed_arm-1) + residual / propensity
               end if
            end if
         end do
         return
      end if

      wr = forest%w(:,1:forest%n_treatments) - forest%w_hat(:,1:forest%n_treatments)
      yr = forest%y(:,1) - forest%y_hat(:,1)
      allocate(gram(forest%n_treatments,forest%n_treatments), gram_inv(forest%n_treatments,forest%n_treatments))
      gram = 0.0_dp
      do j = 1, forest%n_treatments
         do k = 1, forest%n_treatments
            gram(j,k) = weighted_mean(wr(:,j) * wr(:,k), forest%sample_weights)
         end do
      end do
      call invert_matrix(gram, gram_inv, info)
      allocate(scores(forest%n_train,forest%n_treatments), influence(forest%n_treatments))
      if (info /= 0) then
         scores = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      do i = 1, forest%n_train
         residual = yr(i) - dot_product(wr(i,:), tau(i,:))
         influence = matmul(gram_inv, wr(i,:))
         scores(i,:) = tau(i,:) + influence * residual
      end do
   end subroutine get_scores_multi_arm_causal_forest

   pure subroutine get_scores_causal_survival_forest(forest, scores)
      type(grf_forest), intent(in) :: forest !! Trained causal-survival forest supplying stored doubly robust numerator and denominator moments.
      real(dp), allocatable, intent(out) :: scores(:) !! Doubly robust horizon-specific treatment-effect score for each training row.
      real(dp), allocatable :: tau(:)
      real(dp), allocatable :: psi(:)
      real(dp), allocatable :: variance_hat(:)
      real(dp) :: tolerance
      integer :: i

      call predict_causal_survival_forest(forest, forest%x, tau, oob=.true.)
      psi = forest%causal_survival_numerator - forest%causal_survival_denominator * tau
      allocate(variance_hat(forest%n_train), scores(forest%n_train))
      if (allocated(forest%treatment_variance_hat)) then
         variance_hat = forest%treatment_variance_hat
      else
         variance_hat = weighted_mean(forest%causal_survival_denominator, forest%sample_weights)
      end if
      tolerance = 100.0_dp * epsilon(1.0_dp)
      do i = 1, forest%n_train
         if (variance_hat(i) <= tolerance .or. ieee_is_nan(tau(i))) then
            scores(i) = ieee_value(0.0_dp, ieee_quiet_nan)
         else
            scores(i) = tau(i) + psi(i) / variance_hat(i)
         end if
      end do
   end subroutine get_scores_causal_survival_forest

   pure subroutine average_treatment_effect(forest, result, info)
      type(grf_forest), intent(in) :: forest !! Causal, instrumental, causal-survival, or multi-treatment forest to summarize.
      type(grf_ate_result), intent(out) :: result !! Weighted average effect estimates and cluster-robust influence-score standard errors.
      integer, intent(out) :: info !! Zero on success or -1 when the forest family has no treatment-effect score implementation.
      real(dp), allocatable :: scores(:)
      real(dp), allocatable :: multi_scores(:,:)
      real(dp), allocatable :: weights(:)
      integer, allocatable :: clusters(:)
      integer :: j

      info = 0
      call forest_observation_weights(forest, weights)
      call forest_inference_clusters(forest, clusters)
      select case (forest%family)
      case (grf_causal)
         call get_scores_causal_forest(forest, scores)
         allocate(result%estimate(1), result%std_err(1))
         result%estimate(1) = weighted_mean(scores, weights)
         result%std_err(1) = cluster_standard_error_mean(scores, weights, clusters)
      case (grf_instrumental)
         call get_scores_instrumental_forest(forest, scores)
         allocate(result%estimate(1), result%std_err(1))
         result%estimate(1) = weighted_mean(scores, weights)
         result%std_err(1) = cluster_standard_error_mean(scores, weights, clusters)
      case (grf_causal_survival)
         call get_scores_causal_survival_forest(forest, scores)
         allocate(result%estimate(1), result%std_err(1))
         result%estimate(1) = weighted_mean(scores, weights)
         result%std_err(1) = cluster_standard_error_mean(scores, weights, clusters)
      case (grf_lm)
         call get_scores_multi_arm_causal_forest(forest, multi_scores)
         allocate(result%estimate(size(multi_scores,2)), result%std_err(size(multi_scores,2)))
         do j = 1, size(multi_scores,2)
            result%estimate(j) = weighted_mean(multi_scores(:,j), weights)
            result%std_err(j) = cluster_standard_error_mean(multi_scores(:,j), weights, clusters)
         end do
      case default
         info = -1
      end select
   end subroutine average_treatment_effect

   pure subroutine test_calibration(forest, result, info)
      type(grf_forest), intent(in) :: forest !! Scalar treatment-effect forest whose out-of-bag predictions are calibration regressors.
      type(grf_linear_result), intent(out) :: result !! Intercept and slope from robust regression of orthogonal scores on predicted effects.
      integer, intent(out) :: info !! Zero on success; nonzero for unsupported families or singular calibration design.
      real(dp), allocatable :: scores(:)
      real(dp), allocatable :: prediction(:)
      real(dp), allocatable :: design(:,:)

      select case (forest%family)
      case (grf_causal)
         call get_scores_causal_forest(forest, scores)
         call predict_causal_forest(forest, forest%x, prediction, oob=.true.)
      case (grf_instrumental)
         call get_scores_instrumental_forest(forest, scores)
         call predict_instrumental_forest(forest, forest%x, prediction, oob=.true.)
      case (grf_causal_survival)
         call get_scores_causal_survival_forest(forest, scores)
         call predict_causal_survival_forest(forest, forest%x, prediction, oob=.true.)
      case default
         info = -1
         return
      end select
      allocate(design(forest%n_train,2), result%coefficient(2), result%std_err(2))
      design(:,1) = 1.0_dp
      design(:,2) = prediction
      call calibration_regression(forest, design, scores, result%coefficient, result%std_err, info)
      result%info = info
   end subroutine test_calibration

   pure subroutine best_linear_projection(forest, a, result, info)
      type(grf_forest), intent(in) :: forest !! Scalar causal/IV/causal-survival forest whose scores are projected on covariates.
      real(dp), intent(in) :: a(:,:) !! Projection covariates with one row per training observation; an intercept is added automatically.
      type(grf_linear_result), intent(out) :: result !! Weighted robust linear-projection coefficients and HC0 standard errors.
      integer, intent(out) :: info !! Zero on success; negative for unsupported family or row mismatch, otherwise solver status.
      real(dp), allocatable :: scores(:)
      real(dp), allocatable :: design(:,:)

      if (size(a,1) /= forest%n_train) then
         info = -2
         return
      end if
      select case (forest%family)
      case (grf_causal)
         call get_scores_causal_forest(forest, scores)
      case (grf_instrumental)
         call get_scores_instrumental_forest(forest, scores)
      case (grf_causal_survival)
         call get_scores_causal_survival_forest(forest, scores)
      case default
         info = -1
         return
      end select
      allocate(design(forest%n_train,size(a,2)+1))
      allocate(result%coefficient(size(a,2)+1), result%std_err(size(a,2)+1))
      design(:,1) = 1.0_dp
      design(:,2:) = a
      call calibration_regression(forest, design, scores, result%coefficient, result%std_err, info)
      result%info = info
   end subroutine best_linear_projection

   pure subroutine forest_observation_weights(forest, weights)
      type(grf_forest), intent(in) :: forest !! Forest whose inference weights respect optional equalized-cluster training.
      real(dp), allocatable, intent(out) :: weights(:) !! Observation weights for influence-function summaries.
      integer, allocatable :: counts(:)
      integer :: i

      weights = forest%sample_weights
      if (.not. forest%equalize_cluster_weights .or. .not. allocated(forest%clusters)) return
      allocate(counts(forest%n_clusters))
      counts = 0
      do i = 1, forest%n_train
         counts(forest%clusters(i)) = counts(forest%clusters(i)) + 1
      end do
      do i = 1, forest%n_train
         weights(i) = weights(i) / real(max(1,counts(forest%clusters(i))),dp)
      end do
   end subroutine forest_observation_weights

   pure subroutine forest_inference_clusters(forest, clusters)
      type(grf_forest), intent(in) :: forest !! Forest whose training clusters define robust inference units when available.
      integer, allocatable, intent(out) :: clusters(:) !! Canonical cluster labels, or one unique label per row for unclustered forests.
      integer :: i

      if (allocated(forest%clusters)) then
         clusters = forest%clusters
      else
         allocate(clusters(forest%n_train))
         clusters = [(i, i = 1, forest%n_train)]
      end if
   end subroutine forest_inference_clusters

   pure subroutine calibration_regression(forest, design, response, coefficient, std_err, info)
      type(grf_forest), intent(in) :: forest !! Forest supplying inference weights and optional cluster labels.
      real(dp), intent(in) :: design(:,:) !! Calibration/projection design matrix with observations in rows.
      real(dp), intent(in) :: response(:) !! Orthogonal score response vector matching design rows.
      real(dp), intent(out) :: coefficient(:) !! Weighted least-squares coefficients.
      real(dp), intent(out) :: std_err(:) !! HC3 cluster-robust standard errors, reducing to observation clusters when unclustered.
      integer, intent(out) :: info !! Zero on success or solver status for singular designs.
      real(dp), allocatable :: weights(:)
      integer, allocatable :: clusters(:)

      call forest_observation_weights(forest, weights)
      call forest_inference_clusters(forest, clusters)
      call cluster_robust_linear_regression(design, response, weights, clusters, coefficient, std_err, info, hc3=.true.)
   end subroutine calibration_regression

   subroutine rank_average_treatment_effect_fit(scores, priorities, result, info, q, target_qini, weights, &
                                                bootstrap_reps, seed, clusters)
      real(dp), intent(in) :: scores(:) !! Doubly robust treatment-effect scores for the evaluation sample.
      real(dp), intent(in) :: priorities(:) !! Larger values denote observations ranked as higher treatment priority.
      type(grf_rate_result), intent(out) :: result !! Targeting-operator curve, pointwise standard errors, and AUTOC/QINI area estimate.
      integer, intent(out) :: info !! Zero on success; negative values indicate dimension, probability-grid, or weight errors.
      real(dp), intent(in), optional :: q(:) !! Optional top-fraction grid in (0,1]; defaults to ten equally spaced values from 0.1 through 1.
      logical, intent(in), optional :: target_qini !! If true integrate q*TOC(q); otherwise integrate TOC(q) for AUTOC.
      real(dp), intent(in), optional :: weights(:) !! Optional nonnegative evaluation weights; defaults to uniform weights.
      integer, intent(in), optional :: bootstrap_reps !! Optional deterministic half-sample bootstrap replications for area standard error; default 100.
      integer, intent(in), optional :: seed !! Deterministic bootstrap seed; default 42.
      integer, intent(in), optional :: clusters(:) !! Optional cluster labels; bootstrap half-samples clusters rather than individual rows when supplied.
      real(dp), allocatable :: local_q(:)
      real(dp), allocatable :: local_weights(:)
      real(dp), allocatable :: boot_estimates(:)
      real(dp), allocatable :: selected_scores(:)
      real(dp), allocatable :: selected_weights(:)
      integer, allocatable :: order(:)
      integer, allocatable :: keep(:)
      integer, allocatable :: cluster_index(:)
      integer, allocatable :: cluster_labels(:)
      logical, allocatable :: keep_cluster(:)
      real(dp), allocatable :: neg_priority(:)
      type(grf_rng_state) :: rng
      real(dp) :: overall
      real(dp) :: previous_q
      real(dp) :: previous_value
      real(dp) :: current_value
      logical :: qini
      integer :: n
      integer :: i
      integer :: k
      integer :: b
      integer :: reps
      integer :: nkeep
      integer :: local_seed
      integer :: ng
      integer :: g
      integer :: found

      info = 0
      n = size(scores)
      if (n < 2 .or. size(priorities) /= n) then
         info = -1
         return
      end if
      if (present(clusters)) then
         if (size(clusters) /= n) then
            info = -4
            return
         end if
         allocate(cluster_index(n), cluster_labels(n))
         ng = 0
         do i = 1, n
            found = 0
            do g = 1, ng
               if (cluster_labels(g) == clusters(i)) then
                  found = g
                  exit
               end if
            end do
            if (found == 0) then
               ng = ng + 1
               cluster_labels(ng) = clusters(i)
               found = ng
            end if
            cluster_index(i) = found
         end do
      else
         ng = n
         allocate(cluster_index(n))
         cluster_index = [(i, i = 1, n)]
      end if
      allocate(local_weights(n))
      local_weights = 1.0_dp
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp) .or. sum(weights) <= 0.0_dp) then
            info = -2
            return
         end if
         local_weights = weights
      end if
      if (present(q)) then
         allocate(local_q(size(q)))
         local_q = q
      else
         allocate(local_q(10))
         local_q = [(0.1_dp * real(i,dp), i = 1, 10)]
      end if
      if (any(local_q <= 0.0_dp) .or. any(local_q > 1.0_dp)) then
         info = -3
         return
      end if
      allocate(neg_priority(n), order(n))
      neg_priority = -priorities
      call sort_indices_by_values(neg_priority, order)
      allocate(result%q(size(local_q)), result%toc(size(local_q)), result%toc_std_err(size(local_q)))
      result%q = local_q
      overall = weighted_mean(scores, local_weights)
      do i = 1, size(local_q)
         k = max(1, min(n, ceiling(local_q(i) * real(n,dp))))
         allocate(selected_scores(k), selected_weights(k))
         selected_scores = scores(order(1:k))
         selected_weights = local_weights(order(1:k))
         result%toc(i) = weighted_mean(selected_scores, selected_weights) - overall
         result%toc_std_err(i) = cluster_standard_error_mean(selected_scores, selected_weights, cluster_index(order(1:k)))
         deallocate(selected_scores, selected_weights)
      end do
      qini = .false.
      if (present(target_qini)) qini = target_qini
      result%estimate = 0.0_dp
      previous_q = 0.0_dp
      previous_value = result%toc(1)
      do i = 1, size(local_q)
         current_value = result%toc(i)
         if (qini) current_value = local_q(i) * current_value
         if (i == 1) then
            if (qini) previous_value = 0.0_dp
         else
            previous_value = result%toc(i-1)
            if (qini) previous_value = local_q(i-1) * previous_value
         end if
         result%estimate = result%estimate + 0.5_dp * (current_value + previous_value) * (local_q(i) - previous_q)
         previous_q = local_q(i)
      end do

      reps = 100
      if (present(bootstrap_reps)) reps = max(0, bootstrap_reps)
      result%std_err = 0.0_dp
      if (reps > 1) then
         local_seed = 42
         if (present(seed)) local_seed = seed
         call rng%seed(local_seed)
         allocate(boot_estimates(reps), keep(n), keep_cluster(ng))
         do b = 1, reps
            do g = 1, ng
               keep_cluster(g) = rng%uniform() < 0.5_dp
            end do
            if (.not. any(keep_cluster)) keep_cluster(1) = .true.
            nkeep = 0
            do i = 1, n
               if (keep_cluster(cluster_index(i))) then
                  nkeep = nkeep + 1
                  keep(nkeep) = i
               end if
            end do
            if (nkeep < 2) then
               boot_estimates(b) = result%estimate
            else
               call rate_area_subset(scores, priorities, local_weights, keep(1:nkeep), local_q, qini, boot_estimates(b))
            end if
         end do
         result%std_err = sqrt(sum((boot_estimates - sum(boot_estimates) / real(reps,dp))**2) / real(reps - 1,dp))
      end if
   end subroutine rank_average_treatment_effect_fit

   subroutine rank_average_treatment_effect(forest, priorities, result, info, q, target_qini, bootstrap_reps, seed)
      type(grf_forest), intent(in) :: forest !! Scalar causal/IV/causal-survival forest supplying doubly robust evaluation scores.
      real(dp), intent(in) :: priorities(:) !! Larger values rank training observations as higher treatment priority.
      type(grf_rate_result), intent(out) :: result !! AUTOC/QINI targeting curve, pointwise uncertainty, and integrated estimate.
      integer, intent(out) :: info !! Zero on success or negative for unsupported forest family or invalid RATE inputs.
      real(dp), intent(in), optional :: q(:) !! Optional top-fraction grid in (0,1]; defaults to 0.1,0.2,...,1.0.
      logical, intent(in), optional :: target_qini !! If true target QINI weighting; otherwise target AUTOC.
      integer, intent(in), optional :: bootstrap_reps !! Optional deterministic half-sample bootstrap replications; default 100.
      integer, intent(in), optional :: seed !! Deterministic bootstrap seed; default 42.
      real(dp), allocatable :: scores(:)

      select case (forest%family)
      case (grf_causal)
         call get_scores_causal_forest(forest, scores)
      case (grf_instrumental)
         call get_scores_instrumental_forest(forest, scores)
      case (grf_causal_survival)
         call get_scores_causal_survival_forest(forest, scores)
      case default
         info = -1
         return
      end select
      if (allocated(forest%clusters)) then
         call rank_average_treatment_effect_fit(scores, priorities, result, info, q, target_qini, &
                                                forest%sample_weights, bootstrap_reps, seed, forest%clusters)
      else
         call rank_average_treatment_effect_fit(scores, priorities, result, info, q, target_qini, &
                                                forest%sample_weights, bootstrap_reps, seed)
      end if
   end subroutine rank_average_treatment_effect

   pure subroutine rate_area_subset(scores, priorities, weights, subset, q, qini, area)
      real(dp), intent(in) :: scores(:) !! Full evaluation scores from which the bootstrap subset is drawn.
      real(dp), intent(in) :: priorities(:) !! Full priority vector used to rank subset observations.
      real(dp), intent(in) :: weights(:) !! Full nonnegative evaluation-weight vector.
      integer, intent(in) :: subset(:) !! One-based indices retained in the current bootstrap half-sample.
      real(dp), intent(in) :: q(:) !! Ordered top-fraction grid used in the area calculation.
      logical, intent(in) :: qini !! Selects q-weighted QINI area rather than unweighted AUTOC area.
      real(dp), intent(out) :: area !! Integrated targeting-operator curve for this subset.
      real(dp), allocatable :: sub_score(:)
      real(dp), allocatable :: sub_priority(:)
      real(dp), allocatable :: sub_weight(:)
      real(dp), allocatable :: neg_priority(:)
      integer, allocatable :: order(:)
      real(dp) :: overall
      real(dp) :: toc
      real(dp) :: previous
      real(dp) :: current
      real(dp) :: previous_q
      integer :: i
      integer :: k
      integer :: n

      n = size(subset)
      sub_score = scores(subset)
      sub_priority = priorities(subset)
      sub_weight = weights(subset)
      allocate(neg_priority(n), order(n))
      neg_priority = -sub_priority
      call sort_indices_by_values(neg_priority, order)
      overall = weighted_mean(sub_score, sub_weight)
      area = 0.0_dp
      previous = 0.0_dp
      previous_q = 0.0_dp
      do i = 1, size(q)
         k = max(1, min(n, ceiling(q(i) * real(n,dp))))
         toc = weighted_mean(sub_score(order(1:k)), sub_weight(order(1:k))) - overall
         current = toc
         if (qini) current = q(i) * current
         area = area + 0.5_dp * (previous + current) * (q(i) - previous_q)
         previous = current
         previous_q = q(i)
      end do
   end subroutine rate_area_subset

end module grf_analysis
