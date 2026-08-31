! Substantive-model-compatible missing-covariate proposal/rebuild driver for jomo.
! The implementation consolidates the reusable numerical loop in jomo1smcC,
! jomo1ransmcC, jomo1ranhrsmcC, jomo2smcC and jomo2hrsmcC.
! Upstream jomo 2.7-6 by Matteo Quartagno and James Carpenter; License: GPL-2.
! Modern Fortran translation, 2026. Distributed under GPL-2.0-only.
module jomo_smc
   use jomo_kinds, only : dp
   use jomo_rng, only : rng_state, rng_uniform, rng_normal
   use jomo_linalg, only : inverse_spd, solve_spd, quadratic_form
   use jomo_distributions, only : mvnormal_sample, invwishart_sample
   use jomo_latent, only : decode_categories
   use jomo_substantive, only : normal_cdf, cox_partial_loglik_ordered
   use jomo_substantive, only : sample_gaussian_coefficients, sample_linear_variance
   use jomo_substantive, only : cox_coordinate_newton, update_ordinal_thresholds
   implicit none
   private

   integer, parameter, public :: smc_linear = 0
   integer, parameter, public :: smc_binary_probit = 1
   integer, parameter, public :: smc_cox = 2
   integer, parameter, public :: smc_ordinal_probit = 3

   integer, parameter, public :: smc_l1_continuous = 1
   integer, parameter, public :: smc_l1_categorical = 2
   integer, parameter, public :: smc_l2_continuous = 3
   integer, parameter, public :: smc_l2_categorical = 4

   type, public :: smc_factor_spec
      integer :: variable = 0
      integer :: source = smc_l1_continuous
      integer :: power = 1
   end type smc_factor_spec

   type, public :: smc_term_spec
      type(smc_factor_spec), allocatable :: factor(:)
   end type smc_term_spec

   type, public :: smc_design_spec
      logical :: intercept = .true.
      type(smc_term_spec), allocatable :: term(:)
   end type smc_design_spec

   type, public :: smc_substantive_model
      integer :: family = smc_linear
      real(dp), allocatable :: response(:)
      logical, allocatable :: response_observed(:)
      integer, allocatable :: category(:)
      logical, allocatable :: event(:)
      real(dp), allocatable :: beta(:)
      real(dp) :: variance = 1.0_dp
      real(dp), allocatable :: thresholds(:)
      real(dp), allocatable :: latent_response(:)
      real(dp), allocatable :: random_effects(:, :)
      real(dp), allocatable :: random_covariance(:, :)
   end type smc_substantive_model

   type, public :: smc_sweep_stats
      integer :: level1_proposals = 0
      integer :: level1_accepted = 0
      integer :: level1_gibbs = 0
      integer :: level2_proposals = 0
      integer :: level2_accepted = 0
      integer :: level2_gibbs = 0
   contains
      procedure :: level1_acceptance => smc_level1_acceptance
      procedure :: level2_acceptance => smc_level2_acceptance
   end type smc_sweep_stats

   public :: smc_design_columns
   public :: smc_build_design
   public :: smc_mark_level1_predictors
   public :: smc_mark_level2_predictors
   public :: smc_substantive_loglik
   public :: smc_level1_sweep
   public :: smc_level2_sweep
   public :: smc_update_substantive
   public :: smc_compatible_iteration

contains

   pure integer function smc_design_columns(spec, n_levels1, n_levels2) result(ncol)
      type(smc_design_spec), intent(in) :: spec !! Substantive fixed- or random-effect term specification.
      integer, intent(in) :: n_levels1(:) !! Category counts for level-1 categorical covariates.
      integer, intent(in) :: n_levels2(:) !! Category counts for level-2 categorical covariates.
      integer :: j
      integer :: k
      integer :: width

      ncol = 0
      if (spec%intercept) ncol = 1
      if (.not. allocated(spec%term)) return
      do j = 1, size(spec%term)
         if (.not. allocated(spec%term(j)%factor)) error stop "smc_design_columns: term factors are not allocated"
         if (size(spec%term(j)%factor) == 0) error stop "smc_design_columns: empty term"
         width = 1
         do k = 1, size(spec%term(j)%factor)
            width = width * factor_width(spec%term(j)%factor(k), n_levels1, n_levels2)
         end do
         ncol = ncol + width
      end do
   end function smc_design_columns

   subroutine smc_build_design(level1_latent, n_con1, n_levels1, cluster, level2_latent, n_con2, n_levels2, spec, design)
      real(dp), intent(in) :: level1_latent(:, :) !! Complete level-1 latent covariates, shape n by p1.
      integer, intent(in) :: n_con1 !! Number of leading continuous level-1 covariates.
      integer, intent(in) :: n_levels1(:) !! Category counts for level-1 categorical covariates.
      integer, intent(in) :: cluster(:) !! One-based cluster labels for rows; use all ones for single-level models.
      real(dp), intent(in) :: level2_latent(:, :) !! Complete cluster-level latent covariates, shape G by p2; may have zero columns.
      integer, intent(in) :: n_con2 !! Number of leading continuous level-2 covariates.
      integer, intent(in) :: n_levels2(:) !! Category counts for level-2 categorical covariates; may be empty.
      type(smc_design_spec), intent(in) :: spec !! Polynomial/interaction/dummy expansion specification.
      real(dp), intent(out) :: design(:, :) !! Expanded design matrix, shape n by smc_design_columns(spec,...).
      integer :: n
      integer :: pcol
      integer :: row
      integer :: term_index
      integer :: combo
      integer :: width
      integer :: column
      integer, allocatable :: cat1(:, :)
      integer, allocatable :: cat2(:, :)

      n = size(level1_latent, 1)
      if (size(cluster) /= n) error stop "smc_build_design: cluster length mismatch"
      if (any(cluster < 1) .or. any(cluster > size(level2_latent, 1))) error stop "smc_build_design: invalid cluster label"
      pcol = smc_design_columns(spec, n_levels1, n_levels2)
      if (size(design, 1) /= n .or. size(design, 2) /= pcol) error stop "smc_build_design: output shape mismatch"
      if (size(level1_latent, 2) /= latent_width(n_con1, n_levels1)) error stop "smc_build_design: level-1 latent width mismatch"
      if (size(level2_latent, 2) /= latent_width(n_con2, n_levels2)) error stop "smc_build_design: level-2 latent width mismatch"

      allocate(cat1(n, size(n_levels1)), cat2(size(level2_latent, 1), size(n_levels2)))
      if (size(n_levels1) > 0) call decode_categories(level1_latent, n_con1, n_levels1, cat1)
      if (size(n_levels2) > 0) call decode_categories(level2_latent, n_con2, n_levels2, cat2)
      design = 0.0_dp
      column = 0
      if (spec%intercept) then
         column = 1
         design(:, column) = 1.0_dp
      end if
      if (.not. allocated(spec%term)) return
      do term_index = 1, size(spec%term)
         width = term_width(spec%term(term_index), n_levels1, n_levels2)
         do combo = 0, width - 1
            column = column + 1
            do row = 1, n
               design(row, column) = term_value(spec%term(term_index), combo, row, cluster(row), &
                  level1_latent, n_con1, n_levels1, cat1, level2_latent, n_con2, n_levels2, cat2)
            end do
         end do
      end do
   end subroutine smc_build_design

   subroutine smc_mark_level1_predictors(spec, random_spec, n_con, n_levels, used)
      type(smc_design_spec), intent(in) :: spec !! Fixed-effect substantive design specification.
      type(smc_design_spec), intent(in) :: random_spec !! Substantive random-effect design specification.
      integer, intent(in) :: n_con !! Number of leading continuous level-1 latent covariates.
      integer, intent(in) :: n_levels(:) !! Category counts for level-1 categorical covariates.
      logical, intent(out) :: used(:) !! True for latent coordinates entering either substantive design.

      if (size(used) /= latent_width(n_con, n_levels)) error stop "smc_mark_level1_predictors: output size mismatch"
      used = .false.
      call mark_predictors_from_spec(spec, smc_l1_continuous, smc_l1_categorical, n_con, n_levels, used)
      call mark_predictors_from_spec(random_spec, smc_l1_continuous, smc_l1_categorical, n_con, n_levels, used)
   end subroutine smc_mark_level1_predictors

   subroutine smc_mark_level2_predictors(spec, n_con, n_levels, used)
      type(smc_design_spec), intent(in) :: spec !! Fixed-effect substantive design specification; may contain level-2 factors.
      integer, intent(in) :: n_con !! Number of leading continuous level-2 latent covariates.
      integer, intent(in) :: n_levels(:) !! Category counts for level-2 categorical covariates.
      logical, intent(out) :: used(:) !! True for cluster-level latent coordinates entering the substantive design.

      if (size(used) /= latent_width(n_con, n_levels)) error stop "smc_mark_level2_predictors: output size mismatch"
      used = .false.
      call mark_predictors_from_spec(spec, smc_l2_continuous, smc_l2_categorical, n_con, n_levels, used)
   end subroutine smc_mark_level2_predictors

   function smc_substantive_loglik(model, fixed_design, random_design, cluster) result(value)
      type(smc_substantive_model), intent(in) :: model !! Current substantive response, coefficients, and optional random effects.
      real(dp), intent(in) :: fixed_design(:, :) !! Expanded substantive fixed-effect design, shape n by p.
      real(dp), intent(in) :: random_design(:, :) !! Expanded substantive random-effect design, shape n by r; may have zero columns.
      integer, intent(in) :: cluster(:) !! One-based cluster labels corresponding to design rows.
      real(dp) :: value
      integer :: i
      integer :: k
      real(dp) :: eta
      real(dp) :: probability
      real(dp), allocatable :: linear_predictor(:)

      if (.not. allocated(model%beta)) error stop "smc_substantive_loglik: beta is not allocated"
      if (size(fixed_design, 2) /= size(model%beta)) error stop "smc_substantive_loglik: beta/design mismatch"
      if (size(fixed_design, 1) /= size(cluster) .or. size(random_design, 1) /= size(cluster)) &
         error stop "smc_substantive_loglik: row mismatch"
      allocate(linear_predictor(size(cluster)))
      linear_predictor = matmul(fixed_design, model%beta)
      if (size(random_design, 2) > 0) then
         if (.not. allocated(model%random_effects)) error stop "smc_substantive_loglik: random effects are not allocated"
         if (size(model%random_effects, 2) /= size(random_design, 2)) error stop "smc_substantive_loglik: random design mismatch"
         if (any(cluster < 1) .or. any(cluster > size(model%random_effects, 1))) &
            error stop "smc_substantive_loglik: invalid random-effect cluster"
         do i = 1, size(cluster)
            linear_predictor(i) = linear_predictor(i) + dot_product(random_design(i, :), model%random_effects(cluster(i), :))
         end do
      end if

      select case (model%family)
      case (smc_linear)
         if (.not. allocated(model%response)) error stop "smc_substantive_loglik: linear response is not allocated"
         if (size(model%response) /= size(cluster)) error stop "smc_substantive_loglik: response length mismatch"
         if (model%variance <= 0.0_dp) then
            value = -huge(1.0_dp)
         else
            value = -0.5_dp * sum((model%response - linear_predictor)**2) / model%variance
            value = value - 0.5_dp * real(size(cluster), dp) * log(model%variance)
         end if
      case (smc_binary_probit)
         if (.not. allocated(model%category)) error stop "smc_substantive_loglik: binary category is not allocated"
         if (size(model%category) /= size(cluster)) error stop "smc_substantive_loglik: category length mismatch"
         value = 0.0_dp
         do i = 1, size(cluster)
            if (model%category(i) == 1) then
               probability = normal_cdf(-linear_predictor(i))
            else if (model%category(i) == 2) then
               probability = normal_cdf(linear_predictor(i))
            else
               value = -huge(1.0_dp)
               return
            end if
            value = value + log(max(probability, tiny(1.0_dp)))
         end do
      case (smc_cox)
         if (size(random_design, 2) > 0) error stop "smc_substantive_loglik: Cox random effects are not supported upstream"
         if (.not. allocated(model%event)) error stop "smc_substantive_loglik: Cox event vector is not allocated"
         if (size(model%event) /= size(cluster)) error stop "smc_substantive_loglik: Cox event length mismatch"
         value = cox_partial_loglik_ordered(model%event, fixed_design, model%beta)
      case (smc_ordinal_probit)
         if (.not. allocated(model%category)) error stop "smc_substantive_loglik: ordinal category is not allocated"
         if (.not. allocated(model%thresholds)) error stop "smc_substantive_loglik: ordinal thresholds are not allocated"
         if (size(model%category) /= size(cluster)) error stop "smc_substantive_loglik: category length mismatch"
         if (size(model%thresholds) < 1) error stop "smc_substantive_loglik: ordinal model needs thresholds"
         value = 0.0_dp
         k = size(model%thresholds) + 1
         do i = 1, size(cluster)
            if (model%category(i) < 1 .or. model%category(i) > k) then
               value = -huge(1.0_dp)
               return
            end if
            eta = linear_predictor(i)
            if (model%category(i) == 1) then
               probability = normal_cdf(model%thresholds(1) - eta)
            else if (model%category(i) == k) then
               probability = 1.0_dp - normal_cdf(model%thresholds(k - 1) - eta)
            else
               probability = normal_cdf(model%thresholds(model%category(i)) - eta) - &
                  normal_cdf(model%thresholds(model%category(i) - 1) - eta)
            end if
            value = value + log(max(probability, tiny(1.0_dp)))
         end do
      case default
         error stop "smc_substantive_loglik: unknown family"
      end select
   end function smc_substantive_loglik

   subroutine smc_level1_sweep(rng, level1_latent, missing, joint_mean, covariance, n_con1, n_levels1, cluster, &
      level2_latent, n_con2, n_levels2, fixed_spec, random_spec, model, stats, covariance_by_cluster)
      type(rng_state), intent(inout) :: rng !! Mutable generator state for Metropolis and conditional-Gaussian proposals.
      real(dp), intent(inout) :: level1_latent(:, :) !! Complete level-1 latent covariates updated in place.
      logical, intent(in) :: missing(:, :) !! True exactly for originally missing level-1 latent coordinates.
      real(dp), intent(in) :: joint_mean(:, :) !! Current level-1 joint-model mean matrix, same shape as level1_latent.
      real(dp), intent(in) :: covariance(:, :) !! Common level-1 covariance; superseded by supplied cluster covariances.
      integer, intent(in) :: n_con1 !! Number of leading continuous level-1 covariates.
      integer, intent(in) :: n_levels1(:) !! Category counts for level-1 categorical covariates.
      integer, intent(in) :: cluster(:) !! One-based cluster labels; all ones are valid for single-level data.
      real(dp), intent(in) :: level2_latent(:, :) !! Complete cluster-level latent covariates used to rebuild terms.
      integer, intent(in) :: n_con2 !! Number of leading continuous cluster-level covariates.
      integer, intent(in) :: n_levels2(:) !! Category counts for cluster-level categorical covariates.
      type(smc_design_spec), intent(in) :: fixed_spec !! Fixed-effect polynomial/interaction/dummy specification.
      type(smc_design_spec), intent(in) :: random_spec !! Random-effect design specification; may contain no terms and no intercept.
      type(smc_substantive_model), intent(in) :: model !! Current substantive-model state used in MH acceptance ratios.
      type(smc_sweep_stats), intent(inout) :: stats !! Proposal, acceptance, and auxiliary Gibbs counters updated cumulatively.
      real(dp), intent(in), optional :: covariance_by_cluster(:, :, :) !! Optional p by p by G cluster-specific level-1 covariances.
      logical, allocatable :: used(:)
      real(dp), allocatable :: fixed_design(:, :)
      real(dp), allocatable :: random_design(:, :)
      real(dp), allocatable :: fixed_proposed(:, :)
      real(dp), allocatable :: random_proposed(:, :)
      real(dp), allocatable :: precision(:, :)
      real(dp), allocatable :: precision_cluster(:, :, :)
      real(dp), allocatable :: residual(:)
      real(dp), allocatable :: residual_proposed(:)
      real(dp) :: current_ll
      real(dp) :: proposed_ll
      real(dp) :: current_joint
      real(dp) :: proposed_joint
      real(dp) :: proposal
      real(dp) :: proposal_sd
      real(dp) :: old_value
      integer :: i
      integer :: k
      integer :: g
      integer :: info

      if (any(shape(level1_latent) /= shape(missing)) .or. any(shape(level1_latent) /= shape(joint_mean))) &
         error stop "smc_level1_sweep: level-1 shape mismatch"
      if (size(cluster) /= size(level1_latent, 1)) error stop "smc_level1_sweep: cluster length mismatch"
      if (size(covariance, 1) /= size(level1_latent, 2) .or. size(covariance, 2) /= size(level1_latent, 2)) &
         error stop "smc_level1_sweep: covariance shape mismatch"
      if (present(covariance_by_cluster)) then
         if (size(covariance_by_cluster, 1) /= size(level1_latent, 2) .or. &
            size(covariance_by_cluster, 2) /= size(level1_latent, 2)) &
            error stop "smc_level1_sweep: cluster covariance shape mismatch"
         if (size(covariance_by_cluster, 3) < maxval(cluster)) error stop "smc_level1_sweep: missing cluster covariance"
      end if
      allocate(used(size(level1_latent, 2)))
      call smc_mark_level1_predictors(fixed_spec, random_spec, n_con1, n_levels1, used)
      allocate(fixed_design(size(cluster), smc_design_columns(fixed_spec, n_levels1, n_levels2)))
      allocate(random_design(size(cluster), smc_design_columns(random_spec, n_levels1, n_levels2)))
      allocate(fixed_proposed(size(fixed_design, 1), size(fixed_design, 2)))
      allocate(random_proposed(size(random_design, 1), size(random_design, 2)))
      call smc_build_design(level1_latent, n_con1, n_levels1, cluster, level2_latent, n_con2, n_levels2, &
         fixed_spec, fixed_design)
      call smc_build_design(level1_latent, n_con1, n_levels1, cluster, level2_latent, n_con2, n_levels2, &
         random_spec, random_design)
      current_ll = smc_substantive_loglik(model, fixed_design, random_design, cluster)
      allocate(residual(size(level1_latent, 2)), residual_proposed(size(level1_latent, 2)))
      if (present(covariance_by_cluster)) then
         allocate(precision_cluster(size(level1_latent, 2), size(level1_latent, 2), size(covariance_by_cluster, 3)))
         do g = 1, size(covariance_by_cluster, 3)
            call inverse_spd(covariance_by_cluster(:, :, g), precision_cluster(:, :, g), info)
            if (info /= 0) error stop "smc_level1_sweep: cluster covariance is not positive definite"
         end do
      else
         allocate(precision(size(covariance, 1), size(covariance, 2)))
         call inverse_spd(covariance, precision, info)
         if (info /= 0) error stop "smc_level1_sweep: covariance is not positive definite"
      end if

      do i = 1, size(level1_latent, 1)
         do k = 1, size(level1_latent, 2)
            if (.not. missing(i, k)) cycle
            if (.not. used(k)) then
               if (present(covariance_by_cluster)) then
                  call conditional_coordinate_draw(rng, level1_latent(i, :), joint_mean(i, :), &
                     covariance_by_cluster(:, :, cluster(i)), k, proposal, info)
               else
                  call conditional_coordinate_draw(rng, level1_latent(i, :), joint_mean(i, :), covariance, k, proposal, info)
               end if
               if (info /= 0) error stop "smc_level1_sweep: conditional auxiliary draw failed"
               level1_latent(i, k) = proposal
               stats%level1_gibbs = stats%level1_gibbs + 1
               cycle
            end if

            residual = level1_latent(i, :) - joint_mean(i, :)
            if (present(covariance_by_cluster)) then
               current_joint = -0.5_dp * quadratic_form(residual, precision_cluster(:, :, cluster(i)))
               proposal_sd = sqrt(covariance_by_cluster(k, k, cluster(i)) / 10.0_dp)
            else
               current_joint = -0.5_dp * quadratic_form(residual, precision)
               proposal_sd = sqrt(covariance(k, k) / 10.0_dp)
            end if
            if (proposal_sd <= 0.0_dp) error stop "smc_level1_sweep: nonpositive proposal variance"
            old_value = level1_latent(i, k)
            proposal = rng_normal(rng, old_value, proposal_sd)
            level1_latent(i, k) = proposal
            residual_proposed = level1_latent(i, :) - joint_mean(i, :)
            if (present(covariance_by_cluster)) then
               proposed_joint = -0.5_dp * quadratic_form(residual_proposed, precision_cluster(:, :, cluster(i)))
            else
               proposed_joint = -0.5_dp * quadratic_form(residual_proposed, precision)
            end if
            call smc_build_design(level1_latent, n_con1, n_levels1, cluster, level2_latent, n_con2, n_levels2, &
               fixed_spec, fixed_proposed)
            call smc_build_design(level1_latent, n_con1, n_levels1, cluster, level2_latent, n_con2, n_levels2, &
               random_spec, random_proposed)
            proposed_ll = smc_substantive_loglik(model, fixed_proposed, random_proposed, cluster)
            stats%level1_proposals = stats%level1_proposals + 1
            if (log(max(rng_uniform(rng), tiny(1.0_dp))) < &
               min(0.0_dp, proposed_ll + proposed_joint - current_ll - current_joint)) then
               fixed_design = fixed_proposed
               random_design = random_proposed
               current_ll = proposed_ll
               stats%level1_accepted = stats%level1_accepted + 1
            else
               level1_latent(i, k) = old_value
            end if
         end do
      end do
   end subroutine smc_level1_sweep

   subroutine smc_level2_sweep(rng, level1_latent, n_con1, n_levels1, cluster, level2_latent, missing, joint_mean, &
      covariance, n_con2, n_levels2, fixed_spec, random_spec, model, stats)
      type(rng_state), intent(inout) :: rng !! Mutable generator state for cluster-level SMC and conditional-Gaussian draws.
      real(dp), intent(in) :: level1_latent(:, :) !! Current complete level-1 latent covariates used to rebuild terms.
      integer, intent(in) :: n_con1 !! Number of leading continuous level-1 covariates.
      integer, intent(in) :: n_levels1(:) !! Category counts for level-1 categorical covariates.
      integer, intent(in) :: cluster(:) !! One-based cluster labels for level-1 observations.
      real(dp), intent(inout) :: level2_latent(:, :) !! Complete cluster-level latent covariates updated in place.
      logical, intent(in) :: missing(:, :) !! True exactly for originally missing cluster-level latent coordinates.
      real(dp), intent(in) :: joint_mean(:, :) !! Conditional joint-model means for cluster-level latent covariates.
      real(dp), intent(in) :: covariance(:, :) !! Common conditional covariance for cluster-level latent covariates.
      integer, intent(in) :: n_con2 !! Number of leading continuous cluster-level covariates.
      integer, intent(in) :: n_levels2(:) !! Category counts for cluster-level categorical covariates.
      type(smc_design_spec), intent(in) :: fixed_spec !! Fixed-effect substantive design specification.
      type(smc_design_spec), intent(in) :: random_spec !! Random-effect substantive design specification.
      type(smc_substantive_model), intent(in) :: model !! Current substantive-model state used in MH acceptance ratios.
      type(smc_sweep_stats), intent(inout) :: stats !! Proposal, acceptance, and auxiliary Gibbs counters updated cumulatively.
      logical, allocatable :: used(:)
      real(dp), allocatable :: fixed_design(:, :)
      real(dp), allocatable :: random_design(:, :)
      real(dp), allocatable :: fixed_proposed(:, :)
      real(dp), allocatable :: random_proposed(:, :)
      real(dp), allocatable :: precision(:, :)
      real(dp), allocatable :: residual(:)
      real(dp), allocatable :: residual_proposed(:)
      real(dp) :: current_ll
      real(dp) :: proposed_ll
      real(dp) :: current_joint
      real(dp) :: proposed_joint
      real(dp) :: old_value
      real(dp) :: proposal
      real(dp) :: proposal_sd
      integer :: g
      integer :: k
      integer :: info

      if (any(shape(level2_latent) /= shape(missing)) .or. any(shape(level2_latent) /= shape(joint_mean))) &
         error stop "smc_level2_sweep: level-2 shape mismatch"
      if (size(covariance, 1) /= size(level2_latent, 2) .or. size(covariance, 2) /= size(level2_latent, 2)) &
         error stop "smc_level2_sweep: covariance shape mismatch"
      if (size(level2_latent, 2) == 0) return
      allocate(used(size(level2_latent, 2)))
      call smc_mark_level2_predictors(fixed_spec, n_con2, n_levels2, used)
      allocate(fixed_design(size(cluster), smc_design_columns(fixed_spec, n_levels1, n_levels2)))
      allocate(random_design(size(cluster), smc_design_columns(random_spec, n_levels1, n_levels2)))
      allocate(fixed_proposed(size(fixed_design, 1), size(fixed_design, 2)))
      allocate(random_proposed(size(random_design, 1), size(random_design, 2)))
      call smc_build_design(level1_latent, n_con1, n_levels1, cluster, level2_latent, n_con2, n_levels2, &
         fixed_spec, fixed_design)
      call smc_build_design(level1_latent, n_con1, n_levels1, cluster, level2_latent, n_con2, n_levels2, &
         random_spec, random_design)
      current_ll = smc_substantive_loglik(model, fixed_design, random_design, cluster)
      allocate(precision(size(covariance, 1), size(covariance, 2)))
      call inverse_spd(covariance, precision, info)
      if (info /= 0) error stop "smc_level2_sweep: covariance is not positive definite"
      allocate(residual(size(level2_latent, 2)), residual_proposed(size(level2_latent, 2)))

      do g = 1, size(level2_latent, 1)
         do k = 1, size(level2_latent, 2)
            if (.not. missing(g, k)) cycle
            if (.not. used(k)) then
               call conditional_coordinate_draw(rng, level2_latent(g, :), joint_mean(g, :), covariance, k, proposal, info)
               if (info /= 0) error stop "smc_level2_sweep: conditional auxiliary draw failed"
               level2_latent(g, k) = proposal
               stats%level2_gibbs = stats%level2_gibbs + 1
               cycle
            end if
            residual = level2_latent(g, :) - joint_mean(g, :)
            current_joint = -0.5_dp * quadratic_form(residual, precision)
            proposal_sd = sqrt(covariance(k, k) / 10.0_dp)
            if (proposal_sd <= 0.0_dp) error stop "smc_level2_sweep: nonpositive proposal variance"
            old_value = level2_latent(g, k)
            proposal = rng_normal(rng, old_value, proposal_sd)
            level2_latent(g, k) = proposal
            residual_proposed = level2_latent(g, :) - joint_mean(g, :)
            proposed_joint = -0.5_dp * quadratic_form(residual_proposed, precision)
            call smc_build_design(level1_latent, n_con1, n_levels1, cluster, level2_latent, n_con2, n_levels2, &
               fixed_spec, fixed_proposed)
            call smc_build_design(level1_latent, n_con1, n_levels1, cluster, level2_latent, n_con2, n_levels2, &
               random_spec, random_proposed)
            proposed_ll = smc_substantive_loglik(model, fixed_proposed, random_proposed, cluster)
            stats%level2_proposals = stats%level2_proposals + 1
            if (log(max(rng_uniform(rng), tiny(1.0_dp))) < &
               min(0.0_dp, proposed_ll + proposed_joint - current_ll - current_joint)) then
               fixed_design = fixed_proposed
               random_design = random_proposed
               current_ll = proposed_ll
               stats%level2_accepted = stats%level2_accepted + 1
            else
               level2_latent(g, k) = old_value
            end if
         end do
      end do
   end subroutine smc_level2_sweep

   subroutine smc_update_substantive(rng, fixed_design, random_design, cluster, model, random_prior_scale, variance_prior_scale)
      type(rng_state), intent(inout) :: rng !! Mutable generator state for substantive fixed/random-effect and variance updates.
      real(dp), intent(in) :: fixed_design(:, :) !! Current expanded substantive fixed-effect design, shape n by p.
      real(dp), intent(in) :: random_design(:, :) !! Current substantive random-effect design, shape n by r; may have zero columns.
      integer, intent(in) :: cluster(:) !! One-based cluster labels for design rows.
      type(smc_substantive_model), intent(inout) :: model !! Substantive parameters and latent outcomes updated in place.
      real(dp), intent(in), optional :: random_prior_scale(:, :) !! Inverse-Wishart scale for substantive random effects.
      real(dp), intent(in), optional :: variance_prior_scale !! Nonnegative Gaussian residual-variance prior scale.
      real(dp), allocatable :: target(:)
      real(dp), allocatable :: linear_predictor(:)
      real(dp), allocatable :: residual(:)
      real(dp) :: prior_variance

      if (size(fixed_design, 1) /= size(cluster) .or. size(random_design, 1) /= size(cluster)) &
         error stop "smc_update_substantive: row mismatch"
      if (.not. allocated(model%beta)) error stop "smc_update_substantive: beta is not allocated"
      if (size(model%beta) /= size(fixed_design, 2)) error stop "smc_update_substantive: beta/design mismatch"
      prior_variance = 0.0_dp
      if (present(variance_prior_scale)) prior_variance = variance_prior_scale

      select case (model%family)
      case (smc_linear)
         if (.not. allocated(model%response)) error stop "smc_update_substantive: linear response is not allocated"
         call ensure_random_state(model, random_design, cluster)
         target = model%response
         call update_gaussian_mixed_state(rng, target, fixed_design, random_design, cluster, model, &
            max(model%variance, tiny(1.0_dp)), random_prior_scale)
         linear_predictor = substantive_linear_predictor(model, fixed_design, random_design, cluster)
         residual = model%response - linear_predictor
         call sample_linear_variance(rng, residual, prior_variance, model%variance)
         call impute_missing_linear_response(rng, linear_predictor, model)
      case (smc_binary_probit)
         if (.not. allocated(model%category)) error stop "smc_update_substantive: binary category is not allocated"
         call ensure_random_state(model, random_design, cluster)
         call update_binary_latent(rng, fixed_design, random_design, cluster, model)
         call update_gaussian_mixed_state(rng, model%latent_response, fixed_design, random_design, cluster, model, &
            1.0_dp, random_prior_scale)
         model%variance = 1.0_dp
      case (smc_cox)
         if (size(random_design, 2) > 0) error stop "smc_update_substantive: upstream Cox SMC has no random effects"
         if (.not. allocated(model%event)) error stop "smc_update_substantive: Cox event vector is not allocated"
         call cox_coordinate_newton(model%event, fixed_design, model%beta)
      case (smc_ordinal_probit)
         if (.not. allocated(model%category)) error stop "smc_update_substantive: ordinal category is not allocated"
         if (.not. allocated(model%thresholds)) error stop "smc_update_substantive: ordinal thresholds are not allocated"
         call ensure_random_state(model, random_design, cluster)
         call update_ordinal_latent(rng, fixed_design, random_design, cluster, model)
         call update_gaussian_mixed_state(rng, model%latent_response, fixed_design, random_design, cluster, model, &
            1.0_dp, random_prior_scale)
         call update_ordinal_thresholds(rng, model%category, model%latent_response, model%thresholds)
         model%variance = 1.0_dp
      case default
         error stop "smc_update_substantive: unknown family"
      end select
   end subroutine smc_update_substantive

   subroutine smc_compatible_iteration(rng, level1_latent, level1_missing, level1_mean, level1_covariance, n_con1, &
      n_levels1, cluster, level2_latent, level2_missing, level2_mean, level2_covariance, n_con2, n_levels2, &
      fixed_spec, random_spec, model, stats, level1_covariance_by_cluster, random_prior_scale, variance_prior_scale)
      type(rng_state), intent(inout) :: rng !! Mutable generator state for the complete SMC-compatible iteration.
      real(dp), intent(inout) :: level1_latent(:, :) !! Current complete level-1 joint-model latent covariates.
      logical, intent(in) :: level1_missing(:, :) !! Original level-1 missingness mask in expanded latent coordinates.
      real(dp), intent(in) :: level1_mean(:, :) !! Current joint-model mean for level-1 latent covariates.
      real(dp), intent(in) :: level1_covariance(:, :) !! Common level-1 covariance for homogeneous models.
      integer, intent(in) :: n_con1 !! Number of leading continuous level-1 covariates.
      integer, intent(in) :: n_levels1(:) !! Category counts for level-1 categorical covariates.
      integer, intent(in) :: cluster(:) !! One-based cluster labels; all ones are valid for single-level data.
      real(dp), intent(inout) :: level2_latent(:, :) !! Current complete level-2 latent covariates; may have zero columns.
      logical, intent(in) :: level2_missing(:, :) !! Original level-2 missingness mask, same shape as level2_latent.
      real(dp), intent(in) :: level2_mean(:, :) !! Current conditional joint-model means for level-2 latent covariates.
      real(dp), intent(in) :: level2_covariance(:, :) !! Conditional covariance for level-2 latent covariates; may be 0 by 0.
      integer, intent(in) :: n_con2 !! Number of leading continuous level-2 covariates.
      integer, intent(in) :: n_levels2(:) !! Category counts for level-2 categorical covariates.
      type(smc_design_spec), intent(in) :: fixed_spec !! Substantive fixed-effect polynomial/interaction specification.
      type(smc_design_spec), intent(in) :: random_spec !! Substantive random-effect specification.
      type(smc_substantive_model), intent(inout) :: model !! Substantive-model state updated after compatible imputations.
      type(smc_sweep_stats), intent(inout) :: stats !! Cumulative SMC acceptance and auxiliary-Gibbs counters.
      real(dp), intent(in), optional :: level1_covariance_by_cluster(:, :, :) !! Cluster-specific level-1 covariances for HR models.
      real(dp), intent(in), optional :: random_prior_scale(:, :) !! Inverse-Wishart scale for substantive random effects.
      real(dp), intent(in), optional :: variance_prior_scale !! Gaussian substantive residual-variance prior scale.
      real(dp), allocatable :: fixed_design(:, :)
      real(dp), allocatable :: random_design(:, :)

      call smc_level1_sweep(rng, level1_latent, level1_missing, level1_mean, level1_covariance, n_con1, n_levels1, &
         cluster, level2_latent, n_con2, n_levels2, fixed_spec, random_spec, model, stats, level1_covariance_by_cluster)
      if (size(level2_latent, 2) > 0) then
         call smc_level2_sweep(rng, level1_latent, n_con1, n_levels1, cluster, level2_latent, level2_missing, &
            level2_mean, level2_covariance, n_con2, n_levels2, fixed_spec, random_spec, model, stats)
      end if
      allocate(fixed_design(size(cluster), smc_design_columns(fixed_spec, n_levels1, n_levels2)))
      allocate(random_design(size(cluster), smc_design_columns(random_spec, n_levels1, n_levels2)))
      call smc_build_design(level1_latent, n_con1, n_levels1, cluster, level2_latent, n_con2, n_levels2, &
         fixed_spec, fixed_design)
      call smc_build_design(level1_latent, n_con1, n_levels1, cluster, level2_latent, n_con2, n_levels2, &
         random_spec, random_design)
      call smc_update_substantive(rng, fixed_design, random_design, cluster, model, random_prior_scale, variance_prior_scale)
   end subroutine smc_compatible_iteration

   pure function smc_level1_acceptance(self) result(rate)
      class(smc_sweep_stats), intent(in) :: self !! Sweep counters whose level-1 Metropolis acceptance fraction is requested.
      real(dp) :: rate

      if (self%level1_proposals > 0) then
         rate = real(self%level1_accepted, dp) / real(self%level1_proposals, dp)
      else
         rate = 1.0_dp
      end if
   end function smc_level1_acceptance

   pure function smc_level2_acceptance(self) result(rate)
      class(smc_sweep_stats), intent(in) :: self !! Sweep counters whose level-2 Metropolis acceptance fraction is requested.
      real(dp) :: rate

      if (self%level2_proposals > 0) then
         rate = real(self%level2_accepted, dp) / real(self%level2_proposals, dp)
      else
         rate = 1.0_dp
      end if
   end function smc_level2_acceptance

   pure integer function latent_width(n_con, n_levels) result(width)
      integer, intent(in) :: n_con !! Number of continuous variables preceding categorical latent blocks.
      integer, intent(in) :: n_levels(:) !! Category counts for categorical variables.

      if (n_con < 0 .or. any(n_levels < 2)) error stop "latent_width: invalid dimensions"
      width = n_con + sum(n_levels - 1)
   end function latent_width

   pure integer function factor_width(factor, n_levels1, n_levels2) result(width)
      type(smc_factor_spec), intent(in) :: factor !! One continuous/polynomial or categorical design factor.
      integer, intent(in) :: n_levels1(:) !! Category counts for level-1 categorical covariates.
      integer, intent(in) :: n_levels2(:) !! Category counts for level-2 categorical covariates.

      if (factor%variable <= 0) error stop "factor_width: variable indices are one-based"
      select case (factor%source)
      case (smc_l1_continuous, smc_l2_continuous)
         if (factor%power <= 0) error stop "factor_width: continuous powers must be positive"
         width = 1
      case (smc_l1_categorical)
         if (factor%variable > size(n_levels1)) error stop "factor_width: level-1 categorical variable out of range"
         width = n_levels1(factor%variable) - 1
      case (smc_l2_categorical)
         if (factor%variable > size(n_levels2)) error stop "factor_width: level-2 categorical variable out of range"
         width = n_levels2(factor%variable) - 1
      case default
         error stop "factor_width: unknown factor source"
      end select
   end function factor_width

   pure integer function term_width(term, n_levels1, n_levels2) result(width)
      type(smc_term_spec), intent(in) :: term !! One possibly interacting substantive design term.
      integer, intent(in) :: n_levels1(:) !! Category counts for level-1 categorical covariates.
      integer, intent(in) :: n_levels2(:) !! Category counts for level-2 categorical covariates.
      integer :: j

      if (.not. allocated(term%factor)) error stop "term_width: factor array is not allocated"
      if (size(term%factor) == 0) error stop "term_width: empty term"
      width = 1
      do j = 1, size(term%factor)
         width = width * factor_width(term%factor(j), n_levels1, n_levels2)
      end do
   end function term_width

   pure function term_value(term, combo, row, group, level1_latent, n_con1, n_levels1, cat1, &
      level2_latent, n_con2, n_levels2, cat2) result(value)
      type(smc_term_spec), intent(in) :: term !! One polynomial/interaction term to evaluate.
      integer, intent(in) :: combo !! Zero-based dummy-column combination within the expanded term.
      integer, intent(in) :: row !! One-based level-1 observation row.
      integer, intent(in) :: group !! One-based cluster row used for level-2 factors.
      real(dp), intent(in) :: level1_latent(:, :) !! Complete level-1 latent covariates.
      integer, intent(in) :: n_con1 !! Number of leading continuous level-1 covariates.
      integer, intent(in) :: n_levels1(:) !! Category counts for level-1 categorical covariates.
      integer, intent(in) :: cat1(:, :) !! Decoded level-1 categories corresponding to level1_latent.
      real(dp), intent(in) :: level2_latent(:, :) !! Complete cluster-level latent covariates.
      integer, intent(in) :: n_con2 !! Number of leading continuous level-2 covariates.
      integer, intent(in) :: n_levels2(:) !! Category counts for level-2 categorical covariates.
      integer, intent(in) :: cat2(:, :) !! Decoded level-2 categories corresponding to level2_latent.
      real(dp) :: value
      integer :: j
      integer :: stride
      integer :: width
      integer :: level

      value = 1.0_dp
      stride = 1
      do j = 1, size(term%factor)
         width = factor_width(term%factor(j), n_levels1, n_levels2)
         select case (term%factor(j)%source)
         case (smc_l1_continuous)
            if (term%factor(j)%variable > n_con1) error stop "term_value: level-1 continuous variable out of range"
            value = value * level1_latent(row, term%factor(j)%variable)**term%factor(j)%power
         case (smc_l1_categorical)
            level = modulo(combo / stride, width) + 1
            if (cat1(row, term%factor(j)%variable) /= level) value = 0.0_dp
         case (smc_l2_continuous)
            if (term%factor(j)%variable > n_con2) error stop "term_value: level-2 continuous variable out of range"
            value = value * level2_latent(group, term%factor(j)%variable)**term%factor(j)%power
         case (smc_l2_categorical)
            level = modulo(combo / stride, width) + 1
            if (cat2(group, term%factor(j)%variable) /= level) value = 0.0_dp
         case default
            error stop "term_value: unknown factor source"
         end select
         stride = stride * width
      end do
   end function term_value

   subroutine mark_predictors_from_spec(spec, continuous_source, categorical_source, n_con, n_levels, used)
      type(smc_design_spec), intent(in) :: spec !! Design specification searched for factors at one data level.
      integer, intent(in) :: continuous_source !! Source code identifying continuous factors at the requested level.
      integer, intent(in) :: categorical_source !! Source code identifying categorical factors at the requested level.
      integer, intent(in) :: n_con !! Number of continuous latent coordinates at the requested level.
      integer, intent(in) :: n_levels(:) !! Category counts for categorical variables at the requested level.
      logical, intent(inout) :: used(:) !! Latent-coordinate usage mask updated in place.
      integer :: j
      integer :: k
      integer :: first
      integer :: last

      if (.not. allocated(spec%term)) return
      do j = 1, size(spec%term)
         if (.not. allocated(spec%term(j)%factor)) cycle
         do k = 1, size(spec%term(j)%factor)
            if (spec%term(j)%factor(k)%source == continuous_source) then
               if (spec%term(j)%factor(k)%variable > n_con) error stop "mark_predictors_from_spec: continuous variable out of range"
               used(spec%term(j)%factor(k)%variable) = .true.
            else if (spec%term(j)%factor(k)%source == categorical_source) then
               call categorical_block(spec%term(j)%factor(k)%variable, n_con, n_levels, first, last)
               used(first:last) = .true.
            end if
         end do
      end do
   end subroutine mark_predictors_from_spec

   pure subroutine categorical_block(variable, n_con, n_levels, first, last)
      integer, intent(in) :: variable !! One-based categorical-variable index whose latent block is requested.
      integer, intent(in) :: n_con !! Number of leading continuous latent coordinates.
      integer, intent(in) :: n_levels(:) !! Category counts for categorical variables.
      integer, intent(out) :: first !! First one-based latent coordinate of the requested categorical block.
      integer, intent(out) :: last !! Last one-based latent coordinate of the requested categorical block.

      if (variable < 1 .or. variable > size(n_levels)) error stop "categorical_block: variable out of range"
      first = n_con + 1
      if (variable > 1) first = first + sum(n_levels(:variable - 1) - 1)
      last = first + n_levels(variable) - 2
   end subroutine categorical_block

   subroutine conditional_coordinate_draw(rng, values, mean, covariance, coordinate, draw, info)
      type(rng_state), intent(inout) :: rng !! Mutable generator state for the univariate conditional-normal draw.
      real(dp), intent(in) :: values(:) !! Current complete multivariate row.
      real(dp), intent(in) :: mean(:) !! Current multivariate-normal mean for the row.
      real(dp), intent(in) :: covariance(:, :) !! Positive-definite covariance of the row.
      integer, intent(in) :: coordinate !! One-based coordinate to draw conditional on all other coordinates.
      real(dp), intent(out) :: draw !! Sampled conditional value for the requested coordinate.
      integer, intent(out) :: info !! Zero on success; nonzero if the conditioning covariance cannot be solved.
      integer :: p
      integer :: j
      integer :: pos
      integer, allocatable :: other(:)
      real(dp), allocatable :: soo(:, :)
      real(dp), allocatable :: sko(:)
      real(dp), allocatable :: weight(:)
      real(dp), allocatable :: diff(:)
      real(dp) :: conditional_mean
      real(dp) :: conditional_variance

      p = size(values)
      if (size(mean) /= p .or. size(covariance, 1) /= p .or. size(covariance, 2) /= p) &
         error stop "conditional_coordinate_draw: shape mismatch"
      if (coordinate < 1 .or. coordinate > p) error stop "conditional_coordinate_draw: coordinate out of range"
      if (p == 1) then
         if (covariance(1, 1) <= 0.0_dp) then
            info = 1
            draw = values(1)
            return
         end if
         draw = rng_normal(rng, mean(1), sqrt(covariance(1, 1)))
         info = 0
         return
      end if
      allocate(other(p - 1), soo(p - 1, p - 1), sko(p - 1), weight(p - 1), diff(p - 1))
      pos = 0
      do j = 1, p
         if (j == coordinate) cycle
         pos = pos + 1
         other(pos) = j
      end do
      soo = covariance(other, other)
      sko = covariance(other, coordinate)
      call solve_spd(soo, sko, weight, info)
      if (info /= 0) then
         draw = values(coordinate)
         return
      end if
      diff = values(other) - mean(other)
      conditional_mean = mean(coordinate) + dot_product(weight, diff)
      conditional_variance = covariance(coordinate, coordinate) - dot_product(sko, weight)
      if (conditional_variance <= tiny(1.0_dp)) then
         info = 2
         draw = values(coordinate)
         return
      end if
      draw = rng_normal(rng, conditional_mean, sqrt(conditional_variance))
      info = 0
   end subroutine conditional_coordinate_draw

   subroutine ensure_random_state(model, random_design, cluster)
      type(smc_substantive_model), intent(inout) :: model !! Substantive state receiving default random-effect arrays when needed.
      real(dp), intent(in) :: random_design(:, :) !! Random-effect design whose column count determines the random-effect dimension.
      integer, intent(in) :: cluster(:) !! One-based cluster labels determining the number of random-effect vectors.
      integer :: r
      integer :: g
      integer :: i

      r = size(random_design, 2)
      if (r == 0) return
      g = maxval(cluster)
      if (.not. allocated(model%random_effects)) then
         allocate(model%random_effects(g, r))
         model%random_effects = 0.0_dp
      else if (size(model%random_effects, 1) /= g .or. size(model%random_effects, 2) /= r) then
         error stop "ensure_random_state: random-effect shape mismatch"
      end if
      if (.not. allocated(model%random_covariance)) then
         allocate(model%random_covariance(r, r))
         model%random_covariance = 0.0_dp
         do i = 1, r
            model%random_covariance(i, i) = 1.0_dp
         end do
      else if (size(model%random_covariance, 1) /= r .or. size(model%random_covariance, 2) /= r) then
         error stop "ensure_random_state: random-covariance shape mismatch"
      end if
   end subroutine ensure_random_state

   subroutine update_gaussian_mixed_state(rng, target, fixed_design, random_design, cluster, model, variance, random_prior_scale)
      type(rng_state), intent(inout) :: rng !! Mutable generator state for Gaussian fixed/random-effect posterior draws.
      real(dp), intent(in) :: target(:) !! Continuous or latent-Gaussian substantive response.
      real(dp), intent(in) :: fixed_design(:, :) !! Fixed-effect design matrix, shape n by p.
      real(dp), intent(in) :: random_design(:, :) !! Random-effect design matrix, shape n by r; may have zero columns.
      integer, intent(in) :: cluster(:) !! One-based cluster labels for target rows.
      type(smc_substantive_model), intent(inout) :: model !! Fixed/random coefficients and random covariance updated in place.
      real(dp), intent(in) :: variance !! Positive residual variance, fixed to one for probit models.
      real(dp), intent(in), optional :: random_prior_scale(:, :) !! Inverse-Wishart scale for random-effect covariance.
      real(dp), allocatable :: adjusted(:)
      real(dp), allocatable :: random_mean(:)
      real(dp), allocatable :: precision(:, :)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: rhs(:)
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: covu_precision(:, :)
      real(dp), allocatable :: scale(:, :)
      integer :: i
      integer :: g
      integer :: r
      integer :: info

      if (variance <= 0.0_dp) error stop "update_gaussian_mixed_state: variance must be positive"
      if (size(target) /= size(cluster)) error stop "update_gaussian_mixed_state: target length mismatch"
      call ensure_random_state(model, random_design, cluster)
      allocate(random_mean(size(target)))
      random_mean = 0.0_dp
      if (size(random_design, 2) > 0) then
         do i = 1, size(target)
            random_mean(i) = dot_product(random_design(i, :), model%random_effects(cluster(i), :))
         end do
      end if
      adjusted = target - random_mean
      call sample_gaussian_coefficients(rng, adjusted, fixed_design, variance, model%beta, info)
      if (info /= 0) error stop "update_gaussian_mixed_state: fixed-effect draw failed"
      r = size(random_design, 2)
      if (r == 0) return
      allocate(covu_precision(r, r))
      call inverse_spd(model%random_covariance, covu_precision, info)
      if (info /= 0) error stop "update_gaussian_mixed_state: random covariance is not positive definite"
      allocate(precision(r, r), covariance(r, r), rhs(r), mean(r))
      do g = 1, size(model%random_effects, 1)
         precision = covu_precision
         rhs = 0.0_dp
         do i = 1, size(target)
            if (cluster(i) /= g) cycle
            precision = precision + outer_product(random_design(i, :), random_design(i, :)) / variance
            rhs = rhs + random_design(i, :) * (target(i) - dot_product(fixed_design(i, :), model%beta)) / variance
         end do
         call inverse_spd(precision, covariance, info)
         if (info /= 0) error stop "update_gaussian_mixed_state: random-effect posterior inversion failed"
         mean = matmul(covariance, rhs)
         call mvnormal_sample(rng, mean, covariance, model%random_effects(g, :), info)
         if (info /= 0) error stop "update_gaussian_mixed_state: random-effect draw failed"
      end do
      allocate(scale(r, r))
      scale = 0.0_dp
      if (present(random_prior_scale)) then
         if (any(shape(random_prior_scale) /= shape(scale))) &
            error stop "update_gaussian_mixed_state: random prior scale mismatch"
         scale = random_prior_scale
      else
         do i = 1, r
            scale(i, i) = 1.0_dp
         end do
      end if
      do g = 1, size(model%random_effects, 1)
         scale = scale + outer_product(model%random_effects(g, :), model%random_effects(g, :))
      end do
      call invwishart_sample(rng, real(size(model%random_effects, 1) + r, dp), scale, model%random_covariance, info)
      if (info /= 0) error stop "update_gaussian_mixed_state: random covariance draw failed"
   end subroutine update_gaussian_mixed_state

   subroutine update_binary_latent(rng, fixed_design, random_design, cluster, model)
      type(rng_state), intent(inout) :: rng !! Mutable generator state for binary-probit latent response draws.
      real(dp), intent(in) :: fixed_design(:, :) !! Current fixed-effect substantive design.
      real(dp), intent(in) :: random_design(:, :) !! Current random-effect substantive design.
      integer, intent(in) :: cluster(:) !! One-based cluster labels for substantive rows.
      type(smc_substantive_model), intent(inout) :: model !! Binary categories and latent response updated in place.
      real(dp), allocatable :: eta(:)
      real(dp) :: draw
      integer :: i
      integer :: tries
      logical :: observed

      if (size(model%category) /= size(cluster)) error stop "update_binary_latent: category length mismatch"
      if (.not. allocated(model%latent_response)) allocate(model%latent_response(size(cluster)))
      eta = substantive_linear_predictor(model, fixed_design, random_design, cluster)
      do i = 1, size(cluster)
         observed = .true.
         if (allocated(model%response_observed)) observed = model%response_observed(i)
         if (.not. observed) then
            draw = rng_normal(rng, eta(i), 1.0_dp)
            model%latent_response(i) = draw
            if (draw > 0.0_dp) then
               model%category(i) = 2
            else
               model%category(i) = 1
            end if
            cycle
         end if
         if (model%category(i) /= 1 .and. model%category(i) /= 2) error stop "update_binary_latent: invalid category"
         do tries = 1, 10000
            draw = rng_normal(rng, eta(i), 1.0_dp)
            if (model%category(i) == 1 .and. draw < 0.0_dp) exit
            if (model%category(i) == 2 .and. draw > 0.0_dp) exit
         end do
         if (tries > 10000) error stop "update_binary_latent: rejection sampler failed"
         model%latent_response(i) = draw
      end do
   end subroutine update_binary_latent

   subroutine update_ordinal_latent(rng, fixed_design, random_design, cluster, model)
      type(rng_state), intent(inout) :: rng !! Mutable generator state for ordinal-probit latent response draws.
      real(dp), intent(in) :: fixed_design(:, :) !! Current fixed-effect substantive design.
      real(dp), intent(in) :: random_design(:, :) !! Current random-effect substantive design.
      integer, intent(in) :: cluster(:) !! One-based cluster labels for substantive rows.
      type(smc_substantive_model), intent(inout) :: model !! Ordinal categories, thresholds, and latent response updated in place.
      real(dp), allocatable :: eta(:)
      real(dp) :: draw
      real(dp) :: lower
      real(dp) :: upper
      integer :: i
      integer :: k
      integer :: tries
      logical :: observed

      k = size(model%thresholds) + 1
      if (size(model%category) /= size(cluster)) error stop "update_ordinal_latent: category length mismatch"
      if (.not. allocated(model%latent_response)) allocate(model%latent_response(size(cluster)))
      eta = substantive_linear_predictor(model, fixed_design, random_design, cluster)
      do i = 1, size(cluster)
         observed = .true.
         if (allocated(model%response_observed)) observed = model%response_observed(i)
         if (.not. observed) then
            draw = rng_normal(rng, eta(i), 1.0_dp)
            model%latent_response(i) = draw
            model%category(i) = ordinal_category(draw, model%thresholds)
            cycle
         end if
         if (model%category(i) < 1 .or. model%category(i) > k) error stop "update_ordinal_latent: invalid category"
         lower = -huge(1.0_dp)
         upper = huge(1.0_dp)
         if (model%category(i) > 1) lower = model%thresholds(model%category(i) - 1)
         if (model%category(i) < k) upper = model%thresholds(model%category(i))
         do tries = 1, 10000
            draw = rng_normal(rng, eta(i), 1.0_dp)
            if (draw > lower .and. draw < upper) exit
         end do
         if (tries > 10000) error stop "update_ordinal_latent: rejection sampler failed"
         model%latent_response(i) = draw
      end do
   end subroutine update_ordinal_latent

   subroutine impute_missing_linear_response(rng, eta, model)
      type(rng_state), intent(inout) :: rng !! Mutable generator state for missing Gaussian substantive responses.
      real(dp), intent(in) :: eta(:) !! Current substantive linear predictor for every row.
      type(smc_substantive_model), intent(inout) :: model !! Gaussian response and observation mask updated for missing rows.
      integer :: i

      if (.not. allocated(model%response_observed)) return
      if (size(model%response_observed) /= size(model%response)) error stop "impute_missing_linear_response: mask mismatch"
      do i = 1, size(model%response)
         if (.not. model%response_observed(i)) model%response(i) = rng_normal(rng, eta(i), sqrt(model%variance))
      end do
   end subroutine impute_missing_linear_response

   function substantive_linear_predictor(model, fixed_design, random_design, cluster) result(eta)
      type(smc_substantive_model), intent(in) :: model !! Substantive fixed and optional random coefficients.
      real(dp), intent(in) :: fixed_design(:, :) !! Fixed-effect design matrix.
      real(dp), intent(in) :: random_design(:, :) !! Random-effect design matrix; may have zero columns.
      integer, intent(in) :: cluster(:) !! One-based cluster labels for rows.
      real(dp), allocatable :: eta(:)
      integer :: i

      eta = matmul(fixed_design, model%beta)
      if (size(random_design, 2) > 0) then
         if (.not. allocated(model%random_effects)) error stop "substantive_linear_predictor: random effects missing"
         do i = 1, size(cluster)
            eta(i) = eta(i) + dot_product(random_design(i, :), model%random_effects(cluster(i), :))
         end do
      end if
   end function substantive_linear_predictor

   pure integer function ordinal_category(latent, thresholds) result(category)
      real(dp), intent(in) :: latent !! Latent-normal ordinal response value to categorize.
      real(dp), intent(in) :: thresholds(:) !! Strictly increasing ordinal cut points.
      integer :: j

      category = size(thresholds) + 1
      do j = 1, size(thresholds)
         if (latent <= thresholds(j)) then
            category = j
            return
         end if
      end do
   end function ordinal_category

   pure function outer_product(x, y) result(value)
      real(dp), intent(in) :: x(:) !! Left vector of the outer product.
      real(dp), intent(in) :: y(:) !! Right vector of the outer product.
      real(dp) :: value(size(x), size(y))

      value = spread(x, 2, size(y)) * spread(y, 1, size(x))
   end function outer_product

end module jomo_smc
