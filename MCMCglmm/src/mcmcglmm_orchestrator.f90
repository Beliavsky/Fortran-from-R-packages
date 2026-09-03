! SPDX-License-Identifier: GPL-2.0-or-later
! Numeric orchestration for translated MCMCglmm engines; see NOTICE.md and upstream/.
module mcmcglmm_orchestrator
    use r_kinds, only : dp
    use mcmcglmm_rng, only : rng_state
    use mcmcglmm_unified_sampler, only : unified_family_mcmc_result, heterogeneous_multi_term_mixed_mcmc
    use mcmcglmm_grouped_multiterm_sampler, only : grouped_multi_term_mcmc_result, &
        two_part_multi_term_mixed_mcmc, multinomial_multi_term_mixed_mcmc
    use mcmcglmm_multiterm_px_sampler, only : multi_term_family_px_mcmc_result, &
        heterogeneous_multi_term_parameter_expanded_mcmc
    use mcmcglmm_family_sampler, only : ordinal_native_mcmc_result, threshold_cutpoint_mcmc_result, &
        ordinal_native_mixed_mcmc, threshold_cutpoint_mixed_mcmc
    use mcmcglmm_theta_sampler, only : theta_scale_mcmc_result, theta_scale_gaussian_mixed_mcmc
    use mcmcglmm_structural_sampler, only : structural_gaussian_mcmc_result, &
        structural_gaussian_multi_term_mcmc
    use mcmcglmm_covu_sampler, only : covu_gaussian_mcmc_result, covu_gaussian_mixed_mcmc
    use mcmcglmm_sparse, only : mcmcglmm_sparse_matrix, sparse_is_initialized, sparse_to_dense, sparse_validate
    implicit none
    private

    integer, parameter, public :: mcmcglmm_engine_scalar = 1
    integer, parameter, public :: mcmcglmm_engine_scalar_px = 2
    integer, parameter, public :: mcmcglmm_engine_two_process = 3
    integer, parameter, public :: mcmcglmm_engine_multinomial = 4
    integer, parameter, public :: mcmcglmm_engine_ordinal = 5
    integer, parameter, public :: mcmcglmm_engine_threshold = 6
    integer, parameter, public :: mcmcglmm_engine_theta_scale = 7
    integer, parameter, public :: mcmcglmm_engine_structural = 8
    integer, parameter, public :: mcmcglmm_engine_covu = 9

    integer, parameter, public :: mcmcglmm_orchestrator_invalid_engine = 101
    integer, parameter, public :: mcmcglmm_orchestrator_missing_input = 102
    integer, parameter, public :: mcmcglmm_orchestrator_invalid_control = 103

    type, public :: mcmcglmm_numeric_model
        !! Numeric response and design data supplied after formula/model-matrix construction.
        integer :: engine = 0
        integer :: grouped_family = 0
        integer, allocatable :: family(:)
        integer, allocatable :: grouped_response(:, :)
        integer, allocatable :: trials(:)
        integer, allocatable :: y_category(:)
        integer, allocatable :: random_term(:)
        real(dp), allocatable :: y(:, :)
        real(dp), allocatable :: y_vector(:)
        real(dp), allocatable :: additional(:, :)
        real(dp), allocatable :: additional2(:, :)
        real(dp), allocatable :: x(:, :)
        real(dp), allocatable :: z(:, :)
        type(mcmcglmm_sparse_matrix) :: sparse_x
        type(mcmcglmm_sparse_matrix) :: sparse_z
        real(dp), allocatable :: a_inverse(:, :)
        real(dp), allocatable :: cutpoints(:)
        real(dp), allocatable :: x_scale(:, :)
        real(dp), allocatable :: z_scale(:, :)
        real(dp), allocatable :: structural_basis(:, :, :)
        real(dp), allocatable :: random_loading(:, :)
        real(dp), allocatable :: measurement_category_effect(:, :, :)
        real(dp), allocatable :: measurement_prior_probability(:, :)
        integer, allocatable :: measurement_group(:)
        logical, allocatable :: observed(:, :)
        logical, allocatable :: observed_rows(:)
    end type mcmcglmm_numeric_model

    type, public :: mcmcglmm_numeric_prior
        !! Fixed-effect, random-effect, residual, and optional expansion priors.
        real(dp), allocatable :: beta_mean(:, :)
        real(dp), allocatable :: beta_precision(:, :)
        real(dp), allocatable :: g_scale(:, :, :)
        real(dp), allocatable :: g_df(:)
        real(dp), allocatable :: r_scale(:, :)
        real(dp) :: r_df = 0.0_dp
        real(dp), allocatable :: alpha_mean(:, :)
        real(dp), allocatable :: alpha_precision(:, :)
        real(dp) :: theta_mean = 1.0_dp
        real(dp) :: theta_precision = 0.0_dp
        real(dp), allocatable :: structural_mean(:)
        real(dp), allocatable :: structural_precision(:, :)
        real(dp), allocatable :: structural_proposal_sd(:)
        real(dp), allocatable :: joint_scale(:, :)
        real(dp) :: joint_df = 0.0_dp
    end type mcmcglmm_numeric_prior

    type, public :: mcmcglmm_control
        !! MCMC schedule and covariance-update routing shared by supported engines.
        integer :: iterations = 0
        integer :: burn = 0
        integer :: thin = 1
        real(dp) :: proposal_scale = 0.25_dp
        logical :: update_r = .true.
        logical :: adapt_cutpoints = .true.
        logical :: slice_sampling = .false.
        real(dp) :: cutpoint_proposal_sd = 0.10_dp
        real(dp) :: slice_limit = 1.0e6_dp
        real(dp) :: initial_theta = 1.0_dp
        real(dp), allocatable :: initial_structural(:)
        real(dp), allocatable :: initial_joint_covariance(:, :)
        integer :: joint_update_mode = 1
        real(dp), allocatable :: joint_fixed_block(:, :)
        logical, allocatable :: update_g(:)
        integer, allocatable :: g_update_mode(:)
        integer :: r_update_mode = 1
        integer, allocatable :: g_split(:)
        integer :: r_split = -1
        real(dp), allocatable :: initial_g(:, :, :)
        real(dp), allocatable :: initial_r(:, :)
        real(dp), allocatable :: g_fixed_covariance(:, :, :)
        real(dp), allocatable :: r_fixed_covariance(:, :)
        integer, allocatable :: g_ante_order(:)
        integer :: r_ante_order = 1
        logical, allocatable :: g_ante_common_beta(:)
        logical :: r_ante_common_beta = .false.
        logical, allocatable :: g_ante_common_variance(:)
        logical :: r_ante_common_variance = .false.
    end type mcmcglmm_control

    type, public :: mcmcglmm_numeric_result
        !! Tagged result whose component matching engine contains the retained chain.
        integer :: engine = 0
        type(unified_family_mcmc_result) :: scalar
        type(multi_term_family_px_mcmc_result) :: scalar_px
        type(grouped_multi_term_mcmc_result) :: grouped
        type(ordinal_native_mcmc_result) :: ordinal
        type(threshold_cutpoint_mcmc_result) :: threshold
        type(theta_scale_mcmc_result) :: theta_scale
        type(structural_gaussian_mcmc_result) :: structural
        type(covu_gaussian_mcmc_result) :: covu
    end type mcmcglmm_numeric_result

    public :: mcmcglmm_fit_numeric
    public :: mcmcglmm_validate_numeric

contains

    pure subroutine mcmcglmm_fit_numeric(model, prior, control, state, result, info)
        !! Validate and run one supported numeric MCMCglmm engine through a common entry point.
        type(mcmcglmm_numeric_model), intent(in) :: model !! Numeric response, design, and engine specification.
        type(mcmcglmm_numeric_prior), intent(in) :: prior !! Priors shared by the selected engine.
        type(mcmcglmm_control), intent(in) :: control !! MCMC schedule and covariance-routing controls.
        type(rng_state), intent(inout) :: state !! Explicit random-number state consumed by the selected engine.
        type(mcmcglmm_numeric_result), intent(out) :: result !! Tagged result populated for the selected engine.
        integer, intent(out) :: info !! Zero on success; orchestrator or underlying engine status otherwise.
        integer, allocatable :: g_ante_order(:)
        integer, allocatable :: g_mode(:)
        integer, allocatable :: g_split(:)
        logical, allocatable :: g_ante_common_beta(:)
        logical, allocatable :: g_ante_common_variance(:)
        logical, allocatable :: observed(:, :)
        logical, allocatable :: observed_rows(:)
        logical, allocatable :: update_g(:)
        real(dp), allocatable :: g_fixed(:, :, :)
        real(dp), allocatable :: initial_g(:, :, :)
        real(dp), allocatable :: initial_r(:, :)
        real(dp), allocatable :: initial_structural(:)
        real(dp), allocatable :: initial_joint(:, :)
        real(dp), allocatable :: r_fixed(:, :)
        real(dp), allocatable :: x_design(:, :)
        real(dp), allocatable :: z_design(:, :)
        integer :: latent_dimension
        integer :: nterm
        integer :: r_split
        integer :: joint_dimension
        logical :: use_sparse_scalar

        call mcmcglmm_validate_numeric(model, prior, control, latent_dimension, nterm, info)
        if (info /= 0) return
        use_sparse_scalar = model%engine == mcmcglmm_engine_scalar .and. &
            sparse_is_initialized(model%sparse_x) .and. sparse_is_initialized(model%sparse_z)
        if (use_sparse_scalar) then
            allocate(x_design(0, 0), z_design(0, 0))
        else
            call resolve_numeric_design(model, x_design, z_design, info)
            if (info /= 0) return
        end if

        result%engine = model%engine
        if (model%engine == mcmcglmm_engine_covu) then
            joint_dimension = size(prior%joint_scale, 1)
            if (allocated(control%initial_joint_covariance)) then
                if (any(shape(control%initial_joint_covariance) /= [joint_dimension, joint_dimension])) then
                    info = mcmcglmm_orchestrator_invalid_control
                    return
                end if
                initial_joint = control%initial_joint_covariance
            else
                initial_joint = prior%joint_scale / &
                    max(prior%joint_df - real(joint_dimension + 1, dp), 1.0_dp)
            end if
            if (any(control%joint_update_mode == [2, 4])) then
                if (.not. allocated(control%joint_fixed_block)) then
                    info = mcmcglmm_orchestrator_missing_input
                    return
                end if
                call covu_gaussian_mixed_mcmc(model%y, x_design, model%random_loading, prior%beta_mean, &
                    prior%beta_precision, prior%joint_scale, prior%joint_df, control%iterations, control%burn, &
                    control%thin, state, result%covu, info, initial_joint_covariance=initial_joint, &
                    update_mode=control%joint_update_mode, fixed_block=control%joint_fixed_block)
            else
                call covu_gaussian_mixed_mcmc(model%y, x_design, model%random_loading, prior%beta_mean, &
                    prior%beta_precision, prior%joint_scale, prior%joint_df, control%iterations, control%burn, &
                    control%thin, state, result%covu, info, initial_joint_covariance=initial_joint, &
                    update_mode=control%joint_update_mode)
            end if
            return
        end if
        if (model%engine == mcmcglmm_engine_ordinal) then
            call row_observed_mask(model, observed_rows)
            call ordinal_native_mixed_mcmc(model%y_category, model%cutpoints, control%cutpoint_proposal_sd, &
                control%adapt_cutpoints, x_design, z_design, model%a_inverse, prior%beta_mean(:, 1), &
                prior%beta_precision, prior%g_scale(1, 1, 1), prior%g_df(1), prior%r_scale(1, 1), prior%r_df, &
                control%update_r, control%proposal_scale, control%iterations, control%burn, control%thin, state, &
                result%ordinal, info, observed=observed_rows, slice_sampling=control%slice_sampling, &
                slice_limit=control%slice_limit)
            return
        end if
        if (model%engine == mcmcglmm_engine_threshold) then
            call row_observed_mask(model, observed_rows)
            call threshold_cutpoint_mixed_mcmc(model%y_category, model%cutpoints, control%cutpoint_proposal_sd, &
                control%adapt_cutpoints, x_design, z_design, model%a_inverse, prior%beta_mean(:, 1), &
                prior%beta_precision, prior%g_scale(1, 1, 1), prior%g_df(1), control%iterations, control%burn, &
                control%thin, state, result%threshold, info, observed=observed_rows)
            return
        end if
        if (model%engine == mcmcglmm_engine_scalar_px) then
            call scalar_observed_mask(model, observed)
            call heterogeneous_multi_term_parameter_expanded_mcmc(model%family, model%y, model%additional, &
                model%additional2, x_design, z_design, model%random_term, model%a_inverse, prior%beta_mean, &
                prior%beta_precision, prior%g_scale, prior%g_df, prior%r_scale, prior%r_df, prior%alpha_mean, &
                prior%alpha_precision, control%update_r, control%proposal_scale, control%iterations, control%burn, &
                control%thin, state, result%scalar_px, info, observed)
            return
        end if

        call normalize_covariance_control(prior, control, latent_dimension, nterm, update_g, initial_g, initial_r, &
            g_mode, g_split, g_fixed, r_split, r_fixed, g_ante_order, g_ante_common_beta, &
            g_ante_common_variance, info)
        if (info /= 0) return

        if (model%engine == mcmcglmm_engine_theta_scale) then
            call theta_scale_gaussian_mixed_mcmc(model%y, x_design, z_design, model%x_scale, model%z_scale, &
                model%random_term, model%a_inverse, prior%beta_mean, prior%beta_precision, prior%g_scale, &
                prior%g_df, prior%r_scale, prior%r_df, prior%theta_mean, prior%theta_precision, control%iterations, &
                control%burn, control%thin, state, result%theta_scale, info, update_g=update_g, &
                update_r=control%update_r, initial_g=initial_g, initial_r=initial_r, &
                initial_theta=control%initial_theta)
            return
        end if
        if (model%engine == mcmcglmm_engine_structural) then
            initial_structural = prior%structural_mean
            if (allocated(control%initial_structural)) then
                if (size(control%initial_structural) /= size(prior%structural_mean)) then
                    info = mcmcglmm_orchestrator_invalid_control
                    return
                end if
                initial_structural = control%initial_structural
            end if
            call structural_gaussian_multi_term_mcmc(model%y, x_design, z_design, model%random_term, model%a_inverse, &
                prior%beta_mean, prior%beta_precision, prior%g_scale, prior%g_df, prior%r_scale, prior%r_df, &
                model%structural_basis, prior%structural_mean, prior%structural_precision, &
                prior%structural_proposal_sd, control%iterations, control%burn, control%thin, state, &
                result%structural, info, update_g=update_g, update_r=control%update_r, initial_g=initial_g, &
                initial_r=initial_r, initial_structural=initial_structural)
            return
        end if

        select case (model%engine)
        case (mcmcglmm_engine_scalar)
            call scalar_observed_mask(model, observed)
            if (allocated(model%measurement_category_effect) .and. use_sparse_scalar) then
                call heterogeneous_multi_term_mixed_mcmc(model%family, model%y, model%additional, model%additional2, &
                    x_design, z_design, model%random_term, model%a_inverse, prior%beta_mean, prior%beta_precision, &
                    prior%g_scale, prior%g_df, prior%r_scale, prior%r_df, control%update_r, control%proposal_scale, &
                    control%iterations, control%burn, control%thin, state, result%scalar, info, observed=observed, &
                    update_g=update_g, initial_g=initial_g, initial_r=initial_r, g_update_mode=g_mode, &
                    r_update_mode=control%r_update_mode, g_split=g_split, r_split=r_split, &
                    g_fixed_covariance=g_fixed, r_fixed_covariance=r_fixed, g_ante_order=g_ante_order, &
                    r_ante_order=control%r_ante_order, g_ante_common_beta=g_ante_common_beta, &
                    r_ante_common_beta=control%r_ante_common_beta, g_ante_common_variance=g_ante_common_variance, &
                    r_ante_common_variance=control%r_ante_common_variance, &
                    measurement_category_effect=model%measurement_category_effect, &
                    measurement_prior_probability=model%measurement_prior_probability, &
                    measurement_group=model%measurement_group, sparse_x=model%sparse_x, sparse_z=model%sparse_z)
            else if (use_sparse_scalar) then
                call heterogeneous_multi_term_mixed_mcmc(model%family, model%y, model%additional, model%additional2, &
                    x_design, z_design, model%random_term, model%a_inverse, prior%beta_mean, prior%beta_precision, &
                    prior%g_scale, prior%g_df, prior%r_scale, prior%r_df, control%update_r, control%proposal_scale, &
                    control%iterations, control%burn, control%thin, state, result%scalar, info, observed=observed, &
                    update_g=update_g, initial_g=initial_g, initial_r=initial_r, g_update_mode=g_mode, &
                    r_update_mode=control%r_update_mode, g_split=g_split, r_split=r_split, &
                    g_fixed_covariance=g_fixed, r_fixed_covariance=r_fixed, g_ante_order=g_ante_order, &
                    r_ante_order=control%r_ante_order, g_ante_common_beta=g_ante_common_beta, &
                    r_ante_common_beta=control%r_ante_common_beta, g_ante_common_variance=g_ante_common_variance, &
                    r_ante_common_variance=control%r_ante_common_variance, sparse_x=model%sparse_x, &
                    sparse_z=model%sparse_z)
            else if (allocated(model%measurement_category_effect)) then
                call heterogeneous_multi_term_mixed_mcmc(model%family, model%y, model%additional, model%additional2, &
                    x_design, z_design, model%random_term, model%a_inverse, prior%beta_mean, prior%beta_precision, &
                    prior%g_scale, prior%g_df, prior%r_scale, prior%r_df, control%update_r, control%proposal_scale, &
                    control%iterations, control%burn, control%thin, state, result%scalar, info, observed=observed, &
                    update_g=update_g, initial_g=initial_g, initial_r=initial_r, g_update_mode=g_mode, &
                    r_update_mode=control%r_update_mode, g_split=g_split, r_split=r_split, &
                    g_fixed_covariance=g_fixed, r_fixed_covariance=r_fixed, g_ante_order=g_ante_order, &
                    r_ante_order=control%r_ante_order, g_ante_common_beta=g_ante_common_beta, &
                    r_ante_common_beta=control%r_ante_common_beta, g_ante_common_variance=g_ante_common_variance, &
                    r_ante_common_variance=control%r_ante_common_variance, &
                    measurement_category_effect=model%measurement_category_effect, &
                    measurement_prior_probability=model%measurement_prior_probability, &
                    measurement_group=model%measurement_group)
            else
                call heterogeneous_multi_term_mixed_mcmc(model%family, model%y, model%additional, model%additional2, &
                    x_design, z_design, model%random_term, model%a_inverse, prior%beta_mean, prior%beta_precision, &
                    prior%g_scale, prior%g_df, prior%r_scale, prior%r_df, control%update_r, control%proposal_scale, &
                    control%iterations, control%burn, control%thin, state, result%scalar, info, observed=observed, &
                    update_g=update_g, initial_g=initial_g, initial_r=initial_r, g_update_mode=g_mode, &
                    r_update_mode=control%r_update_mode, g_split=g_split, r_split=r_split, &
                    g_fixed_covariance=g_fixed, r_fixed_covariance=r_fixed, g_ante_order=g_ante_order, &
                    r_ante_order=control%r_ante_order, g_ante_common_beta=g_ante_common_beta, &
                    r_ante_common_beta=control%r_ante_common_beta, g_ante_common_variance=g_ante_common_variance, &
                    r_ante_common_variance=control%r_ante_common_variance)
            end if
        case (mcmcglmm_engine_two_process)
            call row_observed_mask(model, observed_rows)
            call two_part_multi_term_mixed_mcmc(model%grouped_family, model%y_vector, model%trials, x_design, &
                z_design, model%random_term, model%a_inverse, prior%beta_mean, prior%beta_precision, prior%g_scale, &
                prior%g_df, prior%r_scale, prior%r_df, control%update_r, control%proposal_scale, &
                control%iterations, control%burn, control%thin, state, result%grouped, info, &
                observed=observed_rows, update_g=update_g, initial_g=initial_g, initial_r=initial_r, &
                g_update_mode=g_mode, r_update_mode=control%r_update_mode, g_split=g_split, r_split=r_split, &
                g_fixed_covariance=g_fixed, r_fixed_covariance=r_fixed, g_ante_order=g_ante_order, &
                r_ante_order=control%r_ante_order, g_ante_common_beta=g_ante_common_beta, &
                r_ante_common_beta=control%r_ante_common_beta, g_ante_common_variance=g_ante_common_variance, &
                r_ante_common_variance=control%r_ante_common_variance)
        case (mcmcglmm_engine_multinomial)
            call row_observed_mask(model, observed_rows)
            call multinomial_multi_term_mixed_mcmc(model%grouped_family, model%grouped_response, x_design, z_design, &
                model%random_term, model%a_inverse, prior%beta_mean, prior%beta_precision, prior%g_scale, &
                prior%g_df, prior%r_scale, prior%r_df, control%update_r, control%proposal_scale, &
                control%iterations, control%burn, control%thin, state, result%grouped, info, &
                observed=observed_rows, update_g=update_g, initial_g=initial_g, initial_r=initial_r, &
                g_update_mode=g_mode, r_update_mode=control%r_update_mode, g_split=g_split, r_split=r_split, &
                g_fixed_covariance=g_fixed, r_fixed_covariance=r_fixed, g_ante_order=g_ante_order, &
                r_ante_order=control%r_ante_order, g_ante_common_beta=g_ante_common_beta, &
                r_ante_common_beta=control%r_ante_common_beta, g_ante_common_variance=g_ante_common_variance, &
                r_ante_common_variance=control%r_ante_common_variance)
        end select
    end subroutine mcmcglmm_fit_numeric

    pure subroutine mcmcglmm_validate_numeric(model, prior, control, latent_dimension, nterm, info)
        !! Check allocations required by the selected numeric engine before entering a sampler.
        type(mcmcglmm_numeric_model), intent(in) :: model !! Numeric response, design, and engine specification.
        type(mcmcglmm_numeric_prior), intent(in) :: prior !! Priors required by all supported multi-term engines.
        type(mcmcglmm_control), intent(in) :: control !! MCMC schedule and covariance-routing controls.
        integer, intent(out) :: latent_dimension !! Number of latent traits implied by the selected engine.
        integer, intent(out) :: nterm !! Number of random-effect covariance structures.
        integer, intent(out) :: info !! Zero when structurally complete; orchestrator status otherwise.
        logical :: x_present
        logical :: z_present

        info = 0
        latent_dimension = 0
        nterm = 0
        if (.not. any(model%engine == [mcmcglmm_engine_scalar, mcmcglmm_engine_scalar_px, &
            mcmcglmm_engine_two_process, mcmcglmm_engine_multinomial, mcmcglmm_engine_ordinal, &
            mcmcglmm_engine_threshold, mcmcglmm_engine_theta_scale, mcmcglmm_engine_structural, &
            mcmcglmm_engine_covu])) then
            info = mcmcglmm_orchestrator_invalid_engine
            return
        end if
        call validate_design_component(model%x, model%sparse_x, x_present, info)
        if (info /= 0) return
        if (.not. x_present .or. .not. allocated(prior%beta_mean) .or. .not. allocated(prior%beta_precision)) then
            info = mcmcglmm_orchestrator_missing_input
            return
        end if
        if (model%engine == mcmcglmm_engine_covu) then
            if (.not. allocated(model%y) .or. .not. allocated(model%random_loading) .or. &
                .not. allocated(prior%joint_scale)) then
                info = mcmcglmm_orchestrator_missing_input
                return
            end if
            latent_dimension = size(model%y, 2)
            if (control%iterations <= control%burn .or. control%burn < 0 .or. control%thin < 1) then
                info = mcmcglmm_orchestrator_invalid_control
            end if
            return
        end if
        call validate_design_component(model%z, model%sparse_z, z_present, info)
        if (info /= 0) return
        if (.not. z_present .or. .not. allocated(model%a_inverse) .or. &
            .not. allocated(prior%g_scale) .or. .not. allocated(prior%g_df) .or. &
            .not. allocated(prior%r_scale)) then
            info = mcmcglmm_orchestrator_missing_input
            return
        end if
        nterm = size(prior%g_scale, 3)
        if (any(model%engine == [mcmcglmm_engine_scalar, mcmcglmm_engine_scalar_px, &
            mcmcglmm_engine_two_process, mcmcglmm_engine_multinomial, mcmcglmm_engine_theta_scale, &
            mcmcglmm_engine_structural])) then
            if (.not. allocated(model%random_term)) then
                info = mcmcglmm_orchestrator_missing_input
                return
            end if
        end if
        select case (model%engine)
        case (mcmcglmm_engine_scalar, mcmcglmm_engine_scalar_px)
            if (.not. allocated(model%family) .or. .not. allocated(model%y) .or. &
                .not. allocated(model%additional) .or. .not. allocated(model%additional2)) then
                info = mcmcglmm_orchestrator_missing_input
                return
            end if
            latent_dimension = size(model%y, 2)
            if (model%engine == mcmcglmm_engine_scalar_px) then
                if (.not. allocated(prior%alpha_mean) .or. .not. allocated(prior%alpha_precision)) then
                    info = mcmcglmm_orchestrator_missing_input
                    return
                end if
            end if
            if (model%engine == mcmcglmm_engine_scalar) then
                if ((allocated(model%measurement_category_effect) .neqv. &
                    allocated(model%measurement_prior_probability)) .or. &
                    (allocated(model%measurement_category_effect) .neqv. allocated(model%measurement_group))) then
                    info = mcmcglmm_orchestrator_missing_input
                    return
                end if
            end if
        case (mcmcglmm_engine_two_process)
            if (.not. allocated(model%y_vector) .or. .not. allocated(model%trials)) then
                info = mcmcglmm_orchestrator_missing_input
                return
            end if
            latent_dimension = 2
        case (mcmcglmm_engine_multinomial)
            if (.not. allocated(model%grouped_response)) then
                info = mcmcglmm_orchestrator_missing_input
                return
            end if
            latent_dimension = size(model%grouped_response, 2) - 1
        case (mcmcglmm_engine_ordinal, mcmcglmm_engine_threshold)
            if (.not. allocated(model%y_category) .or. .not. allocated(model%cutpoints)) then
                info = mcmcglmm_orchestrator_missing_input
                return
            end if
            latent_dimension = 1
        case (mcmcglmm_engine_theta_scale)
            if (.not. allocated(model%y) .or. .not. allocated(model%x_scale) .or. &
                .not. allocated(model%z_scale)) then
                info = mcmcglmm_orchestrator_missing_input
                return
            end if
            latent_dimension = size(model%y, 2)
        case (mcmcglmm_engine_structural)
            if (.not. allocated(model%y) .or. .not. allocated(model%structural_basis) .or. &
                .not. allocated(prior%structural_mean) .or. .not. allocated(prior%structural_precision) .or. &
                .not. allocated(prior%structural_proposal_sd)) then
                info = mcmcglmm_orchestrator_missing_input
                return
            end if
            latent_dimension = size(model%y, 2)
        end select
        if (control%iterations <= control%burn .or. control%burn < 0 .or. control%thin < 1 .or. &
            control%proposal_scale <= 0.0_dp .or. nterm < 1 .or. latent_dimension < 1) then
            info = mcmcglmm_orchestrator_invalid_control
        end if
    end subroutine mcmcglmm_validate_numeric

    pure subroutine validate_design_component(dense, sparse, present_design, info)
        !! Validate one dense-or-CSR design input and reject ambiguous dual representations.
        real(dp), allocatable, intent(in) :: dense(:, :) !! Optional dense design representation.
        type(mcmcglmm_sparse_matrix), intent(in) :: sparse !! Optional CSR design representation.
        logical, intent(out) :: present_design !! True when exactly one valid representation is supplied.
        integer, intent(out) :: info !! Zero on success; invalid-control status for ambiguity or malformed CSR.
        logical :: sparse_present
        integer :: sparse_info

        info = 0
        present_design = .false.
        sparse_present = sparse_storage_present(sparse)
        if (allocated(dense) .and. sparse_present) then
            info = mcmcglmm_orchestrator_invalid_control
            return
        end if
        if (sparse_present) then
            call sparse_validate(sparse, sparse_info)
            if (sparse_info /= 0) then
                info = mcmcglmm_orchestrator_invalid_control
                return
            end if
        end if
        present_design = allocated(dense) .or. sparse_present
    end subroutine validate_design_component

    pure logical function sparse_storage_present(matrix) result(is_present)
        !! Detect an attempted sparse input, including malformed partial storage.
        type(mcmcglmm_sparse_matrix), intent(in) :: matrix !! Sparse design component being queried.

        is_present = matrix%nrow /= 0 .or. matrix%ncol /= 0 .or. allocated(matrix%row_pointer) .or. &
            allocated(matrix%column) .or. allocated(matrix%value)
    end function sparse_storage_present

    pure subroutine resolve_numeric_design(model, x, z, info)
        !! Materialize dense designs for current sampler kernels from either dense or CSR model inputs.
        type(mcmcglmm_numeric_model), intent(in) :: model !! Validated model containing dense or sparse designs.
        real(dp), allocatable, intent(out) :: x(:, :) !! Materialized fixed-effect design matrix.
        real(dp), allocatable, intent(out) :: z(:, :) !! Materialized random-effect design, empty for covu.
        integer, intent(out) :: info !! Zero on success; invalid-control status if conversion fails.
        integer :: sparse_info

        info = 0
        if (allocated(model%x)) then
            x = model%x
        else
            call sparse_to_dense(model%sparse_x, x, sparse_info)
            if (sparse_info /= 0) then
                info = mcmcglmm_orchestrator_invalid_control
                return
            end if
        end if
        if (allocated(model%z)) then
            z = model%z
        else if (sparse_is_initialized(model%sparse_z)) then
            call sparse_to_dense(model%sparse_z, z, sparse_info)
            if (sparse_info /= 0) then
                info = mcmcglmm_orchestrator_invalid_control
                return
            end if
        else
            allocate(z(0, 0))
        end if
    end subroutine resolve_numeric_design

    pure subroutine normalize_covariance_control(prior, control, t, nterm, update_g, initial_g, initial_r, &
                                                   g_mode, g_split, g_fixed, r_split, r_fixed, g_ante_order, &
                                                   g_ante_common_beta, g_ante_common_variance, info)
        !! Materialize optional covariance settings so every routed sampler sees one normalized configuration.
        type(mcmcglmm_numeric_prior), intent(in) :: prior !! Covariance priors used to construct default initial values.
        type(mcmcglmm_control), intent(in) :: control !! User covariance routing with allocatable optional components.
        integer, intent(in) :: t !! Latent trait dimension.
        integer, intent(in) :: nterm !! Number of G covariance structures.
        logical, allocatable, intent(out) :: update_g(:) !! Normalized per-G update mask.
        real(dp), allocatable, intent(out) :: initial_g(:, :, :) !! Normalized initial G matrices.
        real(dp), allocatable, intent(out) :: initial_r(:, :) !! Normalized initial residual covariance.
        integer, allocatable, intent(out) :: g_mode(:) !! Normalized per-G update codes.
        integer, allocatable, intent(out) :: g_split(:) !! Normalized per-G split dimensions.
        real(dp), allocatable, intent(out) :: g_fixed(:, :, :) !! Normalized fixed covariance blocks.
        integer, intent(out) :: r_split !! Normalized residual split dimension.
        real(dp), allocatable, intent(out) :: r_fixed(:, :) !! Normalized fixed residual covariance block.
        integer, allocatable, intent(out) :: g_ante_order(:) !! Normalized antedependence orders.
        logical, allocatable, intent(out) :: g_ante_common_beta(:) !! Normalized shared-coefficient flags.
        logical, allocatable, intent(out) :: g_ante_common_variance(:) !! Normalized shared-variance flags.
        integer, intent(out) :: info !! Zero on success; invalid-control status for inconsistent option shapes.
        integer :: term

        info = 0
        allocate(update_g(nterm), g_mode(nterm), g_split(nterm), g_ante_order(nterm))
        allocate(g_ante_common_beta(nterm), g_ante_common_variance(nterm))
        allocate(initial_g(t, t, nterm), initial_r(t, t), g_fixed(t, t, nterm), r_fixed(t, t))
        update_g = .true.
        g_mode = 1
        g_split = t
        g_ante_order = 1
        g_ante_common_beta = .false.
        g_ante_common_variance = .false.
        g_fixed = 0.0_dp
        r_fixed = 0.0_dp
        r_split = t
        do term = 1, nterm
            initial_g(:, :, term) = prior%g_scale(:, :, term) / &
                max(prior%g_df(term) - real(t + 1, dp), 1.0_dp)
        end do
        initial_r = prior%r_scale / max(prior%r_df - real(t + 1, dp), 1.0_dp)

        if (allocated(control%update_g)) then
            if (size(control%update_g) /= nterm) then
                info = mcmcglmm_orchestrator_invalid_control
                return
            end if
            update_g = control%update_g
        end if
        if (allocated(control%g_update_mode)) then
            if (size(control%g_update_mode) /= nterm) then
                info = mcmcglmm_orchestrator_invalid_control
                return
            end if
            g_mode = control%g_update_mode
        end if
        where (.not. update_g) g_mode = 0
        if (allocated(control%g_split)) then
            if (size(control%g_split) /= nterm) then
                info = mcmcglmm_orchestrator_invalid_control
                return
            end if
            g_split = control%g_split
        end if
        if (control%r_split >= 0) r_split = control%r_split
        if (allocated(control%initial_g)) then
            if (any(shape(control%initial_g) /= [t, t, nterm])) then
                info = mcmcglmm_orchestrator_invalid_control
                return
            end if
            initial_g = control%initial_g
        end if
        if (allocated(control%initial_r)) then
            if (any(shape(control%initial_r) /= [t, t])) then
                info = mcmcglmm_orchestrator_invalid_control
                return
            end if
            initial_r = control%initial_r
        end if
        if (allocated(control%g_fixed_covariance)) then
            if (any(shape(control%g_fixed_covariance) /= [t, t, nterm])) then
                info = mcmcglmm_orchestrator_invalid_control
                return
            end if
            g_fixed = control%g_fixed_covariance
        end if
        if (allocated(control%r_fixed_covariance)) then
            if (any(shape(control%r_fixed_covariance) /= [t, t])) then
                info = mcmcglmm_orchestrator_invalid_control
                return
            end if
            r_fixed = control%r_fixed_covariance
        end if
        if (allocated(control%g_ante_order)) then
            if (size(control%g_ante_order) /= nterm) then
                info = mcmcglmm_orchestrator_invalid_control
                return
            end if
            g_ante_order = control%g_ante_order
        end if
        if (allocated(control%g_ante_common_beta)) then
            if (size(control%g_ante_common_beta) /= nterm) then
                info = mcmcglmm_orchestrator_invalid_control
                return
            end if
            g_ante_common_beta = control%g_ante_common_beta
        end if
        if (allocated(control%g_ante_common_variance)) then
            if (size(control%g_ante_common_variance) /= nterm) then
                info = mcmcglmm_orchestrator_invalid_control
                return
            end if
            g_ante_common_variance = control%g_ante_common_variance
        end if
    end subroutine normalize_covariance_control

    pure subroutine scalar_observed_mask(model, observed)
        !! Return the supplied scalar-response mask or an all-observed default.
        type(mcmcglmm_numeric_model), intent(in) :: model !! Model containing the optional scalar mask.
        logical, allocatable, intent(out) :: observed(:, :) !! Materialized n by trait observation mask.

        allocate(observed(size(model%y, 1), size(model%y, 2)), source=.true.)
        if (allocated(model%observed)) observed = model%observed
    end subroutine scalar_observed_mask

    pure subroutine row_observed_mask(model, observed)
        !! Return the supplied grouped-response mask or an all-observed default.
        type(mcmcglmm_numeric_model), intent(in) :: model !! Model containing the optional grouped-row mask.
        logical, allocatable, intent(out) :: observed(:) !! Materialized row observation mask.
        integer :: n

        if (model%engine == mcmcglmm_engine_two_process) then
            n = size(model%y_vector)
        else if (any(model%engine == [mcmcglmm_engine_ordinal, mcmcglmm_engine_threshold])) then
            n = size(model%y_category)
        else
            n = size(model%grouped_response, 1)
        end if
        allocate(observed(n), source=.true.)
        if (allocated(model%observed_rows)) observed = model%observed_rows
    end subroutine row_observed_mask

end module mcmcglmm_orchestrator
