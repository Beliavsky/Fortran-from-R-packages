! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from MCMCglmm 2.36 computational code; see NOTICE.md and upstream/.
program test_mcmcglmm
    use ieee_arithmetic, only : ieee_is_finite
    use mcmcglmm
    use ape_types, only : phylo_tree
    implicit none

    integer :: failures

    failures = 0
    call test_rng(failures)
    call test_pedigree(failures)
    call test_conditional_normal(failures)
    call test_sparse_foundation(failures)
    call test_distributions(failures)
    call test_covariance_structures(failures)
    call test_covu_sampler(failures)
    call test_matrix_and_posterior(failures)
    call test_utilities(failures)
    call test_gaussian_sampler(failures)
    call test_ante_sampler(failures)
    call test_phylo(failures)
    call test_design_helpers(failures)
    call test_ordinal_sampler(failures)
    call test_spline(failures)
    call test_parameter_expansion(failures)
    call test_engine_features(failures)
    call test_theta_scale_sampler(failures)
    call test_structural_sampler(failures)
    call test_multi_term_sampler(failures)
    call test_multi_term_px_sampler(failures)
    call test_unified_sampler(failures)
    call test_grouped_multi_term_sampler(failures)
    call test_numeric_orchestrator(failures)
    call test_prediction(failures)
    call test_simulation(failures)
    call test_family_likelihoods(failures)
    call test_family_samplers(failures)

    if (failures /= 0) then
        print '(a,i0)', 'MCMCglmm tests failed: ', failures
        error stop 1
    end if
    print '(a)', 'All MCMCglmm tests passed.'

contains

    subroutine check(condition, message, failures)
        logical, intent(in) :: condition !! Assertion condition that must be true for the test to pass.
        character(len=*), intent(in) :: message !! Human-readable assertion description printed when the condition is false.
        integer, intent(inout) :: failures !! Running count of failed assertions.

        if (.not. condition) then
            failures = failures + 1
            print '(a)', 'FAIL: ' // trim(message)
        end if
    end subroutine check

    subroutine check_close(actual, expected, tolerance, message, failures)
        real(dp), intent(in) :: actual !! Computed scalar value under test.
        real(dp), intent(in) :: expected !! Reference scalar value.
        real(dp), intent(in) :: tolerance !! Maximum permitted absolute difference.
        character(len=*), intent(in) :: message !! Human-readable assertion description.
        integer, intent(inout) :: failures !! Running count of failed assertions.

        call check(abs(actual - expected) <= tolerance, message, failures)
    end subroutine check_close

    subroutine test_rng(failures)
        integer, intent(inout) :: failures !! Running count of failed assertions.
        type(rng_state) :: state
        real(dp) :: u

        call rng_seed(state, 1_8)
        call rng_uniform(state, u)
        call check_close(u, 16807.0_dp / 2147483647.0_dp, 1.0e-15_dp, 'Park-Miller first draw', failures)
        call rng_uniform(state, u)
        call check_close(u, 282475249.0_dp / 2147483647.0_dp, 1.0e-15_dp, 'Park-Miller second draw', failures)
    end subroutine test_rng

    subroutine test_pedigree(failures)
        integer, intent(inout) :: failures !! Running count of failed assertions.
        integer :: info
        integer :: dam(4)
        integer :: sire(4)
        logical, allocatable :: mask(:)
        real(dp), allocatable :: a(:, :)
        real(dp), allocatable :: ainv(:, :)
        real(dp), allocatable :: d(:)
        real(dp), allocatable :: f(:)
        real(dp), parameter :: expected_a(4, 4) = reshape([ &
            1.0_dp, 0.0_dp, 0.5_dp, 0.75_dp, &
            0.0_dp, 1.0_dp, 0.5_dp, 0.25_dp, &
            0.5_dp, 0.5_dp, 1.0_dp, 0.75_dp, &
            0.75_dp, 0.25_dp, 0.75_dp, 1.25_dp], [4, 4])
        real(dp), parameter :: expected_inverse(4, 4) = reshape([ &
            2.0_dp, 0.5_dp, -0.5_dp, -1.0_dp, &
            0.5_dp, 1.5_dp, -1.0_dp, 0.0_dp, &
            -0.5_dp, -1.0_dp, 2.5_dp, -1.0_dp, &
            -1.0_dp, 0.0_dp, -1.0_dp, 2.0_dp], [4, 4])

        dam = [0, 0, 1, 1]
        sire = [0, 0, 2, 3]
        call pedigree_relationship(dam, sire, a, f, d, info)
        call check(info == 0, 'pedigree relationship status', failures)
        call check(maxval(abs(a - expected_a)) < 1.0e-12_dp, 'pedigree relationship matrix', failures)
        call check_close(f(4), 0.25_dp, 1.0e-12_dp, 'pedigree inbreeding', failures)
        call check(maxval(abs(d - [1.0_dp, 1.0_dp, 0.5_dp, 0.5_dp])) < 1.0e-12_dp, &
            'Mendelian variances', failures)
        call pedigree_inverse(dam, sire, ainv, f, d, info)
        call check(info == 0, 'pedigree inverse status', failures)
        call check(maxval(abs(ainv - expected_inverse)) < 1.0e-11_dp, 'pedigree inverse matrix', failures)
        call prune_pedigree_mask(dam, sire, [4], mask, info)
        call check(info == 0 .and. all(mask), 'pedigree ancestor pruning mask', failures)
    end subroutine test_pedigree

    subroutine test_conditional_normal(failures)
        integer, intent(inout) :: failures !! Running count of failed assertions.
        integer :: info
        real(dp) :: density
        real(dp) :: probability
        real(dp), allocatable :: ccov(:, :)
        real(dp), allocatable :: cmean(:)
        real(dp) :: covariance(2, 2)
        real(dp) :: mean_value(2)
        real(dp) :: x(2)

        covariance = reshape([2.0_dp, 1.0_dp, 1.0_dp, 3.0_dp], [2, 2])
        mean_value = [0.0_dp, 1.0_dp]
        x = [0.5_dp, 2.0_dp]
        call conditional_mvn_parameters(mean_value, covariance, x, [1], [2], cmean, ccov, info)
        call check(info == 0, 'conditional MVN status', failures)
        call check_close(cmean(1), 1.0_dp / 3.0_dp, 1.0e-12_dp, 'conditional MVN mean', failures)
        call check_close(ccov(1, 1), 5.0_dp / 3.0_dp, 1.0e-12_dp, 'conditional MVN variance', failures)
        call conditional_mvn_log_density(x, mean_value, covariance, [1], [2], density, info)
        call check(info == 0 .and. ieee_is_finite(density), 'conditional MVN log density finite', failures)
        call truncated_conditional_mvn_log_probability(mean_value, covariance, x, 1, -1.0_dp, 1.0_dp, probability, info)
        call check(info == 0 .and. probability < 0.0_dp, 'conditional interval log probability', failures)
    end subroutine test_conditional_normal

    subroutine test_sparse_foundation(failures)
        !! Verify canonical CSR construction, products, conversion, and validation.
        integer, intent(inout) :: failures !! Running count of failed assertions.
        integer :: i
        integer :: info
        integer, parameter :: row(8) = [1, 1, 1, 2, 2, 2, 3, 3]
        integer, parameter :: column(8) = [1, 3, 3, 2, 1, 1, 1, 3]
        real(dp), parameter :: coordinate_value(8) = [1.0_dp, 0.5_dp, 1.5_dp, 3.0_dp, &
            1.0_dp, -1.0_dp, 4.0_dp, 5.0_dp]
        real(dp), parameter :: expected(3, 3) = reshape([1.0_dp, 0.0_dp, 4.0_dp, &
            0.0_dp, 3.0_dp, 0.0_dp, 2.0_dp, 0.0_dp, 5.0_dp], [3, 3])
        real(dp), parameter :: rhs(3, 2) = reshape([1.0_dp, 2.0_dp, 3.0_dp, &
            -2.0_dp, 0.5_dp, 1.0_dp], [3, 2])
        real(dp), parameter :: spd(4, 4) = reshape([4.0_dp, 0.0_dp, 2.0_dp, 2.0_dp, &
            0.0_dp, 4.0_dp, 2.0_dp, -1.0_dp, 2.0_dp, 2.0_dp, 6.0_dp, 0.0_dp, &
            2.0_dp, -1.0_dp, 0.0_dp, 5.3125_dp], [4, 4])
        real(dp), parameter :: precision_rhs(4) = [1.0_dp, -2.0_dp, 0.5_dp, 3.0_dp]
        real(dp), allocatable :: crossproduct(:, :)
        real(dp), allocatable :: dense(:, :)
        real(dp), allocatable :: innovation(:)
        real(dp), allocatable :: mean_value(:)
        real(dp), allocatable :: matrix_product(:, :)
        real(dp), allocatable :: normal_draw(:)
        real(dp), allocatable :: ordered_mean(:)
        real(dp), allocatable :: ordered_sample(:)
        real(dp), allocatable :: sample(:)
        real(dp), allocatable :: solution(:)
        real(dp), allocatable :: stacked_dense(:, :)
        real(dp), allocatable :: stress_design(:, :)
        real(dp), allocatable :: stress_left_dense(:, :)
        real(dp), allocatable :: stress_right_dense(:, :)
        real(dp), allocatable :: expected_stacked(:, :)
        real(dp) :: diagonal_matrix(4, 4)
        real(dp), allocatable :: transpose_product(:, :)
        real(dp), allocatable :: vector_product(:)
        type(mcmcglmm_sparse_matrix) :: factor
        type(mcmcglmm_sparse_matrix) :: analyzed_factor
        type(mcmcglmm_sparse_matrix) :: diagonal_precision
        type(mcmcglmm_sparse_matrix) :: from_dense
        type(mcmcglmm_sparse_matrix) :: invalid
        type(mcmcglmm_sparse_matrix) :: matrix
        type(mcmcglmm_sparse_matrix) :: precision
        type(mcmcglmm_sparse_matrix) :: permuted_precision
        type(mcmcglmm_sparse_matrix) :: stacked
        type(mcmcglmm_sparse_matrix) :: stress_left
        type(mcmcglmm_sparse_matrix) :: stress_right
        type(rng_state) :: reference_state
        type(rng_state) :: sample_state
        integer, allocatable :: permutation(:)
        type(sparse_cholesky_analysis) :: analysis
        type(sparse_precision_cache) :: factor_cache

        call sparse_from_coo(3, 3, row, column, coordinate_value, matrix, info)
        call check(info == 0, 'sparse COO construction status', failures)
        call sparse_validate(matrix, info)
        call check(info == 0 .and. size(matrix%value) == 5, &
            'sparse canonical storage and duplicate cancellation', failures)
        call sparse_to_dense(matrix, dense, info)
        call check(info == 0 .and. maxval(abs(dense - expected)) < 1.0e-14_dp, &
            'sparse-to-dense conversion', failures)
        call sparse_matmul_vector(matrix, [1.0_dp, 2.0_dp, 3.0_dp], vector_product, info)
        call check(info == 0 .and. maxval(abs(vector_product - [7.0_dp, 6.0_dp, 19.0_dp])) < 1.0e-14_dp, &
            'sparse matrix-vector product', failures)
        call sparse_matmul_matrix(matrix, rhs, matrix_product, info)
        call check(info == 0 .and. maxval(abs(matrix_product - matmul(expected, rhs))) < 1.0e-14_dp, &
            'sparse matrix-matrix product', failures)
        call sparse_transpose_matmul_matrix(matrix, rhs, transpose_product, info)
        call check(info == 0 .and. maxval(abs(transpose_product - matmul(transpose(expected), rhs))) < 1.0e-14_dp, &
            'sparse transpose matrix-matrix product', failures)
        call sparse_crossproduct(matrix, crossproduct, info)
        call check(info == 0 .and. maxval(abs(crossproduct - matmul(transpose(expected), expected))) < 1.0e-14_dp, &
            'sparse crossproduct', failures)
        call sparse_from_dense(expected, from_dense, info=info)
        call sparse_transpose_matmul_sparse(matrix, from_dense, crossproduct, info)
        call check(info == 0 .and. maxval(abs(crossproduct - matmul(transpose(expected), expected))) < 1.0e-14_dp, &
            'sparse transpose-times-sparse product', failures)
        call sparse_to_dense(from_dense, dense, info)
        call check(info == 0 .and. maxval(abs(dense - expected)) < 1.0e-14_dp, &
            'dense-to-sparse round trip', failures)
        call sparse_stacked_crossproduct(matrix, from_dense, stacked, info)
        call sparse_to_dense(stacked, stacked_dense, info)
        allocate(expected_stacked(6, 6))
        expected_stacked(:3, :3) = matmul(transpose(expected), expected)
        expected_stacked(:3, 4:) = expected_stacked(:3, :3)
        expected_stacked(4:, :3) = expected_stacked(:3, :3)
        expected_stacked(4:, 4:) = expected_stacked(:3, :3)
        call check(info == 0 .and. maxval(abs(stacked_dense - expected_stacked)) < 1.0e-14_dp, &
            'sparse stacked-design crossproduct', failures)
        allocate(stress_left_dense(200, 8), stress_right_dense(200, 7))
        do i = 1, 200
            stress_left_dense(i, :) = [real(i, dp) / 200.0_dp, 1.0_dp, real(mod(i, 7), dp), &
                real(mod(i, 11), dp), sin(real(i, dp)), cos(real(i, dp)), real(mod(i, 3), dp), 0.0_dp]
            stress_right_dense(i, :) = [real(mod(i, 5), dp), real(i * i, dp) / 40000.0_dp, &
                real(mod(i, 13), dp), 0.25_dp, sin(real(2 * i, dp)), cos(real(3 * i, dp)), 1.0_dp]
        end do
        call sparse_from_dense(stress_left_dense, stress_left, info=info)
        call sparse_from_dense(stress_right_dense, stress_right, info=info)
        call sparse_stacked_crossproduct(stress_left, stress_right, stacked, info)
        call sparse_to_dense(stacked, stacked_dense, info)
        allocate(stress_design(200, 15))
        stress_design(:, :8) = stress_left_dense
        stress_design(:, 9:) = stress_right_dense
        call check(info == 0 .and. maxval(abs(stacked_dense - &
            matmul(transpose(stress_design), stress_design))) < 1.0e-10_dp, &
            'sparse stacked-design crossproduct with repeated dense rows', failures)
        invalid = matrix
        invalid%row_pointer(1) = 0
        call sparse_validate(invalid, info)
        call check(info /= 0, 'sparse validator rejects malformed row pointer', failures)

        call sparse_from_dense(spd, precision, info=info)
        call sparse_cholesky_factor(precision, factor, info)
        call check(info == 0, 'sparse Cholesky factorization status', failures)
        call sparse_cholesky_analyze(precision, analysis, info)
        call sparse_cholesky_factor_analyzed(precision, analysis, analyzed_factor, info)
        call check(info == 0 .and. all(analyzed_factor%row_pointer == factor%row_pointer) .and. &
            all(analyzed_factor%column == factor%column) .and. &
            maxval(abs(analyzed_factor%value - factor%value)) < 1.0e-14_dp, &
            'reused symbolic Cholesky matches combined factorization', failures)
        call sparse_to_dense(factor, dense, info)
        call check(info == 0 .and. maxval(abs(matmul(dense, transpose(dense)) - spd)) < 1.0e-13_dp, &
            'sparse Cholesky reconstructs precision matrix', failures)
        call check(size(factor%value) == 9, 'sparse Cholesky symbolic fill pattern', failures)
        call check_close(factor%value(factor%row_pointer(4) + 2), -0.25_dp, 1.0e-14_dp, &
            'sparse Cholesky numerical fill value', failures)
        call sparse_cholesky_solve(factor, precision_rhs, solution, info)
        call check(info == 0 .and. maxval(abs(matmul(spd, solution) - precision_rhs)) < 1.0e-13_dp, &
            'sparse Cholesky precision solve', failures)
        call rng_seed(sample_state, 830271_8)
        call sample_mvn_sparse_precision(sample_state, precision_rhs, precision, sample, mean_value, info)
        call check(info == 0 .and. maxval(abs(matmul(spd, mean_value) - precision_rhs)) < 1.0e-13_dp, &
            'sparse precision Gaussian conditional mean', failures)
        allocate(normal_draw(4))
        call rng_seed(reference_state, 830271_8)
        call rng_normal(reference_state, normal_draw(1))
        call rng_normal(reference_state, normal_draw(2))
        call rng_normal(reference_state, normal_draw(3))
        call rng_normal(reference_state, normal_draw(4))
        call sparse_cholesky_transpose_solve(factor, normal_draw, innovation, info)
        call check(info == 0 .and. maxval(abs(sample - mean_value - innovation)) < 1.0e-14_dp, &
            'sparse precision Gaussian draw uses inverse-transpose innovation', failures)
        call sparse_reverse_cuthill_mckee(precision, permutation, info)
        call check(info == 0 .and. all(permutation == [2, 4, 3, 1]), &
            'deterministic reverse Cuthill-McKee ordering', failures)
        call sparse_symmetric_permute(precision, permutation, permuted_precision, info)
        call sparse_to_dense(permuted_precision, dense, info)
        call check(info == 0 .and. maxval(abs(dense - spd(permutation, permutation))) < 1.0e-14_dp, &
            'symmetric CSR permutation', failures)
        call rng_seed(sample_state, 830272_8)
        call sample_mvn_sparse_precision(sample_state, precision_rhs, precision, ordered_sample, ordered_mean, info, &
            use_rcm=.true.)
        call check(info == 0 .and. maxval(abs(matmul(spd, ordered_mean) - precision_rhs)) < 1.0e-13_dp, &
            'ordered sparse precision Gaussian conditional mean', failures)
        diagonal_matrix = 0.0_dp
        do i = 1, 4
            diagonal_matrix(i, i) = real(i + 1, dp)
        end do
        call sparse_from_dense(diagonal_matrix, diagonal_precision, info=info)
        call rng_seed(sample_state, 830273_8)
        call sample_mvn_sparse_precision(sample_state, precision_rhs, diagonal_precision, ordered_sample, &
            ordered_mean, info, use_rcm=.true., cache=factor_cache)
        call rng_seed(sample_state, 830274_8)
        call sample_mvn_sparse_precision(sample_state, precision_rhs, diagonal_precision, ordered_sample, &
            ordered_mean, info, use_rcm=.true., cache=factor_cache)
        call check(info == 0 .and. factor_cache%analysis_count == 1, &
            'sparse precision cache reuses symbolic analysis', failures)
        call rng_seed(sample_state, 830275_8)
        call sample_mvn_sparse_precision(sample_state, precision_rhs, precision, ordered_sample, ordered_mean, info, &
            use_rcm=.true., cache=factor_cache)
        call check(info == 0 .and. factor_cache%analysis_count == 2, &
            'sparse precision cache rebuilds after structural expansion', failures)
    end subroutine test_sparse_foundation

    subroutine test_distributions(failures)
        integer, intent(inout) :: failures !! Running count of failed assertions.
        integer :: info
        real(dp) :: draw
        real(dp), allocatable :: conditioned_iw(:, :)
        real(dp), allocatable :: iw(:, :)
        real(dp) :: v(2, 2)
        real(dp) :: cm(1, 1)
        type(rng_state) :: state

        call check_close(pkk_probability([0.5_dp, 0.5_dp], 2.0_dp), 0.5_dp, 1.0e-12_dp, &
            'pkk equal two-category probability', failures)
        call rng_seed(state, 2345_8)
        call truncated_normal_sample(state, 0.0_dp, 1.0_dp, -0.25_dp, 0.5_dp, draw, info)
        call check(info == 0 .and. draw > -0.25_dp .and. draw < 0.5_dp, 'truncated normal bounds', failures)
        v = reshape([1.0_dp, 0.2_dp, 0.2_dp, 2.0_dp], [2, 2])
        call riw_mcmcglmm(state, v, 5.0_dp, iw, info)
        call check(info == 0 .and. all(ieee_is_finite(iw)) .and. iw(1, 1) > 0.0_dp .and. iw(2, 2) > 0.0_dp, &
            'inverse-Wishart sample is finite positive on diagonal', failures)
        cm(1, 1) = 1.75_dp
        call riw_mcmcglmm_conditioned(state, v, 5.0_dp, 2, cm, conditioned_iw, info)
        call check(info == 0 .and. abs(conditioned_iw(2, 2) - cm(1, 1)) < 1.0e-12_dp, &
            'conditional inverse-Wishart fixed block', failures)
        call check(abs(conditioned_iw(1, 2) - conditioned_iw(2, 1)) < 1.0e-12_dp, &
            'conditional inverse-Wishart symmetry', failures)
    end subroutine test_distributions

    subroutine test_covariance_structures(failures)
        integer, intent(inout) :: failures !! Running count of failed assertions.
        logical :: accepted
        integer :: info
        real(dp), allocatable :: covariance(:, :)
        real(dp), allocatable :: fixed_block(:, :)
        real(dp), allocatable :: conditional_covariance(:, :)
        real(dp), allocatable :: joint_covariance(:, :)
        real(dp), allocatable :: leading_covariance(:, :)
        real(dp), allocatable :: reconstructed_covariance(:, :)
        real(dp), allocatable :: regression(:, :)
        real(dp) :: posterior_sum(3, 3)
        real(dp) :: prior_scale(3, 3)
        real(dp) :: old_covariance(3, 3)
        real(dp) :: fixed_one(1, 1)
        real(dp) :: joint_reference(3, 3)
        type(rng_state) :: state

        posterior_sum = reshape([5.0_dp, 0.8_dp, 0.3_dp, 0.8_dp, 4.0_dp, 0.5_dp, &
            0.3_dp, 0.5_dp, 3.0_dp], [3, 3])
        prior_scale = 0.0_dp
        prior_scale(1, 1) = 1.0_dp
        prior_scale(2, 2) = 1.5_dp
        prior_scale(3, 3) = 2.0_dp
        old_covariance = prior_scale

        call rng_seed(state, 44021_8)
        fixed_one(1, 1) = 2.0_dp
        call conditioned_covariance_update(state, posterior_sum, 7.0_dp, 3, fixed_one, covariance, info)
        call check(info == 0, 'conditioned covariance update status', failures)
        if (info == 0) then
            call check_close(covariance(3, 3), 2.0_dp, 1.0e-14_dp, &
                'conditioned covariance fixed lower block', failures)
            call check(maxval(abs(covariance - transpose(covariance))) < 1.0e-12_dp, &
                'conditioned covariance symmetry', failures)
        end if

        call rng_seed(state, 44022_8)
        call identity_direct_sum_update(state, posterior_sum, 7.0_dp, 2, covariance, info)
        call check(info == 0, 'identity direct-sum covariance update status', failures)
        if (info == 0) then
            call check_close(covariance(3, 3), 1.0_dp, 1.0e-14_dp, &
                'identity direct-sum fixed identity coordinate', failures)
            call check(maxval(abs(covariance(1:2, 3))) < 1.0e-14_dp .and. &
                maxval(abs(covariance(3, 1:2))) < 1.0e-14_dp, &
                'identity direct-sum cross block is zero', failures)
        end if

        call rng_seed(state, 44023_8)
        call correlation_structure_update(state, posterior_sum, 8.0_dp, 4.0_dp, prior_scale, old_covariance, &
            covariance, accepted, info)
        call check(info == 0, 'correlation covariance update status', failures)
        if (info == 0) then
            call check_close(covariance(1, 1), prior_scale(1, 1), 1.0e-14_dp, &
                'correlation update fixes first marginal variance', failures)
            call check_close(covariance(2, 2), prior_scale(2, 2), 1.0e-14_dp, &
                'correlation update fixes second marginal variance', failures)
            call check_close(covariance(3, 3), prior_scale(3, 3), 1.0e-14_dp, &
                'correlation update fixes third marginal variance', failures)
            call check(maxval(abs(covariance - transpose(covariance))) < 1.0e-12_dp, &
                'correlation update symmetry', failures)
        end if

        call rng_seed(state, 44024_8)
        call correlation_submatrix_update(state, posterior_sum, 8.0_dp, 4.0_dp, prior_scale, 1, &
            old_covariance(2:3, 2:3), covariance, fixed_block, accepted, info)
        call check(info == 0, 'correlation submatrix covariance update status', failures)
        if (info == 0) then
            call check_close(fixed_block(1, 1), prior_scale(2, 2), 1.0e-14_dp, &
                'correlation submatrix fixes first marginal variance', failures)
            call check_close(fixed_block(2, 2), prior_scale(3, 3), 1.0e-14_dp, &
                'correlation submatrix fixes second marginal variance', failures)
            call check(maxval(abs(covariance(2:3, 2:3) - fixed_block)) < 1.0e-12_dp, &
                'correlation submatrix retained as conditioned block', failures)
        end if

        call rng_seed(state, 44025_8)
        call covariance_update_dispatch(state, 2, posterior_sum, 3.0_dp, 4.0_dp, prior_scale, old_covariance, &
            2, fixed_one, covariance, fixed_block, accepted, info)
        call check(info == 0, 'covariance dispatcher conditioned-IW status', failures)
        if (info == 0) then
            call check_close(covariance(3, 3), 2.0_dp, 1.0e-14_dp, &
                'covariance dispatcher preserves fixed conditioned block', failures)
        end if
        call rng_seed(state, 44026_8)
        call covariance_update_dispatch(state, 6, posterior_sum, 3.0_dp, 4.0_dp, prior_scale, old_covariance, &
            2, reshape([0.0_dp], [0, 0]), covariance, fixed_block, accepted, info)
        call check(info == 0, 'covariance dispatcher direct-sum status', failures)
        if (info == 0) then
            call check_close(covariance(3, 3), 1.0_dp, 1.0e-14_dp, &
                'covariance dispatcher direct-sum identity block', failures)
        end if
        call covariance_update_dispatch(state, 5, posterior_sum, 3.0_dp, 4.0_dp, prior_scale, old_covariance, &
            1, reshape([0.0_dp], [0, 0]), covariance, fixed_block, accepted, info)
        call check(info == 5, 'covariance dispatcher directs code 5 to antedependence sampler', failures)

        joint_reference = reshape([2.0_dp, 1.0_dp, 0.5_dp, 1.0_dp, 3.0_dp, 0.2_dp, &
            0.5_dp, 0.2_dp, 4.0_dp], [3, 3])
        call joint_gr_decompose(joint_reference, 1, leading_covariance, regression, conditional_covariance, info)
        call check(info == 0, 'joint G-R Schur decomposition status', failures)
        if (info == 0) then
            call check_close(regression(1, 1), 0.5_dp, 1.0e-14_dp, &
                'joint G-R first regression coefficient', failures)
            call check_close(regression(2, 1), 0.25_dp, 1.0e-14_dp, &
                'joint G-R second regression coefficient', failures)
            call check_close(conditional_covariance(1, 1), 2.5_dp, 1.0e-14_dp, &
                'joint G-R Schur first variance', failures)
            call check_close(conditional_covariance(2, 2), 3.875_dp, 1.0e-14_dp, &
                'joint G-R Schur second variance', failures)
            call joint_gr_compose(leading_covariance, regression, conditional_covariance, &
                reconstructed_covariance, info)
            call check(info == 0 .and. maxval(abs(reconstructed_covariance - joint_reference)) < 1.0e-13_dp, &
                'joint G-R compose/decompose round trip', failures)
        end if
        call joint_gr_covariance_update(state, 0, posterior_sum, 3.0_dp, 4.0_dp, prior_scale, joint_reference, &
            1, fixed_block, joint_covariance, conditional_covariance, regression, accepted, info)
        call check(info == 0 .and. maxval(abs(joint_covariance - joint_reference)) < 1.0e-14_dp, &
            'joint G-R fixed covariance routing', failures)
    end subroutine test_covariance_structures

    subroutine test_covu_sampler(failures)
        integer, intent(inout) :: failures !! Running count of failed assertions.
        integer :: info
        integer :: j
        real(dp) :: beta_mean(1, 1)
        real(dp) :: beta_precision(1, 1)
        real(dp) :: initial_joint(2, 2)
        real(dp) :: joint_scale(2, 2)
        real(dp) :: loading(1, 1)
        real(dp) :: x(6, 1)
        real(dp) :: y(6, 1)
        type(covu_gaussian_mcmc_result) :: result
        type(rng_state) :: state

        x(:, 1) = 1.0_dp
        y(:, 1) = [-0.8_dp, -0.2_dp, 0.0_dp, 0.4_dp, 0.9_dp, 1.3_dp]
        loading(1, 1) = 1.0_dp
        beta_mean = 0.0_dp
        beta_precision(1, 1) = 0.2_dp
        joint_scale = reshape([1.0_dp, 0.15_dp, 0.15_dp, 1.0_dp], [2, 2])
        initial_joint = reshape([1.0_dp, 0.25_dp, 0.25_dp, 0.8_dp], [2, 2])
        call rng_seed(state, 817246_8)
        call covu_gaussian_mixed_mcmc(y, x, loading, beta_mean, beta_precision, joint_scale, 4.0_dp, &
            64, 24, 4, state, result, info, initial_joint)
        call check(info == 0, 'covu Gaussian sampler status', failures)
        if (info /= 0) return
        call check(all(result%g(1, 1, :) > 0.0_dp), 'covu positive marginal G draws', failures)
        call check(all(result%r(1, 1, :) > 0.0_dp), 'covu positive conditional R draws', failures)
        call check(all(ieee_is_finite(result%regression)), 'covu finite random-residual regression draws', failures)
        call check(all(ieee_is_finite(result%log_likelihood)), 'covu finite conditional likelihoods', failures)
        do j = 1, size(result%joint_covariance, 3)
            call check_close(result%joint_covariance(2, 1, j), &
                result%regression(1, 1, j) * result%g(1, 1, j), 2.0e-12_dp, &
                'covu cross-covariance Schur identity', failures)
            call check_close(result%joint_covariance(2, 2, j), &
                result%r(1, 1, j) + result%regression(1, 1, j)**2 * result%g(1, 1, j), &
                2.0e-12_dp, 'covu trailing-covariance Schur identity', failures)
        end do
    end subroutine test_covu_sampler

    subroutine test_matrix_and_posterior(failures)
        integer, intent(inout) :: failures !! Running count of failed assertions.
        integer :: info
        real(dp), allocatable :: ante_coef(:, :)
        real(dp), allocatable :: ante_var(:)
        real(dp), allocatable :: correlations(:, :, :)
        real(dp), allocatable :: eigenvalues(:, :)
        real(dp), allocatable :: inverses(:, :, :)
        real(dp), allocatable :: modes(:)
        real(dp) :: mode_samples(10, 2)
        real(dp) :: samples(2, 2, 1)

        samples(:, :, 1) = reshape([4.0_dp, 2.0_dp, 2.0_dp, 9.0_dp], [2, 2])
        call posterior_correlations(samples, correlations, info)
        call check(info == 0, 'posterior correlation status', failures)
        call check_close(correlations(1, 2, 1), 1.0_dp / 3.0_dp, 1.0e-12_dp, &
            'posterior correlation value', failures)
        call posterior_inverses(samples, inverses, info)
        call check(info == 0, 'posterior inverse status', failures)
        call check(maxval(abs(matmul(samples(:, :, 1), inverses(:, :, 1)) - &
            reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]))) < 1.0e-11_dp, 'posterior inverse product', failures)
        call posterior_eigenvalues(samples, eigenvalues, info)
        call check(info == 0 .and. eigenvalues(1, 1) >= eigenvalues(2, 1), 'posterior eigenvalue ordering', failures)
        mode_samples(:, 1) = [-0.2_dp, -0.1_dp, -0.05_dp, 0.0_dp, 0.0_dp, &
            0.0_dp, 0.0_dp, 0.05_dp, 0.1_dp, 1.0_dp]
        mode_samples(:, 2) = 3.0_dp
        call posterior_modes(mode_samples, modes=modes, info=info)
        call check(info == 0, 'posterior mode KDE status', failures)
        if (info == 0) then
            call check(abs(modes(1)) < 0.02_dp, 'posterior mode concentrated sample', failures)
            call check(abs(modes(2) - 3.0_dp) < 0.01_dp, 'posterior mode constant sample', failures)
        end if
        call ante_parameters(reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]), 1, ante_var, ante_coef, info)
        call check(info == 0 .and. maxval(abs(ante_var - 1.0_dp)) < 1.0e-12_dp .and. &
            maxval(abs(ante_coef)) < 1.0e-12_dp, 'posterior ante identity case', failures)
    end subroutine test_matrix_and_posterior

    subroutine test_utilities(failures)
        integer, intent(inout) :: failures !! Running count of failed assertions.
        integer :: info
        real(dp) :: angles_expected
        real(dp), allocatable :: angles(:)
        real(dp), allocatable :: bisectors(:, :)
        real(dp), allocatable :: matrix_value(:, :)
        real(dp), allocatable :: moments(:)
        real(dp), allocatable :: normal_moments(:, :)
        real(dp), allocatable :: symmetrizer(:, :)
        real(dp), allocatable :: packed(:)
        real(dp) :: sum_s
        real(dp) :: covariance(3, 3)
        real(dp) :: data(3, 2)

        call triangle_to_matrix([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp], .true., .true., .true., &
            matrix_value, info)
        call check(info == 0 .and. abs(matrix_value(3, 2) - 5.0_dp) < 1.0e-12_dp .and. &
            abs(matrix_value(2, 3) - 5.0_dp) < 1.0e-12_dp, 'Tri2M lower mirrored reconstruction', failures)
        call matrix_to_triangle(matrix_value, .true., .true., packed, info)
        call check(info == 0 .and. maxval(abs(packed - [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp])) &
            < 1.0e-12_dp, 'Tri2M reverse packing', failures)
        call check_close(uniform_central_moment(-1.0_dp, 1.0_dp, 2), 1.0_dp / 3.0_dp, 1.0e-12_dp, &
            'uniform second central moment', failures)
        data = reshape([-1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp, 2.0_dp, 2.0_dp], [3, 2])
        call central_moment_tensor(data, 2, moments, info)
        call check(info == 0 .and. abs(moments(1) - 2.0_dp / 3.0_dp) < 1.0e-12_dp, &
            'Ptensor variance component', failures)
        covariance = 0.0_dp
        covariance(1, 1) = 3.0_dp
        covariance(2, 2) = 2.0_dp
        covariance(3, 3) = 1.0_dp
        call krzanowski_compare(covariance, covariance, [1], [1], .false., sum_s, angles, bisectors, info)
        angles_expected = 0.0_dp
        call check(info == 0, 'Krzanowski status', failures)
        call check_close(sum_s, 1.0_dp, 1.0e-10_dp, 'Krzanowski identical subspace sum', failures)
        call check_close(angles(1), angles_expected, 1.0e-8_dp, 'Krzanowski identical subspace angle', failures)
        call symmetrizer_matrix(2, 2, symmetrizer, info)
        call check(info == 0 .and. size(symmetrizer, 1) == 4, 'KPPM symmetrizer dimensions', failures)
        call check(maxval(abs(matmul(symmetrizer, symmetrizer) - symmetrizer)) < 1.0e-12_dp, &
            'KPPM symmetrizer idempotence', failures)
        call normal_moment_matrix(reshape([2.0_dp], [1, 1]), 4, normal_moments, info)
        call check(info == 0, 'knorm fourth-moment status', failures)
        call check_close(normal_moments(1, 1), 12.0_dp, 1.0e-12_dp, 'knorm univariate fourth moment', failures)
    end subroutine test_utilities

    subroutine test_gaussian_sampler(failures)
        integer, intent(inout) :: failures !! Running count of failed assertions.
        integer :: info
        integer :: i
        real(dp) :: beta_mean
        real(dp) :: ainv(4, 4)
        real(dp) :: beta_prior_mean(1, 1)
        real(dp) :: beta_prior_precision(1, 1)
        real(dp) :: gscale(1, 1)
        real(dp) :: rscale(1, 1)
        real(dp) :: x(8, 1)
        real(dp) :: y(8, 1)
        real(dp) :: z(8, 4)
        type(gaussian_mcmc_result) :: result
        type(rng_state) :: state

        x(:, 1) = 1.0_dp
        z = 0.0_dp
        do i = 1, 4
            z(2 * i - 1:2 * i, i) = 1.0_dp
        end do
        y(:, 1) = [1.0_dp, 1.2_dp, 1.7_dp, 1.9_dp, 2.2_dp, 2.4_dp, 2.8_dp, 3.0_dp]
        ainv = 0.0_dp
        do i = 1, 4
            ainv(i, i) = 1.0_dp
        end do
        beta_prior_mean = 0.0_dp
        beta_prior_precision = 0.01_dp
        gscale = 1.0_dp
        rscale = 1.0_dp
        call rng_seed(state, 987654_8)
        call gaussian_mixed_mcmc(y, x, z, ainv, beta_prior_mean, beta_prior_precision, gscale, 4.0_dp, rscale, 4.0_dp, &
            1500, 500, 5, state, result, info)
        call check(info == 0, 'Gaussian mixed-model Gibbs sampler status', failures)
        if (info /= 0) return
        call check(size(result%beta, 3) == 200, 'Gaussian sampler saved iteration count', failures)
        beta_mean = sum(result%beta(1, 1, :)) / real(size(result%beta, 3), dp)
        call check(beta_mean > 1.0_dp .and. beta_mean < 3.0_dp, 'Gaussian sampler posterior beta location', failures)
        call check(all(result%g(1, 1, :) > 0.0_dp) .and. all(result%r(1, 1, :) > 0.0_dp), &
            'Gaussian sampler positive covariance draws', failures)
        call check(all(ieee_is_finite(result%log_likelihood)), 'Gaussian sampler finite log likelihoods', failures)
    end subroutine test_gaussian_sampler


    subroutine test_ante_sampler(failures)
        integer, intent(inout) :: failures !! Running count of failed assertions.
        integer :: info
        real(dp), allocatable :: samples(:, :, :)
        real(dp) :: location(5, 3)
        type(rng_state) :: state

        location(:, 1) = [-1.0_dp, -0.5_dp, 0.0_dp, 0.5_dp, 1.0_dp]
        location(:, 2) = [-0.7_dp, -0.2_dp, 0.1_dp, 0.4_dp, 0.8_dp]
        location(:, 3) = [-0.5_dp, -0.1_dp, 0.2_dp, 0.5_dp, 0.7_dp]
        call rng_seed(state, 24680_8)
        call ante_covariance_samples(state, location, 1, 4, .true., .false., samples, info)
        call check(info == 0 .and. size(samples, 3) == 4, 'rante antedependence sample count', failures)
        if (info /= 0) return
        call check(all(ieee_is_finite(samples)) .and. all(samples(1, 1, :) > 0.0_dp), &
            'rante antedependence covariance samples finite', failures)
        call check(maxval(abs(samples(:, :, 1) - transpose(samples(:, :, 1)))) < 1.0e-10_dp, &
            'rante antedependence covariance symmetry', failures)
    end subroutine test_ante_sampler


    subroutine test_phylo(failures)
        integer, intent(inout) :: failures !! Running count of failed assertions.
        integer :: info
        integer :: i
        real(dp), allocatable :: inverse_relationship(:, :)
        real(dp), allocatable :: relationship(:, :)
        real(dp), allocatable :: values(:, :)
        real(dp) :: g(1, 1)
        type(rng_state) :: state
        type(phylo_tree) :: tree

        tree%n_tip = 4
        tree%n_node = 3
        allocate(tree%edge(6, 2), tree%edge_length(6))
        tree%edge = reshape([5, 5, 6, 6, 7, 7, 6, 7, 1, 2, 3, 4], [6, 2])
        tree%edge_length = 0.5_dp
        call phylogenetic_precision(tree, [1, 2, 3, 4], .true., inverse_relationship, relationship, info)
        call check(info == 0, 'phylogenetic precision status', failures)
        call check_close(relationship(1, 1), 1.0_dp, 1.0e-12_dp, 'phylogenetic precision tip variance', failures)
        call check_close(relationship(1, 2), 0.5_dp, 1.0e-12_dp, 'phylogenetic precision sister covariance', failures)
        call check_close(relationship(1, 3), 0.0_dp, 1.0e-12_dp, 'phylogenetic precision root covariance', failures)
        call check(maxval(abs(matmul(relationship, inverse_relationship) - &
            reshape([(merge(1.0_dp, 0.0_dp, modulo(i - 1, 5) == 0), i=1, 16)], [4, 4]))) < 1.0e-10_dp, &
            'phylogenetic precision inverse product', failures)
        g(1, 1) = 1.0_dp
        call rng_seed(state, 13579_8)
        call breeding_values_phylo(state, tree, g, .true., values, info)
        call check(info == 0 .and. size(values, 1) == 7, 'phylogenetic breeding value status', failures)
        call check(all(ieee_is_finite(values)) .and. abs(values(5, 1)) < 1.0e-15_dp, &
            'phylogenetic breeding values finite with zero root', failures)
    end subroutine test_phylo


    subroutine test_design_helpers(failures)
        integer, intent(inout) :: failures !! Running count of failed assertions.
        integer :: info
        integer, allocatable :: active(:)
        real(dp), allocatable :: combined(:, :)
        real(dp), allocatable :: covariance(:, :)
        real(dp), allocatable :: marginal(:, :)
        real(dp), allocatable :: matrix_value(:, :)
        real(dp), allocatable :: psi(:, :)
        real(dp) :: parts(2, 3, 2)
        real(dp) :: x1(4, 2)
        real(dp) :: x2(4, 2)
        real(dp) :: z(2, 2)
        real(dp) :: tcov(1, 1)
        real(dp) :: lcov(2, 2)
        real(dp) :: divergence
        type(rng_state) :: state

        call path_matrix([1, 2], [2, 3], 3, psi, info)
        call check(info == 0 .and. abs(psi(2, 1) - 1.0_dp) < 1.0e-12_dp .and. &
            abs(psi(3, 2) - 1.0_dp) < 1.0e-12_dp, 'path matrix construction', failures)
        parts = 0.0_dp
        parts(1, 1, 1) = 1.0_dp
        parts(2, 2, 2) = 1.0_dp
        call multiple_membership_design(parts, combined, active)
        call check(size(active) == 2 .and. all(active == [1, 2]), 'multiple-membership active columns', failures)
        x1(:, 1) = 1.0_dp
        x1(:, 2) = [-1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp]
        x2 = x1
        call gelman_prior_design(x1, x2, 2.0_dp, 1.0_dp, covariance, info)
        call check(info == 0 .and. covariance(1, 1) > 0.0_dp .and. covariance(2, 2) > 0.0_dp, &
            'Gelman prior design covariance', failures)
        z = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
        tcov(1, 1) = 2.0_dp
        lcov = reshape([1.0_dp, 0.25_dp, 0.25_dp, 1.0_dp], [2, 2])
        call random_effect_covariance(z, tcov, lcov, marginal, info)
        call check(info == 0 .and. maxval(abs(marginal - 2.0_dp * lcov)) < 1.0e-12_dp, &
            'buildV random-effect covariance kernel', failures)
        call rng_seed(state, 24680_8)
        call d_divergence_mc(state, lcov, lcov, 8, divergence, info)
        call check(info == 0, 'Ddivergence identical-covariance status', failures)
        call check_close(divergence, 0.0_dp, 1.0e-12_dp, 'Ddivergence identical covariances', failures)
        call sir_matrix(reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]), &
            reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], [2, 2]), .true., matrix_value, info)
        call check(info == 0, 'sir interaction matrix status', failures)
        if (info == 0) then
            call check_close(matrix_value(1, 1), 0.0_dp, 1.0e-14_dp, 'sir zero diagonal first entry', failures)
            call check_close(matrix_value(1, 2), 2.0_dp, 1.0e-14_dp, 'sir off-diagonal interaction', failures)
            call check_close(matrix_value(2, 1), 3.0_dp, 1.0e-14_dp, 'sir reverse interaction', failures)
        end if
    end subroutine test_design_helpers


    subroutine test_ordinal_sampler(failures)
        integer, intent(inout) :: failures !! Running count of failed assertions.
        integer :: info
        integer :: i
        integer :: y(8)
        real(dp) :: ainv(4, 4)
        real(dp) :: beta_mean(1)
        real(dp) :: beta_precision(1, 1)
        real(dp) :: cutpoints(3)
        real(dp) :: x(8, 1)
        real(dp) :: z(8, 4)
        type(ordinal_mcmc_result) :: result
        type(rng_state) :: state

        x(:, 1) = 1.0_dp
        z = 0.0_dp
        do i = 1, 4
            z(2 * i - 1:2 * i, i) = 1.0_dp
        end do
        ainv = 0.0_dp
        do i = 1, 4
            ainv(i, i) = 1.0_dp
        end do
        y = [1, 1, 1, 1, 2, 2, 2, 2]
        cutpoints = [-1.0e33_dp, 0.0_dp, 1.0e33_dp]
        beta_mean = 0.0_dp
        beta_precision(1, 1) = 0.1_dp
        call rng_seed(state, 112233_8)
        call ordinal_probit_mixed_mcmc(y, cutpoints, x, z, ainv, beta_mean, beta_precision, 1.0_dp, 4.0_dp, &
            800, 200, 4, state, result, info)
        call check(info == 0 .and. size(result%beta, 2) == 150, 'ordinal-probit sampler status/count', failures)
        if (info /= 0) return
        call check(all(result%last_liability(1:4) < 0.0_dp) .and. all(result%last_liability(5:8) > 0.0_dp), &
            'ordinal-probit liabilities respect threshold', failures)
        call check(all(result%g > 0.0_dp) .and. all(ieee_is_finite(result%log_likelihood)), &
            'ordinal-probit covariance and likelihood finite', failures)
    end subroutine test_ordinal_sampler

    subroutine test_spline(failures)
        integer, intent(inout) :: failures !! Running count of failed assertions.
        integer :: info
        real(dp), allocatable :: basis(:, :)
        real(dp), parameter :: expected(3, 2) = reshape([ &
            1.0_dp, 0.125_dp, 0.0_dp, 0.0_dp, 0.125_dp, 1.0_dp], [3, 2])

        call spline_lrtp([0.0_dp, 0.5_dp, 1.0_dp], basis, info, knots=[0.0_dp, 1.0_dp])
        call check(info == 0 .and. all(shape(basis) == [3, 2]), 'LRTP spline dimensions', failures)
        call check(maxval(abs(basis - expected)) < 1.0e-11_dp, 'LRTP spline two-knot regression', failures)
    end subroutine test_spline


    subroutine test_parameter_expansion(failures)
        integer, intent(inout) :: failures !! Running count of failed assertions.
        integer :: info
        real(dp) :: alpha_prior_precision(1, 1)
        real(dp) :: design(2, 1)
        real(dp) :: observation_precision(2, 2)
        real(dp) :: px_ainv(2, 2)
        real(dp) :: px_beta_mean(1, 1)
        real(dp) :: px_beta_precision(1, 1)
        real(dp) :: px_gscale(1, 1)
        real(dp) :: px_random_prediction(2, 1)
        real(dp) :: px_residual_base(2, 1)
        real(dp) :: px_r(1, 1)
        real(dp) :: px_rscale(1, 1)
        real(dp) :: px_x(4, 1)
        real(dp) :: px_y(4, 1)
        real(dp) :: px_z(4, 2)
        real(dp) :: prior_precision(1, 1)
        real(dp) :: working_covariance(2, 2)
        real(dp), allocatable :: alpha(:)
        real(dp), allocatable :: covariance(:, :)
        real(dp), allocatable :: mean_alpha(:)
        real(dp), allocatable :: scaled_effects(:, :)
        real(dp) :: effects(2, 2)
        type(gaussian_px_mcmc_result) :: px_result
        type(rng_state) :: state

        design(:, 1) = [1.0_dp, 2.0_dp]
        observation_precision = 0.0_dp
        observation_precision(1, 1) = 1.0_dp
        observation_precision(2, 2) = 1.0_dp
        prior_precision(1, 1) = 1.0_dp
        call rng_seed(state, 76121_8)
        call parameter_expansion_conditional(state, [2.0_dp, 4.0_dp], design, observation_precision, &
            [0.0_dp], prior_precision, alpha, mean_alpha, info)
        call check(info == 0, 'parameter-expansion conditional status', failures)
        if (info == 0) then
            call check_close(mean_alpha(1), 5.0_dp / 3.0_dp, 2.0e-13_dp, &
                'parameter-expansion conditional mean', failures)
            call check(ieee_is_finite(alpha(1)), 'parameter-expansion Gaussian draw finite', failures)
        end if

        effects = reshape([1.0_dp, 3.0_dp, 2.0_dp, 4.0_dp], [2, 2])
        call apply_parameter_expansion(effects, [2.0_dp, -1.0_dp], scaled_effects, info)
        call check(info == 0, 'parameter-expansion effect scaling status', failures)
        if (info == 0) then
            call check(maxval(abs(scaled_effects - reshape([2.0_dp, 6.0_dp, -2.0_dp, -4.0_dp], [2, 2]))) < &
                2.0e-13_dp, 'parameter-expanded effect transformation', failures)
        end if

        working_covariance = reshape([2.0_dp, 0.5_dp, 0.5_dp, 3.0_dp], [2, 2])
        call expanded_covariance(working_covariance, [2.0_dp, -1.0_dp], covariance, info)
        call check(info == 0, 'parameter-expansion covariance status', failures)
        if (info == 0) then
            call check(maxval(abs(covariance - reshape([8.0_dp, -1.0_dp, -1.0_dp, 3.0_dp], [2, 2]))) < &
                2.0e-13_dp, 'parameter-expanded covariance transformation', failures)
        end if

        px_residual_base(:, 1) = [3.0_dp, 4.0_dp]
        px_random_prediction(:, 1) = [1.0_dp, 2.0_dp]
        px_r(1, 1) = 2.0_dp
        alpha_prior_precision(1, 1) = 1.0_dp
        call rng_seed(state, 76122_8)
        call px_alpha_conditional(px_residual_base, px_random_prediction, px_r, [1.0_dp], &
            alpha_prior_precision, state, alpha, mean_alpha, info)
        call check(info == 0, 'integrated PX alpha conditional status', failures)
        if (info == 0) then
            call check_close(mean_alpha(1), 13.0_dp / 7.0_dp, 2.0e-13_dp, &
                'integrated PX alpha conditional mean', failures)
        end if

        px_x(:, 1) = 1.0_dp
        px_z = 0.0_dp
        px_z(1:2, 1) = 1.0_dp
        px_z(3:4, 2) = 1.0_dp
        px_y(:, 1) = [-0.4_dp, 0.2_dp, 0.8_dp, 0.5_dp]
        px_ainv = 0.0_dp
        px_ainv(1, 1) = 1.0_dp
        px_ainv(2, 2) = 1.0_dp
        px_beta_mean = 0.0_dp
        px_beta_precision(1, 1) = 0.1_dp
        px_gscale(1, 1) = 1.0_dp
        px_rscale(1, 1) = 1.0_dp
        call rng_seed(state, 76123_8)
        call gaussian_parameter_expanded_mcmc(px_y, px_x, px_z, px_ainv, px_beta_mean, px_beta_precision, &
            px_gscale, 4.0_dp, px_rscale, 4.0_dp, [1.0_dp], alpha_prior_precision, .true., &
            100, 40, 4, state, px_result, info)
        call check(info == 0, 'integrated Gaussian parameter-expansion sampler status', failures)
        if (info == 0) then
            call check(size(px_result%alpha, 2) == 15, 'integrated PX retained draw count', failures)
            call check(all(ieee_is_finite(px_result%alpha)), 'integrated PX alpha draws finite', failures)
            call check(all(px_result%g(1, 1, :) >= 0.0_dp), 'integrated PX expanded variances nonnegative', failures)
            call check(all(ieee_is_finite(px_result%log_likelihood)), &
                'integrated PX Gaussian finite likelihoods', failures)
        end if
    end subroutine test_parameter_expansion


    subroutine test_engine_features(failures)
        integer, intent(inout) :: failures !! Running count of failed assertions.
        integer :: category_info
        integer :: determinant_sign
        integer :: info
        integer, allocatable :: category(:)
        real(dp) :: accepted_weight
        real(dp) :: attempted_weight
        real(dp) :: conditional_mean
        real(dp) :: conditional_variance
        real(dp) :: log_jacobian
        real(dp) :: log_likelihood
        real(dp) :: lower
        real(dp) :: sample
        real(dp) :: theta_draw
        real(dp) :: upper
        real(dp) :: basis(2, 2, 1)
        real(dp) :: category_effect(2, 1, 2)
        real(dp) :: covariance(1, 1)
        real(dp) :: liability(2, 1)
        real(dp) :: base_mean(2, 1)
        real(dp) :: parameter(1)
        real(dp) :: prior_probability(1, 2)
        real(dp), allocatable :: adjusted_mean(:, :)
        real(dp), allocatable :: posterior_probability(:, :)
        real(dp), allocatable :: transformed(:, :)
        real(dp) :: y_struct(1, 2)
        type(rng_state) :: state

        call check_close(normal_quantile(0.975_dp), 1.959963984540054_dp, 2.0e-12_dp, &
            'standard-normal quantile regression', failures)
        call check_close(optimal_acceptance_ratio(1), 0.4384285714285714_dp, 1.0e-15_dp, &
            'MCMCglmm scalar adaptive-MH target', failures)
        accepted_weight = 2.0_dp
        attempted_weight = 4.0_dp
        call adaptive_mh_decay(.true., accepted_weight, attempted_weight)
        call check_close(accepted_weight, 2.8_dp, 1.0e-15_dp, 'adaptive accepted-weight decay', failures)
        call check_close(attempted_weight, 4.6_dp, 1.0e-15_dp, 'adaptive attempted-weight decay', failures)

        call rng_seed(state, 913001_8)
        call binary_slice_liability_update(state, 3, 1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
            20.0_dp, sample, log_likelihood, lower, upper, info)
        call check(info == 0 .and. sample > lower .and. sample < upper, 'binary multinomial slice interval', failures)
        call check_close(log_likelihood, -log(2.0_dp), 1.0e-14_dp, 'binary multinomial slice likelihood', failures)

        call rng_seed(state, 913002_8)
        call binary_slice_liability_update(state, 14, 2.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
            20.0_dp, sample, log_likelihood, lower, upper, info)
        call check(info == 0 .and. sample > lower .and. sample < upper, 'ordered-probit slice interval', failures)
        call check_close(log_likelihood, -log(2.0_dp), 1.0e-14_dp, 'ordered-probit slice likelihood', failures)

        call rng_seed(state, 913003_8)
        call binary_slice_liability_update(state, 22, 1.0_dp, 5.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
            20.0_dp, sample, log_likelihood, lower, upper, info)
        call check(info == 0 .and. sample > lower .and. sample < upper, 'nonzero-binomial slice interval', failures)
        call check_close(log_likelihood, log(31.0_dp / 32.0_dp), 1.0e-14_dp, &
            'nonzero-binomial slice likelihood', failures)

        call rng_seed(state, 913004_8)
        call theta_scale_conditional(state, 1.0_dp, [1.0_dp, -1.0_dp], [2.0_dp, 1.0_dp], 4.0_dp, &
            0.5_dp, 3.0_dp, theta_draw, conditional_mean, conditional_variance, info)
        call check(info == 0, 'theta-scale scalar conditional status', failures)
        call check_close(conditional_mean, 0.9375_dp, 1.0e-14_dp, 'theta-scale conditional mean', failures)
        call check_close(conditional_variance, 0.5_dp, 1.0e-14_dp, 'theta-scale conditional variance', failures)
        call check(ieee_is_finite(theta_draw), 'theta-scale conditional draw finite', failures)

        liability(:, 1) = 0.0_dp
        base_mean(:, 1) = 0.0_dp
        category_effect(:, 1, 1) = 0.0_dp
        category_effect(:, 1, 2) = 5.0_dp
        covariance(1, 1) = 1.0_dp
        prior_probability(1, :) = 0.5_dp
        call rng_seed(state, 913005_8)
        call categorical_measurement_error_update(state, liability, base_mean, category_effect, covariance, &
            prior_probability, [1, 1], category, posterior_probability, adjusted_mean, category_info)
        call check(category_info == 0, 'discrete Berkson measurement-error update status', failures)
        if (category_info == 0) then
            call check(posterior_probability(1, 1) > 0.999999_dp, 'Berkson posterior favors matching category', failures)
            call check(category(1) == 1, 'Berkson category draw follows concentrated posterior', failures)
            call check_close(adjusted_mean(1, 1), 0.0_dp, 1.0e-14_dp, 'Berkson adjusted predictor', failures)
        end if

        basis = 0.0_dp
        basis(1, 2, 1) = 1.0_dp
        parameter(1) = 0.2_dp
        y_struct(1, :) = [1.0_dp, 2.0_dp]
        call structural_transform(y_struct, basis, parameter, transformed, log_jacobian, determinant_sign, info)
        call check(info == 0 .and. determinant_sign == 1, 'structural Lambda transform status', failures)
        if (info == 0) then
            call check_close(transformed(1, 1), 0.6_dp, 1.0e-14_dp, 'structural transformed first trait', failures)
            call check_close(transformed(1, 2), 2.0_dp, 1.0e-14_dp, 'structural transformed second trait', failures)
            call check_close(log_jacobian, 0.0_dp, 1.0e-14_dp, 'unit-determinant structural Jacobian', failures)
        end if
    end subroutine test_engine_features

    subroutine test_theta_scale_sampler(failures)
        integer, intent(inout) :: failures !! Running count of failed assertions.
        integer :: info
        integer :: random_term(2)
        real(dp) :: ainv(2, 2)
        real(dp) :: beta_mean(1, 1)
        real(dp) :: beta_precision(1, 1)
        real(dp) :: g_df(1)
        real(dp) :: g_scale(1, 1, 1)
        real(dp) :: r_scale(1, 1)
        real(dp) :: x(6, 1)
        real(dp) :: x_scale(6, 1)
        real(dp) :: y(6, 1)
        real(dp) :: z(6, 2)
        real(dp) :: z_scale(6, 2)
        type(rng_state) :: state
        type(theta_scale_mcmc_result) :: result

        x(:, 1) = 1.0_dp
        x_scale(:, 1) = 1.0_dp
        z = 0.0_dp
        z(1:3, 1) = 1.0_dp
        z(4:6, 2) = 1.0_dp
        z_scale = 0.0_dp
        random_term = [1, 1]
        ainv = 0.0_dp
        ainv(1, 1) = 1.0_dp
        ainv(2, 2) = 1.0_dp
        beta_mean = 0.0_dp
        beta_precision(1, 1) = 0.2_dp
        g_scale(1, 1, 1) = 1.0_dp
        g_df(1) = 4.0_dp
        r_scale(1, 1) = 1.0_dp
        y(:, 1) = [-0.5_dp, -0.2_dp, 0.1_dp, 0.4_dp, 0.7_dp, 1.0_dp]
        call rng_seed(state, 913006_8)
        call theta_scale_gaussian_mixed_mcmc(y, x, z, x_scale, z_scale, random_term, ainv, beta_mean, &
            beta_precision, g_scale, g_df, r_scale, 4.0_dp, 1.0_dp, 1.0_dp, 60, 20, 4, state, result, info)
        call check(info == 0, 'theta-scale Gaussian MCMC status', failures)
        if (info == 0) then
            call check(size(result%theta_scale) == 10, 'theta-scale retained draw count', failures)
            call check(all(ieee_is_finite(result%theta_scale)), 'theta-scale retained values finite', failures)
            call check(all(result%g(1, 1, 1, :) > 0.0_dp), 'theta-scale positive G draws', failures)
            call check(all(result%r(1, 1, :) > 0.0_dp), 'theta-scale positive R draws', failures)
            call check(all(ieee_is_finite(result%log_likelihood)), 'theta-scale likelihoods finite', failures)
        end if
    end subroutine test_theta_scale_sampler

    subroutine test_structural_sampler(failures)
        integer, intent(inout) :: failures !! Running count of failed assertions.
        integer :: i
        integer :: info
        integer :: random_term(2)
        real(dp) :: ainv(2, 2)
        real(dp) :: basis(2, 2, 1)
        real(dp) :: beta_mean(1, 2)
        real(dp) :: beta_precision(2, 2)
        real(dp) :: g_df(1)
        real(dp) :: g_scale(2, 2, 1)
        real(dp) :: r_scale(2, 2)
        real(dp) :: structural_mean(1)
        real(dp) :: structural_precision(1, 1)
        real(dp) :: structural_sd(1)
        real(dp) :: x(8, 1)
        real(dp) :: y(8, 2)
        real(dp) :: z(8, 2)
        type(rng_state) :: state
        type(structural_gaussian_mcmc_result) :: result

        x(:, 1) = 1.0_dp
        z = 0.0_dp
        z(1:4, 1) = 1.0_dp
        z(5:8, 2) = 1.0_dp
        random_term = [1, 1]
        ainv = 0.0_dp
        do i = 1, 2
            ainv(i, i) = 1.0_dp
        end do
        y(:, 1) = [-0.8_dp, -0.4_dp, -0.1_dp, 0.2_dp, 0.5_dp, 0.7_dp, 1.0_dp, 1.3_dp]
        y(:, 2) = [-0.2_dp, -0.1_dp, 0.1_dp, 0.3_dp, 0.5_dp, 0.8_dp, 1.0_dp, 1.2_dp]
        beta_mean = 0.0_dp
        beta_precision = 0.0_dp
        beta_precision(1, 1) = 0.2_dp
        beta_precision(2, 2) = 0.2_dp
        g_scale = 0.0_dp
        g_scale(1, 1, 1) = 1.0_dp
        g_scale(2, 2, 1) = 1.0_dp
        g_df(1) = 5.0_dp
        r_scale = 0.0_dp
        r_scale(1, 1) = 1.0_dp
        r_scale(2, 2) = 1.0_dp
        basis = 0.0_dp
        basis(2, 1, 1) = 1.0_dp
        structural_mean = 0.0_dp
        structural_precision(1, 1) = 0.5_dp
        structural_sd(1) = 0.08_dp
        call rng_seed(state, 913007_8)
        call structural_gaussian_multi_term_mcmc(y, x, z, random_term, ainv, beta_mean, beta_precision, &
            g_scale, g_df, r_scale, 5.0_dp, basis, structural_mean, structural_precision, structural_sd, &
            80, 30, 5, state, result, info)
        call check(info == 0, 'structural path Gaussian sampler status', failures)
        if (info == 0) then
            call check(size(result%structural_parameter, 2) == 10, 'structural path retained draw count', failures)
            call check(all(ieee_is_finite(result%structural_parameter)), 'structural path parameters finite', failures)
            call check(result%structural_acceptance_rate > 0.0_dp .and. &
                result%structural_acceptance_rate < 1.0_dp, 'structural path nondegenerate MH acceptance', failures)
            call check(all(result%g(1, 1, 1, :) > 0.0_dp) .and. all(result%g(2, 2, 1, :) > 0.0_dp), &
                'structural path positive G diagonals', failures)
            call check(all(ieee_is_finite(result%log_likelihood)), 'structural path finite likelihoods', failures)
        end if
    end subroutine test_structural_sampler

    subroutine test_multi_term_sampler(failures)
        integer, intent(inout) :: failures !! Running count of failed assertions.
        integer :: info
        integer :: i
        integer :: random_term(4)
        real(dp) :: ainv(4, 4)
        real(dp) :: beta_mean(1, 1)
        real(dp) :: beta_precision(1, 1)
        real(dp) :: g_df(2)
        real(dp) :: g_scale(1, 1, 2)
        real(dp) :: r_scale(1, 1)
        real(dp) :: x(8, 1)
        real(dp) :: y(8, 1)
        real(dp) :: y2(8, 2)
        real(dp) :: z(8, 4)
        real(dp) :: beta_mean2(1, 2)
        real(dp) :: beta_precision2(2, 2)
        real(dp) :: g_scale2(2, 2, 2)
        real(dp) :: r_scale2(2, 2)
        type(multi_term_gaussian_mcmc_result) :: result
        type(multi_term_gaussian_mcmc_result) :: routed_result
        type(multi_term_gaussian_mcmc_result) :: fixed_result
        type(multi_term_gaussian_mcmc_result) :: ante_result
        type(rng_state) :: state
        logical :: update_g(2)
        real(dp) :: initial_g(1, 1, 2)
        real(dp), allocatable :: beta_dense(:, :)
        real(dp), allocatable :: beta_sparse(:, :)
        real(dp), allocatable :: beta_sparse_repeat(:, :)
        real(dp), allocatable :: random_dense(:, :)
        real(dp), allocatable :: random_sparse(:, :)
        real(dp), allocatable :: random_sparse_repeat(:, :)
        real(dp), allocatable :: mean_dense(:)
        real(dp), allocatable :: mean_sparse(:)
        type(mcmcglmm_sparse_matrix) :: sparse_x
        type(mcmcglmm_sparse_matrix) :: sparse_z
        type(rng_state) :: sparse_state

        x(:, 1) = 1.0_dp
        z = 0.0_dp
        z(1:4, 1) = 1.0_dp
        z(5:8, 2) = 1.0_dp
        z([1, 3, 5, 7], 3) = 1.0_dp
        z([2, 4, 6, 8], 4) = 1.0_dp
        random_term = [1, 1, 2, 2]
        ainv = 0.0_dp
        do i = 1, 4
            ainv(i, i) = 1.0_dp
        end do
        y(:, 1) = [-0.6_dp, -0.3_dp, -0.2_dp, 0.0_dp, 0.4_dp, 0.7_dp, 0.8_dp, 1.0_dp]
        beta_mean = 0.0_dp
        beta_precision(1, 1) = 0.1_dp
        g_scale = 1.0_dp
        g_df = 4.0_dp
        r_scale(1, 1) = 1.0_dp
        call sparse_from_dense(x, sparse_x, info=info)
        call sparse_from_dense(z, sparse_z, info=info)
        call rng_seed(state, 817232_8)
        call rng_seed(sparse_state, 817232_8)
        call multi_term_coefficient_conditional(y, x, z, random_term, ainv, g_scale, r_scale, beta_mean, &
            beta_precision, state, beta_dense, random_dense, info, mean_dense)
        call check(info == 0, 'dense multi-term coefficient conditional status', failures)
        call multi_term_coefficient_conditional_sparse(y, sparse_x, sparse_z, random_term, ainv, g_scale, r_scale, &
            beta_mean, beta_precision, sparse_state, beta_sparse, random_sparse, info, mean_sparse)
        call check(info == 0, 'sparse multi-term coefficient conditional status', failures)
        call rng_seed(sparse_state, 817232_8)
        call multi_term_coefficient_conditional_sparse(y, sparse_x, sparse_z, random_term, ainv, g_scale, r_scale, &
            beta_mean, beta_precision, sparse_state, beta_sparse_repeat, random_sparse_repeat, info)
        if (info == 0) then
            call check(maxval(abs(mean_sparse - mean_dense)) < 1.0e-12_dp, &
                'sparse coefficient precision preserves conditional mean', failures)
            call check(all(ieee_is_finite(beta_sparse)) .and. all(ieee_is_finite(random_sparse)), &
                'sparse coefficient conditional returns finite draws', failures)
            call check(maxval(abs(beta_sparse - beta_sparse_repeat)) < 1.0e-14_dp, &
                'sparse coefficient conditional fixed-effect reproducibility', failures)
            call check(maxval(abs(random_sparse - random_sparse_repeat)) < 1.0e-14_dp, &
                'sparse coefficient conditional random-effect reproducibility', failures)
        end if
        call rng_seed(state, 817233_8)
        call multi_term_gaussian_mixed_mcmc(y, x, z, random_term, ainv, beta_mean, beta_precision, &
            g_scale, g_df, r_scale, 4.0_dp, 100, 40, 4, state, result, info)
        call check(info == 0, 'multi-term Gaussian sampler status', failures)
        if (info /= 0) return
        call check(size(result%beta, 3) == 15, 'multi-term Gaussian retained draw count', failures)
        call check(size(result%g, 3) == 2, 'multi-term Gaussian independent G block count', failures)
        call check(all(result%g(1, 1, :, :) > 0.0_dp), 'multi-term Gaussian positive G draws', failures)
        call check(all(result%r(1, 1, :) > 0.0_dp), 'multi-term Gaussian positive R draws', failures)
        call check(all(ieee_is_finite(result%random_effects)), 'multi-term Gaussian random effects finite', failures)
        call check(all(ieee_is_finite(result%log_likelihood)), 'multi-term Gaussian likelihoods finite', failures)

        update_g = [.true., .false.]
        initial_g(1, 1, 1) = 1.0_dp
        initial_g(1, 1, 2) = 2.0_dp
        call rng_seed(state, 817235_8)
        call multi_term_gaussian_mixed_mcmc(y, x, z, random_term, ainv, beta_mean, beta_precision, &
            g_scale, g_df, r_scale, 4.0_dp, 60, 20, 4, state, fixed_result, info, update_g=update_g, &
            initial_g=initial_g)
        call check(info == 0, 'multi-term fixed-G sampler status', failures)
        if (info == 0) then
            call check(maxval(abs(fixed_result%g(1, 1, 2, :) - 2.0_dp)) < 1.0e-14_dp, &
                'fixed G block remains unchanged for mev-style term', failures)
            call check(all(fixed_result%g(1, 1, 1, :) > 0.0_dp), 'updated G block remains positive', failures)
        end if

        y2(:, 1) = y(:, 1)
        y2(:, 2) = [0.8_dp, 0.6_dp, 0.4_dp, 0.2_dp, 0.0_dp, -0.2_dp, -0.4_dp, -0.6_dp]
        beta_mean2 = 0.0_dp
        beta_precision2 = 0.0_dp
        beta_precision2(1, 1) = 0.1_dp
        beta_precision2(2, 2) = 0.1_dp
        g_scale2 = 0.0_dp
        g_scale2(1, 1, :) = 1.0_dp
        g_scale2(2, 2, :) = 1.0_dp
        r_scale2 = 0.0_dp
        r_scale2(1, 1) = 1.0_dp
        r_scale2(2, 2) = 1.0_dp
        call rng_seed(state, 817237_8)
        call rng_seed(sparse_state, 817237_8)
        call multi_term_coefficient_conditional(y2, x, z, random_term, ainv, g_scale2, r_scale2, beta_mean2, &
            beta_precision2, state, beta_dense, random_dense, info, mean_dense)
        call multi_term_coefficient_conditional_sparse(y2, sparse_x, sparse_z, random_term, ainv, g_scale2, &
            r_scale2, beta_mean2, beta_precision2, sparse_state, beta_sparse, random_sparse, info, mean_sparse)
        call check(info == 0 .and. maxval(abs(mean_sparse - mean_dense)) < 1.0e-12_dp, &
            'sparse multivariate coefficient precision preserves conditional mean', failures)
        call rng_seed(state, 817238_8)
        call multi_term_gaussian_mixed_mcmc(y2, x, z, random_term, ainv, beta_mean2, beta_precision2, &
            g_scale2, [5.0_dp, 5.0_dp], r_scale2, 5.0_dp, 60, 20, 4, state, routed_result, info, &
            r_update_mode=6, r_split=1)
        call check(info == 0, 'multi-term routed direct-sum R sampler status', failures)
        if (info == 0) then
            call check(maxval(abs(routed_result%r(2, 2, :) - 1.0_dp)) < 1.0e-14_dp, &
                'routed direct-sum R fixes identity residual coordinate', failures)
            call check(maxval(abs(routed_result%r(1, 2, :))) < 1.0e-14_dp, &
                'routed direct-sum R has zero cross covariance', failures)
        end if
        call rng_seed(state, 817242_8)
        call multi_term_gaussian_mixed_mcmc(y2, x, z, random_term, ainv, beta_mean2, beta_precision2, &
            g_scale2, [5.0_dp, 5.0_dp], r_scale2, 5.0_dp, 48, 16, 4, state, ante_result, info, &
            r_update_mode=5, r_ante_order=1)
        call check(info == 0, 'multi-term antedependence R sampler status', failures)
        if (info == 0) then
            call check(all(ante_result%r(1, 1, :) > 0.0_dp) .and. all(ante_result%r(2, 2, :) > 0.0_dp), &
                'multi-term antedependence positive residual variances', failures)
            call check(all(ieee_is_finite(ante_result%r)), 'multi-term antedependence finite covariance draws', failures)
        end if
    end subroutine test_multi_term_sampler


    subroutine test_multi_term_px_sampler(failures)
        integer, intent(inout) :: failures !! Running count of failed assertions.
        integer :: family(1)
        integer :: info
        integer :: i
        integer :: random_term(4)
        logical :: observed(8, 1)
        real(dp), allocatable :: alpha_draw(:, :)
        real(dp), allocatable :: alpha_mean(:, :)
        real(dp) :: additional(8, 1)
        real(dp) :: additional2(8, 1)
        real(dp) :: ainv(4, 4)
        real(dp) :: alpha_prior_mean(1, 2)
        real(dp) :: alpha_prior_precision(2, 2)
        real(dp) :: beta_mean(1, 1)
        real(dp) :: beta_precision(1, 1)
        real(dp) :: g_df(2)
        real(dp) :: g_scale(1, 1, 2)
        real(dp) :: r_scale(1, 1)
        real(dp) :: residual_base(2, 1)
        real(dp) :: term_predictions(2, 1, 2)
        real(dp) :: x(8, 1)
        real(dp) :: y(8, 1)
        real(dp) :: y_family(8, 1)
        real(dp) :: z(8, 4)
        type(multi_term_family_px_mcmc_result) :: family_result
        type(multi_term_gaussian_px_mcmc_result) :: result
        type(rng_state) :: state

        residual_base(:, 1) = [3.0_dp, 3.0_dp]
        term_predictions = 0.0_dp
        term_predictions(1, 1, 1) = 1.0_dp
        term_predictions(2, 1, 2) = 1.0_dp
        alpha_prior_mean = 1.0_dp
        alpha_prior_precision = 0.0_dp
        alpha_prior_precision(1, 1) = 1.0_dp
        alpha_prior_precision(2, 2) = 1.0_dp
        call rng_seed(state, 817239_8)
        call multi_term_px_alpha_conditional(residual_base, term_predictions, reshape([1.0_dp], [1, 1]), &
            alpha_prior_mean, alpha_prior_precision, state, alpha_draw, alpha_mean, info)
        call check(info == 0, 'multi-term PX alpha conditional status', failures)
        if (info == 0) then
            call check(maxval(abs(alpha_mean - 2.0_dp)) < 1.0e-12_dp, &
                'multi-term PX alpha conditional mean', failures)
        end if

        x(:, 1) = 1.0_dp
        z = 0.0_dp
        z(1:4, 1) = 1.0_dp
        z(5:8, 2) = 1.0_dp
        z([1, 3, 5, 7], 3) = 1.0_dp
        z([2, 4, 6, 8], 4) = 1.0_dp
        random_term = [1, 1, 2, 2]
        ainv = 0.0_dp
        do i = 1, 4
            ainv(i, i) = 1.0_dp
        end do
        y(:, 1) = [-0.6_dp, -0.3_dp, -0.2_dp, 0.0_dp, 0.4_dp, 0.7_dp, 0.8_dp, 1.0_dp]
        beta_mean = 0.0_dp
        beta_precision(1, 1) = 0.1_dp
        g_scale = 1.0_dp
        g_df = 4.0_dp
        r_scale(1, 1) = 1.0_dp
        alpha_prior_mean = 1.0_dp
        alpha_prior_precision = 0.0_dp
        alpha_prior_precision(1, 1) = 2.0_dp
        alpha_prior_precision(2, 2) = 2.0_dp
        call rng_seed(state, 817240_8)
        call multi_term_gaussian_parameter_expanded_mcmc(y, x, z, random_term, ainv, beta_mean, &
            beta_precision, g_scale, g_df, r_scale, 4.0_dp, alpha_prior_mean, alpha_prior_precision, &
            .true., 60, 20, 4, state, result, info)
        call check(info == 0, 'multi-term Gaussian parameter-expanded sampler status', failures)
        if (info == 0) then
            call check(size(result%g, 3) == 2, 'multi-term PX independent G block count', failures)
            call check(all(result%g(1, 1, :, :) > 0.0_dp), 'multi-term PX positive expanded G draws', failures)
            call check(all(result%r(1, 1, :) > 0.0_dp), 'multi-term PX positive R draws', failures)
            call check(all(ieee_is_finite(result%alpha)), 'multi-term PX alpha draws finite', failures)
            call check(all(ieee_is_finite(result%log_likelihood)), 'multi-term PX likelihoods finite', failures)
        end if

        family = [2]
        y_family(:, 1) = [0.0_dp, 1.0_dp, 0.0_dp, 2.0_dp, 1.0_dp, 3.0_dp, 2.0_dp, 4.0_dp]
        additional = 0.0_dp
        additional2 = 0.0_dp
        observed = .true.
        observed(2, 1) = .false.
        y_family(2, 1) = -999.0_dp
        call rng_seed(state, 817241_8)
        call heterogeneous_multi_term_parameter_expanded_mcmc(family, y_family, additional, additional2, x, z, &
            random_term, ainv, beta_mean, beta_precision, g_scale, g_df, r_scale, 4.0_dp, alpha_prior_mean, &
            alpha_prior_precision, .true., 0.30_dp, 56, 20, 4, state, family_result, info, observed)
        call check(info == 0, 'heterogeneous multi-term parameter-expanded sampler status', failures)
        if (info == 0) then
            call check(size(family_result%g, 3) == 2, 'heterogeneous PX multiple G block count', failures)
            call check(all(family_result%g(1, 1, :, :) > 0.0_dp), 'heterogeneous PX positive G draws', failures)
            call check(all(ieee_is_finite(family_result%alpha)), 'heterogeneous PX alpha draws finite', failures)
            call check(ieee_is_finite(family_result%last_liability(2, 1)), &
                'heterogeneous PX missing response imputed', failures)
            call check(family_result%acceptance_rate > 0.0_dp .and. family_result%acceptance_rate < 1.0_dp, &
                'heterogeneous PX nondegenerate liability acceptance', failures)
        end if
    end subroutine test_multi_term_px_sampler

    subroutine test_unified_sampler(failures)
        integer, intent(inout) :: failures !! Running count of failed assertions.
        integer :: family(2)
        integer :: info
        integer :: i
        integer :: random_term(4)
        logical :: observed(8, 2)
        real(dp) :: additional(8, 2)
        real(dp) :: additional2(8, 2)
        real(dp) :: ainv(4, 4)
        real(dp) :: beta_mean(1, 2)
        real(dp) :: beta_precision(2, 2)
        real(dp) :: g_df(2)
        real(dp) :: g_scale(2, 2, 2)
        real(dp) :: r_scale(2, 2)
        real(dp) :: x(8, 1)
        real(dp) :: y(8, 2)
        real(dp) :: z(8, 4)
        type(unified_family_mcmc_result) :: result
        type(unified_family_mcmc_result) :: routed_result
        type(unified_family_mcmc_result) :: ante_result
        type(rng_state) :: state

        family = [1, 2]
        x(:, 1) = 1.0_dp
        z = 0.0_dp
        z(1:4, 1) = 1.0_dp
        z(5:8, 2) = 1.0_dp
        z([1, 3, 5, 7], 3) = 1.0_dp
        z([2, 4, 6, 8], 4) = 1.0_dp
        random_term = [1, 1, 2, 2]
        ainv = 0.0_dp
        do i = 1, 4
            ainv(i, i) = 1.0_dp
        end do
        y(:, 1) = [-0.4_dp, -0.1_dp, 0.1_dp, 0.2_dp, 0.5_dp, 0.7_dp, 0.8_dp, 1.1_dp]
        y(:, 2) = [0.0_dp, 1.0_dp, 0.0_dp, 2.0_dp, 1.0_dp, 3.0_dp, 2.0_dp, 4.0_dp]
        additional = 0.0_dp
        additional2 = 0.0_dp
        observed = .true.
        observed(2, 2) = .false.
        y(2, 2) = -999.0_dp
        beta_mean = 0.0_dp
        beta_precision = 0.0_dp
        beta_precision(1, 1) = 0.1_dp
        beta_precision(2, 2) = 0.1_dp
        g_scale = 0.0_dp
        g_scale(1, 1, :) = 1.0_dp
        g_scale(2, 2, :) = 1.0_dp
        g_df = 5.0_dp
        r_scale = 0.0_dp
        r_scale(1, 1) = 1.0_dp
        r_scale(2, 2) = 1.0_dp
        call rng_seed(state, 817234_8)
        call heterogeneous_multi_term_mixed_mcmc(family, y, additional, additional2, x, z, random_term, ainv, &
            beta_mean, beta_precision, g_scale, g_df, r_scale, 5.0_dp, .true., 0.35_dp, &
            80, 32, 4, state, result, info, observed)
        call check(info == 0, 'unified mixed-family multi-G sampler status', failures)
        if (info /= 0) return
        call check(size(result%beta, 3) == 12, 'unified sampler retained draw count', failures)
        call check(size(result%g, 3) == 2, 'unified sampler multiple G block count', failures)
        call check(all(result%g(1, 1, :, :) > 0.0_dp) .and. all(result%g(2, 2, :, :) > 0.0_dp), &
            'unified sampler positive G diagonals', failures)
        call check(all(result%r(1, 1, :) > 0.0_dp) .and. all(result%r(2, 2, :) > 0.0_dp), &
            'unified sampler positive R diagonals', failures)
        call check(ieee_is_finite(result%last_liability(2, 2)), 'unified sampler missing response imputed', failures)
        call check(result%acceptance_rate > 0.0_dp .and. result%acceptance_rate < 1.0_dp, &
            'unified sampler nondegenerate latent acceptance', failures)
        call check(all(ieee_is_finite(result%log_likelihood)), 'unified sampler likelihoods finite', failures)

        call rng_seed(state, 817235_8)
        call heterogeneous_multi_term_mixed_mcmc(family, y, additional, additional2, x, z, random_term, ainv, &
            beta_mean, beta_precision, g_scale, g_df, r_scale, 5.0_dp, .true., 0.35_dp, &
            40, 16, 4, state, routed_result, info, observed, r_update_mode=6, r_split=1)
        call check(info == 0, 'unified routed covariance sampler status', failures)
        if (info == 0) then
            call check(all(abs(routed_result%r(2, 2, :) - 1.0_dp) < 1.0e-12_dp), &
                'unified mode-6 residual identity block', failures)
            call check(all(abs(routed_result%r(1, 2, :)) < 1.0e-12_dp), &
                'unified mode-6 residual cross covariance', failures)
        end if
        call rng_seed(state, 817243_8)
        call heterogeneous_multi_term_mixed_mcmc(family, y, additional, additional2, x, z, random_term, ainv, &
            beta_mean, beta_precision, g_scale, g_df, r_scale, 5.0_dp, .true., 0.35_dp, &
            36, 12, 4, state, ante_result, info, observed, r_update_mode=5, r_ante_order=1)
        call check(info == 0, 'unified antedependence residual sampler status', failures)
        if (info == 0) then
            call check(all(ante_result%r(1, 1, :) > 0.0_dp) .and. all(ante_result%r(2, 2, :) > 0.0_dp), &
                'unified antedependence positive residual variances', failures)
        end if
    end subroutine test_unified_sampler

    subroutine test_grouped_multi_term_sampler(failures)
        integer, intent(inout) :: failures !! Running count of failed assertions.
        integer :: info
        integer :: i
        integer :: random_term(4)
        integer :: response(8, 3)
        integer :: trials(8)
        logical :: observed(8)
        logical :: update_g(2)
        real(dp) :: ainv(4, 4)
        real(dp) :: beta_mean(1, 2)
        real(dp) :: beta_precision(2, 2)
        real(dp) :: g_df(2)
        real(dp) :: g_scale(2, 2, 2)
        real(dp) :: initial_g(2, 2, 2)
        real(dp) :: r_scale(2, 2)
        real(dp) :: x(8, 1)
        real(dp) :: y(8)
        real(dp) :: z(8, 4)
        type(grouped_multi_term_mcmc_result) :: multinomial_result
        type(grouped_multi_term_mcmc_result) :: routed_multinomial_result
        type(grouped_multi_term_mcmc_result) :: routed_two_part_result
        type(grouped_multi_term_mcmc_result) :: ante_multinomial_result
        type(grouped_multi_term_mcmc_result) :: ante_two_part_result
        type(grouped_multi_term_mcmc_result) :: two_part_result
        type(rng_state) :: state

        x(:, 1) = 1.0_dp
        z = 0.0_dp
        z(1:4, 1) = 1.0_dp
        z(5:8, 2) = 1.0_dp
        z([1, 3, 5, 7], 3) = 1.0_dp
        z([2, 4, 6, 8], 4) = 1.0_dp
        random_term = [1, 1, 2, 2]
        ainv = 0.0_dp
        do i = 1, 4
            ainv(i, i) = 1.0_dp
        end do
        beta_mean = 0.0_dp
        beta_precision = 0.0_dp
        beta_precision(1, 1) = 0.1_dp
        beta_precision(2, 2) = 0.1_dp
        g_scale = 0.0_dp
        g_scale(1, 1, :) = 1.0_dp
        g_scale(2, 2, :) = 1.0_dp
        g_df = 5.0_dp
        r_scale = 0.0_dp
        r_scale(1, 1) = 1.0_dp
        r_scale(2, 2) = 1.0_dp
        y = [0.0_dp, 1.0_dp, 0.0_dp, 2.0_dp, 3.0_dp, 0.0_dp, 1.0_dp, 4.0_dp]
        trials = 1
        observed = .true.
        observed(3) = .false.
        call rng_seed(state, 817236_8)
        call two_part_multi_term_mixed_mcmc(11, y, trials, x, z, random_term, ainv, beta_mean, &
            beta_precision, g_scale, g_df, r_scale, 5.0_dp, .true., 0.30_dp, 60, 20, 4, state, &
            two_part_result, info, observed=observed)
        call check(info == 0, 'two-process multi-G sampler status', failures)
        if (info == 0) then
            call check(size(two_part_result%g, 3) == 2, 'two-process multi-G block count', failures)
            call check(all(two_part_result%g(1, 1, :, :) > 0.0_dp), 'two-process multi-G positive G draws', failures)
            call check(all(ieee_is_finite(two_part_result%last_liability)), &
                'two-process multi-G missing row imputed', failures)
            call check(all(ieee_is_finite(two_part_result%log_likelihood)), &
                'two-process multi-G likelihood finite', failures)
        end if
        call rng_seed(state, 817238_8)
        call two_part_multi_term_mixed_mcmc(11, y, trials, x, z, random_term, ainv, beta_mean, &
            beta_precision, g_scale, g_df, r_scale, 5.0_dp, .true., 0.30_dp, 36, 12, 4, state, &
            routed_two_part_result, info, observed=observed, r_update_mode=6, r_split=1)
        call check(info == 0, 'two-process routed covariance sampler status', failures)
        if (info == 0) then
            call check(all(abs(routed_two_part_result%r(2, 2, :) - 1.0_dp) < 1.0e-12_dp), &
                'two-process mode-6 residual identity block', failures)
        end if
        call rng_seed(state, 817244_8)
        call two_part_multi_term_mixed_mcmc(11, y, trials, x, z, random_term, ainv, beta_mean, &
            beta_precision, g_scale, g_df, r_scale, 5.0_dp, .true., 0.30_dp, 32, 12, 4, state, &
            ante_two_part_result, info, observed=observed, r_update_mode=5, r_ante_order=1)
        call check(info == 0, 'two-process antedependence residual sampler status', failures)
        if (info == 0) then
            call check(all(ante_two_part_result%r(1, 1, :) > 0.0_dp) .and. &
                all(ante_two_part_result%r(2, 2, :) > 0.0_dp), &
                'two-process antedependence positive residual variances', failures)
        end if

        response = 0
        response(1, 1) = 1
        response(2, 2) = 1
        response(3, 3) = 1
        response(4, 1) = 1
        response(5, 2) = 1
        response(6, 3) = 1
        response(7, 1) = 1
        response(8, 2) = 1
        observed = .true.
        observed(4) = .false.
        update_g = [.true., .false.]
        initial_g = 0.0_dp
        initial_g(1, 1, 1) = 1.0_dp
        initial_g(2, 2, 1) = 1.0_dp
        initial_g(1, 1, 2) = 2.0_dp
        initial_g(2, 2, 2) = 2.0_dp
        call rng_seed(state, 817237_8)
        call multinomial_multi_term_mixed_mcmc(3, response, x, z, random_term, ainv, beta_mean, &
            beta_precision, g_scale, g_df, r_scale, 5.0_dp, .true., 0.25_dp, 60, 20, 4, state, &
            multinomial_result, info, observed=observed, update_g=update_g, initial_g=initial_g)
        call check(info == 0, 'multinomial multi-G sampler status', failures)
        if (info == 0) then
            call check(size(multinomial_result%g, 3) == 2, 'multinomial multi-G block count', failures)
            call check(maxval(abs(multinomial_result%g(:, :, 2, :) - &
                spread(initial_g(:, :, 2), 3, size(multinomial_result%g, 4)))) < 1.0e-14_dp, &
                'multinomial fixed G block remains unchanged', failures)
            call check(all(ieee_is_finite(multinomial_result%last_liability)), &
                'multinomial multi-G missing group imputed', failures)
            call check(all(ieee_is_finite(multinomial_result%log_likelihood)), &
                'multinomial multi-G likelihood finite', failures)
        end if
        call rng_seed(state, 817239_8)
        call multinomial_multi_term_mixed_mcmc(3, response, x, z, random_term, ainv, beta_mean, &
            beta_precision, g_scale, g_df, r_scale, 5.0_dp, .true., 0.25_dp, 36, 12, 4, state, &
            routed_multinomial_result, info, observed=observed, r_update_mode=6, r_split=1)
        call check(info == 0, 'multinomial routed covariance sampler status', failures)
        if (info == 0) then
            call check(all(abs(routed_multinomial_result%r(2, 2, :) - 1.0_dp) < 1.0e-12_dp), &
                'multinomial mode-6 residual identity block', failures)
        end if
        call rng_seed(state, 817245_8)
        call multinomial_multi_term_mixed_mcmc(3, response, x, z, random_term, ainv, beta_mean, &
            beta_precision, g_scale, g_df, r_scale, 5.0_dp, .true., 0.25_dp, 32, 12, 4, state, &
            ante_multinomial_result, info, observed=observed, r_update_mode=5, r_ante_order=1)
        call check(info == 0, 'multinomial antedependence residual sampler status', failures)
        if (info == 0) then
            call check(all(ante_multinomial_result%r(1, 1, :) > 0.0_dp) .and. &
                all(ante_multinomial_result%r(2, 2, :) > 0.0_dp), &
                'multinomial antedependence positive residual variances', failures)
        end if
    end subroutine test_grouped_multi_term_sampler

    subroutine test_numeric_orchestrator(failures)
        integer, intent(inout) :: failures !! Running count of failed assertions.
        integer :: i
        integer :: info
        integer :: response(6, 3)
        integer :: direct_info
        integer :: family(2)
        integer :: random_term(2)
        integer :: trials(6)
        integer :: y_category(6)
        integer :: measurement_group(6)
        logical :: observed(6, 2)
        logical :: observed_rows(6)
        real(dp) :: additional(6, 2)
        real(dp) :: additional2(6, 2)
        real(dp) :: a_inverse(2, 2)
        real(dp) :: basis(2, 2, 1)
        real(dp) :: beta_mean(1, 2)
        real(dp) :: beta_precision(2, 2)
        real(dp) :: covu_y(6, 1)
        real(dp) :: g_df(2)
        real(dp) :: g_scale(2, 2, 2)
        real(dp) :: r_scale(2, 2)
        real(dp) :: joint_scale(2, 2)
        real(dp) :: measurement_effect(6, 2, 2)
        real(dp) :: measurement_prior(1, 2)
        real(dp) :: random_loading(1, 1)
        real(dp) :: x(6, 1)
        real(dp) :: x_scale(6, 1)
        real(dp) :: y(6, 2)
        real(dp) :: y_vector(6)
        real(dp) :: z(6, 2)
        real(dp) :: z_scale(6, 2)
        type(covu_gaussian_mcmc_result) :: direct_covu
        type(unified_family_mcmc_result) :: direct_result
        type(mcmcglmm_control) :: control
        type(mcmcglmm_numeric_model) :: model
        type(mcmcglmm_numeric_model) :: sparse_model
        type(mcmcglmm_numeric_prior) :: prior
        type(mcmcglmm_numeric_result) :: result
        type(mcmcglmm_numeric_result) :: sparse_result
        type(mcmcglmm_numeric_result) :: sparse_repeat_result
        type(ordinal_native_mcmc_result) :: direct_ordinal
        type(rng_state) :: direct_state
        type(rng_state) :: state
        type(structural_gaussian_mcmc_result) :: direct_structural
        type(threshold_cutpoint_mcmc_result) :: direct_threshold
        type(theta_scale_mcmc_result) :: direct_theta

        family = [1, 2]
        x(:, 1) = 1.0_dp
        z = 0.0_dp
        z(1:3, 1) = 1.0_dp
        z(4:6, 2) = 1.0_dp
        random_term = [1, 2]
        a_inverse = 0.0_dp
        do i = 1, 2
            a_inverse(i, i) = 1.0_dp
        end do
        y(:, 1) = [-0.5_dp, -0.1_dp, 0.2_dp, 0.4_dp, 0.8_dp, 1.0_dp]
        y(:, 2) = [0.0_dp, 1.0_dp, 0.0_dp, 2.0_dp, 1.0_dp, 3.0_dp]
        additional = 0.0_dp
        additional2 = 0.0_dp
        observed = .true.
        observed(2, 2) = .false.
        beta_mean = 0.0_dp
        beta_precision = 0.0_dp
        beta_precision(1, 1) = 0.1_dp
        beta_precision(2, 2) = 0.1_dp
        g_scale = 0.0_dp
        g_scale(1, 1, :) = 1.0_dp
        g_scale(2, 2, :) = 1.0_dp
        g_df = 5.0_dp
        r_scale = 0.0_dp
        r_scale(1, 1) = 1.0_dp
        r_scale(2, 2) = 1.0_dp

        model%engine = mcmcglmm_engine_scalar
        model%family = family
        model%y = y
        model%additional = additional
        model%additional2 = additional2
        model%x = x
        model%z = z
        model%random_term = random_term
        model%a_inverse = a_inverse
        model%observed = observed
        prior%beta_mean = beta_mean
        prior%beta_precision = beta_precision
        prior%g_scale = g_scale
        prior%g_df = g_df
        prior%r_scale = r_scale
        prior%r_df = 5.0_dp
        control%iterations = 24
        control%burn = 8
        control%thin = 4
        control%proposal_scale = 0.30_dp
        call rng_seed(state, 917301_8)
        call mcmcglmm_fit_numeric(model, prior, control, state, result, info)
        call check(info == 0, 'numeric orchestrator scalar route status', failures)
        if (info == 0) then
            call check(result%engine == mcmcglmm_engine_scalar, 'numeric orchestrator scalar result tag', failures)
            call check(size(result%scalar%g, 3) == 2, 'numeric orchestrator scalar multi-G result', failures)
        end if

        call rng_seed(direct_state, 917301_8)
        call heterogeneous_multi_term_mixed_mcmc(family, y, additional, additional2, x, z, random_term, a_inverse, &
            beta_mean, beta_precision, g_scale, g_df, r_scale, 5.0_dp, .true., 0.30_dp, 24, 8, 4, direct_state, &
            direct_result, direct_info, observed)
        call check(direct_info == 0, 'numeric orchestrator direct scalar reference status', failures)
        if (info == 0 .and. direct_info == 0) then
            call check(maxval(abs(result%scalar%beta - direct_result%beta)) < 1.0e-14_dp, &
                'numeric orchestrator scalar route preserves draws', failures)
            call check(maxval(abs(result%scalar%log_likelihood - direct_result%log_likelihood)) < 1.0e-14_dp, &
                'numeric orchestrator scalar route preserves likelihoods', failures)
        end if

        sparse_model = model
        call sparse_from_dense(x, sparse_model%sparse_x, info=info)
        call check(info == 0, 'numeric orchestrator sparse X construction', failures)
        call sparse_from_dense(z, sparse_model%sparse_z, info=info)
        call check(info == 0, 'numeric orchestrator sparse Z construction', failures)
        deallocate(sparse_model%x, sparse_model%z)
        call rng_seed(state, 917301_8)
        call mcmcglmm_fit_numeric(sparse_model, prior, control, state, sparse_result, info)
        call check(info == 0, 'numeric orchestrator sparse-design route status', failures)
        call rng_seed(state, 917301_8)
        call mcmcglmm_fit_numeric(sparse_model, prior, control, state, sparse_repeat_result, info)
        if (info == 0) then
            call check(all(ieee_is_finite(sparse_result%scalar%beta)) .and. &
                all(ieee_is_finite(sparse_result%scalar%log_likelihood)), &
                'numeric orchestrator sparse design returns finite chain', failures)
            call check(sparse_result%scalar%symbolic_analyses >= 1 .and. &
                sparse_result%scalar%symbolic_analyses < control%iterations, &
                'numeric orchestrator sparse design reuses symbolic analysis', failures)
            call check(maxval(abs(sparse_result%scalar%beta - sparse_repeat_result%scalar%beta)) < 1.0e-14_dp, &
                'numeric orchestrator sparse design reproducible draws', failures)
            call check(maxval(abs(sparse_result%scalar%log_likelihood - &
                sparse_repeat_result%scalar%log_likelihood)) < 1.0e-14_dp, &
                'numeric orchestrator sparse design reproducible likelihoods', failures)
        end if
        sparse_model%x = x
        call rng_seed(state, 917301_8)
        call mcmcglmm_fit_numeric(sparse_model, prior, control, state, sparse_result, info)
        call check(info == mcmcglmm_orchestrator_invalid_control, &
            'numeric orchestrator rejects simultaneous dense and sparse X', failures)

        allocate(prior%alpha_mean(2, 2), source=1.0_dp)
        allocate(prior%alpha_precision(4, 4), source=0.0_dp)
        do i = 1, 4
            prior%alpha_precision(i, i) = 0.2_dp
        end do
        model%engine = mcmcglmm_engine_scalar_px
        call rng_seed(state, 917302_8)
        call mcmcglmm_fit_numeric(model, prior, control, state, result, info)
        call check(info == 0, 'numeric orchestrator parameter-expanded scalar route status', failures)
        if (info == 0) then
            call check(result%engine == mcmcglmm_engine_scalar_px, &
                'numeric orchestrator parameter-expanded result tag', failures)
            call check(all(ieee_is_finite(result%scalar_px%alpha)), &
                'numeric orchestrator parameter-expanded alpha draws', failures)
        end if

        y_vector = [0.0_dp, 1.0_dp, 0.0_dp, 2.0_dp, 3.0_dp, 1.0_dp]
        trials = 1
        observed_rows = .true.
        observed_rows(3) = .false.
        model%engine = mcmcglmm_engine_two_process
        model%grouped_family = 11
        model%y_vector = y_vector
        model%trials = trials
        model%observed_rows = observed_rows
        call rng_seed(state, 917303_8)
        call mcmcglmm_fit_numeric(model, prior, control, state, result, info)
        call check(info == 0, 'numeric orchestrator two-process route status', failures)
        if (info == 0) then
            call check(result%engine == mcmcglmm_engine_two_process, &
                'numeric orchestrator two-process result tag', failures)
            call check(all(ieee_is_finite(result%grouped%last_liability)), &
                'numeric orchestrator two-process liabilities', failures)
        end if

        response = 0
        response(1, 1) = 1
        response(2, 2) = 1
        response(3, 3) = 1
        response(4, 1) = 1
        response(5, 2) = 1
        response(6, 3) = 1
        model%engine = mcmcglmm_engine_multinomial
        model%grouped_family = 3
        model%grouped_response = response
        call rng_seed(state, 917304_8)
        call mcmcglmm_fit_numeric(model, prior, control, state, result, info)
        call check(info == 0, 'numeric orchestrator multinomial route status', failures)
        if (info == 0) then
            call check(result%engine == mcmcglmm_engine_multinomial, &
                'numeric orchestrator multinomial result tag', failures)
            call check(all(ieee_is_finite(result%grouped%log_likelihood)), &
                'numeric orchestrator multinomial likelihoods', failures)
        end if

        y_category = [1, 1, 1, 2, 2, 2]
        model%engine = mcmcglmm_engine_ordinal
        model%y_category = y_category
        model%cutpoints = [-1.0e40_dp, 0.0_dp, 1.0e40_dp]
        observed_rows = .true.
        model%observed_rows = observed_rows
        prior%beta_mean = reshape([0.0_dp], [1, 1])
        prior%beta_precision = reshape([0.1_dp], [1, 1])
        prior%g_scale = reshape([1.0_dp], [1, 1, 1])
        prior%g_df = [3.0_dp]
        prior%r_scale = reshape([1.0_dp], [1, 1])
        prior%r_df = 3.0_dp
        control%proposal_scale = 0.25_dp
        control%cutpoint_proposal_sd = 0.10_dp
        control%slice_sampling = .true.
        call rng_seed(state, 917305_8)
        call mcmcglmm_fit_numeric(model, prior, control, state, result, info)
        call check(info == 0, 'numeric orchestrator ordinal route status', failures)
        call rng_seed(direct_state, 917305_8)
        call ordinal_native_mixed_mcmc(y_category, model%cutpoints, 0.10_dp, .true., x, z, a_inverse, &
            prior%beta_mean(:, 1), prior%beta_precision, 1.0_dp, 3.0_dp, 1.0_dp, 3.0_dp, .true., 0.25_dp, &
            24, 8, 4, direct_state, direct_ordinal, direct_info, observed_rows, .true., 1.0e6_dp)
        call check(direct_info == 0, 'numeric orchestrator direct ordinal reference status', failures)
        if (info == 0 .and. direct_info == 0) then
            call check(maxval(abs(result%ordinal%beta - direct_ordinal%beta)) < 1.0e-14_dp, &
                'numeric orchestrator ordinal route preserves draws', failures)
            call check(maxval(abs(result%ordinal%last_liability - direct_ordinal%last_liability)) < 1.0e-14_dp, &
                'numeric orchestrator ordinal route preserves liabilities', failures)
        end if

        y_category = [1, 1, 2, 2, 3, 3]
        model%engine = mcmcglmm_engine_threshold
        model%y_category = y_category
        model%cutpoints = [-1.0e40_dp, 0.0_dp, 1.0_dp, 1.0e40_dp]
        control%slice_sampling = .false.
        call rng_seed(state, 917306_8)
        call mcmcglmm_fit_numeric(model, prior, control, state, result, info)
        call check(info == 0, 'numeric orchestrator threshold route status', failures)
        call rng_seed(direct_state, 917306_8)
        call threshold_cutpoint_mixed_mcmc(y_category, model%cutpoints, 0.10_dp, .true., x, z, a_inverse, &
            prior%beta_mean(:, 1), prior%beta_precision, 1.0_dp, 3.0_dp, 24, 8, 4, direct_state, &
            direct_threshold, direct_info, observed_rows)
        call check(direct_info == 0, 'numeric orchestrator direct threshold reference status', failures)
        if (info == 0 .and. direct_info == 0) then
            call check(maxval(abs(result%threshold%beta - direct_threshold%beta)) < 1.0e-14_dp, &
                'numeric orchestrator threshold route preserves draws', failures)
            call check(maxval(abs(result%threshold%cutpoints - direct_threshold%cutpoints)) < 1.0e-14_dp, &
                'numeric orchestrator threshold route preserves cutpoints', failures)
        end if

        prior%beta_mean = beta_mean
        prior%beta_precision = beta_precision
        prior%g_scale = g_scale
        prior%g_df = g_df
        prior%r_scale = r_scale
        prior%r_df = 5.0_dp
        prior%theta_mean = 1.0_dp
        prior%theta_precision = 0.5_dp
        model%y = y
        x_scale = 0.0_dp
        x_scale(:, 1) = 1.0_dp
        z_scale = 0.0_dp
        z_scale(:, 1) = z(:, 1)
        model%x_scale = x_scale
        model%z_scale = z_scale
        model%engine = mcmcglmm_engine_theta_scale
        control%initial_theta = 1.0_dp
        call rng_seed(state, 917307_8)
        call mcmcglmm_fit_numeric(model, prior, control, state, result, info)
        call check(info == 0, 'numeric orchestrator theta-scale route status', failures)
        call rng_seed(direct_state, 917307_8)
        call theta_scale_gaussian_mixed_mcmc(y, x, z, x_scale, z_scale, random_term, a_inverse, beta_mean, &
            beta_precision, g_scale, g_df, r_scale, 5.0_dp, 1.0_dp, 0.5_dp, 24, 8, 4, direct_state, &
            direct_theta, direct_info, update_r=.true., initial_theta=1.0_dp)
        call check(direct_info == 0, 'numeric orchestrator direct theta-scale reference status', failures)
        if (info == 0 .and. direct_info == 0) then
            call check(maxval(abs(result%theta_scale%beta - direct_theta%beta)) < 1.0e-14_dp, &
                'numeric orchestrator theta-scale route preserves draws', failures)
            call check(maxval(abs(result%theta_scale%theta_scale - direct_theta%theta_scale)) < 1.0e-14_dp, &
                'numeric orchestrator theta-scale route preserves scales', failures)
        end if

        basis = 0.0_dp
        basis(2, 1, 1) = 1.0_dp
        model%structural_basis = basis
        prior%structural_mean = [0.0_dp]
        prior%structural_precision = reshape([0.5_dp], [1, 1])
        prior%structural_proposal_sd = [0.05_dp]
        model%engine = mcmcglmm_engine_structural
        call rng_seed(state, 917308_8)
        call mcmcglmm_fit_numeric(model, prior, control, state, result, info)
        call check(info == 0, 'numeric orchestrator structural route status', failures)
        call rng_seed(direct_state, 917308_8)
        call structural_gaussian_multi_term_mcmc(y, x, z, random_term, a_inverse, beta_mean, beta_precision, &
            g_scale, g_df, r_scale, 5.0_dp, basis, prior%structural_mean, prior%structural_precision, &
            prior%structural_proposal_sd, 24, 8, 4, direct_state, direct_structural, direct_info, &
            update_r=.true., initial_structural=prior%structural_mean)
        call check(direct_info == 0, 'numeric orchestrator direct structural reference status', failures)
        if (info == 0 .and. direct_info == 0) then
            call check(maxval(abs(result%structural%beta - direct_structural%beta)) < 1.0e-14_dp, &
                'numeric orchestrator structural route preserves draws', failures)
            call check(maxval(abs(result%structural%structural_parameter - &
                direct_structural%structural_parameter)) < 1.0e-14_dp, &
                'numeric orchestrator structural route preserves path parameters', failures)
        end if

        covu_y(:, 1) = y(:, 1)
        random_loading(1, 1) = 1.0_dp
        joint_scale = 0.0_dp
        joint_scale(1, 1) = 1.0_dp
        joint_scale(2, 2) = 1.0_dp
        model%engine = mcmcglmm_engine_covu
        model%y = covu_y
        model%random_loading = random_loading
        prior%beta_mean = reshape([0.0_dp], [1, 1])
        prior%beta_precision = reshape([0.1_dp], [1, 1])
        prior%joint_scale = joint_scale
        prior%joint_df = 4.0_dp
        control%joint_update_mode = 1
        call rng_seed(state, 917309_8)
        call mcmcglmm_fit_numeric(model, prior, control, state, result, info)
        call check(info == 0, 'numeric orchestrator covu route status', failures)
        call rng_seed(direct_state, 917309_8)
        call covu_gaussian_mixed_mcmc(covu_y, x, random_loading, prior%beta_mean, prior%beta_precision, &
            joint_scale, 4.0_dp, 24, 8, 4, direct_state, direct_covu, direct_info)
        call check(direct_info == 0, 'numeric orchestrator direct covu reference status', failures)
        if (info == 0 .and. direct_info == 0) then
            call check(maxval(abs(result%covu%beta - direct_covu%beta)) < 1.0e-14_dp, &
                'numeric orchestrator covu route preserves draws', failures)
            call check(maxval(abs(result%covu%joint_covariance - direct_covu%joint_covariance)) < 1.0e-14_dp, &
                'numeric orchestrator covu route preserves joint covariance', failures)
        end if

        family = [1, 1]
        model%engine = mcmcglmm_engine_scalar
        model%family = family
        model%y = y
        model%y = 0.0_dp
        model%additional = 0.0_dp
        model%additional2 = 0.0_dp
        model%observed = .true.
        model%random_term = random_term
        prior%beta_mean = beta_mean
        prior%beta_precision = beta_precision
        prior%g_scale = g_scale
        prior%g_df = g_df
        prior%r_scale = r_scale
        prior%r_df = 5.0_dp
        measurement_effect = 0.0_dp
        measurement_effect(:, 1, 2) = 5.0_dp
        measurement_prior = 0.5_dp
        measurement_group = 1
        model%measurement_category_effect = measurement_effect
        model%measurement_prior_probability = measurement_prior
        model%measurement_group = measurement_group
        call rng_seed(state, 917310_8)
        call mcmcglmm_fit_numeric(model, prior, control, state, result, info)
        call check(info == 0, 'numeric orchestrator measurement-error route status', failures)
        call rng_seed(direct_state, 917310_8)
        call heterogeneous_multi_term_mixed_mcmc(family, model%y, model%additional, model%additional2, x, z, &
            random_term, a_inverse, beta_mean, beta_precision, g_scale, g_df, r_scale, 5.0_dp, .true., &
            control%proposal_scale, 24, 8, 4, direct_state, direct_result, direct_info, observed=model%observed, &
            measurement_category_effect=measurement_effect, measurement_prior_probability=measurement_prior, &
            measurement_group=measurement_group)
        call check(direct_info == 0, 'numeric orchestrator direct measurement-error reference status', failures)
        if (info == 0 .and. direct_info == 0) then
            call check(all(result%scalar%measurement_category == direct_result%measurement_category), &
                'numeric orchestrator measurement-error route preserves categories', failures)
            call check(maxval(abs(result%scalar%measurement_probability - &
                direct_result%measurement_probability)) < 1.0e-14_dp, &
                'numeric orchestrator measurement-error route preserves probabilities', failures)
            call check(result%scalar%measurement_probability(1, 1) > 0.99_dp, &
                'integrated measurement-error posterior favors matching category', failures)
        end if
        sparse_model = model
        call sparse_from_dense(x, sparse_model%sparse_x, info=info)
        call sparse_from_dense(z, sparse_model%sparse_z, info=info)
        deallocate(sparse_model%x, sparse_model%z)
        call rng_seed(state, 917310_8)
        call mcmcglmm_fit_numeric(sparse_model, prior, control, state, sparse_result, info)
        call check(info == 0, 'sparse measurement-error orchestration status', failures)
        if (info == 0 .and. direct_info == 0) then
            call check(all(sparse_result%scalar%measurement_category == direct_result%measurement_category), &
                'sparse measurement-error route preserves categories', failures)
            call check(maxval(abs(sparse_result%scalar%measurement_probability - &
                direct_result%measurement_probability)) < 1.0e-12_dp, &
                'sparse measurement-error route preserves probabilities', failures)
        end if

        model%engine = 999
        call mcmcglmm_fit_numeric(model, prior, control, state, result, info)
        call check(info == mcmcglmm_orchestrator_invalid_engine, &
            'numeric orchestrator rejects unknown engine', failures)
    end subroutine test_numeric_orchestrator

    subroutine test_prediction(failures)
        integer, intent(inout) :: failures !! Running count of failed assertions.
        integer :: info
        integer :: random_term(2)
        logical :: include_g(2)
        logical :: marginalize_g(2)
        real(dp) :: a_inverse(2, 2)
        real(dp) :: beta_draws(1, 2, 2)
        real(dp) :: expectation
        real(dp), allocatable :: covariance(:, :)
        real(dp), allocatable :: diagonal(:, :)
        real(dp) :: g_matrix(2, 2, 2)
        real(dp), allocatable :: predictor(:, :, :)
        real(dp) :: r_matrix(2, 2)
        real(dp) :: random_draws(2, 2, 2)
        real(dp) :: x(2, 1)
        real(dp) :: z(2, 2)

        z = 0.0_dp
        z(1, 1) = 1.0_dp
        z(2, 2) = 1.0_dp
        random_term = [1, 2]
        a_inverse = 0.0_dp
        a_inverse(1, 1) = 1.0_dp
        a_inverse(2, 2) = 1.0_dp
        g_matrix = 0.0_dp
        g_matrix(:, :, 1) = reshape([1.0_dp, 0.2_dp, 0.2_dp, 2.0_dp], [2, 2])
        g_matrix(:, :, 2) = reshape([0.5_dp, 0.1_dp, 0.1_dp, 0.7_dp], [2, 2])
        r_matrix = reshape([0.3_dp, 0.05_dp, 0.05_dp, 0.4_dp], [2, 2])
        include_g = .true.
        call multi_term_build_v(z, random_term, a_inverse, g_matrix, r_matrix, include_g, &
            covariance, diagonal, info)
        call check(info == 0, 'buildV multi-term covariance status', failures)
        if (info == 0) then
            call check_close(diagonal(1, 1), 1.3_dp, 1.0e-13_dp, 'buildV first trait/level variance', failures)
            call check_close(diagonal(2, 1), 0.8_dp, 1.0e-13_dp, 'buildV second level variance', failures)
            call check_close(diagonal(1, 2), 2.4_dp, 1.0e-13_dp, 'buildV second trait first variance', failures)
            call check_close(covariance(1, 3), 0.25_dp, 1.0e-13_dp, 'buildV first cross-trait covariance', failures)
            call check_close(covariance(2, 4), 0.15_dp, 1.0e-13_dp, 'buildV second cross-trait covariance', failures)
            call check_close(covariance(1, 2), 0.0_dp, 1.0e-13_dp, 'buildV independent-level covariance', failures)
        end if

        x(:, 1) = 1.0_dp
        beta_draws(1, :, 1) = [1.0_dp, 10.0_dp]
        beta_draws(1, :, 2) = [2.0_dp, 20.0_dp]
        random_draws(:, :, 1) = reshape([0.5_dp, 2.0_dp, 1.5_dp, 3.0_dp], [2, 2])
        random_draws(:, :, 2) = reshape([1.0_dp, 4.0_dp, 2.0_dp, 6.0_dp], [2, 2])
        marginalize_g = [.false., .true.]
        call posterior_linear_predictor(x, z, random_term, beta_draws, random_draws, marginalize_g, predictor, info)
        call check(info == 0, 'posterior linear predictor status', failures)
        if (info == 0) then
            call check(maxval(abs(predictor(:, :, 1) - &
                reshape([1.5_dp, 1.0_dp, 11.5_dp, 10.0_dp], [2, 2]))) < 1.0e-13_dp, &
                'posterior predictor marginal random-term selection', failures)
        end if
        call scalar_response_expectation(2, 0.3_dp, 0.4_dp, 0.0_dp, 0.0_dp, expectation, info)
        call check(info == 0, 'Poisson posterior response expectation status', failures)
        call check_close(expectation, exp(0.5_dp), 2.0e-13_dp, 'Poisson latent-variance response mean', failures)
        call scalar_response_expectation(5, 0.3_dp, 0.4_dp, 0.0_dp, 0.0_dp, expectation, info)
        call check_close(expectation, exp(-0.1_dp), 2.0e-13_dp, 'exponential latent-variance response mean', failures)
        call scalar_response_expectation(3, 0.0_dp, 0.8_dp, 5.0_dp, 0.0_dp, expectation, info)
        call check(info == 0, 'binomial posterior response expectation status', failures)
        call check_close(expectation, 2.5_dp, 2.0e-13_dp, 'symmetric logistic-normal binomial mean', failures)
        call scalar_response_expectation(16, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, expectation, info)
        call check_close(expectation, 1.0_dp / (1.0_dp - exp(-1.0_dp)), 2.0e-13_dp, &
            'zero-truncated Poisson zero-variance response mean', failures)
        call scalar_response_expectation(22, 0.2_dp, 0.0_dp, 5.0_dp, 0.0_dp, expectation, info)
        call check_close(expectation, 1.0_dp - (1.0_dp - logistic(0.2_dp))**5, 2.0e-13_dp, &
            'nonzero-binomial zero-variance response mean', failures)
    end subroutine test_prediction

    subroutine test_simulation(failures)
        integer, intent(inout) :: failures !! Running count of failed assertions.
        integer :: category
        integer :: info
        integer, allocatable :: grouped_response(:)
        integer :: random_term(2)
        real(dp), allocatable :: latent(:, :)
        real(dp), allocatable :: latent_repeat(:, :)
        real(dp), allocatable :: random_effects(:, :)
        real(dp), allocatable :: random_effects_repeat(:, :)
        real(dp) :: a_inverse(2, 2)
        real(dp) :: beta(1, 2)
        real(dp) :: g_matrix(2, 2, 2)
        real(dp) :: r_matrix(2, 2)
        real(dp) :: response
        real(dp) :: x(4, 1)
        real(dp) :: z(4, 2)
        type(rng_state) :: state

        call rng_seed(state, 817241_8)
        call simulate_scalar_response(state, 1, 0.7_dp, 0.0_dp, 0.0_dp, response, info)
        call check(info == 0, 'Gaussian response simulation status', failures)
        call check_close(response, 0.7_dp, 1.0e-14_dp, 'Gaussian response keeps latent value', failures)
        call simulate_scalar_response(state, 2, 0.2_dp, 0.0_dp, 0.0_dp, response, info)
        call check(info == 0 .and. response >= 0.0_dp .and. &
            abs(response - real(nint(response), dp)) < 1.0e-14_dp, 'Poisson response simulation', failures)
        call simulate_scalar_response(state, 16, -1.0_dp, 0.0_dp, 0.0_dp, response, info)
        call check(info == 0 .and. response >= 1.0_dp, 'zero-truncated Poisson response simulation', failures)
        call simulate_scalar_response(state, 17, 0.3_dp, 0.0_dp, 0.0_dp, response, info)
        call check(info == 0 .and. response >= 0.0_dp, 'geometric response simulation', failures)
        call simulate_scalar_response(state, 22, 0.1_dp, 5.0_dp, 0.0_dp, response, info)
        call check(info == 0 .and. (nint(response) == 0 .or. nint(response) == 1), &
            'nonzero-binomial response simulation', failures)
        call simulate_scalar_response(state, 23, 0.2_dp, 0.5_dp, 5.0_dp, response, info)
        call check(info == 0 .and. ieee_is_finite(response), 'noncentral-t response simulation', failures)
        call simulate_scalar_response(state, 24, 0.2_dp, 0.5_dp, 5.0_dp, response, info)
        call check(info == 0 .and. ieee_is_finite(response), 'mean-shifted t response simulation', failures)

        call simulate_ordinal_response(state, 0.0_dp, [-1.0e40_dp, 0.0_dp, 1.0e40_dp], category, info)
        call check(info == 0 .and. (category == 1 .or. category == 2), 'ordinal response simulation', failures)
        call simulate_threshold_response(0.5_dp, [-1.0e40_dp, 0.0_dp, 1.0e40_dp], category, info)
        call check(info == 0 .and. category == 2, 'threshold response category transformation', failures)
        call simulate_two_part_response(state, 15, 0.2_dp, -0.3_dp, 1, response, info)
        call check(info == 0 .and. response >= 0.0_dp, 'hurdle-Poisson response simulation', failures)

        call simulate_multinomial_response(state, 3, [0.2_dp, -0.1_dp], 5, grouped_response, info)
        call check(info == 0 .and. size(grouped_response) == 3 .and. sum(grouped_response) == 5, &
            'multinomial response simulation', failures)
        call simulate_multinomial_response(state, 26, [0.1_dp, -0.2_dp, 0.3_dp], 1, grouped_response, info)
        call check(info == 0 .and. sum(grouped_response) > 0, 'ztmb response simulation', failures)

        x(:, 1) = 1.0_dp
        z = 0.0_dp
        z(1:2, 1) = 1.0_dp
        z(3:4, 2) = 1.0_dp
        random_term = [1, 2]
        a_inverse = 0.0_dp
        a_inverse(1, 1) = 1.0_dp
        a_inverse(2, 2) = 1.0_dp
        beta(1, :) = [0.5_dp, -0.5_dp]
        g_matrix = 0.0_dp
        g_matrix(1, 1, 1) = 1.0_dp
        g_matrix(2, 2, 1) = 0.5_dp
        g_matrix(1, 1, 2) = 0.7_dp
        g_matrix(2, 2, 2) = 1.2_dp
        r_matrix = 0.0_dp
        r_matrix(1, 1) = 0.3_dp
        r_matrix(2, 2) = 0.4_dp
        call rng_seed(state, 817242_8)
        call simulate_multi_term_gaussian_latent(state, x, z, random_term, a_inverse, beta, g_matrix, r_matrix, &
            latent, random_effects, info)
        call check(info == 0 .and. all(ieee_is_finite(latent)) .and. all(ieee_is_finite(random_effects)), &
            'multi-term Gaussian latent predictive simulation', failures)
        call rng_seed(state, 817242_8)
        call simulate_multi_term_gaussian_latent(state, x, z, random_term, a_inverse, beta, g_matrix, r_matrix, &
            latent_repeat, random_effects_repeat, info)
        call check(info == 0 .and. maxval(abs(latent_repeat - latent)) < 1.0e-14_dp .and. &
            maxval(abs(random_effects_repeat - random_effects)) < 1.0e-14_dp, &
            'multi-term predictive simulation reproducibility', failures)
    end subroutine test_simulation

    subroutine test_family_likelihoods(failures)
        integer, intent(inout) :: failures !! Running count of failed assertions.
        integer :: counts(3)
        integer :: binary_response(3)
        real(dp) :: eta(2)
        real(dp) :: eta_binary(3)
        real(dp) :: value

        call check_close(poisson_cdf(2, 1.5_dp), 0.8088468305380581_dp, 2.0e-13_dp, &
            'Poisson CDF regression', failures)
        call check_close(scalar_family_loglik(2, 3.0_dp, 0.4_dp, 0.0_dp, 0.0_dp), &
            -2.0835841668693256_dp, 2.0e-13_dp, 'native Poisson log likelihood', failures)
        call check_close(scalar_family_loglik(5, 1.2_dp, 0.4_dp, 0.0_dp, 0.0_dp), &
            -1.3901896371695242_dp, 2.0e-13_dp, 'native exponential parameterization', failures)
        call check_close(scalar_family_loglik(4, 1.2_dp, 0.4_dp, 0.0_dp, 0.0_dp), &
            -1.3901896371695242_dp, 2.0e-13_dp, 'shape-one Weibull native kernel', failures)
        call check_close(scalar_family_loglik(8, 0.4_dp, 0.2_dp, 1.6_dp, 0.0_dp), &
            scalar_family_loglik(9, 0.4_dp, 0.2_dp, 1.6_dp, 0.0_dp), 2.0e-13_dp, &
            'shape-one censored Weibull kernel equals censored exponential', failures)
        call check_close(scalar_family_loglik(16, 2.0_dp, -0.2_dp, 0.0_dp, 0.0_dp), &
            -1.3302881814512904_dp, 2.0e-13_dp, 'zero-truncated Poisson likelihood', failures)
        call check_close(scalar_family_loglik(17, 4.0_dp, 0.3_dp, 0.0_dp, 0.0_dp), &
            -3.971776222342636_dp, 2.0e-13_dp, 'native geometric link parameterization', failures)
        call check_close(scalar_family_loglik(22, 1.0_dp, 0.3_dp, 5.0_dp, 0.0_dp), &
            -0.014055284805030393_dp, 2.0e-13_dp, 'nonzero-binomial positive likelihood', failures)
        call check_close(scalar_family_loglik(24, 1.2_dp, 0.4_dp, 0.5_dp, 5.0_dp), &
            -1.5157722417668031_dp, 2.0e-13_dp, 'mean-shifted scaled-t likelihood', failures)
        call check_close(noncentral_t_logpdf(1.2_dp, 5.0_dp, 0.7_dp), -1.1682336333823042_dp, &
            2.0e-10_dp, 'noncentral scaled-t density versus SciPy reference', failures)
        call check_close(ordinal_probit_loglik(1, 0.0_dp, [-1.0e40_dp, 0.0_dp, 1.0_dp, 1.0e40_dp]), &
            -0.6931471805599453_dp, 2.0e-13_dp, 'ordered-probit first category likelihood', failures)
        call check_close(ordinal_probit_loglik(2, 0.0_dp, [-1.0e40_dp, 0.0_dp, 1.0_dp, 1.0e40_dp]), &
            -1.0748623268620714_dp, 2.0e-13_dp, 'ordered-probit interior category likelihood', failures)
        call check_close(two_part_family_loglik(11, 0.0_dp, 0.2_dp, -0.3_dp, 1), &
            -0.5193409221213295_dp, 2.0e-13_dp, 'zero-inflated Poisson zero likelihood', failures)
        call check_close(two_part_family_loglik(15, 2.0_dp, 0.2_dp, -0.3_dp, 1), &
            -1.7196082110574753_dp, 2.0e-13_dp, 'hurdle Poisson positive likelihood', failures)
        call check_close(two_part_family_loglik(18, 0.0_dp, 0.2_dp, -0.3_dp, 1), &
            -0.7408182206817179_dp, 2.0e-13_dp, 'zero-altered Poisson zero likelihood', failures)
        call check_close(two_part_family_loglik(19, 2.0_dp, 0.2_dp, -0.3_dp, 5), &
            -1.8424644983824408_dp, 2.0e-13_dp, 'zero-inflated binomial likelihood', failures)
        call check_close(two_part_family_loglik(25, 2.0_dp, 0.2_dp, -0.3_dp, 5), &
            -1.823804607088509_dp, 2.0e-13_dp, 'hurdle binomial likelihood', failures)
        counts = [2, 3, 5]
        eta = [0.2_dp, -0.1_dp]
        call check_close(multinomial_log_kernel(counts, eta), -11.2983106084446_dp, 2.0e-13_dp, &
            'multinomial native log kernel', failures)
        binary_response = [1, 0, 1]
        eta_binary = [0.1_dp, -0.2_dp, 0.3_dp]
        value = ztmb_log_kernel(binary_response, eta_binary)
        call check(ieee_is_finite(value), 'ztmb likelihood finite', failures)
        value = ztmultinomial_log_kernel(counts, eta)
        call check(ieee_is_finite(value), 'zero-truncated multinomial likelihood finite', failures)
    end subroutine test_family_likelihoods

    subroutine test_family_samplers(failures)
        integer, intent(inout) :: failures !! Running count of failed assertions.
        integer :: i
        integer :: info
        integer :: family_codes(2)
        integer :: counts(6, 3)
        integer :: threshold_y(9)
        integer :: ordinal_binary_y(6)
        integer :: trials(6)
        logical :: observed_group(6)
        logical :: observed_multi(6, 2)
        logical :: observed_threshold(9)
        logical :: observed_two(6)
        real(dp) :: additional(6)
        real(dp) :: additional2(6)
        real(dp) :: ainv(2, 2)
        real(dp) :: beta_mean_scalar(1)
        real(dp) :: beta_precision_scalar(1, 1)
        real(dp) :: beta_mean_two(1, 2)
        real(dp) :: beta_precision_two(2, 2)
        real(dp) :: gscale_two(2, 2)
        real(dp) :: rscale_two(2, 2)
        real(dp) :: x(6, 1)
        real(dp) :: x_threshold(9, 1)
        real(dp) :: y(6)
        real(dp) :: y_multi(6, 2)
        real(dp) :: auxiliary_multi(6, 2)
        real(dp) :: auxiliary2_multi(6, 2)
        real(dp) :: z(6, 2)
        real(dp) :: z_threshold(9, 2)
        type(family_mcmc_result) :: scalar_result
        type(multivariate_family_mcmc_result) :: group_result
        type(multivariate_family_mcmc_result) :: heterogeneous_result
        type(multivariate_family_mcmc_result) :: two_result
        type(ordinal_native_mcmc_result) :: ordinal_native_result
        type(rng_state) :: state
        type(threshold_cutpoint_mcmc_result) :: cutpoint_result

        x(:, 1) = 1.0_dp
        z = 0.0_dp
        z(1:3, 1) = 1.0_dp
        z(4:6, 2) = 1.0_dp
        ainv = 0.0_dp
        ainv(1, 1) = 1.0_dp
        ainv(2, 2) = 1.0_dp
        beta_mean_scalar = 0.0_dp
        beta_precision_scalar(1, 1) = 0.1_dp
        y = [0.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, 3.0_dp, 2.0_dp]
        additional = 1.0_dp
        additional2 = 1.0_dp
        call rng_seed(state, 81421_8)
        call latent_family_mixed_mcmc(2, y, additional, additional2, x, z, ainv, beta_mean_scalar, &
            beta_precision_scalar, 1.0_dp, 4.0_dp, 1.0_dp, 4.0_dp, .false., 0.7_dp, &
            150, 50, 5, state, scalar_result, info)
        call check(info == 0, 'Poisson latent-family sampler status', failures)
        if (info == 0) then
            call check(size(scalar_result%beta, 2) == 20, 'Poisson latent-family retained draws', failures)
            call check(scalar_result%acceptance_rate > 0.0_dp .and. scalar_result%acceptance_rate < 1.0_dp, &
                'Poisson latent-family MH acceptance', failures)
            call check(all(ieee_is_finite(scalar_result%log_likelihood)), &
                'Poisson latent-family finite likelihoods', failures)
        end if

        y = [0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp]
        additional = 1.0_dp
        call rng_seed(state, 81423_8)
        call latent_family_mixed_mcmc(3, y, additional, additional2, x, z, ainv, beta_mean_scalar, &
            beta_precision_scalar, 1.0_dp, 4.0_dp, 1.0_dp, 4.0_dp, .false., 0.7_dp, &
            80, 30, 5, state, scalar_result, info, slice_sampling=.true.)
        call check(info == 0, 'binary multinomial native-slice sampler status', failures)
        if (info == 0) then
            call check(abs(scalar_result%acceptance_rate) < tiny(1.0_dp), 'slice sampler bypasses MH proposals', failures)
            call check(all(ieee_is_finite(scalar_result%last_liability)), 'slice liabilities finite', failures)
        end if

        y = [0.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, 3.0_dp, 2.0_dp]
        call rng_seed(state, 81424_8)
        call latent_family_mixed_mcmc(2, y, additional, additional2, x, z, ainv, beta_mean_scalar, &
            beta_precision_scalar, 1.0_dp, 4.0_dp, 1.0_dp, 4.0_dp, .false., 0.7_dp, &
            80, 30, 5, state, scalar_result, info, adaptive_mh=.true.)
        call check(info == 0, 'adaptive Poisson latent-family sampler status', failures)
        if (info == 0) then
            call check(scalar_result%final_proposal_sd > 0.0_dp .and. &
                ieee_is_finite(scalar_result%final_proposal_sd), 'adaptive sampler finite proposal scale', failures)
        end if

        beta_mean_two = 0.0_dp
        beta_precision_two = 0.0_dp
        beta_precision_two(1, 1) = 0.1_dp
        beta_precision_two(2, 2) = 0.1_dp
        gscale_two = 0.0_dp
        rscale_two = 0.0_dp
        do i = 1, 2
            gscale_two(i, i) = 1.0_dp
            rscale_two(i, i) = 1.0_dp
        end do
        trials = 1
        y = [0.0_dp, 999.0_dp, 0.0_dp, 2.0_dp, 3.0_dp, 0.0_dp]
        observed_two = .true.
        observed_two(2) = .false.
        call rng_seed(state, 91422_8)
        call two_part_mixed_mcmc(11, y, trials, x, z, ainv, beta_mean_two, beta_precision_two, &
            gscale_two, 5.0_dp, rscale_two, 5.0_dp, .false., 0.55_dp, 150, 50, 5, state, two_result, info, &
            observed=observed_two)
        call check(info == 0, 'zero-inflated Poisson two-process sampler status', failures)
        if (info == 0) then
            call check(two_result%acceptance_rate > 0.0_dp .and. two_result%acceptance_rate < 1.0_dp, &
                'zero-inflated Poisson MH acceptance', failures)
            call check(all(ieee_is_finite(two_result%log_likelihood)), &
                'zero-inflated Poisson finite likelihoods', failures)
        end if

        counts = reshape([2, 1, 0, 3, 1, 2, 1, 2, 2, 0, 3, 1, 2, 1, 3, 1, 1, 2], [6, 3])
        observed_group = .true.
        observed_group(3) = .false.
        counts(3, :) = -999
        call rng_seed(state, 101423_8)
        call multinomial_family_mixed_mcmc(3, counts, x, z, ainv, beta_mean_two, beta_precision_two, &
            gscale_two, 5.0_dp, rscale_two, 5.0_dp, .false., 0.5_dp, 150, 50, 5, state, group_result, info, &
            observed=observed_group)
        call check(info == 0, 'multinomial latent-family sampler status', failures)
        if (info == 0) then
            call check(group_result%acceptance_rate > 0.0_dp .and. group_result%acceptance_rate < 1.0_dp, &
                'multinomial latent-family MH acceptance', failures)
            call check(all(ieee_is_finite(group_result%log_likelihood)), &
                'multinomial latent-family finite likelihoods', failures)
        end if

        family_codes = [1, 2]
        y_multi(:, 1) = [-0.4_dp, 999.0_dp, 0.5_dp, -0.1_dp, 0.7_dp, 0.3_dp]
        y_multi(:, 2) = [0.0_dp, 1.0_dp, 2.0_dp, 999.0_dp, 3.0_dp, 2.0_dp]
        auxiliary_multi = 1.0_dp
        auxiliary2_multi = 1.0_dp
        observed_multi = .true.
        observed_multi(2, 1) = .false.
        observed_multi(4, 2) = .false.
        call rng_seed(state, 106424_8)
        call heterogeneous_family_mixed_mcmc(family_codes, y_multi, auxiliary_multi, auxiliary2_multi, x, z, ainv, &
            beta_mean_two, beta_precision_two, gscale_two, 5.0_dp, rscale_two, 5.0_dp, .true., 0.55_dp, &
            120, 40, 4, state, heterogeneous_result, info, observed=observed_multi)
        call check(info == 0, 'heterogeneous Gaussian-Poisson sampler status', failures)
        if (info == 0) then
            call check(heterogeneous_result%acceptance_rate > 0.0_dp .and. &
                heterogeneous_result%acceptance_rate < 1.0_dp, 'heterogeneous non-Gaussian MH acceptance', failures)
            call check(abs(heterogeneous_result%last_liability(1, 1) - y_multi(1, 1)) < 1.0e-14_dp .and. &
                abs(heterogeneous_result%last_liability(3, 1) - y_multi(3, 1)) < 1.0e-14_dp .and. &
                abs(heterogeneous_result%last_liability(6, 1) - y_multi(6, 1)) < 1.0e-14_dp, &
                'heterogeneous observed Gaussian liabilities remain fixed', failures)
            call check(abs(heterogeneous_result%last_liability(2, 1) - 999.0_dp) > 1.0_dp, &
                'heterogeneous missing Gaussian response is imputed', failures)
            call check(all(ieee_is_finite(heterogeneous_result%log_likelihood)), &
                'heterogeneous sampler finite complete-data likelihoods', failures)
        end if

        x_threshold(:, 1) = 1.0_dp
        z_threshold = 0.0_dp
        z_threshold(1:5, 1) = 1.0_dp
        z_threshold(6:9, 2) = 1.0_dp
        threshold_y = [1, 1, 0, 2, 2, 2, 3, 3, 3]
        observed_threshold = .true.
        observed_threshold(3) = .false.
        call rng_seed(state, 111424_8)
        call threshold_cutpoint_mixed_mcmc(threshold_y, [-1.0e40_dp, 0.0_dp, 1.0_dp, 1.0e40_dp], &
            0.25_dp, .true., x_threshold, z_threshold, ainv, beta_mean_scalar, beta_precision_scalar, &
            1.0_dp, 4.0_dp, 160, 60, 5, state, cutpoint_result, info, observed=observed_threshold)
        call check(info == 0, 'adaptive threshold cutpoint sampler status', failures)
        if (info == 0) then
            call check(cutpoint_result%cutpoint_acceptance_rate > 0.0_dp .and. &
                cutpoint_result%cutpoint_acceptance_rate < 1.0_dp, 'threshold cutpoint MH acceptance', failures)
            call check(all(cutpoint_result%cutpoints(3, :) > 0.0_dp), &
                'threshold sampled free cutpoint stays ordered', failures)
            call check(maxval(abs(cutpoint_result%cutpoints(1, :) + 1.0e40_dp)) < tiny(1.0_dp) .and. &
                maxval(abs(cutpoint_result%cutpoints(2, :))) < tiny(1.0_dp) .and. &
                maxval(abs(cutpoint_result%cutpoints(4, :) - 1.0e40_dp)) < tiny(1.0_dp), &
                'threshold fixed cutpoints remain fixed', failures)
            call check(ieee_is_finite(cutpoint_result%last_liability(3)), &
                'threshold missing response liability is imputed', failures)
        end if

        call rng_seed(state, 121425_8)
        call ordinal_native_mixed_mcmc(threshold_y, [-1.0e40_dp, 0.0_dp, 1.0_dp, 1.0e40_dp], &
            0.25_dp, .true., x_threshold, z_threshold, ainv, beta_mean_scalar, beta_precision_scalar, &
            1.0_dp, 4.0_dp, 1.0_dp, 4.0_dp, .true., 0.65_dp, 160, 60, 5, state, ordinal_native_result, info, &
            observed=observed_threshold)
        call check(info == 0, 'native family-14 ordered-probit sampler status', failures)
        if (info == 0) then
            call check(ordinal_native_result%liability_acceptance_rate > 0.0_dp .and. &
                ordinal_native_result%liability_acceptance_rate < 1.0_dp, &
                'native family-14 liability MH acceptance', failures)
            call check(ordinal_native_result%cutpoint_acceptance_rate > 0.0_dp .and. &
                ordinal_native_result%cutpoint_acceptance_rate < 1.0_dp, &
                'native family-14 cutpoint MH acceptance', failures)
            call check(all(ieee_is_finite(ordinal_native_result%log_likelihood)), &
                'native family-14 finite likelihoods', failures)
            call check(all(ordinal_native_result%r > 0.0_dp), &
                'native family-14 positive latent residual variance', failures)
            call check(ieee_is_finite(ordinal_native_result%last_liability(3)), &
                'native family-14 missing response liability is imputed', failures)
        end if

        ordinal_binary_y = [1, 2, 1, 2, 2, 1]
        call rng_seed(state, 121426_8)
        call ordinal_native_mixed_mcmc(ordinal_binary_y, [-1.0e40_dp, 0.0_dp, 1.0e40_dp], &
            0.25_dp, .false., x, z, ainv, beta_mean_scalar, beta_precision_scalar, &
            1.0_dp, 4.0_dp, 1.0_dp, 4.0_dp, .false., 0.65_dp, 60, 20, 4, state, ordinal_native_result, info, &
            slice_sampling=.true.)
        call check(info == 0, 'native binary family-14 slice sampler status', failures)
        if (info == 0) then
            call check_close(ordinal_native_result%liability_acceptance_rate, 1.0_dp, 1.0e-14_dp, &
                'native binary family-14 slice updates bypass MH rejection', failures)
            call check(all(ieee_is_finite(ordinal_native_result%last_liability)), &
                'native binary family-14 slice liabilities finite', failures)
        end if
    end subroutine test_family_samplers

end program test_mcmcglmm
