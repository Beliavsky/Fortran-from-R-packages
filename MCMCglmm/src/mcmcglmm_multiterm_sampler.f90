! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from MCMCglmm 2.36 multi-G mixed-model conditionals; see NOTICE.md and upstream/.
module mcmcglmm_multiterm_sampler
    use r_kinds, only : dp
    use r_linalg, only : inverse_matrix
    use mcmcglmm_rng, only : rng_state
    use mcmcglmm_matrix, only : mvn_log_density, sample_mvn_precision
    use mcmcglmm_covariance, only : covariance_update_dispatch
    use mcmcglmm_sparse, only : mcmcglmm_sparse_matrix, sparse_from_coo, sparse_stacked_crossproduct, &
        sparse_transpose_matmul_matrix, sparse_validate
    use mcmcglmm_sparse_factorization, only : sample_mvn_sparse_precision, sparse_precision_cache
    implicit none
    private

    type, public :: multi_term_gaussian_mcmc_result
        real(dp), allocatable :: beta(:, :, :)
        real(dp), allocatable :: random_effects(:, :, :)
        real(dp), allocatable :: g(:, :, :, :)
        real(dp), allocatable :: r(:, :, :)
        real(dp), allocatable :: log_likelihood(:)
    end type multi_term_gaussian_mcmc_result

    public :: multi_term_coefficient_conditional
    public :: multi_term_coefficient_conditional_sparse
    public :: multi_term_gaussian_mixed_mcmc
    public :: multi_term_gaussian_loglik

contains

    pure integer function mt_joint_index(trait, effect, effects_per_trait) result(index_value)
        integer, intent(in) :: trait !! One-based trait index.
        integer, intent(in) :: effect !! One-based fixed-then-random coefficient index within a trait block.
        integer, intent(in) :: effects_per_trait !! Number of fixed plus random coefficients per trait.

        index_value = (trait - 1) * effects_per_trait + effect
    end function mt_joint_index

    pure integer function mt_beta_index(trait, effect, fixed_effects) result(index_value)
        integer, intent(in) :: trait !! One-based trait index.
        integer, intent(in) :: effect !! One-based fixed-effect index within a trait block.
        integer, intent(in) :: fixed_effects !! Number of fixed-effect coefficients per trait.

        index_value = (trait - 1) * fixed_effects + effect
    end function mt_beta_index

    pure subroutine mt_unpack_coefficients(packed, fixed_effects, random_levels, traits, beta, random_effects)
        real(dp), intent(in) :: packed(:) !! Trait-major vector containing fixed then random coefficients.
        integer, intent(in) :: fixed_effects !! Number of fixed-effect coefficients per trait.
        integer, intent(in) :: random_levels !! Total random-effect coefficients per trait across all terms.
        integer, intent(in) :: traits !! Number of response traits.
        real(dp), allocatable, intent(out) :: beta(:, :) !! Allocated p by traits fixed-effect matrix.
        real(dp), allocatable, intent(out) :: random_effects(:, :) !! Allocated q by traits random-effect matrix.
        integer :: a
        integer :: k
        integer :: m

        m = fixed_effects + random_levels
        allocate(beta(fixed_effects, traits), random_effects(random_levels, traits))
        do a = 1, traits
            do k = 1, fixed_effects
                beta(k, a) = packed(mt_joint_index(a, k, m))
            end do
            do k = 1, random_levels
                random_effects(k, a) = packed(mt_joint_index(a, fixed_effects + k, m))
            end do
        end do
    end subroutine mt_unpack_coefficients

    pure subroutine multi_term_coefficient_conditional(y, x, z, random_term, a_inverse, g_matrix, r_matrix, &
                                                       beta_prior_mean, beta_prior_precision, state, beta, &
                                                       random_effects, info, mean_value)
        real(dp), intent(in) :: y(:, :) !! n by traits Gaussian response matrix.
        real(dp), intent(in) :: x(:, :) !! n by p fixed-effect design shared across traits.
        real(dp), intent(in) :: z(:, :) !! n by q concatenated random-effect design shared across traits.
        integer, intent(in) :: random_term(:) !! Length-q one-based G-structure label for each random-effect column.
        real(dp), intent(in) :: a_inverse(:, :) !! q by q block row precision for all random effects.
        real(dp), intent(in) :: g_matrix(:, :, :) !! traits by traits by n_term random-effect covariance blocks.
        real(dp), intent(in) :: r_matrix(:, :) !! traits by traits residual covariance matrix.
        real(dp), intent(in) :: beta_prior_mean(:, :) !! p by traits Gaussian fixed-effect prior mean.
        real(dp), intent(in) :: beta_prior_precision(:, :) !! p*traits square Gaussian fixed-effect prior precision.
        type(rng_state), intent(inout) :: state !! Generator state consumed by the joint coefficient Gaussian draw.
        real(dp), allocatable, intent(out) :: beta(:, :) !! Allocated sampled fixed-effect matrix.
        real(dp), allocatable, intent(out) :: random_effects(:, :) !! Allocated sampled concatenated random effects.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid shape or SPD failure.
        real(dp), allocatable, optional, intent(out) :: mean_value(:) !! Optional packed conditional mean.
        real(dp), allocatable :: w(:, :)
        real(dp), allocatable :: wtw(:, :)
        real(dp), allocatable :: wty(:, :)
        integer :: m
        integer :: p
        integer :: q
        integer :: traits

        info = 0
        p = size(x, 2)
        q = size(z, 2)
        traits = size(y, 2)
        m = p + q
        if (size(y, 1) /= size(x, 1) .or. size(y, 1) /= size(z, 1)) then
            allocate(beta(0, 0), random_effects(0, 0))
            info = 1
            return
        end if
        allocate(w(size(y, 1), m))
        w(:, 1:p) = x
        w(:, p + 1:m) = z
        wtw = matmul(transpose(w), w)
        wty = matmul(transpose(w), y)

        call coefficient_conditional_from_products(wtw, wty, p, q, traits, random_term, a_inverse, g_matrix, &
            r_matrix, beta_prior_mean, beta_prior_precision, state, beta, random_effects, info, &
            packed_mean=mean_value)
    end subroutine multi_term_coefficient_conditional

    pure subroutine multi_term_coefficient_conditional_sparse(y, x, z, random_term, a_inverse, g_matrix, r_matrix, &
                                                              beta_prior_mean, beta_prior_precision, state, beta, &
                                                              random_effects, info, mean_value, factor_cache)
        !! Draw coefficients using normal-equation products formed directly from CSR designs.
        real(dp), intent(in) :: y(:, :) !! n by traits Gaussian response matrix.
        type(mcmcglmm_sparse_matrix), intent(in) :: x !! n by p fixed-effect design in canonical CSR storage.
        type(mcmcglmm_sparse_matrix), intent(in) :: z !! n by q random-effect design in canonical CSR storage.
        integer, intent(in) :: random_term(:) !! Length-q one-based G-structure label for each random-effect column.
        real(dp), intent(in) :: a_inverse(:, :) !! q by q block row precision for all random effects.
        real(dp), intent(in) :: g_matrix(:, :, :) !! traits by traits by n_term random-effect covariance blocks.
        real(dp), intent(in) :: r_matrix(:, :) !! traits by traits residual covariance matrix.
        real(dp), intent(in) :: beta_prior_mean(:, :) !! p by traits Gaussian fixed-effect prior mean.
        real(dp), intent(in) :: beta_prior_precision(:, :) !! p*traits square Gaussian fixed-effect prior precision.
        type(rng_state), intent(inout) :: state !! Generator state consumed by the joint coefficient Gaussian draw.
        real(dp), allocatable, intent(out) :: beta(:, :) !! Allocated sampled fixed-effect matrix.
        real(dp), allocatable, intent(out) :: random_effects(:, :) !! Allocated sampled concatenated random effects.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid shape or SPD failure.
        real(dp), allocatable, optional, intent(out) :: mean_value(:) !! Optional packed conditional mean.
        type(sparse_precision_cache), optional, intent(inout) :: factor_cache !! Optional reusable symbolic cache.
        real(dp), allocatable :: xty(:, :)
        real(dp), allocatable :: zty(:, :)
        real(dp), allocatable :: wtw(:, :)
        real(dp), allocatable :: wty(:, :)
        type(mcmcglmm_sparse_matrix) :: sparse_wtw
        integer :: p
        integer :: q
        integer :: traits

        call sparse_validate(x, info)
        if (info /= 0) then
            allocate(beta(0, 0), random_effects(0, 0))
            return
        end if
        call sparse_validate(z, info)
        if (info /= 0 .or. x%nrow /= z%nrow .or. size(y, 1) /= x%nrow) then
            allocate(beta(0, 0), random_effects(0, 0))
            if (info == 0) info = 1
            return
        end if
        p = x%ncol
        q = z%ncol
        traits = size(y, 2)
        call sparse_stacked_crossproduct(x, z, sparse_wtw, info)
        if (info /= 0) return
        call sparse_transpose_matmul_matrix(x, y, xty, info)
        if (info /= 0) return
        call sparse_transpose_matmul_matrix(z, y, zty, info)
        if (info /= 0) return
        allocate(wtw(0, 0), wty(p + q, traits))
        wty(:p, :) = xty
        wty(p + 1:, :) = zty

        call coefficient_conditional_from_products(wtw, wty, p, q, traits, random_term, a_inverse, g_matrix, &
            r_matrix, beta_prior_mean, beta_prior_precision, state, beta, random_effects, info, &
            use_sparse_factorization=.true., packed_mean=mean_value, sparse_wtw=sparse_wtw, &
            factor_cache=factor_cache)
    end subroutine multi_term_coefficient_conditional_sparse

    pure subroutine coefficient_conditional_from_products(wtw, wty, p, q, traits, random_term, a_inverse, &
                                                           g_matrix, r_matrix, beta_prior_mean, &
                                                           beta_prior_precision, state, beta, random_effects, info, &
                                                           use_sparse_factorization, packed_mean, sparse_wtw, &
                                                           factor_cache)
        !! Draw joint coefficients from precomputed design crossproducts and cross-response products.
        real(dp), intent(in) :: wtw(:, :) !! (p+q)-square transpose(W)*W matrix.
        real(dp), intent(in) :: wty(:, :) !! (p+q) by traits transpose(W)*Y matrix.
        integer, intent(in) :: p !! Number of fixed-effect columns.
        integer, intent(in) :: q !! Number of random-effect columns.
        integer, intent(in) :: traits !! Number of response traits.
        integer, intent(in) :: random_term(:) !! Length-q one-based covariance-term labels.
        real(dp), intent(in) :: a_inverse(:, :) !! q by q block row precision for random effects.
        real(dp), intent(in) :: g_matrix(:, :, :) !! traits by traits by n_term random-effect covariance blocks.
        real(dp), intent(in) :: r_matrix(:, :) !! traits by traits residual covariance matrix.
        real(dp), intent(in) :: beta_prior_mean(:, :) !! p by traits fixed-effect prior mean.
        real(dp), intent(in) :: beta_prior_precision(:, :) !! p*traits square fixed-effect prior precision.
        type(rng_state), intent(inout) :: state !! Generator state consumed by the Gaussian draw.
        real(dp), allocatable, intent(out) :: beta(:, :) !! Allocated sampled fixed effects.
        real(dp), allocatable, intent(out) :: random_effects(:, :) !! Allocated sampled random effects.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid shape or SPD failure.
        logical, optional, intent(in) :: use_sparse_factorization !! Use sparse posterior assembly and Cholesky.
        real(dp), allocatable, optional, intent(out) :: packed_mean(:) !! Optional packed conditional mean.
        type(mcmcglmm_sparse_matrix), optional, intent(in) :: sparse_wtw !! Optional CSR design crossproduct.
        type(sparse_precision_cache), optional, intent(inout) :: factor_cache !! Optional reusable symbolic cache.
        real(dp), allocatable :: conditional_mean(:)
        real(dp), allocatable :: g_inverse(:, :, :)
        real(dp), allocatable :: inverse_block(:, :)
        real(dp), allocatable :: packed(:)
        real(dp), allocatable :: precision(:, :)
        type(mcmcglmm_sparse_matrix) :: sparse_precision
        real(dp), allocatable :: r_inverse(:, :)
        real(dp), allocatable :: rhs(:)
        integer :: m
        integer :: nterm
        integer :: term
        logical :: sparse_factorization

        info = 0
        sparse_factorization = .false.
        if (present(use_sparse_factorization)) sparse_factorization = use_sparse_factorization
        m = p + q
        nterm = size(g_matrix, 3)
        if (any(shape(wty) /= [m, traits]) .or. size(random_term) /= q .or. &
            any(shape(a_inverse) /= [q, q]) .or. size(g_matrix, 1) /= traits .or. &
            size(g_matrix, 2) /= traits .or. any(shape(r_matrix) /= [traits, traits]) .or. &
            any(shape(beta_prior_mean) /= [p, traits]) .or. &
            any(shape(beta_prior_precision) /= [p * traits, p * traits]) .or. &
            any(random_term < 1) .or. any(random_term > nterm)) then
            allocate(beta(0, 0), random_effects(0, 0))
            info = 1
            return
        end if
        if (sparse_factorization) then
            if (.not. present(sparse_wtw)) then
                allocate(beta(0, 0), random_effects(0, 0))
                info = 1
                return
            end if
            call sparse_validate(sparse_wtw, info)
            if (info /= 0 .or. sparse_wtw%nrow /= m .or. sparse_wtw%ncol /= m) then
                allocate(beta(0, 0), random_effects(0, 0))
                info = 1
                return
            end if
        else if (any(shape(wtw) /= [m, m])) then
            allocate(beta(0, 0), random_effects(0, 0))
            info = 1
            return
        end if

        allocate(g_inverse(traits, traits, nterm))
        do term = 1, nterm
            call inverse_matrix(g_matrix(:, :, term), inverse_block, info)
            if (info /= 0) then
                allocate(beta(0, 0), random_effects(0, 0))
                return
            end if
            g_inverse(:, :, term) = inverse_block
        end do
        call inverse_matrix(r_matrix, r_inverse, info)
        if (info /= 0) then
            allocate(beta(0, 0), random_effects(0, 0))
            return
        end if
        call assemble_coefficient_rhs(wty, p, q, traits, r_inverse, beta_prior_mean, beta_prior_precision, rhs)
        if (sparse_factorization) then
            call assemble_coefficient_precision_sparse(sparse_wtw, p, q, traits, random_term, a_inverse, g_inverse, &
                r_inverse, beta_prior_precision, sparse_precision, info)
            if (info /= 0) then
                allocate(beta(0, 0), random_effects(0, 0))
                return
            end if
            call sample_mvn_sparse_precision(state, rhs, sparse_precision, packed, conditional_mean, info, &
                use_rcm=.true., cache=factor_cache)
        else
            call assemble_coefficient_precision_dense(wtw, p, q, traits, random_term, a_inverse, g_inverse, &
                r_inverse, beta_prior_precision, precision)
            call sample_mvn_precision(state, rhs, precision, packed, conditional_mean, info)
        end if
        if (info /= 0) then
            allocate(beta(0, 0), random_effects(0, 0))
            return
        end if
        if (present(packed_mean)) packed_mean = conditional_mean
        call mt_unpack_coefficients(packed, p, q, traits, beta, random_effects)
    end subroutine coefficient_conditional_from_products

    pure subroutine assemble_coefficient_rhs(wty, p, q, traits, r_inverse, beta_prior_mean, &
                                               beta_prior_precision, rhs)
        !! Assemble the precision-weighted coefficient mean independently of precision storage.
        real(dp), intent(in) :: wty(:, :) !! (p+q) by traits transpose(W)*Y matrix.
        integer, intent(in) :: p !! Number of fixed-effect columns.
        integer, intent(in) :: q !! Number of random-effect columns.
        integer, intent(in) :: traits !! Number of response traits.
        real(dp), intent(in) :: r_inverse(:, :) !! traits-square residual precision.
        real(dp), intent(in) :: beta_prior_mean(:, :) !! p by traits fixed-effect prior mean.
        real(dp), intent(in) :: beta_prior_precision(:, :) !! p*traits square prior precision.
        real(dp), allocatable, intent(out) :: rhs(:) !! Allocated joint precision-weighted mean.
        real(dp), allocatable :: beta_mean_vector(:)
        real(dp), allocatable :: beta_rhs(:)
        integer :: a
        integer :: b
        integer :: i
        integer :: k
        integer :: m

        m = p + q
        allocate(rhs(m * traits), beta_mean_vector(p * traits))
        rhs = 0.0_dp
        do a = 1, traits
            do k = 1, p
                beta_mean_vector(mt_beta_index(a, k, p)) = beta_prior_mean(k, a)
            end do
        end do
        beta_rhs = matmul(beta_prior_precision, beta_mean_vector)
        do a = 1, traits
            do i = 1, m
                do b = 1, traits
                    rhs(mt_joint_index(a, i, m)) = rhs(mt_joint_index(a, i, m)) + r_inverse(a, b) * wty(i, b)
                end do
            end do
            do i = 1, p
                rhs(mt_joint_index(a, i, m)) = rhs(mt_joint_index(a, i, m)) + beta_rhs(mt_beta_index(a, i, p))
            end do
        end do
    end subroutine assemble_coefficient_rhs

    pure subroutine assemble_coefficient_precision_dense(wtw, p, q, traits, random_term, a_inverse, g_inverse, &
                                                          r_inverse, beta_prior_precision, precision)
        !! Assemble the joint coefficient precision in dense storage for the legacy solver path.
        real(dp), intent(in) :: wtw(:, :) !! (p+q)-square design crossproduct.
        integer, intent(in) :: p !! Number of fixed-effect columns.
        integer, intent(in) :: q !! Number of random-effect columns.
        integer, intent(in) :: traits !! Number of response traits.
        integer, intent(in) :: random_term(:) !! Length-q covariance-term labels.
        real(dp), intent(in) :: a_inverse(:, :) !! q-square random-level precision.
        real(dp), intent(in) :: g_inverse(:, :, :) !! Per-term trait precision matrices.
        real(dp), intent(in) :: r_inverse(:, :) !! Residual trait precision matrix.
        real(dp), intent(in) :: beta_prior_precision(:, :) !! Fixed-effect prior precision.
        real(dp), allocatable, intent(out) :: precision(:, :) !! Allocated dense joint precision.
        integer :: a
        integer :: b
        integer :: i
        integer :: j
        integer :: m
        integer :: term

        m = p + q
        allocate(precision(m * traits, m * traits), source=0.0_dp)
        do a = 1, traits
            do b = 1, traits
                do j = 1, m
                    do i = 1, m
                        precision(mt_joint_index(a, i, m), mt_joint_index(b, j, m)) = &
                            precision(mt_joint_index(a, i, m), mt_joint_index(b, j, m)) + r_inverse(a, b) * wtw(i, j)
                    end do
                end do
                do j = 1, q
                    do i = 1, q
                        if (random_term(i) /= random_term(j)) cycle
                        term = random_term(i)
                        precision(mt_joint_index(a, p + i, m), mt_joint_index(b, p + j, m)) = &
                            precision(mt_joint_index(a, p + i, m), mt_joint_index(b, p + j, m)) + &
                            g_inverse(a, b, term) * a_inverse(i, j)
                    end do
                end do
                do j = 1, p
                    do i = 1, p
                        precision(mt_joint_index(a, i, m), mt_joint_index(b, j, m)) = &
                            precision(mt_joint_index(a, i, m), mt_joint_index(b, j, m)) + &
                            beta_prior_precision(mt_beta_index(a, i, p), mt_beta_index(b, j, p))
                    end do
                end do
            end do
        end do
    end subroutine assemble_coefficient_precision_dense

    pure subroutine assemble_coefficient_precision_sparse(wtw, p, q, traits, random_term, a_inverse, g_inverse, &
                                                           r_inverse, beta_prior_precision, precision, info)
        !! Assemble the joint coefficient precision as canonical CSR coordinate sums.
        type(mcmcglmm_sparse_matrix), intent(in) :: wtw !! (p+q)-square CSR design crossproduct.
        integer, intent(in) :: p !! Number of fixed-effect columns.
        integer, intent(in) :: q !! Number of random-effect columns.
        integer, intent(in) :: traits !! Number of response traits.
        integer, intent(in) :: random_term(:) !! Length-q covariance-term labels.
        real(dp), intent(in) :: a_inverse(:, :) !! q-square random-level precision.
        real(dp), intent(in) :: g_inverse(:, :, :) !! Per-term trait precision matrices.
        real(dp), intent(in) :: r_inverse(:, :) !! Residual trait precision matrix.
        real(dp), intent(in) :: beta_prior_precision(:, :) !! Fixed-effect prior precision.
        type(mcmcglmm_sparse_matrix), intent(out) :: precision !! Canonical CSR joint precision.
        integer, intent(out) :: info !! Zero on success; nonzero if CSR construction fails.
        integer, allocatable :: coordinate_column(:)
        integer, allocatable :: coordinate_row(:)
        real(dp), allocatable :: coordinate_value(:)
        integer :: a
        integer :: b
        integer :: i
        integer :: j
        integer :: k
        integer :: m
        integer :: maximum_entries
        integer :: next_entry
        integer :: term
        real(dp) :: contribution

        m = p + q
        maximum_entries = size(wtw%value) * traits * traits + &
            count(beta_prior_precision /= 0.0_dp)
        do j = 1, q
            do i = 1, q
                if (random_term(i) /= random_term(j) .or. a_inverse(i, j) == 0.0_dp) cycle
                term = random_term(i)
                maximum_entries = maximum_entries + count(g_inverse(:, :, term) /= 0.0_dp)
            end do
        end do
        allocate(coordinate_row(maximum_entries), coordinate_column(maximum_entries), &
            coordinate_value(maximum_entries))
        next_entry = 0
        do a = 1, traits
            do b = 1, traits
                if (r_inverse(a, b) /= 0.0_dp) then
                    do i = 1, m
                        do k = wtw%row_pointer(i), wtw%row_pointer(i + 1) - 1
                            j = wtw%column(k)
                            contribution = r_inverse(a, b) * wtw%value(k)
                            if (contribution == 0.0_dp) cycle
                            next_entry = next_entry + 1
                            coordinate_row(next_entry) = mt_joint_index(a, i, m)
                            coordinate_column(next_entry) = mt_joint_index(b, j, m)
                            coordinate_value(next_entry) = contribution
                        end do
                    end do
                end if
                do j = 1, q
                    do i = 1, q
                        if (random_term(i) /= random_term(j)) cycle
                        term = random_term(i)
                        contribution = g_inverse(a, b, term) * a_inverse(i, j)
                        if (contribution == 0.0_dp) cycle
                        next_entry = next_entry + 1
                        coordinate_row(next_entry) = mt_joint_index(a, p + i, m)
                        coordinate_column(next_entry) = mt_joint_index(b, p + j, m)
                        coordinate_value(next_entry) = contribution
                    end do
                end do
                do j = 1, p
                    do i = 1, p
                        contribution = beta_prior_precision(mt_beta_index(a, i, p), mt_beta_index(b, j, p))
                        if (contribution == 0.0_dp) cycle
                        next_entry = next_entry + 1
                        coordinate_row(next_entry) = mt_joint_index(a, i, m)
                        coordinate_column(next_entry) = mt_joint_index(b, j, m)
                        coordinate_value(next_entry) = contribution
                    end do
                end do
            end do
        end do
        call sparse_from_coo(m * traits, m * traits, coordinate_row(:next_entry), coordinate_column(:next_entry), &
            coordinate_value(:next_entry), precision, info)
    end subroutine assemble_coefficient_precision_sparse

    pure subroutine multi_term_gaussian_loglik(y, x, z, beta, random_effects, r_matrix, log_likelihood, info)
        real(dp), intent(in) :: y(:, :) !! n by traits Gaussian response matrix.
        real(dp), intent(in) :: x(:, :) !! n by p fixed-effect design matrix.
        real(dp), intent(in) :: z(:, :) !! n by q concatenated random-effect design matrix.
        real(dp), intent(in) :: beta(:, :) !! p by traits fixed-effect matrix.
        real(dp), intent(in) :: random_effects(:, :) !! q by traits concatenated random-effect coefficients.
        real(dp), intent(in) :: r_matrix(:, :) !! traits by traits residual covariance matrix.
        real(dp), intent(out) :: log_likelihood !! Sum of row-wise multivariate-normal log densities.
        integer, intent(out) :: info !! Zero on success; nonzero for a covariance or shape failure.
        real(dp), allocatable :: mean_row(:)
        real(dp), allocatable :: zero_mean(:)
        real(dp) :: row_log_density
        integer :: i
        integer :: traits

        traits = size(y, 2)
        if (size(x, 1) /= size(y, 1) .or. size(z, 1) /= size(y, 1) .or. size(beta, 1) /= size(x, 2) .or. &
            size(beta, 2) /= traits .or. size(random_effects, 1) /= size(z, 2) .or. &
            size(random_effects, 2) /= traits .or. size(r_matrix, 1) /= traits .or. size(r_matrix, 2) /= traits) then
            log_likelihood = -huge(1.0_dp)
            info = 1
            return
        end if
        allocate(mean_row(traits), zero_mean(traits))
        zero_mean = 0.0_dp
        log_likelihood = 0.0_dp
        do i = 1, size(y, 1)
            mean_row = matmul(x(i, :), beta) + matmul(z(i, :), random_effects)
            call mvn_log_density(y(i, :) - mean_row, zero_mean, r_matrix, row_log_density, info)
            if (info /= 0) then
                log_likelihood = -huge(1.0_dp)
                return
            end if
            log_likelihood = log_likelihood + row_log_density
        end do
    end subroutine multi_term_gaussian_loglik

    pure subroutine multi_term_gaussian_mixed_mcmc(y, x, z, random_term, a_inverse, beta_prior_mean, &
                                                   beta_prior_precision, g_prior_scale, g_prior_df, &
                                                   r_prior_scale, r_prior_df, iterations, burn, thin, state, result, info, &
                                                   update_g, update_r, initial_g, initial_r, g_update_mode, &
                                                   r_update_mode, g_split, r_split, g_fixed_covariance, &
                                                   r_fixed_covariance, g_ante_order, r_ante_order, &
                                                   g_ante_common_beta, r_ante_common_beta, &
                                                   g_ante_common_variance, r_ante_common_variance)
        real(dp), intent(in) :: y(:, :) !! n by traits Gaussian response matrix.
        real(dp), intent(in) :: x(:, :) !! n by p fixed-effect design shared across traits.
        real(dp), intent(in) :: z(:, :) !! n by q concatenated random-effect design shared across traits.
        integer, intent(in) :: random_term(:) !! Length-q one-based G-structure label for each random-effect column.
        real(dp), intent(in) :: a_inverse(:, :) !! q by q block precision for all random-effect levels.
        real(dp), intent(in) :: beta_prior_mean(:, :) !! p by traits Gaussian fixed-effect prior mean.
        real(dp), intent(in) :: beta_prior_precision(:, :) !! p*traits square Gaussian fixed-effect prior precision.
        real(dp), intent(in) :: g_prior_scale(:, :, :) !! traits by traits by n_term inverse-Wishart prior scales.
        real(dp), intent(in) :: g_prior_df(:) !! Length-n_term inverse-Wishart prior degrees of freedom.
        real(dp), intent(in) :: r_prior_scale(:, :) !! traits by traits inverse-Wishart residual prior scale.
        real(dp), intent(in) :: r_prior_df !! Residual inverse-Wishart prior degrees of freedom.
        integer, intent(in) :: iterations !! Total Gibbs iterations including burn-in.
        integer, intent(in) :: burn !! Number of initial iterations discarded.
        integer, intent(in) :: thin !! Positive thinning interval after burn-in.
        type(rng_state), intent(inout) :: state !! Generator state consumed by coefficient and covariance draws.
        type(multi_term_gaussian_mcmc_result), intent(out) :: result !! Retained multi-G posterior samples.
        integer, intent(out) :: info !! Zero on success; nonzero for invalid settings or a numerical failure.
        logical, optional, intent(in) :: update_g(:) !! Optional length-n_term mask; false keeps that G covariance fixed.
        logical, optional, intent(in) :: update_r !! Optional residual-covariance update flag; defaults to true.
        real(dp), optional, intent(in) :: initial_g(:, :, :) !! Optional traits by traits by n_term starting/fixed G matrices.
        real(dp), optional, intent(in) :: initial_r(:, :) !! Optional traits by traits starting/fixed residual covariance.
        integer, optional, intent(in) :: g_update_mode(:) !! Optional per-G covariance mode: 0, 1, 2, 3, 4, or 6.
        integer, optional, intent(in) :: r_update_mode !! Optional residual covariance mode: 0, 1, 2, 3, 4, or 6.
        integer, optional, intent(in) :: g_split(:) !! Optional per-G leading dimensions used by covariance modes 2, 4, and 6.
        integer, optional, intent(in) :: r_split !! Optional residual leading dimension used by covariance modes 2, 4, and 6.
        real(dp), optional, intent(in) :: g_fixed_covariance(:, :, :) !! Full matrices supplying fixed lower-right G blocks.
        real(dp), optional, intent(in) :: r_fixed_covariance(:, :) !! Full matrix supplying a fixed lower-right residual block.
        integer, optional, intent(in) :: g_ante_order(:) !! Optional per-G antedependence lag orders; defaults to one.
        integer, optional, intent(in) :: r_ante_order !! Optional residual antedependence lag order; defaults to one.
        logical, optional, intent(in) :: g_ante_common_beta(:) !! Optional per-G flags sharing one coefficient per ante lag.
        logical, optional, intent(in) :: r_ante_common_beta !! Optional residual flag sharing one coefficient per ante lag.
        logical, optional, intent(in) :: g_ante_common_variance(:) !! Optional per-G flags sharing one ante innovation variance.
        logical, optional, intent(in) :: r_ante_common_variance !! Optional residual flag sharing one ante innovation variance.
        real(dp), allocatable :: beta(:, :)
        real(dp), allocatable :: g_matrix(:, :, :)
        real(dp), allocatable :: fixed_block_in(:, :)
        real(dp), allocatable :: fixed_block_out(:, :)
        integer, allocatable :: g_ante_order_value(:)
        logical, allocatable :: g_ante_common_beta_value(:)
        logical, allocatable :: g_ante_common_variance_value(:)
        integer, allocatable :: g_mode(:)
        integer, allocatable :: g_split_value(:)
        real(dp), allocatable :: g_scale_post(:, :)
        real(dp), allocatable :: random_effects(:, :)
        real(dp), allocatable :: r_matrix(:, :)
        real(dp), allocatable :: sampled_g(:, :)
        real(dp), allocatable :: r_scale_post(:, :)
        real(dp), allocatable :: residual(:, :)
        real(dp), allocatable :: term_a_inverse(:, :)
        real(dp), allocatable :: term_effects(:, :)
        integer, allocatable :: term_indices(:)
        real(dp) :: log_likelihood
        integer :: i
        integer :: iteration
        integer :: n
        integer :: nsave
        integer :: nterm
        integer :: p
        integer :: q
        integer :: qterm
        integer :: random_index
        integer :: save_index
        integer :: term
        integer :: traits
        integer :: r_ante_order_value
        logical :: r_ante_common_beta_value
        logical :: r_ante_common_variance_value
        integer :: r_mode
        integer :: r_split_value
        logical :: accepted
        logical, allocatable :: update_g_mask(:)
        logical :: update_r_value

        info = 0
        n = size(y, 1)
        traits = size(y, 2)
        p = size(x, 2)
        q = size(z, 2)
        nterm = size(g_prior_scale, 3)
        if (n < 1 .or. traits < 1 .or. p < 1 .or. q < 1 .or. nterm < 1 .or. size(x, 1) /= n .or. &
            size(z, 1) /= n .or. size(random_term) /= q .or. size(a_inverse, 1) /= q .or. &
            size(a_inverse, 2) /= q .or. size(beta_prior_mean, 1) /= p .or. &
            size(beta_prior_mean, 2) /= traits .or. size(beta_prior_precision, 1) /= p * traits .or. &
            size(beta_prior_precision, 2) /= p * traits .or. size(g_prior_scale, 1) /= traits .or. &
            size(g_prior_scale, 2) /= traits .or. size(g_prior_df) /= nterm .or. &
            size(r_prior_scale, 1) /= traits .or. size(r_prior_scale, 2) /= traits .or. &
            any(random_term < 1) .or. any(random_term > nterm)) then
            info = 1
            return
        end if
        if (present(update_g)) then
            if (size(update_g) /= nterm) then
                info = 1
                return
            end if
        end if
        if (present(initial_g)) then
            if (size(initial_g, 1) /= traits .or. size(initial_g, 2) /= traits .or. size(initial_g, 3) /= nterm) then
                info = 1
                return
            end if
        end if
        if (present(initial_r)) then
            if (size(initial_r, 1) /= traits .or. size(initial_r, 2) /= traits) then
                info = 1
                return
            end if
        end if
        if (present(g_update_mode)) then
            if (size(g_update_mode) /= nterm .or. any(.not. valid_covariance_mode(g_update_mode))) then
                info = 1
                return
            end if
        end if
        if (present(r_update_mode)) then
            if (.not. valid_covariance_mode(r_update_mode)) then
                info = 1
                return
            end if
        end if
        if (present(g_split)) then
            if (size(g_split) /= nterm .or. any(g_split < 0) .or. any(g_split > traits)) then
                info = 1
                return
            end if
        end if
        if (present(r_split)) then
            if (r_split < 0 .or. r_split > traits) then
                info = 1
                return
            end if
        end if
        if (present(g_fixed_covariance)) then
            if (size(g_fixed_covariance, 1) /= traits .or. size(g_fixed_covariance, 2) /= traits .or. &
                size(g_fixed_covariance, 3) /= nterm) then
                info = 1
                return
            end if
        end if
        if (present(r_fixed_covariance)) then
            if (size(r_fixed_covariance, 1) /= traits .or. size(r_fixed_covariance, 2) /= traits) then
                info = 1
                return
            end if
        end if
        if (present(g_ante_order)) then
            if (size(g_ante_order) /= nterm) then
                info = 1
                return
            end if
        end if
        if (present(g_ante_common_beta)) then
            if (size(g_ante_common_beta) /= nterm) then
                info = 1
                return
            end if
        end if
        if (present(g_ante_common_variance)) then
            if (size(g_ante_common_variance) /= nterm) then
                info = 1
                return
            end if
        end if
        do term = 1, nterm
            if (count(random_term == term) < 1 .or. g_prior_df(term) <= real(traits - 1, dp)) then
                info = 2
                return
            end if
        end do
        if (r_prior_df <= real(traits - 1, dp) .or. iterations <= burn .or. burn < 0 .or. thin < 1) then
            info = 2
            return
        end if
        nsave = (iterations - burn) / thin
        if (nsave < 1) then
            info = 3
            return
        end if

        allocate(update_g_mask(nterm), g_mode(nterm), g_split_value(nterm), g_ante_order_value(nterm))
        allocate(g_ante_common_beta_value(nterm), g_ante_common_variance_value(nterm))
        update_g_mask = .true.
        if (present(update_g)) update_g_mask = update_g
        g_mode = 1
        if (present(g_update_mode)) g_mode = g_update_mode
        where (.not. update_g_mask) g_mode = 0
        g_split_value = traits
        if (present(g_split)) g_split_value = g_split
        update_r_value = .true.
        if (present(update_r)) update_r_value = update_r
        r_mode = 1
        if (present(r_update_mode)) r_mode = r_update_mode
        if (.not. update_r_value) r_mode = 0
        r_split_value = traits
        if (present(r_split)) r_split_value = r_split
        g_ante_order_value = 1
        if (present(g_ante_order)) g_ante_order_value = g_ante_order
        g_ante_common_beta_value = .false.
        if (present(g_ante_common_beta)) g_ante_common_beta_value = g_ante_common_beta
        g_ante_common_variance_value = .false.
        if (present(g_ante_common_variance)) g_ante_common_variance_value = g_ante_common_variance
        r_ante_order_value = 1
        if (present(r_ante_order)) r_ante_order_value = r_ante_order
        r_ante_common_beta_value = .false.
        if (present(r_ante_common_beta)) r_ante_common_beta_value = r_ante_common_beta
        r_ante_common_variance_value = .false.
        if (present(r_ante_common_variance)) r_ante_common_variance_value = r_ante_common_variance
        do term = 1, nterm
            if (g_mode(term) == 5) then
                if (traits < 2 .or. g_ante_order_value(term) < 1 .or. g_ante_order_value(term) >= traits) then
                    info = 5
                    return
                end if
            end if
        end do
        if (r_mode == 5) then
            if (traits < 2 .or. r_ante_order_value < 1 .or. r_ante_order_value >= traits) then
                info = 5
                return
            end if
        end if
        do term = 1, nterm
            if ((g_mode(term) == 2 .or. g_mode(term) == 4) .and. .not. present(g_fixed_covariance)) then
                info = 1
                return
            end if
        end do
        if ((r_mode == 2 .or. r_mode == 4) .and. .not. present(r_fixed_covariance)) then
            info = 1
            return
        end if

        allocate(result%beta(p, traits, nsave), result%random_effects(q, traits, nsave))
        allocate(result%g(traits, traits, nterm, nsave), result%r(traits, traits, nsave))
        allocate(result%log_likelihood(nsave), g_matrix(traits, traits, nterm))
        beta = beta_prior_mean
        allocate(random_effects(q, traits))
        random_effects = 0.0_dp
        do term = 1, nterm
            if (present(initial_g)) then
                g_matrix(:, :, term) = initial_g(:, :, term)
            else
                g_matrix(:, :, term) = g_prior_scale(:, :, term) / &
                    max(g_prior_df(term) - real(traits + 1, dp), 1.0_dp)
            end if
        end do
        if (present(initial_r)) then
            r_matrix = initial_r
        else
            r_matrix = r_prior_scale / max(r_prior_df - real(traits + 1, dp), 1.0_dp)
        end if
        save_index = 0

        do iteration = 1, iterations
            call multi_term_coefficient_conditional(y, x, z, random_term, a_inverse, g_matrix, r_matrix, &
                beta_prior_mean, beta_prior_precision, state, beta, random_effects, info)
            if (info /= 0) return

            do term = 1, nterm
                if (.not. update_g_mask(term)) cycle
                qterm = count(random_term == term)
                allocate(term_indices(qterm))
                i = 0
                do random_index = 1, q
                    if (random_term(random_index) /= term) cycle
                    i = i + 1
                    term_indices(i) = random_index
                end do
                qterm = size(term_indices)
                term_effects = random_effects(term_indices, :)
                term_a_inverse = a_inverse(term_indices, term_indices)
                g_scale_post = g_prior_scale(:, :, term) + &
                    matmul(transpose(term_effects), matmul(term_a_inverse, term_effects))
                if (allocated(fixed_block_in)) deallocate(fixed_block_in)
                if ((g_mode(term) == 2 .or. g_mode(term) == 4) .and. present(g_fixed_covariance)) then
                    fixed_block_in = g_fixed_covariance(g_split_value(term) + 1:traits, &
                        g_split_value(term) + 1:traits, term)
                else
                    allocate(fixed_block_in(0, 0))
                end if
                if (g_mode(term) == 5) then
                    call covariance_update_dispatch(state, g_mode(term), g_scale_post, real(qterm, dp), &
                        g_prior_df(term), g_prior_scale(:, :, term), g_matrix(:, :, term), g_split_value(term), &
                        fixed_block_in, sampled_g, fixed_block_out, accepted, info, ante_location=term_effects, &
                        ante_lag_order=g_ante_order_value(term), &
                        ante_common_beta=g_ante_common_beta_value(term), &
                        ante_common_variance=g_ante_common_variance_value(term), ante_a_inverse=term_a_inverse)
                else
                    call covariance_update_dispatch(state, g_mode(term), g_scale_post, real(qterm, dp), &
                        g_prior_df(term), g_prior_scale(:, :, term), g_matrix(:, :, term), g_split_value(term), &
                        fixed_block_in, sampled_g, fixed_block_out, accepted, info)
                end if
                deallocate(term_indices)
                if (info /= 0) return
                g_matrix(:, :, term) = sampled_g
            end do

            if (r_mode /= 0) then
                residual = y - matmul(x, beta) - matmul(z, random_effects)
                r_scale_post = r_prior_scale + matmul(transpose(residual), residual)
                if (allocated(fixed_block_in)) deallocate(fixed_block_in)
                if ((r_mode == 2 .or. r_mode == 4) .and. present(r_fixed_covariance)) then
                    fixed_block_in = r_fixed_covariance(r_split_value + 1:traits, r_split_value + 1:traits)
                else
                    allocate(fixed_block_in(0, 0))
                end if
                if (r_mode == 5) then
                    call covariance_update_dispatch(state, r_mode, r_scale_post, real(n, dp), r_prior_df, &
                        r_prior_scale, r_matrix, r_split_value, fixed_block_in, sampled_g, fixed_block_out, accepted, &
                        info, ante_location=residual, ante_lag_order=r_ante_order_value, &
                        ante_common_beta=r_ante_common_beta_value, &
                        ante_common_variance=r_ante_common_variance_value)
                else
                    call covariance_update_dispatch(state, r_mode, r_scale_post, real(n, dp), r_prior_df, &
                        r_prior_scale, r_matrix, r_split_value, fixed_block_in, sampled_g, fixed_block_out, accepted, info)
                end if
                if (info /= 0) return
                r_matrix = sampled_g
            end if

            if (iteration > burn .and. modulo(iteration - burn, thin) == 0) then
                save_index = save_index + 1
                result%beta(:, :, save_index) = beta
                result%random_effects(:, :, save_index) = random_effects
                result%g(:, :, :, save_index) = g_matrix
                result%r(:, :, save_index) = r_matrix
                call multi_term_gaussian_loglik(y, x, z, beta, random_effects, r_matrix, log_likelihood, info)
                if (info /= 0) return
                result%log_likelihood(save_index) = log_likelihood
            end if
        end do
    end subroutine multi_term_gaussian_mixed_mcmc


    pure elemental logical function valid_covariance_mode(mode) result(valid)
        integer, intent(in) :: mode !! Candidate MCMCglmm covariance update code.

        valid = any(mode == [0, 1, 2, 3, 4, 5, 6])
    end function valid_covariance_mode

end module mcmcglmm_multiterm_sampler
