! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from MCMCglmm 2.36 buildV/predict numerical covariance assembly; see NOTICE.md and upstream/.
module mcmcglmm_prediction
    use r_kinds, only : dp
    use r_linalg, only : inverse_matrix
    use mcmcglmm_math, only : logistic
    implicit none
    private

    public :: multi_term_build_v
    public :: posterior_linear_predictor
    public :: scalar_response_expectation

contains

    pure integer function prediction_vector_index(row, trait, observations) result(index_value)
        integer, intent(in) :: row !! One-based observation row within a response trait.
        integer, intent(in) :: trait !! One-based response-trait index.
        integer, intent(in) :: observations !! Number of observation rows in each response trait.

        index_value = (trait - 1) * observations + row
    end function prediction_vector_index

    pure subroutine multi_term_build_v(z, random_term, a_inverse, g_matrix, r_matrix, include_g, &
                                       covariance, diagonal, info)
        real(dp), intent(in) :: z(:, :) !! n by q concatenated random-effect design shared across response traits.
        integer, intent(in) :: random_term(:) !! Length-q one-based G-structure label for each random-effect column.
        real(dp), intent(in) :: a_inverse(:, :) !! q by q random-level precision, block-partitioned by random_term.
        real(dp), intent(in) :: g_matrix(:, :, :) !! traits by traits by n_term random-effect covariance draws.
        real(dp), intent(in) :: r_matrix(:, :) !! traits by traits residual covariance for independent observation rows.
        logical, intent(in) :: include_g(:) !! Length-n_term mask; true includes that random structure in marginal V.
        real(dp), allocatable, intent(out) :: covariance(:, :) !! Allocated nt by nt covariance in trait-major row order.
        real(dp), allocatable, intent(out) :: diagonal(:, :) !! Allocated n by traits marginal variance matrix.
        integer, intent(out) :: info !! Zero on success; nonzero for shapes, missing terms, or singular level precision.
        real(dp), allocatable :: a_covariance(:, :)
        real(dp), allocatable :: a_term_inverse(:, :)
        real(dp), allocatable :: k_term(:, :)
        real(dp), allocatable :: z_term(:, :)
        integer, allocatable :: term_indices(:)
        integer :: a
        integer :: b
        integer :: i
        integer :: j
        integer :: n
        integer :: nterm
        integer :: q
        integer :: qterm
        integer :: random_index
        integer :: t
        integer :: term

        info = 0
        n = size(z, 1)
        q = size(z, 2)
        t = size(g_matrix, 1)
        nterm = size(g_matrix, 3)
        if (n < 1 .or. q < 1 .or. t < 1 .or. nterm < 1 .or. size(random_term) /= q .or. &
            size(a_inverse, 1) /= q .or. size(a_inverse, 2) /= q .or. size(g_matrix, 2) /= t .or. &
            size(r_matrix, 1) /= t .or. size(r_matrix, 2) /= t .or. size(include_g) /= nterm .or. &
            any(random_term < 1) .or. any(random_term > nterm)) then
            allocate(covariance(0, 0), diagonal(0, 0))
            info = 1
            return
        end if
        do term = 1, nterm
            if (count(random_term == term) < 1) then
                allocate(covariance(0, 0), diagonal(0, 0))
                info = 1
                return
            end if
        end do

        allocate(covariance(n * t, n * t), diagonal(n, t))
        covariance = 0.0_dp
        diagonal = 0.0_dp
        do a = 1, t
            do b = 1, t
                do i = 1, n
                    covariance(prediction_vector_index(i, a, n), prediction_vector_index(i, b, n)) = r_matrix(a, b)
                end do
            end do
            diagonal(:, a) = r_matrix(a, a)
        end do

        do term = 1, nterm
            if (.not. include_g(term)) cycle
            qterm = count(random_term == term)
            allocate(term_indices(qterm))
            i = 0
            do random_index = 1, q
                if (random_term(random_index) /= term) cycle
                i = i + 1
                term_indices(i) = random_index
            end do
            a_term_inverse = a_inverse(term_indices, term_indices)
            call inverse_matrix(a_term_inverse, a_covariance, info)
            if (info /= 0) then
                deallocate(term_indices)
                deallocate(covariance, diagonal)
                allocate(covariance(0, 0), diagonal(0, 0))
                return
            end if
            z_term = z(:, term_indices)
            k_term = matmul(z_term, matmul(a_covariance, transpose(z_term)))
            do a = 1, t
                do b = 1, t
                    do j = 1, n
                        do i = 1, n
                            covariance(prediction_vector_index(i, a, n), prediction_vector_index(j, b, n)) = &
                                covariance(prediction_vector_index(i, a, n), prediction_vector_index(j, b, n)) + &
                                g_matrix(a, b, term) * k_term(i, j)
                        end do
                    end do
                end do
                diagonal(:, a) = diagonal(:, a) + g_matrix(a, a, term) * [(k_term(i, i), i = 1, n)]
            end do
            deallocate(term_indices)
        end do
    end subroutine multi_term_build_v

    pure real(dp) function standard_normal_density(value) result(density)
        real(dp), intent(in) :: value !! Standard-normal argument used by deterministic response integration.
        real(dp), parameter :: inv_sqrt_two_pi = 0.398942280401432677939946059934381868_dp

        density = inv_sqrt_two_pi * exp(-0.5_dp * value * value)
    end function standard_normal_density

    pure real(dp) function scalar_response_at_link(family, link_value, additional) result(response_mean)
        integer, intent(in) :: family !! Native scalar family code requiring numerical Gaussian-link integration.
        real(dp), intent(in) :: link_value !! Conditional scalar link value at one Gaussian quadrature point.
        real(dp), intent(in) :: additional !! Family auxiliary value, currently binomial trial count when required.
        integer :: trials
        real(dp) :: lambda
        real(dp) :: probability

        select case (family)
        case (3)
            trials = nint(additional)
            response_mean = real(trials, dp) * logistic(link_value)
        case (16)
            if (link_value > log(huge(1.0_dp)) - 2.0_dp) then
                response_mean = huge(1.0_dp)
            else
                lambda = exp(link_value)
                if (lambda < sqrt(epsilon(1.0_dp))) then
                    response_mean = 1.0_dp + 0.5_dp * lambda
                else if (lambda > 40.0_dp) then
                    response_mean = lambda
                else
                    response_mean = lambda / (1.0_dp - exp(-lambda))
                end if
            end if
        case (22)
            trials = nint(additional)
            probability = logistic(link_value)
            response_mean = 1.0_dp - (1.0_dp - probability) ** trials
        case default
            response_mean = 0.0_dp
        end select
    end function scalar_response_at_link

    pure subroutine scalar_response_expectation(family, mean_link, variance_link, additional, additional2, &
                                                expectation, info)
        integer, intent(in) :: family !! Native scalar family code using the maintained simulation/prediction convention.
        real(dp), intent(in) :: mean_link !! Posterior Gaussian mean of the latent/link-scale predictor.
        real(dp), intent(in) :: variance_link !! Nonnegative marginal Gaussian variance of the latent/link predictor.
        real(dp), intent(in) :: additional !! Family auxiliary value such as binomial trials or Student-t scale.
        real(dp), intent(in) :: additional2 !! Second auxiliary value, used as Student-t degrees of freedom.
        real(dp), intent(out) :: expectation !! Posterior response-scale mean after integrating latent Gaussian variance.
        integer, intent(out) :: info !! Zero on success; nonzero for unsupported family or invalid auxiliary parameters.
        integer, parameter :: intervals = 512
        integer :: i
        real(dp) :: h
        real(dp) :: link_value
        real(dp) :: scale
        real(dp) :: total
        real(dp) :: value
        real(dp) :: weight
        real(dp) :: z

        info = 0
        if (variance_link < 0.0_dp) then
            expectation = 0.0_dp
            info = 1
            return
        end if
        select case (family)
        case (1, 6, 24)
            expectation = mean_link
        case (2, 7)
            expectation = exp(mean_link + 0.5_dp * variance_link)
        case (4, 5, 8, 9, 17)
            expectation = exp(-mean_link + 0.5_dp * variance_link)
        case (23)
            if (additional <= 0.0_dp .or. additional2 <= 1.0_dp) then
                expectation = 0.0_dp
                info = 1
                return
            end if
            scale = exp(0.5_dp * log(additional2 / 2.0_dp) + log_gamma(0.5_dp * (additional2 - 1.0_dp)) - &
                log_gamma(0.5_dp * additional2))
            expectation = scale * mean_link
        case (3, 16, 22)
            if ((family == 3 .or. family == 22) .and. nint(additional) < 1) then
                expectation = 0.0_dp
                info = 1
                return
            end if
            if (variance_link <= epsilon(1.0_dp)) then
                expectation = scalar_response_at_link(family, mean_link, additional)
                return
            end if
            h = 16.0_dp / real(intervals, dp)
            total = 0.0_dp
            do i = 0, intervals
                z = -8.0_dp + h * real(i, dp)
                link_value = mean_link + sqrt(variance_link) * z
                value = scalar_response_at_link(family, link_value, additional)
                if (i == 0 .or. i == intervals) then
                    weight = 1.0_dp
                else if (modulo(i, 2) == 0) then
                    weight = 2.0_dp
                else
                    weight = 4.0_dp
                end if
                total = total + weight * value * standard_normal_density(z)
            end do
            expectation = h * total / 3.0_dp
        case default
            expectation = 0.0_dp
            info = 2
        end select
    end subroutine scalar_response_expectation

    pure subroutine posterior_linear_predictor(x, z, random_term, beta_draws, random_effect_draws, &
                                               marginalize_g, predictor, info)
        real(dp), intent(in) :: x(:, :) !! n by p fixed-effect design matrix used for posterior prediction.
        real(dp), intent(in) :: z(:, :) !! n by q concatenated random-effect design matrix used for prediction.
        integer, intent(in) :: random_term(:) !! Length-q one-based covariance-structure label for each random column.
        real(dp), intent(in) :: beta_draws(:, :, :) !! p by traits by samples retained fixed-effect posterior draws.
        real(dp), intent(in) :: random_effect_draws(:, :, :) !! q by traits by samples retained random-effect draws.
        logical, intent(in) :: marginalize_g(:) !! Length-n_term mask; true omits that random term from the predictor.
        real(dp), allocatable, intent(out) :: predictor(:, :, :) !! Allocated n by traits by samples linear predictors.
        integer, intent(out) :: info !! Zero on success; nonzero for incompatible posterior/design dimensions.
        real(dp), allocatable :: retained_random(:, :)
        integer :: draw
        integer :: n
        integer :: nterm
        integer :: p
        integer :: q
        integer :: random_index
        integer :: samples
        integer :: traits

        info = 0
        n = size(x, 1)
        p = size(x, 2)
        q = size(z, 2)
        traits = size(beta_draws, 2)
        samples = size(beta_draws, 3)
        nterm = size(marginalize_g)
        if (n < 1 .or. p < 1 .or. q < 1 .or. traits < 1 .or. samples < 1 .or. size(z, 1) /= n .or. &
            size(beta_draws, 1) /= p .or. size(random_effect_draws, 1) /= q .or. &
            size(random_effect_draws, 2) /= traits .or. size(random_effect_draws, 3) /= samples .or. &
            size(random_term) /= q .or. any(random_term < 1) .or. any(random_term > nterm)) then
            allocate(predictor(0, 0, 0))
            info = 1
            return
        end if
        allocate(predictor(n, traits, samples), retained_random(q, traits))
        do draw = 1, samples
            retained_random = random_effect_draws(:, :, draw)
            do random_index = 1, q
                if (marginalize_g(random_term(random_index))) retained_random(random_index, :) = 0.0_dp
            end do
            predictor(:, :, draw) = matmul(x, beta_draws(:, :, draw)) + matmul(z, retained_random)
        end do
    end subroutine posterior_linear_predictor

end module mcmcglmm_prediction
