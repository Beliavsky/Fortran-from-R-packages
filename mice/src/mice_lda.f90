! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from R package mice 3.19.0 by Stef van Buuren, Karin Groothuis-Oudshoorn,
! and mice contributors; see NOTICE.md and PROVENANCE.md for attribution.
! Linear-discriminant imputation derived from mice.impute.lda.R and standard LDA likelihoods.
module mice_lda
    use r_kinds, only : dp
    use r_linalg, only : inverse_matrix
    use mice_categorical, only : draw_category
    use mice_rng, only : mice_rng_state
    use mice_status, only : mice_ok, mice_invalid_argument, mice_invalid_shape, mice_no_observed, mice_singular
    implicit none
    private

    public :: impute_lda

contains

    subroutine impute_lda(category, observed, x, where, ncat, rng, imputed, info, ridge)
        integer, intent(in) :: category(:) !! Incomplete one-based categorical response codes.
        logical, intent(in) :: observed(:) !! True for rows used to estimate category means and pooled covariance.
        real(dp), intent(in) :: x(:, :) !! Complete numeric predictor matrix.
        logical, intent(in) :: where(:) !! True for rows where LDA posterior draws are requested.
        integer, intent(in), value :: ncat !! Number of response categories; every category must be observed at least once.
        type(mice_rng_state), intent(inout) :: rng !! RNG state used for posterior category draws.
        integer, allocatable, intent(out) :: imputed(:) !! One-based category draws ordered by true entries of `where`.
        integer, intent(out) :: info !! `mice_ok` on success or a package status code on failure.
        real(dp), intent(in), optional :: ridge !! Relative ridge added to pooled covariance; default `1e-8`.

        real(dp), allocatable :: means(:, :), covariance(:, :), inv_cov(:, :), linear(:, :), probs(:), center(:)
        real(dp) :: max_score, prior, ridge_value, score, scale
        integer, allocatable :: counts(:)
        integer :: c, i, j, nobs, p, status

        if (size(observed) /= size(category) .or. size(where) /= size(category) .or. size(x, 1) /= size(category)) then
            info = mice_invalid_shape
            return
        end if
        nobs = count(observed)
        p = size(x, 2)
        if (nobs < 2) then
            info = mice_no_observed
            return
        end if
        if (ncat < 2 .or. p < 1) then
            info = mice_invalid_argument
            return
        end if
        ridge_value = 1.0e-8_dp
        if (present(ridge)) ridge_value = ridge
        if (ridge_value < 0.0_dp) then
            info = mice_invalid_argument
            return
        end if
        allocate(means(ncat, p), covariance(p, p), counts(ncat), center(p))
        means = 0.0_dp
        counts = 0
        do i = 1, size(category)
            if (.not. observed(i)) cycle
            c = category(i)
            if (c < 1 .or. c > ncat) then
                info = mice_invalid_argument
                return
            end if
            counts(c) = counts(c) + 1
            means(c, :) = means(c, :) + x(i, :)
        end do
        if (any(counts < 1)) then
            info = mice_invalid_argument
            return
        end if
        do c = 1, ncat
            means(c, :) = means(c, :) / real(counts(c), dp)
        end do
        covariance = 0.0_dp
        do i = 1, size(category)
            if (.not. observed(i)) cycle
            c = category(i)
            center = x(i, :) - means(c, :)
            covariance = covariance + outer_product(center, center)
        end do
        covariance = covariance / real(max(nobs - ncat, 1), dp)
        scale = max(maxval(abs(covariance)), 1.0_dp)
        do j = 1, p
            covariance(j, j) = covariance(j, j) + max(ridge_value * max(abs(covariance(j, j)), scale), &
                                                        sqrt(epsilon(1.0_dp)) * scale)
        end do
        call inverse_matrix(covariance, inv_cov, status)
        if (status /= 0) then
            info = mice_singular
            return
        end if
        allocate(linear(ncat, p))
        do c = 1, ncat
            linear(c, :) = matmul(inv_cov, means(c, :))
        end do
        allocate(imputed(count(where)), probs(ncat))
        j = 0
        do i = 1, size(category)
            if (.not. where(i)) cycle
            j = j + 1
            do c = 1, ncat
                prior = real(counts(c), dp) / real(nobs, dp)
                score = dot_product(x(i, :), linear(c, :)) - &
                        0.5_dp * dot_product(means(c, :), linear(c, :)) + log(prior)
                probs(c) = score
            end do
            max_score = maxval(probs)
            probs = exp(probs - max_score)
            probs = probs / sum(probs)
            imputed(j) = draw_category(probs, rng)
        end do
        info = mice_ok
    end subroutine impute_lda

    pure function outer_product(a, b) result(matrix)
        real(dp), intent(in) :: a(:) !! Left vector of the outer product.
        real(dp), intent(in) :: b(:) !! Right vector of the outer product.
        real(dp) :: matrix(size(a), size(b))
        integer :: i, j

        do j = 1, size(b)
            do i = 1, size(a)
                matrix(i, j) = a(i) * b(j)
            end do
        end do
    end function outer_product

end module mice_lda
