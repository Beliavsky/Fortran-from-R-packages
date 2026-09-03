! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from R package mice 3.19.0 by Stef van Buuren, Karin Groothuis-Oudshoorn,
! and mice contributors; see NOTICE.md and PROVENANCE.md for attribution.
! Computational translation derived from mice 3.19.0 categorical imputers.
module mice_categorical
    use r_kinds, only : dp
    use r_linalg, only : inverse_matrix, solve_spd
    use mice_regression, only : chol_lower
    use mice_rng, only : mice_rng_state, rng_normal, rng_uniform, rng_integer
    use mice_status, only : mice_ok, mice_invalid_argument, mice_invalid_shape, mice_no_observed, mice_not_converged
    implicit none
    private

    public :: logistic_fit
    public :: impute_logreg
    public :: impute_logreg_boot
    public :: multinomial_fit
    public :: impute_polyreg
    public :: augment_multinomial
    public :: draw_category

contains

    subroutine logistic_fit(x, y, weights, coef, vcov, info, maxit, tol)
        real(dp), intent(in) :: x(:, :) !! Logistic design matrix, normally including an intercept column.
        real(dp), intent(in) :: y(:) !! Binary response coded zero or one.
        real(dp), intent(in) :: weights(:) !! Nonnegative case weights corresponding to rows of `x`.
        real(dp), allocatable, intent(out) :: coef(:) !! Fitted logistic regression coefficients.
        real(dp), allocatable, intent(out) :: vcov(:, :) !! Inverse observed Fisher information at the fitted coefficients.
        integer, intent(out) :: info !! `mice_ok`, `mice_not_converged`, or an argument/linear-algebra status code.
        integer, intent(in), optional :: maxit !! Maximum IRLS iterations; default 50.
        real(dp), intent(in), optional :: tol !! Relative coefficient convergence tolerance; default `1e-8`.

        real(dp), allocatable :: eta(:), mu(:), w(:), z(:), xtwx(:, :), rhs(:), next(:)
        real(dp) :: eps, tolerance
        integer :: i, iter, iterations, la_info, n, p

        n = size(x, 1)
        p = size(x, 2)
        if (size(y) /= n .or. size(weights) /= n .or. n < 1 .or. p < 1) then
            info = mice_invalid_shape
            return
        end if
        if (any(weights < 0.0_dp) .or. any(y < 0.0_dp) .or. any(y > 1.0_dp)) then
            info = mice_invalid_argument
            return
        end if
        iterations = 50
        if (present(maxit)) iterations = maxit
        tolerance = 1.0e-8_dp
        if (present(tol)) tolerance = tol
        if (iterations < 1 .or. tolerance <= 0.0_dp) then
            info = mice_invalid_argument
            return
        end if
        allocate(coef(p), eta(n), mu(n), w(n), z(n), xtwx(p, p), rhs(p), next(p))
        coef = 0.0_dp
        eps = sqrt(epsilon(1.0_dp))
        do iter = 1, iterations
            eta = matmul(x, coef)
            do i = 1, n
                mu(i) = logistic(eta(i))
                mu(i) = min(1.0_dp - eps, max(eps, mu(i)))
                w(i) = weights(i) * mu(i) * (1.0_dp - mu(i))
                z(i) = eta(i) + (y(i) - mu(i)) / (mu(i) * (1.0_dp - mu(i)))
            end do
            call weighted_crossprod(x, w, xtwx)
            call weighted_rhs(x, w, z, rhs)
            call solve_spd(xtwx, rhs, next, la_info)
            if (la_info /= 0) then
                call add_diagonal_ridge(xtwx, 1.0e-8_dp)
                call solve_spd(xtwx, rhs, next, la_info)
            end if
            if (la_info /= 0) then
                info = mice_not_converged
                return
            end if
            if (maxval(abs(next - coef)) <= tolerance * (1.0_dp + maxval(abs(coef)))) then
                coef = next
                call inverse_matrix(xtwx, vcov, la_info)
                if (la_info /= 0) then
                    info = mice_not_converged
                else
                    info = mice_ok
                end if
                return
            end if
            coef = next
        end do
        eta = matmul(x, coef)
        do i = 1, n
            mu(i) = logistic(eta(i))
            mu(i) = min(1.0_dp - eps, max(eps, mu(i)))
            w(i) = weights(i) * mu(i) * (1.0_dp - mu(i))
        end do
        call weighted_crossprod(x, w, xtwx)
        call inverse_matrix(xtwx, vcov, la_info)
        if (la_info /= 0) then
            info = mice_not_converged
        else
            info = mice_not_converged
        end if
    end subroutine logistic_fit

    subroutine impute_logreg(y, observed, x, where, rng, imputed, info, prediction_offset)
        real(dp), intent(in) :: y(:) !! Binary incomplete response coded zero or one.
        logical, intent(in) :: observed(:) !! True for observed binary responses used to fit the model.
        real(dp), intent(in) :: x(:, :) !! Complete numeric predictors without an intercept column.
        logical, intent(in) :: where(:) !! True for rows where binary imputations are requested.
        type(mice_rng_state), intent(inout) :: rng !! RNG state used for posterior coefficient and Bernoulli draws.
        real(dp), allocatable, intent(out) :: imputed(:) !! Binary imputations ordered by true entries of `where`.
        integer, intent(out) :: info !! `mice_ok` on success or a package status code on failure.
        real(dp), intent(in), optional :: prediction_offset(:) !! Additive row offsets applied only to imputation linear predictors.

        real(dp), allocatable :: xa(:, :), ya(:), wa(:), design(:, :), coef(:), vcov(:, :), factor(:, :), beta(:), z(:)
        real(dp) :: eta
        integer :: i, j, p

        if (size(observed) /= size(y) .or. size(where) /= size(y) .or. size(x, 1) /= size(y)) then
            info = mice_invalid_shape
            return
        end if
        if (present(prediction_offset)) then
            if (size(prediction_offset) /= size(y)) then
                info = mice_invalid_shape
                return
            end if
        end if
        if (count(observed) < 1) then
            info = mice_no_observed
            return
        end if
        call augment_binary(y, observed, x, xa, ya, wa)
        allocate(design(size(ya), size(xa, 2) + 1))
        design(:, 1) = 1.0_dp
        if (size(xa, 2) > 0) design(:, 2:) = xa
        call logistic_fit(design, ya, wa, coef, vcov, info)
        if (info /= mice_ok .and. info /= mice_not_converged) return
        p = size(coef)
        allocate(factor(p, p), beta(p), z(p))
        call chol_lower(vcov, factor, info)
        if (info /= mice_ok) return
        do j = 1, p
            z(j) = rng_normal(rng)
        end do
        beta = coef + matmul(factor, z)
        allocate(imputed(count(where)))
        j = 0
        do i = 1, size(y)
            if (.not. where(i)) cycle
            j = j + 1
            eta = beta(1) + dot_product(x(i, :), beta(2:))
            if (present(prediction_offset)) eta = eta + prediction_offset(i)
            if (rng_uniform(rng) <= logistic(eta)) then
                imputed(j) = 1.0_dp
            else
                imputed(j) = 0.0_dp
            end if
        end do
        info = mice_ok
    end subroutine impute_logreg

    subroutine impute_logreg_boot(y, observed, x, where, rng, imputed, info)
        real(dp), intent(in) :: y(:) !! Binary incomplete response coded zero or one.
        logical, intent(in) :: observed(:) !! True for rows eligible for bootstrap resampling.
        real(dp), intent(in) :: x(:, :) !! Complete numeric predictors without an intercept column.
        logical, intent(in) :: where(:) !! True for rows where imputations are requested.
        type(mice_rng_state), intent(inout) :: rng !! RNG state used for bootstrap and Bernoulli draws.
        real(dp), allocatable, intent(out) :: imputed(:) !! Bootstrap-logistic imputations ordered by true entries of `where`.
        integer, intent(out) :: info !! `mice_ok` on success or a package status code on failure.

        real(dp), allocatable :: xo(:, :), yo(:), weights(:), design(:, :), coef(:), vcov(:, :)
        integer, allocatable :: obs_index(:)
        integer :: i, j, nobs, pick

        if (size(observed) /= size(y) .or. size(where) /= size(y) .or. size(x, 1) /= size(y)) then
            info = mice_invalid_shape
            return
        end if
        nobs = count(observed)
        if (nobs < 1) then
            info = mice_no_observed
            return
        end if
        allocate(obs_index(nobs), xo(nobs, size(x, 2)), yo(nobs), weights(nobs), design(nobs, size(x, 2) + 1))
        j = 0
        do i = 1, size(y)
            if (.not. observed(i)) cycle
            j = j + 1
            obs_index(j) = i
        end do
        do i = 1, nobs
            pick = obs_index(rng_integer(rng, 1, nobs))
            xo(i, :) = x(pick, :)
            yo(i) = y(pick)
        end do
        design(:, 1) = 1.0_dp
        if (size(x, 2) > 0) design(:, 2:) = xo
        weights = 1.0_dp
        call logistic_fit(design, yo, weights, coef, vcov, info)
        if (info /= mice_ok .and. info /= mice_not_converged) return
        allocate(imputed(count(where)))
        j = 0
        do i = 1, size(y)
            if (.not. where(i)) cycle
            j = j + 1
            if (rng_uniform(rng) <= logistic(coef(1) + dot_product(x(i, :), coef(2:)))) then
                imputed(j) = 1.0_dp
            else
                imputed(j) = 0.0_dp
            end if
        end do
        info = mice_ok
    end subroutine impute_logreg_boot

    subroutine multinomial_fit(x, category, weights, ncat, coef, info, maxit, tol)
        real(dp), intent(in) :: x(:, :) !! Multinomial design matrix, normally including an intercept.
        integer, intent(in) :: category(:) !! One-based category codes in the range `1:ncat`.
        real(dp), intent(in) :: weights(:) !! Nonnegative case weights.
        integer, intent(in), value :: ncat !! Number of response categories; must be at least two.
        real(dp), allocatable, intent(out) :: coef(:, :) !! Predictor-by-nonbaseline-category coefficient matrix.
        integer, intent(out) :: info !! `mice_ok`, `mice_not_converged`, or an argument status code.
        integer, intent(in), optional :: maxit !! Maximum Newton iterations; default 100.
        real(dp), intent(in), optional :: tol !! Relative Newton convergence tolerance; default `1e-8`.

        real(dp), allocatable :: prob(:, :), grad(:), hess(:, :), step(:), flat(:)
        real(dp) :: tolerance, wi
        integer :: a, b, c, d, i, ia, ib, iter, iterations, j, k, n, p, q, la_info

        n = size(x, 1)
        p = size(x, 2)
        q = p * (ncat - 1)
        if (size(category) /= n .or. size(weights) /= n .or. ncat < 2 .or. p < 1) then
            info = mice_invalid_shape
            return
        end if
        if (any(category < 1) .or. any(category > ncat) .or. any(weights < 0.0_dp)) then
            info = mice_invalid_argument
            return
        end if
        iterations = 100
        if (present(maxit)) iterations = maxit
        tolerance = 1.0e-8_dp
        if (present(tol)) tolerance = tol
        allocate(coef(p, ncat - 1), prob(n, ncat), grad(q), hess(q, q), step(q), flat(q))
        coef = 0.0_dp
        do iter = 1, iterations
            call multinomial_probabilities(x, coef, prob)
            grad = 0.0_dp
            hess = 0.0_dp
            do i = 1, n
                wi = weights(i)
                do a = 1, ncat - 1
                    do j = 1, p
                        ia = (a - 1) * p + j
                        if (category(i) == a + 1) grad(ia) = grad(ia) + wi * x(i, j)
                        grad(ia) = grad(ia) - wi * prob(i, a + 1) * x(i, j)
                        do b = 1, ncat - 1
                            do k = 1, p
                                ib = (b - 1) * p + k
                                c = 0
                                if (a == b) c = 1
                                hess(ia, ib) = hess(ia, ib) + wi * prob(i, a + 1) * &
                                    (real(c, dp) - prob(i, b + 1)) * x(i, j) * x(i, k)
                            end do
                        end do
                    end do
                end do
            end do
            do d = 1, q
                hess(d, d) = hess(d, d) + 1.0e-8_dp
            end do
            call solve_spd(hess, grad, step, la_info)
            if (la_info /= 0) then
                info = mice_not_converged
                return
            end if
            do a = 1, ncat - 1
                flat((a - 1) * p + 1:a * p) = coef(:, a)
            end do
            flat = flat + step
            do a = 1, ncat - 1
                coef(:, a) = flat((a - 1) * p + 1:a * p)
            end do
            if (maxval(abs(step)) <= tolerance * (1.0_dp + maxval(abs(flat)))) then
                info = mice_ok
                return
            end if
        end do
        info = mice_not_converged
    end subroutine multinomial_fit

    subroutine impute_polyreg(category, observed, x, where, ncat, rng, imputed, info)
        integer, intent(in) :: category(:) !! Incomplete one-based nominal category codes; observed entries must lie in `1:ncat`.
        logical, intent(in) :: observed(:) !! True for observed response categories used for fitting.
        real(dp), intent(in) :: x(:, :) !! Complete numeric predictors without an intercept column.
        logical, intent(in) :: where(:) !! True for rows where categorical imputations are requested.
        integer, intent(in), value :: ncat !! Number of nominal categories.
        type(mice_rng_state), intent(inout) :: rng !! RNG state used for categorical probability draws.
        integer, allocatable, intent(out) :: imputed(:) !! One-based category draws ordered by true entries of `where`.
        integer, intent(out) :: info !! `mice_ok` on success or a package status code on failure.

        real(dp), allocatable :: xa(:, :), weights(:), design(:, :), coef(:, :), probs(:, :), one(:, :)
        integer, allocatable :: ya(:)
        integer :: i, j, naug

        if (size(observed) /= size(category) .or. size(where) /= size(category) .or. size(x, 1) /= size(category)) then
            info = mice_invalid_shape
            return
        end if
        if (count(observed) < 1 .or. ncat < 2) then
            info = mice_no_observed
            return
        end if
        call augment_multinomial(category, observed, x, ncat, xa, ya, weights)
        naug = size(ya)
        allocate(design(naug, size(x, 2) + 1))
        design(:, 1) = 1.0_dp
        if (size(x, 2) > 0) design(:, 2:) = xa
        call multinomial_fit(design, ya, weights, ncat, coef, info)
        if (info /= mice_ok .and. info /= mice_not_converged) return
        allocate(imputed(count(where)), one(1, size(design, 2)), probs(1, ncat))
        j = 0
        do i = 1, size(category)
            if (.not. where(i)) cycle
            j = j + 1
            one(1, 1) = 1.0_dp
            if (size(x, 2) > 0) one(1, 2:) = x(i, :)
            call multinomial_probabilities(one, coef, probs)
            imputed(j) = draw_category(probs(1, :), rng)
        end do
        info = mice_ok
    end subroutine impute_polyreg

    pure real(dp) function logistic(eta) result(p)
        real(dp), intent(in), value :: eta !! Linear predictor whose inverse-logit is returned.

        if (eta >= 0.0_dp) then
            p = 1.0_dp / (1.0_dp + exp(-eta))
        else
            p = exp(eta) / (1.0_dp + exp(eta))
        end if
    end function logistic

    pure subroutine weighted_crossprod(x, weights, result)
        real(dp), intent(in) :: x(:, :) !! Design matrix.
        real(dp), intent(in) :: weights(:) !! Row weights.
        real(dp), intent(out) :: result(:, :) !! Weighted cross-product `X' W X`.
        integer :: i, j, k

        result = 0.0_dp
        do i = 1, size(x, 1)
            do j = 1, size(x, 2)
                do k = 1, size(x, 2)
                    result(j, k) = result(j, k) + weights(i) * x(i, j) * x(i, k)
                end do
            end do
        end do
    end subroutine weighted_crossprod

    pure subroutine weighted_rhs(x, weights, y, result)
        real(dp), intent(in) :: x(:, :) !! Design matrix.
        real(dp), intent(in) :: weights(:) !! Row weights.
        real(dp), intent(in) :: y(:) !! Working response.
        real(dp), intent(out) :: result(:) !! Weighted right-hand side `X' W y`.
        integer :: i, j

        result = 0.0_dp
        do i = 1, size(x, 1)
            do j = 1, size(x, 2)
                result(j) = result(j) + weights(i) * x(i, j) * y(i)
            end do
        end do
    end subroutine weighted_rhs

    pure subroutine add_diagonal_ridge(matrix, ridge)
        real(dp), intent(inout) :: matrix(:, :) !! Symmetric matrix whose diagonal is regularized in place.
        real(dp), intent(in), value :: ridge !! Positive additive diagonal ridge.
        integer :: j

        do j = 1, min(size(matrix, 1), size(matrix, 2))
            matrix(j, j) = matrix(j, j) + ridge
        end do
    end subroutine add_diagonal_ridge

    subroutine augment_binary(y, observed, x, xa, ya, weights)
        real(dp), intent(in) :: y(:) !! Original binary response values.
        logical, intent(in) :: observed(:) !! Original observed-response mask.
        real(dp), intent(in) :: x(:, :) !! Original complete predictor matrix.
        real(dp), allocatable, intent(out) :: xa(:, :) !! Observed rows plus White-Daniel-Royston pseudo predictor rows.
        real(dp), allocatable, intent(out) :: ya(:) !! Binary responses for original observed and pseudo rows.
        real(dp), allocatable, intent(out) :: weights(:) !! Case weights, with small weights assigned to pseudo rows.

        real(dp), allocatable :: meanx(:), sdx(:), minx(:), maxx(:)
        integer :: i, j, k, nobs, p, row, nr

        nobs = count(observed)
        p = size(x, 2)
        if (p == 0 .or. count(.not. observed) == 1) then
            allocate(xa(nobs, p), ya(nobs), weights(nobs))
            row = 0
            do i = 1, size(y)
                if (.not. observed(i)) cycle
                row = row + 1
                if (p > 0) xa(row, :) = x(i, :)
                ya(row) = y(i)
            end do
            weights = 1.0_dp
            return
        end if
        nr = 4 * p
        allocate(meanx(p), sdx(p), minx(p), maxx(p))
        do j = 1, p
            meanx(j) = sum(x(:, j)) / real(size(x, 1), dp)
            if (size(x, 1) > 1) then
                sdx(j) = sqrt(sum((x(:, j) - meanx(j))**2) / real(size(x, 1) - 1, dp))
            else
                sdx(j) = 0.0_dp
            end if
            minx(j) = minval(x(:, j))
            maxx(j) = maxval(x(:, j))
        end do
        allocate(xa(nobs + nr, p), ya(nobs + nr), weights(nobs + nr))
        row = 0
        do i = 1, size(y)
            if (.not. observed(i)) cycle
            row = row + 1
            xa(row, :) = x(i, :)
            ya(row) = y(i)
            weights(row) = 1.0_dp
        end do
        do j = 1, p
            do k = 0, 1
                row = row + 1
                xa(row, :) = meanx
                xa(row, j) = min(maxx(j), max(minx(j), meanx(j) + 0.5_dp * sdx(j)))
                ya(row) = real(k, dp)
                weights(row) = real(p + 1, dp) / real(nr, dp)
                row = row + 1
                xa(row, :) = meanx
                xa(row, j) = min(maxx(j), max(minx(j), meanx(j) - 0.5_dp * sdx(j)))
                ya(row) = real(k, dp)
                weights(row) = real(p + 1, dp) / real(nr, dp)
            end do
        end do
    end subroutine augment_binary

    subroutine augment_multinomial(category, observed, x, ncat, xa, ya, weights)
        integer, intent(in) :: category(:) !! Original one-based nominal category codes.
        logical, intent(in) :: observed(:) !! Original observed-response mask.
        real(dp), intent(in) :: x(:, :) !! Original complete predictor matrix.
        integer, intent(in), value :: ncat !! Number of nominal categories.
        real(dp), allocatable, intent(out) :: xa(:, :) !! Observed and pseudo predictor rows.
        integer, allocatable, intent(out) :: ya(:) !! Category codes for observed and pseudo rows.
        real(dp), allocatable, intent(out) :: weights(:) !! Original unit and pseudo-row case weights.

        real(dp), allocatable :: meanx(:), sdx(:), minx(:), maxx(:)
        integer :: cat, i, j, nobs, nr, p, row

        nobs = count(observed)
        p = size(x, 2)
        if (p == 0 .or. count(.not. observed) == 1) then
            allocate(xa(nobs, p), ya(nobs), weights(nobs))
            row = 0
            do i = 1, size(category)
                if (.not. observed(i)) cycle
                row = row + 1
                if (p > 0) xa(row, :) = x(i, :)
                ya(row) = category(i)
            end do
            weights = 1.0_dp
            return
        end if
        nr = 2 * p * ncat
        allocate(meanx(p), sdx(p), minx(p), maxx(p))
        do j = 1, p
            meanx(j) = sum(x(:, j)) / real(size(x, 1), dp)
            if (size(x, 1) > 1) then
                sdx(j) = sqrt(sum((x(:, j) - meanx(j))**2) / real(size(x, 1) - 1, dp))
            else
                sdx(j) = 0.0_dp
            end if
            minx(j) = minval(x(:, j))
            maxx(j) = maxval(x(:, j))
        end do
        allocate(xa(nobs + nr, p), ya(nobs + nr), weights(nobs + nr))
        row = 0
        do i = 1, size(category)
            if (.not. observed(i)) cycle
            row = row + 1
            xa(row, :) = x(i, :)
            ya(row) = category(i)
            weights(row) = 1.0_dp
        end do
        do j = 1, p
            do cat = 1, ncat
                row = row + 1
                xa(row, :) = meanx
                xa(row, j) = min(maxx(j), max(minx(j), meanx(j) + 0.5_dp * sdx(j)))
                ya(row) = cat
                weights(row) = real(p + 1, dp) / real(nr, dp)
                row = row + 1
                xa(row, :) = meanx
                xa(row, j) = min(maxx(j), max(minx(j), meanx(j) - 0.5_dp * sdx(j)))
                ya(row) = cat
                weights(row) = real(p + 1, dp) / real(nr, dp)
            end do
        end do
    end subroutine augment_multinomial

    pure subroutine multinomial_probabilities(x, coef, prob)
        real(dp), intent(in) :: x(:, :) !! Multinomial design matrix.
        real(dp), intent(in) :: coef(:, :) !! Coefficients for nonbaseline categories.
        real(dp), intent(out) :: prob(:, :) !! Rowwise softmax probabilities including the baseline category in column one.
        real(dp) :: denom, maxeta
        integer :: i, k

        do i = 1, size(x, 1)
            maxeta = 0.0_dp
            do k = 1, size(coef, 2)
                maxeta = max(maxeta, dot_product(x(i, :), coef(:, k)))
            end do
            prob(i, 1) = exp(-maxeta)
            denom = prob(i, 1)
            do k = 1, size(coef, 2)
                prob(i, k + 1) = exp(dot_product(x(i, :), coef(:, k)) - maxeta)
                denom = denom + prob(i, k + 1)
            end do
            prob(i, :) = prob(i, :) / denom
        end do
    end subroutine multinomial_probabilities

    integer function draw_category(prob, rng) result(category)
        real(dp), intent(in) :: prob(:) !! Category probabilities expected to sum approximately to one.
        type(mice_rng_state), intent(inout) :: rng !! RNG state used for one uniform category draw.
        real(dp) :: cumulative, u
        integer :: k

        u = rng_uniform(rng)
        cumulative = 0.0_dp
        category = size(prob)
        do k = 1, size(prob)
            cumulative = cumulative + prob(k)
            if (u <= cumulative) then
                category = k
                return
            end if
        end do
    end function draw_category

end module mice_categorical
