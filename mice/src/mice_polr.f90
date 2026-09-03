! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from R package mice 3.19.0 by Stef van Buuren, Karin Groothuis-Oudshoorn,
! and mice contributors; see NOTICE.md and PROVENANCE.md for attribution.
! Proportional-odds imputation derived from mice.impute.polr.R and MASS::polr semantics.
module mice_polr
    use r_kinds, only : dp
    use mice_categorical, only : augment_multinomial, draw_category
    use mice_rng, only : mice_rng_state
    use mice_status, only : mice_ok, mice_invalid_argument, mice_invalid_shape, mice_no_observed, &
                            mice_not_converged
    implicit none
    private

    public :: proportional_odds_fit
    public :: proportional_odds_probabilities
    public :: impute_polr

contains

    subroutine proportional_odds_fit(x, category, weights, ncat, beta, thresholds, info, maxit, tol)
        real(dp), intent(in) :: x(:, :) !! Predictor matrix without an intercept; ordered cut points provide the intercepts.
        integer, intent(in) :: category(:) !! One-based ordered response codes in the range `1:ncat`.
        real(dp), intent(in) :: weights(:) !! Nonnegative case weights for the ordered likelihood.
        integer, intent(in), value :: ncat !! Number of ordered categories; must be at least two.
        real(dp), allocatable, intent(out) :: beta(:) !! Fitted common-slope proportional-odds coefficients.
        real(dp), allocatable, intent(out) :: thresholds(:) !! Increasing fitted cumulative-logit cut points.
        integer, intent(out) :: info !! `mice_ok`, `mice_not_converged`, or an argument status code.
        integer, intent(in), optional :: maxit !! Maximum BFGS iterations; default 200.
        real(dp), intent(in), optional :: tol !! Relative gradient tolerance; default `1e-8`.

        real(dp), allocatable :: theta(:), gradient(:), new_gradient(:), direction(:), candidate(:)
        real(dp), allocatable :: h_inv(:, :), s(:), ydiff(:), identity(:, :), work(:, :)
        real(dp) :: alpha_step, f, f_new, gnorm, rho, sy, tolerance, yy
        integer :: i, iter, iterations, p, q

        if (size(category) /= size(x, 1) .or. size(weights) /= size(x, 1)) then
            info = mice_invalid_shape
            return
        end if
        if (ncat < 2 .or. size(x, 1) < 1 .or. any(weights < 0.0_dp)) then
            info = mice_invalid_argument
            return
        end if
        if (any(category < 1) .or. any(category > ncat) .or. sum(weights) <= 0.0_dp) then
            info = mice_invalid_argument
            return
        end if
        iterations = 200
        if (present(maxit)) iterations = maxit
        tolerance = 1.0e-8_dp
        if (present(tol)) tolerance = tol
        if (iterations < 1 .or. tolerance <= 0.0_dp) then
            info = mice_invalid_argument
            return
        end if

        p = size(x, 2)
        q = p + ncat - 1
        allocate(theta(q), gradient(q), new_gradient(q), direction(q), candidate(q))
        allocate(h_inv(q, q), s(q), ydiff(q), identity(q, q), work(q, q))
        call initial_parameters(category, weights, ncat, p, theta)
        identity = 0.0_dp
        do i = 1, q
            identity(i, i) = 1.0_dp
        end do
        h_inv = identity
        call ordered_objective_gradient(theta, x, category, weights, ncat, f, gradient)

        do iter = 1, iterations
            gnorm = maxval(abs(gradient))
            if (gnorm <= tolerance * (1.0_dp + abs(f))) then
                call unpack_parameters(theta, p, beta, thresholds)
                info = mice_ok
                return
            end if
            direction = -matmul(h_inv, gradient)
            if (dot_product(direction, gradient) >= 0.0_dp) direction = -gradient
            alpha_step = 1.0_dp
            do
                candidate = theta + alpha_step * direction
                call ordered_objective_gradient(candidate, x, category, weights, ncat, f_new, new_gradient)
                if (f_new <= f + 1.0e-4_dp * alpha_step * dot_product(gradient, direction)) exit
                alpha_step = 0.5_dp * alpha_step
                if (alpha_step < 1.0e-10_dp) exit
            end do
            if (alpha_step < 1.0e-10_dp) then
                call unpack_parameters(theta, p, beta, thresholds)
                info = mice_not_converged
                return
            end if
            s = candidate - theta
            ydiff = new_gradient - gradient
            sy = dot_product(s, ydiff)
            yy = dot_product(ydiff, ydiff)
            if (sy > sqrt(epsilon(1.0_dp)) * sqrt(max(dot_product(s, s) * yy, tiny(1.0_dp)))) then
                rho = 1.0_dp / sy
                work = identity - rho * outer_product(s, ydiff)
                h_inv = matmul(matmul(work, h_inv), transpose(work)) + rho * outer_product(s, s)
            else
                h_inv = identity
            end if
            theta = candidate
            gradient = new_gradient
            f = f_new
        end do
        call unpack_parameters(theta, p, beta, thresholds)
        info = mice_not_converged
    end subroutine proportional_odds_fit

    pure subroutine proportional_odds_probabilities(x, beta, thresholds, probabilities)
        real(dp), intent(in) :: x(:, :) !! Predictor matrix without an intercept.
        real(dp), intent(in) :: beta(:) !! Common-slope proportional-odds coefficients.
        real(dp), intent(in) :: thresholds(:) !! Increasing cumulative-logit cut points.
        real(dp), intent(out) :: probabilities(:, :) !! Row-by-category probabilities; second extent is `size(thresholds)+1`.

        real(dp) :: eta, lower, upper, total
        integer :: c, i, ncat

        ncat = size(thresholds) + 1
        probabilities = 0.0_dp
        do i = 1, size(x, 1)
            eta = dot_product(x(i, :), beta)
            lower = 0.0_dp
            do c = 1, ncat - 1
                upper = inverse_logit(thresholds(c) - eta)
                probabilities(i, c) = max(upper - lower, 0.0_dp)
                lower = upper
            end do
            probabilities(i, ncat) = max(1.0_dp - lower, 0.0_dp)
            total = sum(probabilities(i, :))
            if (total > 0.0_dp) then
                probabilities(i, :) = probabilities(i, :) / total
            else
                probabilities(i, :) = 1.0_dp / real(ncat, dp)
            end if
        end do
    end subroutine proportional_odds_probabilities

    subroutine impute_polr(category, observed, x, where, ncat, rng, imputed, info)
        integer, intent(in) :: category(:) !! Incomplete ordered one-based category codes.
        logical, intent(in) :: observed(:) !! True for observed categories used for fitting.
        real(dp), intent(in) :: x(:, :) !! Complete numeric predictors without an intercept column.
        logical, intent(in) :: where(:) !! True for rows where ordered-category imputations are requested.
        integer, intent(in), value :: ncat !! Number of ordered categories.
        type(mice_rng_state), intent(inout) :: rng !! RNG state used for category draws.
        integer, allocatable, intent(out) :: imputed(:) !! One-based ordered-category draws for true entries of `where`.
        integer, intent(out) :: info !! `mice_ok`, `mice_not_converged`, or a package status code.

        real(dp), allocatable :: xa(:, :), weights(:), beta(:), thresholds(:), one(:, :), probs(:, :)
        integer, allocatable :: ya(:)
        integer :: fit_info, i, j

        if (size(observed) /= size(category) .or. size(where) /= size(category) .or. &
            size(x, 1) /= size(category)) then
            info = mice_invalid_shape
            return
        end if
        if (ncat < 2) then
            info = mice_invalid_argument
            return
        end if
        if (count(observed) < 1) then
            info = mice_no_observed
            return
        end if
        call augment_multinomial(category, observed, x, ncat, xa, ya, weights)
        call proportional_odds_fit(xa, ya, weights, ncat, beta, thresholds, fit_info)
        if (fit_info /= mice_ok .and. fit_info /= mice_not_converged) then
            info = fit_info
            return
        end if
        allocate(imputed(count(where)), one(1, size(x, 2)), probs(1, ncat))
        j = 0
        do i = 1, size(category)
            if (.not. where(i)) cycle
            j = j + 1
            if (size(x, 2) > 0) one(1, :) = x(i, :)
            call proportional_odds_probabilities(one, beta, thresholds, probs)
            imputed(j) = draw_category(probs(1, :), rng)
        end do
        info = fit_info
    end subroutine impute_polr

    subroutine initial_parameters(category, weights, ncat, p, theta)
        integer, intent(in) :: category(:) !! Ordered response codes used to initialize cumulative proportions.
        real(dp), intent(in) :: weights(:) !! Case weights corresponding to `category`.
        integer, intent(in), value :: ncat !! Number of ordered response categories.
        integer, intent(in), value :: p !! Number of common-slope coefficients preceding threshold parameters.
        real(dp), intent(out) :: theta(:) !! Unconstrained initial coefficient/threshold parameter vector.

        real(dp), allocatable :: alpha(:)
        real(dp) :: cumulative, eps, total
        integer :: c

        theta = 0.0_dp
        allocate(alpha(ncat - 1))
        total = sum(weights)
        eps = 1.0e-6_dp
        cumulative = 0.0_dp
        do c = 1, ncat - 1
            cumulative = cumulative + sum(weights, mask=category == c)
            alpha(c) = logit(min(1.0_dp - eps, max(eps, cumulative / total)))
            if (c > 1) alpha(c) = max(alpha(c), alpha(c - 1) + 1.0e-3_dp)
        end do
        theta(p + 1) = alpha(1)
        do c = 2, ncat - 1
            theta(p + c) = log(max(alpha(c) - alpha(c - 1), 1.0e-3_dp))
        end do
    end subroutine initial_parameters

    subroutine ordered_objective_gradient(theta, x, category, weights, ncat, objective, gradient)
        real(dp), intent(in) :: theta(:) !! Unconstrained common-slope and transformed-threshold parameters.
        real(dp), intent(in) :: x(:, :) !! Predictor matrix without an intercept.
        integer, intent(in) :: category(:) !! One-based ordered category codes.
        real(dp), intent(in) :: weights(:) !! Nonnegative likelihood weights.
        integer, intent(in), value :: ncat !! Number of ordered categories.
        real(dp), intent(out) :: objective !! Weighted negative log likelihood.
        real(dp), intent(out) :: gradient(:) !! Gradient of `objective` with respect to `theta`.

        real(dp), allocatable :: alpha(:), g_alpha(:), g_beta(:)
        real(dp) :: eta, fl, fu, lower, pcat, upper, wi
        integer :: c, i, j, p

        p = size(x, 2)
        allocate(alpha(ncat - 1), g_alpha(ncat - 1), g_beta(p))
        call transformed_thresholds(theta(p + 1:), alpha)
        objective = 0.0_dp
        g_alpha = 0.0_dp
        g_beta = 0.0_dp
        do i = 1, size(x, 1)
            wi = weights(i)
            if (wi <= 0.0_dp) cycle
            c = category(i)
            eta = dot_product(x(i, :), theta(1:p))
            if (c == 1) then
                lower = 0.0_dp
                fl = 0.0_dp
            else
                lower = inverse_logit(alpha(c - 1) - eta)
                fl = lower * (1.0_dp - lower)
            end if
            if (c == ncat) then
                upper = 1.0_dp
                fu = 0.0_dp
            else
                upper = inverse_logit(alpha(c) - eta)
                fu = upper * (1.0_dp - upper)
            end if
            pcat = max(upper - lower, tiny(1.0_dp))
            objective = objective - wi * log(pcat)
            if (p > 0) g_beta = g_beta + wi * (fu - fl) * x(i, :) / pcat
            if (c < ncat) g_alpha(c) = g_alpha(c) - wi * fu / pcat
            if (c > 1) g_alpha(c - 1) = g_alpha(c - 1) + wi * fl / pcat
        end do
        if (p > 0) gradient(1:p) = g_beta
        gradient(p + 1) = sum(g_alpha)
        do j = 2, ncat - 1
            gradient(p + j) = exp(theta(p + j)) * sum(g_alpha(j:))
        end do
    end subroutine ordered_objective_gradient

    subroutine unpack_parameters(theta, p, beta, thresholds)
        real(dp), intent(in) :: theta(:) !! Unconstrained fitted parameter vector.
        integer, intent(in), value :: p !! Number of common-slope coefficients.
        real(dp), allocatable, intent(out) :: beta(:) !! Extracted common-slope coefficients.
        real(dp), allocatable, intent(out) :: thresholds(:) !! Extracted ordered thresholds.

        allocate(beta(p), thresholds(size(theta) - p))
        if (p > 0) beta = theta(1:p)
        call transformed_thresholds(theta(p + 1:), thresholds)
    end subroutine unpack_parameters

    pure subroutine transformed_thresholds(gamma, alpha)
        real(dp), intent(in) :: gamma(:) !! Unconstrained threshold parameters; later elements encode log gaps.
        real(dp), intent(out) :: alpha(:) !! Strictly increasing threshold sequence.
        integer :: j

        if (size(gamma) < 1) return
        alpha(1) = gamma(1)
        do j = 2, size(gamma)
            alpha(j) = alpha(j - 1) + exp(min(gamma(j), log(huge(1.0_dp)) / 4.0_dp))
        end do
    end subroutine transformed_thresholds

    pure real(dp) function inverse_logit(value) result(probability)
        real(dp), intent(in), value :: value !! Real-valued logit argument.

        if (value >= 0.0_dp) then
            probability = 1.0_dp / (1.0_dp + exp(-value))
        else
            probability = exp(value) / (1.0_dp + exp(value))
        end if
    end function inverse_logit

    pure real(dp) function logit(probability) result(value)
        real(dp), intent(in), value :: probability !! Probability strictly between zero and one.

        value = log(probability / (1.0_dp - probability))
    end function logit

    pure function outer_product(a, b) result(matrix)
        real(dp), intent(in) :: a(:) !! Left vector in the outer product.
        real(dp), intent(in) :: b(:) !! Right vector in the outer product.
        real(dp) :: matrix(size(a), size(b))
        integer :: i, j

        do j = 1, size(b)
            do i = 1, size(a)
                matrix(i, j) = a(i) * b(j)
            end do
        end do
    end function outer_product

end module mice_polr
