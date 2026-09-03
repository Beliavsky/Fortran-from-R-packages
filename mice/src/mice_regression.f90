! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from R package mice 3.19.0 by Stef van Buuren, Karin Groothuis-Oudshoorn,
! and mice contributors; see NOTICE.md and PROVENANCE.md for attribution.
! Computational translation derived from mice 3.19.0, especially R/mice.impute.norm*.R.
module mice_regression
    use r_kinds, only : dp
    use r_linalg, only : inverse_matrix, solve_spd
    use mice_rng, only : mice_rng_state, rng_normal, rng_chisq, rng_integer
    use mice_status, only : mice_ok, mice_invalid_argument, mice_invalid_shape
    use mice_status, only : mice_singular, mice_no_observed
    implicit none
    private

    type, public :: mice_norm_draw
        real(dp), allocatable :: coef(:)
        real(dp), allocatable :: beta(:)
        real(dp), allocatable :: vcov(:, :)
        real(dp) :: sigma = 0.0_dp
        integer :: df = 0
        logical :: ridge_used = .false.
    end type mice_norm_draw

    public :: estimice
    public :: norm_draw
    public :: impute_norm
    public :: impute_norm_nob
    public :: impute_norm_predict
    public :: impute_norm_boot
    public :: chol_lower

contains

    subroutine estimice(x, y, ridge, coef, residuals, vcov, df, ridge_used, info)
        real(dp), intent(in) :: x(:, :) !! Complete design matrix with observations in rows and predictors in columns.
        real(dp), intent(in) :: y(:) !! Complete response vector with one element for each row of `x`.
        real(dp), intent(in), value :: ridge !! Relative ridge penalty multiplied by each diagonal of `X'X` when inversion fails.
        real(dp), allocatable, intent(out) :: coef(:) !! Estimated least-squares or ridge coefficients.
        real(dp), allocatable, intent(out) :: residuals(:) !! Residual vector `y - X*coef`.
        real(dp), allocatable, intent(out) :: vcov(:, :) !! Unscaled coefficient covariance matrix, approximately `(X'X)^-1`.
        integer, intent(out) :: df !! Residual degrees of freedom, lower-bounded at one as in upstream `estimice`.
        logical, intent(out) :: ridge_used !! True when a ridge penalty was needed to obtain the covariance inverse.
        integer, intent(out) :: info !! `mice_ok` on success or a package status code on failure.

        real(dp), allocatable :: xtx(:, :), rhs(:), xtx_reg(:, :)
        real(dp) :: scale
        integer :: j, n, p, la_info

        n = size(x, 1)
        p = size(x, 2)
        if (size(y) /= n) then
            info = mice_invalid_shape
            return
        end if
        if (n < 1 .or. p < 1 .or. ridge < 0.0_dp) then
            info = mice_invalid_argument
            return
        end if

        allocate(coef(p), residuals(n), vcov(p, p), xtx(p, p), rhs(p), xtx_reg(p, p))
        xtx = matmul(transpose(x), x)
        rhs = matmul(transpose(x), y)
        call inverse_matrix(xtx, vcov, la_info)
        ridge_used = la_info /= 0

        if (ridge_used) then
            xtx_reg = xtx
            scale = max(maxval(abs(diagonal_values(xtx))), 1.0_dp)
            do j = 1, p
                xtx_reg(j, j) = xtx_reg(j, j) + ridge * max(abs(xtx(j, j)), scale * epsilon(1.0_dp))
            end do
            call inverse_matrix(xtx_reg, vcov, la_info)
            if (la_info /= 0) then
                info = mice_singular
                return
            end if
            coef = matmul(vcov, rhs)
        else
            call solve_spd(xtx, rhs, coef, la_info)
            if (la_info /= 0) coef = matmul(vcov, rhs)
        end if

        residuals = y - matmul(x, coef)
        df = max(n - p, 1)
        info = mice_ok
    end subroutine estimice

    subroutine norm_draw(y, observed, x, rng, ridge, draw, info)
        real(dp), intent(in) :: y(:) !! Response vector containing observed and potentially placeholder missing entries.
        logical, intent(in) :: observed(:) !! True for rows used to estimate the regression model.
        real(dp), intent(in) :: x(:, :) !! Complete design matrix including any desired intercept column.
        type(mice_rng_state), intent(inout) :: rng !! RNG state used for chi-square and Gaussian posterior draws.
        real(dp), intent(in), value :: ridge !! Relative ridge penalty passed to `estimice`.
        type(mice_norm_draw), intent(out) :: draw !! Least-squares estimates and one Bayesian draw of coefficients and sigma.
        integer, intent(out) :: info !! `mice_ok` on success or a package status code on failure.

        real(dp), allocatable :: xo(:, :), yo(:), residuals(:), factor(:, :), z(:)
        real(dp) :: chisq, rss
        integer :: i, j, nobs, p, row

        if (size(observed) /= size(y) .or. size(x, 1) /= size(y)) then
            info = mice_invalid_shape
            return
        end if
        nobs = count(observed)
        p = size(x, 2)
        if (nobs < 1) then
            info = mice_no_observed
            return
        end if

        allocate(xo(nobs, p), yo(nobs))
        row = 0
        do i = 1, size(y)
            if (.not. observed(i)) cycle
            row = row + 1
            xo(row, :) = x(i, :)
            yo(row) = y(i)
        end do

        call estimice(xo, yo, ridge, draw%coef, residuals, draw%vcov, draw%df, draw%ridge_used, info)
        if (info /= mice_ok) return
        rss = dot_product(residuals, residuals)
        chisq = rng_chisq(rng, real(draw%df, dp))
        if (chisq <= tiny(1.0_dp)) chisq = tiny(1.0_dp)
        draw%sigma = sqrt(max(rss, 0.0_dp) / chisq)

        allocate(factor(p, p), z(p), draw%beta(p))
        call chol_lower(draw%vcov, factor, info)
        if (info /= mice_ok) return
        do j = 1, p
            z(j) = rng_normal(rng)
        end do
        draw%beta = draw%coef + draw%sigma * matmul(factor, z)
        info = mice_ok
    end subroutine norm_draw

    subroutine impute_norm(y, observed, x, where, rng, imputed, info, ridge)
        real(dp), intent(in) :: y(:) !! Incomplete numeric response vector.
        logical, intent(in) :: observed(:) !! True for rows used to estimate the normal imputation model.
        real(dp), intent(in) :: x(:, :) !! Complete predictor matrix without an intercept column.
        logical, intent(in) :: where(:) !! True for rows where an imputed value should be generated.
        type(mice_rng_state), intent(inout) :: rng !! RNG state used for posterior and residual draws.
        real(dp), allocatable, intent(out) :: imputed(:) !! Generated values in the order of true entries of `where`.
        integer, intent(out) :: info !! `mice_ok` on success or a package status code on failure.
        real(dp), intent(in), optional :: ridge !! Relative ridge penalty; default is `1e-5` as in upstream mice.

        real(dp), allocatable :: design(:, :)
        real(dp) :: ridge_value
        type(mice_norm_draw) :: draw
        integer :: i, k, n, p

        call make_design_with_intercept(x, design)
        ridge_value = 1.0e-5_dp
        if (present(ridge)) ridge_value = ridge
        call norm_draw(y, observed, design, rng, ridge_value, draw, info)
        if (info /= mice_ok) return
        n = size(y)
        p = size(design, 2)
        allocate(imputed(count(where)))
        k = 0
        do i = 1, n
            if (.not. where(i)) cycle
            k = k + 1
            imputed(k) = dot_product(design(i, 1:p), draw%beta) + draw%sigma * rng_normal(rng)
        end do
    end subroutine impute_norm

    subroutine impute_norm_nob(y, observed, x, where, rng, imputed, info, ridge)
        real(dp), intent(in) :: y(:) !! Incomplete numeric response vector.
        logical, intent(in) :: observed(:) !! True for rows used to fit the regression.
        real(dp), intent(in) :: x(:, :) !! Complete predictor matrix without an intercept column.
        logical, intent(in) :: where(:) !! True for rows where an imputed value should be generated.
        type(mice_rng_state), intent(inout) :: rng !! RNG state used only for residual Gaussian draws.
        real(dp), allocatable, intent(out) :: imputed(:) !! Generated values in the order of true entries of `where`.
        integer, intent(out) :: info !! `mice_ok` on success or a package status code on failure.
        real(dp), intent(in), optional :: ridge !! Relative ridge penalty; default is `1e-5`.

        real(dp), allocatable :: design(:, :), xo(:, :), yo(:), coef(:), residuals(:), vcov(:, :)
        real(dp) :: ridge_value, sigma
        integer :: df, i, k, nobs, p, row
        logical :: ridge_used

        call make_design_with_intercept(x, design)
        ridge_value = 1.0e-5_dp
        if (present(ridge)) ridge_value = ridge
        nobs = count(observed)
        p = size(design, 2)
        if (nobs < 1) then
            info = mice_no_observed
            return
        end if
        allocate(xo(nobs, p), yo(nobs))
        row = 0
        do i = 1, size(y)
            if (.not. observed(i)) cycle
            row = row + 1
            xo(row, :) = design(i, :)
            yo(row) = y(i)
        end do
        call estimice(xo, yo, ridge_value, coef, residuals, vcov, df, ridge_used, info)
        if (info /= mice_ok) return
        sigma = sqrt(dot_product(residuals, residuals) / real(max(nobs - p - 1, 1), dp))
        allocate(imputed(count(where)))
        k = 0
        do i = 1, size(y)
            if (.not. where(i)) cycle
            k = k + 1
            imputed(k) = dot_product(design(i, :), coef) + sigma * rng_normal(rng)
        end do
    end subroutine impute_norm_nob

    subroutine impute_norm_predict(y, observed, x, where, imputed, info, ridge)
        real(dp), intent(in) :: y(:) !! Incomplete numeric response vector.
        logical, intent(in) :: observed(:) !! True for rows used to fit the regression.
        real(dp), intent(in) :: x(:, :) !! Complete predictor matrix without an intercept column.
        logical, intent(in) :: where(:) !! True for rows where a deterministic prediction should be returned.
        real(dp), allocatable, intent(out) :: imputed(:) !! Predicted values in the order of true entries of `where`.
        integer, intent(out) :: info !! `mice_ok` on success or a package status code on failure.
        real(dp), intent(in), optional :: ridge !! Relative ridge penalty; default is `1e-5`.

        real(dp), allocatable :: design(:, :), xo(:, :), yo(:), coef(:), residuals(:), vcov(:, :)
        real(dp) :: ridge_value
        integer :: df, i, k, nobs, p, row
        logical :: ridge_used

        call make_design_with_intercept(x, design)
        ridge_value = 1.0e-5_dp
        if (present(ridge)) ridge_value = ridge
        nobs = count(observed)
        p = size(design, 2)
        if (nobs < 1) then
            info = mice_no_observed
            return
        end if
        allocate(xo(nobs, p), yo(nobs))
        row = 0
        do i = 1, size(y)
            if (.not. observed(i)) cycle
            row = row + 1
            xo(row, :) = design(i, :)
            yo(row) = y(i)
        end do
        call estimice(xo, yo, ridge_value, coef, residuals, vcov, df, ridge_used, info)
        if (info /= mice_ok) return
        allocate(imputed(count(where)))
        k = 0
        do i = 1, size(y)
            if (.not. where(i)) cycle
            k = k + 1
            imputed(k) = dot_product(design(i, :), coef)
        end do
    end subroutine impute_norm_predict

    subroutine impute_norm_boot(y, observed, x, where, rng, imputed, info, ridge)
        real(dp), intent(in) :: y(:) !! Incomplete numeric response vector.
        logical, intent(in) :: observed(:) !! True for rows eligible for bootstrap fitting.
        real(dp), intent(in) :: x(:, :) !! Complete predictor matrix without an intercept column.
        logical, intent(in) :: where(:) !! True for rows where an imputed value should be generated.
        type(mice_rng_state), intent(inout) :: rng !! RNG state used for bootstrap and residual draws.
        real(dp), allocatable, intent(out) :: imputed(:) !! Bootstrap-regression imputations ordered by true entries of `where`.
        integer, intent(out) :: info !! `mice_ok` on success or a package status code on failure.
        real(dp), intent(in), optional :: ridge !! Relative ridge penalty; default is `1e-5`.

        real(dp), allocatable :: design(:, :), xo(:, :), yo(:), bx(:, :), by(:)
        real(dp), allocatable :: coef(:), residuals(:), vcov(:, :)
        real(dp) :: ridge_value, sigma
        integer, allocatable :: observed_rows(:)
        integer :: df, i, j, k, nobs, p, pick, row
        logical :: ridge_used

        call make_design_with_intercept(x, design)
        ridge_value = 1.0e-5_dp
        if (present(ridge)) ridge_value = ridge
        nobs = count(observed)
        p = size(design, 2)
        if (nobs < 1) then
            info = mice_no_observed
            return
        end if
        allocate(observed_rows(nobs), xo(nobs, p), yo(nobs), bx(nobs, p), by(nobs))
        row = 0
        do i = 1, size(y)
            if (.not. observed(i)) cycle
            row = row + 1
            observed_rows(row) = i
            xo(row, :) = design(i, :)
            yo(row) = y(i)
        end do
        do j = 1, nobs
            pick = rng_integer(rng, 1, nobs)
            bx(j, :) = xo(pick, :)
            by(j) = yo(pick)
        end do
        call estimice(bx, by, ridge_value, coef, residuals, vcov, df, ridge_used, info)
        if (info /= mice_ok) return
        sigma = sqrt(dot_product(residuals, residuals) / real(max(nobs - p - 1, 1), dp))
        allocate(imputed(count(where)))
        k = 0
        do i = 1, size(y)
            if (.not. where(i)) cycle
            k = k + 1
            imputed(k) = dot_product(design(i, :), coef) + sigma * rng_normal(rng)
        end do
    end subroutine impute_norm_boot

    pure subroutine chol_lower(a, l, info)
        real(dp), intent(in) :: a(:, :) !! Symmetric positive-definite input matrix.
        real(dp), intent(out) :: l(:, :) !! Lower-triangular Cholesky factor satisfying approximately `A=L*L'`.
        integer, intent(out) :: info !! `mice_ok` on success or `mice_singular` if a nonpositive pivot is found.

        real(dp) :: s
        integer :: i, j, k, n

        n = size(a, 1)
        if (size(a, 2) /= n .or. size(l, 1) /= n .or. size(l, 2) /= n) then
            info = mice_invalid_shape
            return
        end if
        l = 0.0_dp
        do i = 1, n
            do j = 1, i
                s = 0.5_dp * (a(i, j) + a(j, i))
                do k = 1, j - 1
                    s = s - l(i, k) * l(j, k)
                end do
                if (i == j) then
                    if (s <= max(tiny(1.0_dp), epsilon(1.0_dp) * max(1.0_dp, abs(a(i, i))))) then
                        info = mice_singular
                        return
                    end if
                    l(i, j) = sqrt(s)
                else
                    l(i, j) = s / l(j, j)
                end if
            end do
        end do
        info = mice_ok
    end subroutine chol_lower

    pure subroutine make_design_with_intercept(x, design)
        real(dp), intent(in) :: x(:, :) !! Predictor matrix without an intercept column.
        real(dp), allocatable, intent(out) :: design(:, :) !! Matrix with a leading column of ones followed by `x`.

        allocate(design(size(x, 1), size(x, 2) + 1))
        design(:, 1) = 1.0_dp
        if (size(x, 2) > 0) design(:, 2:) = x
    end subroutine make_design_with_intercept

    pure function diagonal_values(a) result(d)
        real(dp), intent(in) :: a(:, :) !! Matrix whose main diagonal is requested.
        real(dp) :: d(min(size(a, 1), size(a, 2)))
        integer :: i

        do i = 1, size(d)
            d(i) = a(i, i)
        end do
    end function diagonal_values

end module mice_regression
