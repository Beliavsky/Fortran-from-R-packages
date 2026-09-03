! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from R package mice 3.19.0 by Stef van Buuren, Karin Groothuis-Oudshoorn,
! and mice contributors; see NOTICE.md and PROVENANCE.md for attribution.
! Computational translation derived from mice 3.19.0 continuous univariate imputers.
module mice_impute_continuous
    use r_kinds, only : dp
    use mice_matching, only : matchindex, matcher
    use mice_regression, only : mice_norm_draw, norm_draw
    use mice_rng, only : mice_rng_state, rng_integer, rng_normal
    use mice_status, only : mice_ok, mice_invalid_argument, mice_invalid_shape, mice_no_observed
    implicit none
    private

    public :: impute_mean
    public :: impute_sample
    public :: impute_pmm
    public :: impute_random_indicator
    public :: impute_quadratic

contains

    pure subroutine impute_mean(y, observed, where, imputed, info)
        real(dp), intent(in) :: y(:) !! Incomplete numeric response vector.
        logical, intent(in) :: observed(:) !! True for response values used to calculate the observed mean.
        logical, intent(in) :: where(:) !! True for rows at which imputations are requested.
        real(dp), allocatable, intent(out) :: imputed(:) !! Mean imputations ordered by true entries of `where`.
        integer, intent(out) :: info !! `mice_ok` on success or a package status code on failure.

        real(dp) :: mu
        integer :: i, nobs

        if (size(observed) /= size(y) .or. size(where) /= size(y)) then
            info = mice_invalid_shape
            return
        end if
        nobs = count(observed)
        if (nobs < 1) then
            info = mice_no_observed
            return
        end if
        mu = 0.0_dp
        do i = 1, size(y)
            if (observed(i)) mu = mu + y(i)
        end do
        mu = mu / real(nobs, dp)
        allocate(imputed(count(where)))
        imputed = mu
        info = mice_ok
    end subroutine impute_mean

    subroutine impute_sample(y, observed, where, rng, imputed, info)
        real(dp), intent(in) :: y(:) !! Incomplete numeric response vector.
        logical, intent(in) :: observed(:) !! True for values eligible to be sampled as donors.
        logical, intent(in) :: where(:) !! True for rows at which imputations are requested.
        type(mice_rng_state), intent(inout) :: rng !! RNG state used to sample observed donor values with replacement.
        real(dp), allocatable, intent(out) :: imputed(:) !! Sampled donor values ordered by true entries of `where`.
        integer, intent(out) :: info !! `mice_ok` on success or a package status code on failure.

        real(dp), allocatable :: donors(:)
        integer :: i, j, nobs

        if (size(observed) /= size(y) .or. size(where) /= size(y)) then
            info = mice_invalid_shape
            return
        end if
        nobs = count(observed)
        allocate(imputed(count(where)))
        if (nobs < 1) then
            do i = 1, size(imputed)
                imputed(i) = rng_normal(rng)
            end do
            info = mice_ok
            return
        end if
        allocate(donors(nobs))
        j = 0
        do i = 1, size(y)
            if (.not. observed(i)) cycle
            j = j + 1
            donors(j) = y(i)
        end do
        do i = 1, size(imputed)
            imputed(i) = donors(rng_integer(rng, 1, nobs))
        end do
        info = mice_ok
    end subroutine impute_sample

    subroutine impute_pmm(y, observed, x, where, rng, imputed, info, donors, matchtype, ridge, use_matcher)
        real(dp), intent(in) :: y(:) !! Incomplete numeric response vector; observed values form the PMM donor pool.
        logical, intent(in) :: observed(:) !! True for rows used to estimate the predictive model and provide donors.
        real(dp), intent(in) :: x(:, :) !! Complete predictor matrix without an intercept column.
        logical, intent(in) :: where(:) !! True for rows at which PMM imputations are requested.
        type(mice_rng_state), intent(inout) :: rng !! RNG state used for posterior regression draws and donor selection.
        real(dp), allocatable, intent(out) :: imputed(:) !! PMM donor values ordered by true entries of `where`.
        integer, intent(out) :: info !! `mice_ok` on success or a package status code on failure.
        integer, intent(in), optional :: donors !! Donor-pool size; defaults to five and is clipped to available donors.
        integer, intent(in), optional :: matchtype !! PMM distance type: 0 predicted/predicted, 1 predicted/drawn, 2 drawn/drawn.
        real(dp), intent(in), optional :: ridge !! Relative ridge penalty for the normal regression draw; default `1e-5`.
        logical, intent(in), optional :: use_matcher !! Use the legacy jitter-based matcher instead of `matchindex` when true.

        real(dp), allocatable :: design(:, :), yobs(:), obs_metric(:), mis_metric(:)
        integer, allocatable :: idx(:)
        type(mice_norm_draw) :: draw
        real(dp) :: ridge_value
        integer :: donor_count, i, j, mt, nobs, nmis
        logical :: old_matcher

        if (size(observed) /= size(y) .or. size(where) /= size(y) .or. size(x, 1) /= size(y)) then
            info = mice_invalid_shape
            return
        end if
        nobs = count(observed)
        nmis = count(where)
        if (nobs < 1) then
            info = mice_no_observed
            return
        end if
        donor_count = 5
        if (present(donors)) donor_count = donors
        if (donor_count < 1) then
            info = mice_invalid_argument
            return
        end if
        mt = 1
        if (present(matchtype)) mt = matchtype
        if (mt < 0 .or. mt > 2) then
            info = mice_invalid_argument
            return
        end if
        ridge_value = 1.0e-5_dp
        if (present(ridge)) ridge_value = ridge
        old_matcher = .false.
        if (present(use_matcher)) old_matcher = use_matcher

        allocate(design(size(y), size(x, 2) + 1))
        design(:, 1) = 1.0_dp
        if (size(x, 2) > 0) design(:, 2:) = x
        call norm_draw(y, observed, design, rng, ridge_value, draw, info)
        if (info /= mice_ok) return

        allocate(yobs(nobs), obs_metric(nobs), mis_metric(nmis), imputed(nmis))
        j = 0
        do i = 1, size(y)
            if (.not. observed(i)) cycle
            j = j + 1
            yobs(j) = y(i)
            select case (mt)
            case (0, 1)
                obs_metric(j) = dot_product(design(i, :), draw%coef)
            case (2)
                obs_metric(j) = dot_product(design(i, :), draw%beta)
            end select
        end do
        j = 0
        do i = 1, size(y)
            if (.not. where(i)) cycle
            j = j + 1
            select case (mt)
            case (0)
                mis_metric(j) = dot_product(design(i, :), draw%coef)
            case (1, 2)
                mis_metric(j) = dot_product(design(i, :), draw%beta)
            end select
        end do

        if (old_matcher) then
            call matcher(obs_metric, mis_metric, donor_count, rng, idx, info)
        else
            call matchindex(obs_metric, mis_metric, donor_count, rng, idx, info)
        end if
        if (info /= mice_ok) return
        do i = 1, nmis
            imputed(i) = yobs(idx(i))
        end do
        info = mice_ok
    end subroutine impute_pmm

    subroutine impute_random_indicator(y, observed, x, where, rng, imputed, info, ridge)
        real(dp), intent(in) :: y(:) !! Incomplete numeric response vector for random-indicator imputation.
        logical, intent(in) :: observed(:) !! True for observed response values used to fit the regression.
        real(dp), intent(in) :: x(:, :) !! Complete predictor matrix without an intercept or response indicator.
        logical, intent(in) :: where(:) !! True for rows at which imputations are requested.
        type(mice_rng_state), intent(inout) :: rng !! RNG state used for generated indicators and posterior regression draws.
        real(dp), allocatable, intent(out) :: imputed(:) !! Random-indicator imputations ordered by true entries of `where`.
        integer, intent(out) :: info !! `mice_ok` on success or a package status code on failure.
        real(dp), intent(in), optional :: ridge !! Relative ridge penalty for the normal regression; default `1e-5`.

        real(dp), allocatable :: extended(:, :), design(:, :)
        type(mice_norm_draw) :: draw
        real(dp) :: ridge_value
        integer :: i, j

        if (size(observed) /= size(y) .or. size(where) /= size(y) .or. size(x, 1) /= size(y)) then
            info = mice_invalid_shape
            return
        end if
        ridge_value = 1.0e-5_dp
        if (present(ridge)) ridge_value = ridge
        allocate(extended(size(y), size(x, 2) + 1))
        if (size(x, 2) > 0) extended(:, 1:size(x, 2)) = x
        do i = 1, size(y)
            if (observed(i)) then
                extended(i, size(extended, 2)) = real(rng_integer(rng, 0, 1), dp)
            else
                extended(i, size(extended, 2)) = real(rng_integer(rng, 0, 1), dp)
            end if
        end do
        allocate(design(size(y), size(extended, 2) + 1))
        design(:, 1) = 1.0_dp
        design(:, 2:) = extended
        call norm_draw(y, observed, design, rng, ridge_value, draw, info)
        if (info /= mice_ok) return
        allocate(imputed(count(where)))
        j = 0
        do i = 1, size(y)
            if (.not. where(i)) cycle
            j = j + 1
            imputed(j) = dot_product(design(i, :), draw%beta) + draw%sigma * rng_normal(rng)
        end do
    end subroutine impute_random_indicator

    subroutine impute_quadratic(y, observed, x, where, rng, imputed, info, donors, ridge)
        real(dp), intent(in) :: y(:) !! Numeric response interpreted as a squared quantity whose positive root is imputed by PMM.
        logical, intent(in) :: observed(:) !! True for observed response values.
        real(dp), intent(in) :: x(:, :) !! Complete predictor matrix without an intercept column.
        logical, intent(in) :: where(:) !! True for rows at which imputations are requested.
        type(mice_rng_state), intent(inout) :: rng !! RNG state used by PMM.
        real(dp), allocatable, intent(out) :: imputed(:) !! Nonnegative quadratic imputations ordered by true entries of `where`.
        integer, intent(out) :: info !! `mice_ok` on success or a package status code on failure.
        integer, intent(in), optional :: donors !! PMM donor-pool size; default five.
        real(dp), intent(in), optional :: ridge !! Relative ridge penalty; default `1e-5`.

        real(dp), allocatable :: roots(:), root_imp(:)
        integer :: i, donor_count
        real(dp) :: ridge_value

        allocate(roots(size(y)))
        do i = 1, size(y)
            roots(i) = sqrt(max(y(i), 0.0_dp))
        end do
        donor_count = 5
        if (present(donors)) donor_count = donors
        ridge_value = 1.0e-5_dp
        if (present(ridge)) ridge_value = ridge
        call impute_pmm(roots, observed, x, where, rng, root_imp, info, donors=donor_count, ridge=ridge_value)
        if (info /= mice_ok) return
        allocate(imputed(size(root_imp)))
        imputed = root_imp * root_imp
    end subroutine impute_quadratic

end module mice_impute_continuous
