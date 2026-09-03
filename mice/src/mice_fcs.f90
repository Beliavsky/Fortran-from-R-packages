! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from R package mice 3.19.0 by Stef van Buuren, Karin Groothuis-Oudshoorn,
! and mice contributors; see NOTICE.md and PROVENANCE.md for attribution.
! Numerical chained-equations sampler derived from mice 3.19.0 sampler.R and univariate imputers.
module mice_fcs
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use iso_fortran_env, only : int64
    use r_kinds, only : dp
    use mice_categorical, only : impute_logreg, impute_logreg_boot, impute_polyreg
    use mice_polr, only : impute_polr
    use mice_lda, only : impute_lda
    use mice_impute_continuous, only : impute_mean, impute_sample, impute_pmm
    use mice_midastouch, only : impute_midastouch
    use mice_regression, only : impute_norm, impute_norm_nob, impute_norm_predict, impute_norm_boot
    use mice_rng, only : mice_rng_state, rng_seed
    use mice_status, only : mice_ok, mice_invalid_argument, mice_invalid_shape, mice_not_converged
    implicit none
    private

    integer, parameter, public :: mice_method_skip = 0
    integer, parameter, public :: mice_method_mean = 1
    integer, parameter, public :: mice_method_sample = 2
    integer, parameter, public :: mice_method_norm = 3
    integer, parameter, public :: mice_method_norm_nob = 4
    integer, parameter, public :: mice_method_norm_predict = 5
    integer, parameter, public :: mice_method_norm_boot = 6
    integer, parameter, public :: mice_method_pmm = 7
    integer, parameter, public :: mice_method_logreg = 8
    integer, parameter, public :: mice_method_logreg_boot = 9
    integer, parameter, public :: mice_method_polyreg = 10
    integer, parameter, public :: mice_method_polr = 11
    integer, parameter, public :: mice_method_midastouch = 12
    integer, parameter, public :: mice_method_lda = 13

    type, public :: mice_fcs_result
        real(dp), allocatable :: completed(:, :, :)
        real(dp), allocatable :: chain_mean(:, :, :)
        real(dp), allocatable :: chain_variance(:, :, :)
    end type mice_fcs_result

    public :: mice_fcs_impute

contains

    subroutine mice_fcs_impute(data, methods, predictor_matrix, maxit, m, seed, result, info, category_count, donors)
        real(dp), intent(in) :: data(:, :) !! Numeric incomplete data matrix with IEEE NaNs marking cells to impute.
        integer, intent(in) :: methods(:) !! Per-column imputation method codes exported by this module.
        integer, intent(in) :: predictor_matrix(:, :) !! Target-by-predictor inclusion matrix; nonzero entries select predictors.
        integer, intent(in), value :: maxit !! Number of Gibbs/FCS iterations for each imputation chain; must be nonnegative.
        integer, intent(in), value :: m !! Number of independently initialized completed datasets to generate; must be positive.
        integer(int64), intent(in), value :: seed !! Base deterministic RNG seed; each chain receives a reproducible derived seed.
        type(mice_fcs_result), intent(out) :: result !! Completed datasets and per-iteration chain means/variances of imputed cells.
        integer, intent(out) :: info !! `mice_ok` on success or the first method/argument status code encountered.
        integer, intent(in), optional :: category_count(:) !! Categories per column for polyreg/polr; inferred when absent.
        integer, intent(in), optional :: donors !! PMM donor-pool size used for every PMM target; default five.

        logical, allocatable :: missing(:, :), observed(:), where(:)
        real(dp), allocatable :: current(:, :), x(:, :), y(:), imputed(:), initial(:)
        integer, allocatable :: predictors(:), category(:), cat_imp(:)
        type(mice_rng_state) :: rng
        integer :: chain, donor_count, i, iter, j, k, n, ncat, p, np, status

        n = size(data, 1)
        p = size(data, 2)
        if (size(methods) /= p .or. size(predictor_matrix, 1) /= p .or. size(predictor_matrix, 2) /= p) then
            info = mice_invalid_shape
            return
        end if
        if (present(category_count)) then
            if (size(category_count) /= p) then
                info = mice_invalid_shape
                return
            end if
        end if
        if (maxit < 0 .or. m < 1 .or. any(methods < mice_method_skip) .or. any(methods > mice_method_lda)) then
            info = mice_invalid_argument
            return
        end if
        donor_count = 5
        if (present(donors)) donor_count = donors
        if (donor_count < 1) then
            info = mice_invalid_argument
            return
        end if

        allocate(missing(n, p), result%completed(n, p, m), &
                 result%chain_mean(maxit, p, m), result%chain_variance(maxit, p, m))
        do j = 1, p
            do i = 1, n
                missing(i, j) = ieee_is_nan(data(i, j))
            end do
        end do
        result%chain_mean = 0.0_dp
        result%chain_variance = 0.0_dp

        do chain = 1, m
            call rng_seed(rng, seed + int(104729 * (chain - 1), int64))
            allocate(current(n, p))
            current = data
            do j = 1, p
                if (.not. any(missing(:, j))) cycle
                allocate(y(n), observed(n), where(n))
                y = current(:, j)
                observed = .not. missing(:, j)
                where = missing(:, j)
                call impute_sample(y, observed, where, rng, initial, status)
                if (status /= mice_ok) then
                    info = status
                    return
                end if
                call fill_where(current(:, j), where, initial)
                deallocate(y, observed, where, initial)
            end do

            do iter = 1, maxit
                do j = 1, p
                    if (methods(j) == mice_method_skip .or. .not. any(missing(:, j))) cycle
                    np = count(predictor_matrix(j, :) /= 0)
                    allocate(predictors(np), x(n, np), y(n), observed(n), where(n))
                    k = 0
                    do i = 1, p
                        if (predictor_matrix(j, i) == 0) cycle
                        k = k + 1
                        predictors(k) = i
                        x(:, k) = current(:, i)
                    end do
                    y = current(:, j)
                    observed = .not. missing(:, j)
                    where = missing(:, j)
                    select case (methods(j))
                    case (mice_method_mean)
                        call impute_mean(y, observed, where, imputed, status)
                    case (mice_method_sample)
                        call impute_sample(y, observed, where, rng, imputed, status)
                    case (mice_method_norm)
                        call impute_norm(y, observed, x, where, rng, imputed, status)
                    case (mice_method_norm_nob)
                        call impute_norm_nob(y, observed, x, where, rng, imputed, status)
                    case (mice_method_norm_predict)
                        call impute_norm_predict(y, observed, x, where, imputed, status)
                    case (mice_method_norm_boot)
                        call impute_norm_boot(y, observed, x, where, rng, imputed, status)
                    case (mice_method_pmm)
                        call impute_pmm(y, observed, x, where, rng, imputed, status, donors=donor_count)
                    case (mice_method_midastouch)
                        call impute_midastouch(y, observed, x, where, rng, imputed, status)
                    case (mice_method_logreg)
                        call impute_logreg(y, observed, x, where, rng, imputed, status)
                    case (mice_method_logreg_boot)
                        call impute_logreg_boot(y, observed, x, where, rng, imputed, status)
                    case (mice_method_polyreg, mice_method_polr, mice_method_lda)
                        allocate(category(n))
                        category = nint(y)
                        if (present(category_count)) then
                            ncat = category_count(j)
                        else
                            ncat = infer_category_count(category, observed)
                        end if
                        if (methods(j) == mice_method_polyreg) then
                            call impute_polyreg(category, observed, x, where, ncat, rng, cat_imp, status)
                        else if (methods(j) == mice_method_polr) then
                            call impute_polr(category, observed, x, where, ncat, rng, cat_imp, status)
                        else
                            call impute_lda(category, observed, x, where, ncat, rng, cat_imp, status)
                        end if
                        if (status == mice_ok .or. status == mice_not_converged) then
                            allocate(imputed(size(cat_imp)))
                            imputed = real(cat_imp, dp)
                        end if
                        deallocate(category)
                        if (allocated(cat_imp)) deallocate(cat_imp)
                    case default
                        status = mice_invalid_argument
                    end select
                    if (status /= mice_ok) then
                        ! A nonconverged categorical fit may still yield usable probabilities.
                        if (status /= mice_not_converged) then
                            info = status
                            return
                        end if
                    end if
                    call fill_where(current(:, j), where, imputed)
                    call record_chain(current(:, j), where, result%chain_mean(iter, j, chain), &
                                      result%chain_variance(iter, j, chain))
                    deallocate(predictors, x, y, observed, where, imputed)
                end do
            end do
            result%completed(:, :, chain) = current
            deallocate(current)
        end do
        info = mice_ok
    end subroutine mice_fcs_impute

    pure subroutine fill_where(column, where, values)
        real(dp), intent(inout) :: column(:) !! Data column whose selected entries are overwritten.
        logical, intent(in) :: where(:) !! True entries identify positions to receive `values`.
        real(dp), intent(in) :: values(:) !! Replacement values ordered by true entries of `where`.
        integer :: i, k

        k = 0
        do i = 1, size(column)
            if (.not. where(i)) cycle
            k = k + 1
            column(i) = values(k)
        end do
    end subroutine fill_where

    pure subroutine record_chain(column, where, mean_value, variance_value)
        real(dp), intent(in) :: column(:) !! Current completed target column.
        logical, intent(in) :: where(:) !! Original missingness mask for the target column.
        real(dp), intent(out) :: mean_value !! Mean of current imputed cells, or zero if none are selected.
        real(dp), intent(out) :: variance_value !! Sample variance of current imputed cells, or zero for fewer than two cells.
        real(dp) :: delta
        integer :: i, n

        n = count(where)
        if (n < 1) then
            mean_value = 0.0_dp
            variance_value = 0.0_dp
            return
        end if
        mean_value = 0.0_dp
        do i = 1, size(column)
            if (where(i)) mean_value = mean_value + column(i)
        end do
        mean_value = mean_value / real(n, dp)
        variance_value = 0.0_dp
        if (n < 2) return
        do i = 1, size(column)
            if (.not. where(i)) cycle
            delta = column(i) - mean_value
            variance_value = variance_value + delta * delta
        end do
        variance_value = variance_value / real(n - 1, dp)
    end subroutine record_chain

    pure integer function infer_category_count(category, observed) result(ncat)
        integer, intent(in) :: category(:) !! One-based category codes with arbitrary placeholders in unobserved rows.
        logical, intent(in) :: observed(:) !! True entries select valid response codes.
        integer :: i

        ncat = 0
        do i = 1, size(category)
            if (observed(i)) ncat = max(ncat, category(i))
        end do
    end function infer_category_count

end module mice_fcs
