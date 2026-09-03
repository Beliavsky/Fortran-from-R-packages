! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from R package mice 3.19.0 by Stef van Buuren, Karin Groothuis-Oudshoorn,
! and mice contributors; see NOTICE.md and PROVENANCE.md for attribution.
! Computational translation derived from mice 3.19.0 missing-data diagnostics.
module mice_diagnostics
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use r_kinds, only : dp
    use mice_status, only : mice_ok, mice_invalid_argument, mice_invalid_shape
    implicit none
    private

    type, public :: md_pairs_result
        integer, allocatable :: rr(:, :)
        integer, allocatable :: rm(:, :)
        integer, allocatable :: mr(:, :)
        integer, allocatable :: mm(:, :)
    end type md_pairs_result

    type, public :: flux_result
        real(dp), allocatable :: pobs(:)
        real(dp), allocatable :: influx(:)
        real(dp), allocatable :: outflux(:)
        real(dp), allocatable :: ainb(:)
        real(dp), allocatable :: aout(:)
        real(dp), allocatable :: fico(:)
    end type flux_result

    type, public :: md_pattern_result
        integer, allocatable :: pattern(:, :)
        integer, allocatable :: frequency(:)
        integer, allocatable :: missing_per_pattern(:)
        integer, allocatable :: missing_per_column(:)
        integer, allocatable :: column_order(:)
    end type md_pattern_result

    public :: missing_mask
    public :: md_pairs
    public :: md_pairs_from_mask
    public :: flux
    public :: md_pattern
    public :: quickpred

contains

    pure subroutine missing_mask(data, missing)
        real(dp), intent(in) :: data(:, :) !! Numeric data matrix in which IEEE NaNs represent missing cells.
        logical, allocatable, intent(out) :: missing(:, :) !! True at cells where `data` is NaN.
        integer :: i, j

        allocate(missing(size(data, 1), size(data, 2)))
        do j = 1, size(data, 2)
            do i = 1, size(data, 1)
                missing(i, j) = ieee_is_nan(data(i, j))
            end do
        end do
    end subroutine missing_mask

    pure subroutine md_pairs(data, result, info)
        real(dp), intent(in) :: data(:, :) !! Numeric incomplete data matrix with IEEE NaNs marking missing cells.
        type(md_pairs_result), intent(out) :: result !! Pairwise observed/missing count matrices from `md.pairs`.
        integer, intent(out) :: info !! `mice_ok` on success or `mice_invalid_argument` when fewer than two columns are supplied.
        logical, allocatable :: missing(:, :)

        call missing_mask(data, missing)
        call md_pairs_from_mask(missing, result, info)
    end subroutine md_pairs

    pure subroutine md_pairs_from_mask(missing, result, info)
        logical, intent(in) :: missing(:, :) !! Missingness mask with rows as cases and columns as variables.
        type(md_pairs_result), intent(out) :: result !! Pairwise observed/missing count matrices.
        integer, intent(out) :: info !! `mice_ok` on success or `mice_invalid_argument` for fewer than two variables.
        integer :: i, j, k, p

        p = size(missing, 2)
        if (p < 2) then
            info = mice_invalid_argument
            return
        end if
        allocate(result%rr(p, p), result%rm(p, p), result%mr(p, p), result%mm(p, p))
        result%rr = 0
        result%rm = 0
        result%mr = 0
        result%mm = 0
        do j = 1, p
            do k = 1, p
                do i = 1, size(missing, 1)
                    if (.not. missing(i, j) .and. .not. missing(i, k)) result%rr(j, k) = result%rr(j, k) + 1
                    if (.not. missing(i, j) .and. missing(i, k)) result%rm(j, k) = result%rm(j, k) + 1
                    if (missing(i, j) .and. .not. missing(i, k)) result%mr(j, k) = result%mr(j, k) + 1
                    if (missing(i, j) .and. missing(i, k)) result%mm(j, k) = result%mm(j, k) + 1
                end do
            end do
        end do
        info = mice_ok
    end subroutine md_pairs_from_mask

    pure subroutine flux(data, result, info)
        real(dp), intent(in) :: data(:, :) !! Numeric incomplete data matrix with IEEE NaNs marking missing cells.
        type(flux_result), intent(out) :: result !! Proportion observed, influx, outflux, average in/out statistics, and FICO.
        integer, intent(out) :: info !! `mice_ok` on success or a package status code on failure.

        logical, allocatable :: missing(:, :)
        type(md_pairs_result) :: pairs
        real(dp) :: denom, sum_ratio
        integer :: i, j, n, p, observed_count, incomplete_observed

        call missing_mask(data, missing)
        call md_pairs_from_mask(missing, pairs, info)
        if (info /= mice_ok) return
        n = size(data, 1)
        p = size(data, 2)
        allocate(result%pobs(p), result%influx(p), result%outflux(p), result%ainb(p), result%aout(p), result%fico(p))
        do j = 1, p
            result%pobs(j) = real(count(.not. missing(:, j)), dp) / real(n, dp)
            denom = real(sum(pairs%mr(j, :) + pairs%rr(j, :)), dp)
            if (denom > 0.0_dp) then
                result%influx(j) = real(sum(pairs%mr(j, :)), dp) / denom
            else
                result%influx(j) = 0.0_dp
            end if
            denom = real(sum(pairs%rm(j, :) + pairs%mm(j, :)), dp)
            if (denom > 0.0_dp) then
                result%outflux(j) = real(sum(pairs%rm(j, :)), dp) / denom
            else
                result%outflux(j) = 0.0_dp
            end if
            sum_ratio = 0.0_dp
            do i = 1, p
                if (i == j) cycle
                denom = real(pairs%mr(j, i) + pairs%mm(j, i), dp)
                if (denom > 0.0_dp) sum_ratio = sum_ratio + real(pairs%mr(j, i), dp) / denom
            end do
            result%ainb(j) = sum_ratio / real(max(p - 1, 1), dp)
            sum_ratio = 0.0_dp
            do i = 1, p
                if (i == j) cycle
                denom = real(pairs%rm(j, i) + pairs%rr(j, i), dp)
                if (denom > 0.0_dp) sum_ratio = sum_ratio + real(pairs%rm(j, i), dp) / denom
            end do
            result%aout(j) = sum_ratio / real(max(p - 1, 1), dp)
            observed_count = 0
            incomplete_observed = 0
            do i = 1, n
                if (missing(i, j)) cycle
                observed_count = observed_count + 1
                if (any(missing(i, :))) incomplete_observed = incomplete_observed + 1
            end do
            if (observed_count > 0) then
                result%fico(j) = real(incomplete_observed, dp) / real(observed_count, dp)
            else
                result%fico(j) = 0.0_dp
            end if
        end do
        info = mice_ok
    end subroutine flux

    pure subroutine quickpred(data, mincor, minpuc, predictor_matrix, info)
        real(dp), intent(in) :: data(:, :) !! Numeric incomplete data matrix with IEEE NaNs marking missing cells.
        real(dp), intent(in), value :: mincor !! Minimum absolute pairwise Pearson or response-indicator correlation.
        real(dp), intent(in), value :: minpuc !! Minimum proportion of usable cases for a target-predictor pair.
        integer, allocatable, intent(out) :: predictor_matrix(:, :) !! Binary target-by-predictor `quickpred` matrix.
        integer, intent(out) :: info !! `mice_ok` on success or an argument status code on failure.

        logical, allocatable :: missing(:, :)
        type(md_pairs_result) :: pairs
        real(dp) :: c_data, c_response, puc
        integer :: j, k, p

        if (mincor < 0.0_dp .or. minpuc < 0.0_dp .or. minpuc > 1.0_dp) then
            info = mice_invalid_argument
            return
        end if
        call missing_mask(data, missing)
        call md_pairs_from_mask(missing, pairs, info)
        if (info /= mice_ok) return
        p = size(data, 2)
        allocate(predictor_matrix(p, p))
        predictor_matrix = 0
        do j = 1, p
            if (count(missing(:, j)) == 0) cycle
            do k = 1, p
                if (j == k) cycle
                c_data = abs(pairwise_corr(data(:, j), data(:, k)))
                c_response = abs(response_corr(missing(:, j), data(:, k)))
                if (max(c_data, c_response) > mincor) predictor_matrix(j, k) = 1
                if (minpuc > 0.0_dp) then
                    if (pairs%mr(j, k) + pairs%mm(j, k) > 0) then
                        puc = real(pairs%mr(j, k), dp) / real(pairs%mr(j, k) + pairs%mm(j, k), dp)
                    else
                        puc = 0.0_dp
                    end if
                    if (puc < minpuc) predictor_matrix(j, k) = 0
                end if
            end do
        end do
        info = mice_ok
    end subroutine quickpred

    pure subroutine md_pattern(data, result, info)
        real(dp), intent(in) :: data(:, :) !! Numeric incomplete data matrix with IEEE NaNs marking missing cells.
        type(md_pattern_result), intent(out) :: result !! Unique patterns, frequencies, margins, and sorted column order.
        integer, intent(out) :: info !! `mice_ok` on success or `mice_invalid_argument` for fewer than two variables.

        logical, allocatable :: missing(:, :)
        integer, allocatable :: order(:), row_pattern(:), temp_pattern(:, :), frequency(:)
        integer :: i, j, k, n, n_unique, p
        logical :: found

        n = size(data, 1)
        p = size(data, 2)
        if (p < 2) then
            info = mice_invalid_argument
            return
        end if
        call missing_mask(data, missing)
        allocate(result%missing_per_column(p), order(p))
        do j = 1, p
            result%missing_per_column(j) = count(missing(:, j))
            order(j) = j
        end do
        call stable_sort_columns(result%missing_per_column, order)
        allocate(temp_pattern(max(n, 1), p), frequency(max(n, 1)), row_pattern(p))
        frequency = 0
        n_unique = 0
        do i = 1, n
            do j = 1, p
                if (missing(i, order(j))) then
                    row_pattern(j) = 0
                else
                    row_pattern(j) = 1
                end if
            end do
            found = .false.
            do k = 1, n_unique
                if (all(temp_pattern(k, :) == row_pattern)) then
                    frequency(k) = frequency(k) + 1
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                n_unique = n_unique + 1
                temp_pattern(n_unique, :) = row_pattern
                frequency(n_unique) = 1
            end if
        end do
        call sort_patterns(temp_pattern, frequency, n_unique)
        allocate(result%pattern(n_unique, p), result%frequency(n_unique), &
                 result%missing_per_pattern(n_unique), result%column_order(p))
        result%pattern = temp_pattern(1:n_unique, :)
        result%frequency = frequency(1:n_unique)
        result%column_order = order
        do i = 1, n_unique
            result%missing_per_pattern(i) = count(result%pattern(i, :) == 0)
        end do
        info = mice_ok
    end subroutine md_pattern

    pure real(dp) function pairwise_corr(x, y) result(correlation)
        real(dp), intent(in) :: x(:) !! First numeric vector with NaNs marking unavailable values.
        real(dp), intent(in) :: y(:) !! Second numeric vector with NaNs marking unavailable values.
        real(dp) :: mx, my, sxx, syy, sxy
        integer :: i, n

        mx = 0.0_dp
        my = 0.0_dp
        n = 0
        do i = 1, size(x)
            if (ieee_is_nan(x(i)) .or. ieee_is_nan(y(i))) cycle
            n = n + 1
            mx = mx + x(i)
            my = my + y(i)
        end do
        if (n < 2) then
            correlation = 0.0_dp
            return
        end if
        mx = mx / real(n, dp)
        my = my / real(n, dp)
        sxx = 0.0_dp
        syy = 0.0_dp
        sxy = 0.0_dp
        do i = 1, size(x)
            if (ieee_is_nan(x(i)) .or. ieee_is_nan(y(i))) cycle
            sxx = sxx + (x(i) - mx)**2
            syy = syy + (y(i) - my)**2
            sxy = sxy + (x(i) - mx) * (y(i) - my)
        end do
        if (sxx <= 0.0_dp .or. syy <= 0.0_dp) then
            correlation = 0.0_dp
        else
            correlation = sxy / sqrt(sxx * syy)
        end if
    end function pairwise_corr

    pure real(dp) function response_corr(missing, predictor) result(correlation)
        logical, intent(in) :: missing(:) !! Missingness indicator for the target response.
        real(dp), intent(in) :: predictor(:) !! Numeric predictor with NaNs marking unavailable values.
        real(dp) :: mr, mx, srr, sxx, srx, response
        integer :: i, n

        mr = 0.0_dp
        mx = 0.0_dp
        n = 0
        do i = 1, size(predictor)
            if (ieee_is_nan(predictor(i))) cycle
            n = n + 1
            if (missing(i)) then
                response = 0.0_dp
            else
                response = 1.0_dp
            end if
            mr = mr + response
            mx = mx + predictor(i)
        end do
        if (n < 2) then
            correlation = 0.0_dp
            return
        end if
        mr = mr / real(n, dp)
        mx = mx / real(n, dp)
        srr = 0.0_dp
        sxx = 0.0_dp
        srx = 0.0_dp
        do i = 1, size(predictor)
            if (ieee_is_nan(predictor(i))) cycle
            if (missing(i)) then
                response = 0.0_dp
            else
                response = 1.0_dp
            end if
            srr = srr + (response - mr)**2
            sxx = sxx + (predictor(i) - mx)**2
            srx = srx + (response - mr) * (predictor(i) - mx)
        end do
        if (srr <= 0.0_dp .or. sxx <= 0.0_dp) then
            correlation = 0.0_dp
        else
            correlation = srx / sqrt(srr * sxx)
        end if
    end function response_corr

    pure subroutine stable_sort_columns(missing_count, order)
        integer, intent(in) :: missing_count(:) !! Missing-cell counts for original columns.
        integer, intent(inout) :: order(:) !! Original column indices, sorted by increasing missing count.
        integer :: i, j, key

        do i = 2, size(order)
            key = order(i)
            j = i - 1
            do while (j >= 1)
                if (missing_count(order(j)) <= missing_count(key)) exit
                order(j + 1) = order(j)
                j = j - 1
            end do
            order(j + 1) = key
        end do
    end subroutine stable_sort_columns

    pure subroutine sort_patterns(pattern, frequency, n_unique)
        integer, intent(inout) :: pattern(:, :) !! Pattern rows sorted in place by increasing binary missingness code.
        integer, intent(inout) :: frequency(:) !! Pattern frequencies moved with their corresponding rows.
        integer, intent(in), value :: n_unique !! Number of initialized rows in `pattern` and `frequency`.
        integer, allocatable :: key(:)
        integer :: i, j, fkey

        allocate(key(size(pattern, 2)))
        do i = 2, n_unique
            key = pattern(i, :)
            fkey = frequency(i)
            j = i - 1
            do while (j >= 1)
                if (.not. pattern_less(pattern(j, :), key)) exit
                pattern(j + 1, :) = pattern(j, :)
                frequency(j + 1) = frequency(j)
                j = j - 1
            end do
            pattern(j + 1, :) = key
            frequency(j + 1) = fkey
        end do
    end subroutine sort_patterns

    pure logical function pattern_less(left, right) result(less)
        integer, intent(in) :: left(:) !! First binary observed/missing pattern.
        integer, intent(in) :: right(:) !! Second binary observed/missing pattern.
        integer :: j

        less = .false.
        do j = 1, size(left)
            if (left(j) == right(j)) cycle
            less = left(j) < right(j)
            return
        end do
    end function pattern_less

end module mice_diagnostics
