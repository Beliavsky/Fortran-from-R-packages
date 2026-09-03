! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from R package mice 3.19.0 by Stef van Buuren, Karin Groothuis-Oudshoorn,
! and mice contributors; see NOTICE.md and PROVENANCE.md for attribution.
! Translation of mice src/matchindex.cpp, src/match.cpp, and PMM matching logic.
module mice_matching
    use r_kinds, only : dp
    use mice_rng, only : mice_rng_state, rng_uniform, rng_integer, rng_shuffle
    use mice_status, only : mice_ok, mice_invalid_argument, mice_invalid_shape
    implicit none
    private

    public :: matchindex
    public :: matcher
    public :: pmm_match_value

contains

    subroutine matchindex(donor_metric, target_metric, donors, rng, indices, info)
        real(dp), intent(in) :: donor_metric(:) !! Predictive metric for observed donor cases.
        real(dp), intent(in) :: target_metric(:) !! Predictive metric for cases requiring imputation.
        integer, intent(in), value :: donors !! Requested donor-pool size; clipped to `[1, size(donor_metric)]`.
        type(mice_rng_state), intent(inout) :: rng !! RNG state used for tie-shuffling and donor-rank draws.
        integer, allocatable, intent(out) :: indices(:) !! One-based indices into `donor_metric`, one per target case.
        integer, intent(out) :: info !! `mice_ok` on success or a package status code on failure.

        integer, allocatable :: shuffled(:), order(:), ids(:)
        real(dp), allocatable :: sorted_values(:)
        integer :: count_found, h, i, k, l, n0, n1, r
        real(dp) :: value

        n1 = size(donor_metric)
        n0 = size(target_metric)
        if (n1 < 1) then
            info = mice_invalid_argument
            return
        end if
        k = max(1, min(donors, n1))
        allocate(shuffled(n1), order(n1), ids(n1), sorted_values(n1), indices(n0))
        shuffled = [(i, i = 1, n1)]
        call rng_shuffle(rng, shuffled)
        order = [(i, i = 1, n1)]
        call stable_order_by_values(donor_metric(shuffled), order)
        do i = 1, n1
            ids(i) = shuffled(order(i))
            sorted_values(i) = donor_metric(ids(i))
        end do

        do i = 1, n0
            value = target_metric(i)
            h = rng_integer(rng, 1, k)
            r = lower_bound(sorted_values, value)
            l = r - 1
            count_found = 0
            indices(i) = ids(max(1, min(n1, r)))
            do while (count_found < h .and. l >= 1 .and. r <= n1)
                if (value - sorted_values(l) < sorted_values(r) - value) then
                    indices(i) = ids(l)
                    l = l - 1
                else
                    indices(i) = ids(r)
                    r = r + 1
                end if
                count_found = count_found + 1
            end do
            do while (count_found < h .and. l >= 1)
                indices(i) = ids(l)
                l = l - 1
                count_found = count_found + 1
            end do
            do while (count_found < h .and. r <= n1)
                indices(i) = ids(r)
                r = r + 1
                count_found = count_found + 1
            end do
        end do
        info = mice_ok
    end subroutine matchindex

    subroutine matcher(observed_metric, missing_metric, donors, rng, indices, info)
        real(dp), intent(in) :: observed_metric(:) !! Predictive metric for observed donor cases.
        real(dp), intent(in) :: missing_metric(:) !! Predictive metric for target cases.
        integer, intent(in), value :: donors !! Requested donor-pool size, clipped to the available donors.
        type(mice_rng_state), intent(inout) :: rng !! RNG state used for tie-breaking noise and donor selection.
        integer, allocatable, intent(out) :: indices(:) !! One-based donor indices corresponding to each target.
        integer, intent(out) :: info !! `mice_ok` on success or a package status code on failure.

        real(dp), allocatable :: distance(:), work(:)
        real(dp) :: maximum, minimum, small, threshold
        integer :: count_found, goal, i, j, k, n0, n1

        n1 = size(observed_metric)
        n0 = size(missing_metric)
        if (n1 < 1) then
            info = mice_invalid_argument
            return
        end if
        k = max(1, min(donors, n1))
        minimum = minval(observed_metric)
        maximum = maxval(observed_metric)
        small = (maximum - minimum) / 65536.0_dp
        allocate(indices(n0), distance(n1), work(n1))
        do i = 1, n0
            do j = 1, n1
                distance(j) = abs(observed_metric(j) - missing_metric(i)) + small * rng_uniform(rng)
            end do
            work = distance
            call sort_real(work)
            threshold = work(k)
            goal = rng_integer(rng, 1, k)
            count_found = 0
            indices(i) = 1
            do j = 1, n1
                if (distance(j) <= threshold) count_found = count_found + 1
                if (count_found == goal) then
                    indices(i) = j
                    exit
                end if
            end do
        end do
        info = mice_ok
    end subroutine matcher

    function pmm_match_value(z, yhat, y, donors, rng, info) result(value)
        real(dp), intent(in), value :: z !! Predicted metric for one target case.
        real(dp), intent(in) :: yhat(:) !! Predicted metric for all observed donor cases.
        real(dp), intent(in) :: y(:) !! Observed response values corresponding one-to-one with `yhat`.
        integer, intent(in), value :: donors !! Requested donor-pool size.
        type(mice_rng_state), intent(inout) :: rng !! RNG state used to break ties and choose a donor.
        integer, intent(out) :: info !! `mice_ok` on success or a package status code on failure.
        real(dp) :: value

        real(dp), allocatable :: d(:), positive(:), work(:)
        real(dp) :: jitter, threshold
        integer, allocatable :: candidates(:)
        integer :: i, k, nc, np

        if (size(yhat) /= size(y) .or. size(y) < 1) then
            info = mice_invalid_shape
            value = 0.0_dp
            return
        end if
        k = max(1, min(donors, size(y)))
        allocate(d(size(y)), positive(size(y)), work(size(y)), candidates(size(y)))
        np = 0
        do i = 1, size(y)
            d(i) = abs(yhat(i) - z)
            if (d(i) > 0.0_dp) then
                np = np + 1
                positive(np) = d(i)
            end if
        end do
        if (np > 0) then
            jitter = minval(positive(1:np)) / 1.0e10_dp
        else
            jitter = 1.0e-10_dp
        end if
        do i = 1, size(y)
            d(i) = d(i) + jitter * rng_uniform(rng)
        end do
        work = d
        call sort_real(work)
        threshold = work(k)
        nc = 0
        do i = 1, size(y)
            if (d(i) <= threshold) then
                nc = nc + 1
                candidates(nc) = i
            end if
        end do
        value = y(candidates(rng_integer(rng, 1, nc)))
        info = mice_ok
    end function pmm_match_value

    pure integer function lower_bound(values, target) result(position)
        real(dp), intent(in) :: values(:) !! Nondecreasing vector to search.
        real(dp), intent(in), value :: target !! Search value.
        integer :: lo, hi, mid

        lo = 1
        hi = size(values) + 1
        do while (lo < hi)
            mid = lo + (hi - lo) / 2
            if (mid <= size(values)) then
                if (values(mid) < target) then
                    lo = mid + 1
                else
                    hi = mid
                end if
            else
                hi = mid
            end if
        end do
        position = lo
    end function lower_bound

    pure subroutine stable_order_by_values(values, order)
        real(dp), intent(in) :: values(:) !! Values whose stable ascending order is requested.
        integer, intent(inout) :: order(:) !! On input indices; on output the stable order by `values(index)`.
        integer :: i, j, key

        do i = 2, size(order)
            key = order(i)
            j = i - 1
            do while (j >= 1)
                if (values(order(j)) <= values(key)) exit
                order(j + 1) = order(j)
                j = j - 1
            end do
            order(j + 1) = key
        end do
    end subroutine stable_order_by_values

    pure subroutine sort_real(values)
        real(dp), intent(inout) :: values(:) !! Real vector sorted in ascending order in place.
        real(dp) :: key
        integer :: i, j

        do i = 2, size(values)
            key = values(i)
            j = i - 1
            do while (j >= 1)
                if (values(j) <= key) exit
                values(j + 1) = values(j)
                j = j - 1
            end do
            values(j + 1) = key
        end do
    end subroutine sort_real

end module mice_matching
