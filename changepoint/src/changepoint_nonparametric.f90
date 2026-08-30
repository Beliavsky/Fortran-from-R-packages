! Computational translation of R package changepoint 2.3.
! Upstream license declaration: GPL (unversioned).
! See NOTICE.md and UPSTREAM.md for authorship and provenance.
module changepoint_nonparametric
use r_kinds, only : dp
use changepoint_types, only : binseg_result, segneigh_result, cp_invalid_argument, cp_invalid_data
implicit none
private
public :: cp_binseg_cusum
public :: cp_binseg_css
public :: cp_segneigh_cusum
public :: cp_segneigh_css

contains

subroutine cp_binseg_cusum(data, penalty, minseglen, max_cpts, result)
    real(dp), intent(in) :: data(:)
    real(dp), intent(in) :: penalty
    integer, intent(in) :: minseglen, max_cpts
    type(binseg_result), intent(out) :: result
    real(dp), allocatable :: sx(:)
    integer :: n, q, split, b, left, right, best_split, nbounds, i
    integer, allocatable :: bounds(:)
    real(dp) :: stat, best_stat, old_stat, total, left_sum

    result = binseg_result()
    n = size(data)
    if (n < 2 * minseglen .or. minseglen < 1 .or. max_cpts < 1 .or. penalty < 0.0_dp) then
        result%status = cp_invalid_argument
        return
    end if
    allocate(sx(0:n), bounds(0:max_cpts + 1), result%candidates(max_cpts), result%scores(max_cpts))
    sx = 0.0_dp
    do i = 1, n
        sx(i) = sx(i - 1) + data(i)
    end do
    bounds = 0
    bounds(0) = 0
    bounds(1) = n
    nbounds = 2
    result%candidates = 0
    result%scores = -huge(1.0_dp)
    old_stat = huge(1.0_dp)
    do q = 1, max_cpts
        best_stat = -1.0_dp
        best_split = 0
        do b = 0, nbounds - 2
            left = bounds(b)
            right = bounds(b + 1)
            if (right - left < 2 * minseglen) cycle
            total = sx(right) - sx(left)
            do split = left + minseglen, right - minseglen
                left_sum = sx(split) - sx(left)
                stat = abs(left_sum - real(split - left, dp) / real(right - left, dp) * total) / real(right - left, dp)
                if (stat > best_stat) then
                    best_stat = stat
                    best_split = split
                end if
            end do
        end do
        if (best_split == 0) exit
        old_stat = min(old_stat, best_stat)
        result%candidates(q) = best_split
        result%scores(q) = old_stat
        call insert_bound(bounds, nbounds, best_split)
    end do
    call select_binseg(result, penalty)
end subroutine cp_binseg_cusum

subroutine cp_binseg_css(data, penalty, minseglen, max_cpts, result)
    real(dp), intent(in) :: data(:)
    real(dp), intent(in) :: penalty
    integer, intent(in) :: minseglen, max_cpts
    type(binseg_result), intent(out) :: result
    real(dp), allocatable :: sx2(:)
    integer :: n, q, split, b, left, right, best_split, nbounds, i
    integer, allocatable :: bounds(:)
    real(dp) :: stat, best_stat, old_stat, total, left_ss

    result = binseg_result()
    n = size(data)
    if (n < 2 * minseglen .or. minseglen < 1 .or. max_cpts < 1 .or. penalty < 0.0_dp) then
        result%status = cp_invalid_argument
        return
    end if
    allocate(sx2(0:n), bounds(0:max_cpts + 1), result%candidates(max_cpts), result%scores(max_cpts))
    sx2 = 0.0_dp
    do i = 1, n
        sx2(i) = sx2(i - 1) + data(i) * data(i)
    end do
    bounds = 0
    bounds(0) = 0
    bounds(1) = n
    nbounds = 2
    result%candidates = 0
    result%scores = -huge(1.0_dp)
    old_stat = huge(1.0_dp)
    do q = 1, max_cpts
        best_stat = -1.0_dp
        best_split = 0
        do b = 0, nbounds - 2
            left = bounds(b)
            right = bounds(b + 1)
            if (right - left < 2 * minseglen) cycle
            total = sx2(right) - sx2(left)
            if (total <= 0.0_dp) cycle
            do split = left + minseglen, right - minseglen
                left_ss = sx2(split) - sx2(left)
                stat = sqrt(real(right - left, dp) / 2.0_dp) * &
                    abs(left_ss / total - real(split - left, dp) / real(right - left, dp))
                if (stat > best_stat) then
                    best_stat = stat
                    best_split = split
                end if
            end do
        end do
        if (best_split == 0) exit
        old_stat = min(old_stat, best_stat)
        result%candidates(q) = best_split
        result%scores(q) = old_stat
        call insert_bound(bounds, nbounds, best_split)
    end do
    call select_binseg(result, penalty)
end subroutine cp_binseg_css

subroutine cp_segneigh_cusum(data, penalty, max_cpts, result)
    real(dp), intent(in) :: data(:)
    real(dp), intent(in) :: penalty
    integer, intent(in) :: max_cpts
    type(segneigh_result), intent(out) :: result
    real(dp), allocatable :: sx(:), score(:, :)
    integer, allocatable :: prev(:, :), rev(:), tmpcpts(:)
    integer :: n, max_segments, q, j, v, best_v, i, op_cps
    real(dp) :: cand, best, local_stat
    logical :: accepted

    result = segneigh_result()
    n = size(data)
    if (n < 2 .or. max_cpts < 0 .or. penalty < 0.0_dp) then
        result%status = cp_invalid_argument
        return
    end if
    max_segments = min(max_cpts + 1, n)
    allocate(sx(0:n), score(0:max_segments, 0:n), prev(0:max_segments, 0:n))
    allocate(result%cost_by_ncpts(0:max_segments - 1))
    sx = 0.0_dp
    do i = 1, n
        sx(i) = sx(i - 1) + data(i)
    end do
    score = -huge(1.0_dp)
    prev = -1
    result%cost_by_ncpts = -huge(1.0_dp)
    do q = 2, max_segments
        do j = q, n
            best = -huge(1.0_dp)
            best_v = -1
            do v = q - 1, j - 1
                if (q == 2) then
                    cand = abs((sx(v) - real(v, dp) / real(j, dp) * sx(j)) / real(j, dp))
                else
                    if (prev(q - 1, v) < 0 .or. score(q - 1, v) <= -huge(1.0_dp) / 4.0_dp) cycle
                    local_stat = abs(((sx(v) - sx(prev(q - 1, v))) - &
                        real(v - prev(q - 1, v), dp) / real(j - prev(q - 1, v), dp) * &
                        (sx(j) - sx(prev(q - 1, v)))) / real(j - prev(q - 1, v), dp))
                    cand = score(q - 1, v) + local_stat
                end if
                if (cand > best) then
                    best = cand
                    best_v = v
                end if
            end do
            score(q, j) = best
            prev(q, j) = best_v
        end do
        result%cost_by_ncpts(q - 1) = score(q, n)
    end do
    result%cost_by_ncpts(0) = 0.0_dp

    op_cps = 0
    do q = 2, max_segments
        call recover_score_cpts(prev, q, n, tmpcpts)
        accepted = .true.
        call cusum_local_check(sx, n, tmpcpts, penalty, accepted)
        if (.not. accepted) exit
        op_cps = q - 1
    end do
    result%ncpts = op_cps
    allocate(result%cpts(op_cps))
    if (op_cps > 0) then
        call recover_score_cpts(prev, op_cps + 1, n, rev)
        result%cpts = rev
    end if
end subroutine cp_segneigh_cusum

subroutine cp_segneigh_css(data, penalty, max_cpts, result)
    real(dp), intent(in) :: data(:)
    real(dp), intent(in) :: penalty
    integer, intent(in) :: max_cpts
    type(segneigh_result), intent(out) :: result
    real(dp), allocatable :: sx2(:), score(:, :)
    integer, allocatable :: prev(:, :), tmpcpts(:), rev(:)
    integer :: n, max_segments, q, j, v, best_v, i, op_cps
    real(dp) :: cand, best, local_stat, denom
    logical :: accepted

    result = segneigh_result()
    n = size(data)
    if (n < 4 .or. max_cpts < 0 .or. penalty < 0.0_dp) then
        result%status = cp_invalid_argument
        return
    end if
    max_segments = min(max_cpts + 1, n)
    allocate(sx2(0:n), score(0:max_segments, 0:n), prev(0:max_segments, 0:n))
    allocate(result%cost_by_ncpts(0:max_segments - 1))
    sx2 = 0.0_dp
    do i = 1, n
        sx2(i) = sx2(i - 1) + data(i) * data(i)
    end do
    if (sx2(n) <= 0.0_dp) then
        result%status = cp_invalid_data
        return
    end if
    score = -huge(1.0_dp)
    prev = -1
    result%cost_by_ncpts = -huge(1.0_dp)
    do q = 2, max_segments
        do j = q, n
            best = -huge(1.0_dp)
            best_v = -1
            do v = q - 1, j - 1
                if (q == 2) then
                    cand = abs(sqrt(real(j, dp) / 2.0_dp) * (sx2(v) / sx2(j) - real(v, dp) / real(j, dp)))
                else
                    if (prev(q - 1, v) < 0 .or. score(q - 1, v) <= -huge(1.0_dp) / 4.0_dp) cycle
                    denom = sx2(j) - sx2(prev(q - 1, v))
                    if (denom <= 0.0_dp) cycle
                    local_stat = abs(sqrt(real(j - prev(q - 1, v), dp) / 2.0_dp) * &
                        ((sx2(v) - sx2(prev(q - 1, v))) / denom - &
                        real(v - prev(q - 1, v), dp) / real(j - prev(q - 1, v), dp)))
                    cand = score(q - 1, v) + local_stat
                end if
                if (cand > best) then
                    best = cand
                    best_v = v
                end if
            end do
            score(q, j) = best
            prev(q, j) = best_v
        end do
        result%cost_by_ncpts(q - 1) = score(q, n)
    end do
    result%cost_by_ncpts(0) = 0.0_dp

    op_cps = 0
    do q = 2, max_segments
        call recover_score_cpts(prev, q, n, tmpcpts)
        accepted = .true.
        call css_local_check(sx2, n, tmpcpts, penalty, accepted)
        if (.not. accepted) exit
        op_cps = q - 1
    end do
    result%ncpts = op_cps
    allocate(result%cpts(op_cps))
    if (op_cps > 0) then
        call recover_score_cpts(prev, op_cps + 1, n, rev)
        result%cpts = rev
    end if
end subroutine cp_segneigh_css

subroutine select_binseg(result, penalty)
    type(binseg_result), intent(inout) :: result
    real(dp), intent(in) :: penalty
    integer :: q

    result%ncpts = 0
    do q = 1, size(result%candidates)
        if (result%candidates(q) == 0) exit
        if (result%scores(q) >= penalty) then
            result%ncpts = q
        else
            exit
        end if
    end do
    allocate(result%cpts(result%ncpts))
    if (result%ncpts > 0) then
        result%cpts = result%candidates(1:result%ncpts)
        call sort_int(result%cpts)
    end if
end subroutine select_binseg

subroutine insert_bound(bounds, nbounds, value)
    integer, intent(inout) :: bounds(0:)
    integer, intent(inout) :: nbounds
    integer, intent(in) :: value
    integer :: pos, i

    pos = 0
    do while (pos < nbounds .and. bounds(pos) < value)
        pos = pos + 1
    end do
    do i = nbounds - 1, pos, -1
        bounds(i + 1) = bounds(i)
    end do
    bounds(pos) = value
    nbounds = nbounds + 1
end subroutine insert_bound

subroutine recover_score_cpts(prev, segments, n, cpts)
    integer, intent(in) :: prev(0:, 0:)
    integer, intent(in) :: segments, n
    integer, allocatable, intent(out) :: cpts(:)
    integer :: m, last, k

    allocate(cpts(max(0, segments - 1)))
    last = n
    do m = segments, 2, -1
        k = m - 1
        cpts(k) = prev(m, last)
        last = prev(m, last)
    end do
    if (size(cpts) > 1) call sort_int(cpts)
end subroutine recover_score_cpts

subroutine cusum_local_check(sx, n, cpts, penalty, accepted)
    real(dp), intent(in) :: sx(0:)
    integer, intent(in) :: n, cpts(:)
    real(dp), intent(in) :: penalty
    logical, intent(out) :: accepted
    integer, allocatable :: b(:)
    integer :: i
    real(dp) :: stat

    allocate(b(0:size(cpts) + 1))
    b(0) = 0
    if (size(cpts) > 0) b(1:size(cpts)) = cpts
    b(size(cpts) + 1) = n
    accepted = .true.
    do i = 1, size(cpts)
        stat = abs(((sx(b(i)) - sx(b(i - 1))) - &
            real(b(i) - b(i - 1), dp) / real(b(i + 1) - b(i - 1), dp) * &
            (sx(b(i + 1)) - sx(b(i - 1)))) / real(b(i + 1) - b(i - 1), dp))
        if (stat < penalty) then
            accepted = .false.
            return
        end if
    end do
end subroutine cusum_local_check

subroutine css_local_check(sx2, n, cpts, penalty, accepted)
    real(dp), intent(in) :: sx2(0:)
    integer, intent(in) :: n, cpts(:)
    real(dp), intent(in) :: penalty
    logical, intent(out) :: accepted
    integer, allocatable :: b(:)
    integer :: i
    real(dp) :: stat, denom

    allocate(b(0:size(cpts) + 1))
    b(0) = 0
    if (size(cpts) > 0) b(1:size(cpts)) = cpts
    b(size(cpts) + 1) = n
    accepted = .true.
    do i = 1, size(cpts)
        denom = sx2(b(i + 1)) - sx2(b(i - 1))
        if (denom <= 0.0_dp) then
            accepted = .false.
            return
        end if
        stat = abs(sqrt(real(b(i + 1) - b(i - 1), dp) / 2.0_dp) * &
            ((sx2(b(i)) - sx2(b(i - 1))) / denom - &
            real(b(i) - b(i - 1), dp) / real(b(i + 1) - b(i - 1), dp)))
        if (stat < penalty) then
            accepted = .false.
            return
        end if
    end do
end subroutine css_local_check

subroutine sort_int(x)
    integer, intent(inout) :: x(:)
    integer :: i, j, key
    do i = 2, size(x)
        key = x(i)
        j = i - 1
        do while (j >= 1)
            if (x(j) <= key) exit
            x(j + 1) = x(j)
            j = j - 1
        end do
        x(j + 1) = key
    end do
end subroutine sort_int

end module changepoint_nonparametric
