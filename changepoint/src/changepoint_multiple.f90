! Computational translation of R package changepoint 2.3.
! Upstream license declaration: GPL (unversioned).
! See NOTICE.md and UPSTREAM.md for authorship and provenance.
module changepoint_multiple
use r_kinds, only : dp
use changepoint_costs, only : cp_cost_model, cp_build_cost_model, cp_segment_cost
use changepoint_types, only : changepoint_result, binseg_result, segneigh_result
use changepoint_types, only : cp_invalid_argument, cp_invalid_data
implicit none
private
public :: cp_pelt
public :: cp_binseg
public :: cp_segneigh

contains

subroutine cp_pelt(data, cost_code, penalty, minseglen, result, shape, known_mean, mbic)
    real(dp), intent(in) :: data(:)
    integer, intent(in) :: cost_code
    real(dp), intent(in) :: penalty
    integer, intent(in) :: minseglen
    type(changepoint_result), intent(out) :: result
    real(dp), intent(in), optional :: shape, known_mean
    logical, intent(in), optional :: mbic
    type(cp_cost_model) :: model
    real(dp), allocatable :: f(:), trial(:)
    integer, allocatable :: previous(:), nsegments(:), active(:), next_active(:), rev(:)
    integer :: n, t, i, s, nactive, nnext, best_i, new_s, last, k, status
    real(dp) :: best, val
    logical :: use_mbic

    result = changepoint_result()
    n = size(data)
    if (minseglen < 1 .or. n < minseglen .or. penalty < 0.0_dp) then
        result%status = cp_invalid_argument
        return
    end if
    use_mbic = .false.
    if (present(mbic)) use_mbic = mbic
    call build_model_optional(data, cost_code, model, status, shape, known_mean, use_mbic)
    if (status /= 0) then
        result%status = cp_invalid_data
        return
    end if
    allocate(f(0:n), previous(0:n), nsegments(0:n), active(0:n), next_active(0:n), trial(0:n))
    f = huge(1.0_dp)
    previous = -1
    nsegments = 0
    f(0) = -penalty
    previous(0) = 0
    active = 0
    active(0) = 0
    nactive = 1

    do t = minseglen, n
        best = huge(1.0_dp)
        best_i = 0
        trial = huge(1.0_dp)
        do i = 0, nactive - 1
            s = active(i)
            if (t - s < minseglen) cycle
            if (f(s) >= huge(1.0_dp) / 4.0_dp) cycle
            val = f(s) + cp_segment_cost(model, s, t) + penalty
            trial(i) = val
            if (val < best) then
                best = val
                best_i = i
            end if
        end do
        if (best >= huge(1.0_dp) / 4.0_dp) cycle
        f(t) = best
        previous(t) = active(best_i)
        nsegments(t) = nsegments(previous(t)) + 1

        nnext = 0
        do i = 0, nactive - 1
            if (trial(i) <= f(t) + penalty) then
                next_active(nnext) = active(i)
                nnext = nnext + 1
            end if
        end do
        new_s = t - minseglen + 1
        if (new_s >= 0 .and. new_s <= n) then
            if (f(new_s) < huge(1.0_dp) / 4.0_dp) then
                if (nnext == 0) then
                    next_active(nnext) = new_s
                    nnext = nnext + 1
                else if (.not. any(next_active(0:nnext - 1) == new_s)) then
                    next_active(nnext) = new_s
                    nnext = nnext + 1
                end if
            end if
        end if
        if (nnext > 0) active(0:nnext - 1) = next_active(0:nnext - 1)
        nactive = nnext
    end do

    if (f(n) >= huge(1.0_dp) / 4.0_dp) then
        result%status = cp_invalid_argument
        return
    end if
    result%ncpts = max(0, nsegments(n) - 1)
    result%objective = f(n)
    result%unpenalized_cost = f(n) - real(result%ncpts, dp) * penalty
    allocate(result%cpts(result%ncpts))
    if (result%ncpts > 0) then
        allocate(rev(result%ncpts))
        last = n
        k = 0
        do while (previous(last) > 0)
            k = k + 1
            rev(k) = previous(last)
            last = previous(last)
        end do
        do i = 1, result%ncpts
            result%cpts(i) = rev(result%ncpts - i + 1)
        end do
    end if
end subroutine cp_pelt

subroutine cp_binseg(data, cost_code, penalty, minseglen, max_cpts, result, shape, known_mean, mbic)
    real(dp), intent(in) :: data(:)
    integer, intent(in) :: cost_code
    real(dp), intent(in) :: penalty
    integer, intent(in) :: minseglen, max_cpts
    type(binseg_result), intent(out) :: result
    real(dp), intent(in), optional :: shape, known_mean
    logical, intent(in), optional :: mbic
    type(cp_cost_model) :: model
    integer, allocatable :: bounds(:)
    integer :: n, q, nbounds, b, left, right, split, best_split, pos, status, i
    real(dp) :: whole, score, best_score, old_score
    logical :: use_mbic

    result = binseg_result()
    n = size(data)
    if (minseglen < 1 .or. max_cpts < 1 .or. n < 2 * minseglen .or. penalty < 0.0_dp) then
        result%status = cp_invalid_argument
        return
    end if
    use_mbic = .false.
    if (present(mbic)) use_mbic = mbic
    call build_model_optional(data, cost_code, model, status, shape, known_mean, use_mbic)
    if (status /= 0) then
        result%status = cp_invalid_data
        return
    end if
    allocate(result%candidates(max_cpts), result%scores(max_cpts), bounds(0:max_cpts + 1))
    result%candidates = 0
    result%scores = -huge(1.0_dp)
    bounds = 0
    bounds(0) = 0
    bounds(1) = n
    nbounds = 2
    old_score = huge(1.0_dp)

    do q = 1, max_cpts
        best_score = -huge(1.0_dp)
        best_split = 0
        do b = 0, nbounds - 2
            left = bounds(b)
            right = bounds(b + 1)
            if (right - left < 2 * minseglen) cycle
            whole = cp_segment_cost(model, left, right)
            do split = left + minseglen, right - minseglen
                score = 0.5_dp * (whole - cp_segment_cost(model, left, split) - cp_segment_cost(model, split, right))
                if (score > best_score) then
                    best_score = score
                    best_split = split
                end if
            end do
        end do
        if (best_split == 0) exit
        old_score = min(old_score, best_score)
        result%candidates(q) = best_split
        result%scores(q) = old_score
        pos = 0
        do while (pos < nbounds .and. bounds(pos) < best_split)
            pos = pos + 1
        end do
        do i = nbounds - 1, pos, -1
            bounds(i + 1) = bounds(i)
        end do
        bounds(pos) = best_split
        nbounds = nbounds + 1
    end do

    result%ncpts = 0
    do q = 1, max_cpts
        if (result%candidates(q) == 0) exit
        if (2.0_dp * result%scores(q) >= penalty) then
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
end subroutine cp_binseg

subroutine cp_segneigh(data, cost_code, penalty, minseglen, max_cpts, result, shape, known_mean, mbic)
    real(dp), intent(in) :: data(:)
    integer, intent(in) :: cost_code
    real(dp), intent(in) :: penalty
    integer, intent(in) :: minseglen, max_cpts
    type(segneigh_result), intent(out) :: result
    real(dp), intent(in), optional :: shape, known_mean
    logical, intent(in), optional :: mbic
    type(cp_cost_model) :: model
    real(dp), allocatable :: dp_cost(:, :)
    integer, allocatable :: prev(:, :), rev(:)
    integer :: n, max_segments, m, t, s, best_s, best_m, status, k, last, i
    real(dp) :: best, value, criterion, best_criterion
    logical :: use_mbic

    result = segneigh_result()
    n = size(data)
    max_segments = max_cpts + 1
    if (minseglen < 1 .or. max_cpts < 0 .or. n < minseglen .or. penalty < 0.0_dp) then
        result%status = cp_invalid_argument
        return
    end if
    use_mbic = .false.
    if (present(mbic)) use_mbic = mbic
    call build_model_optional(data, cost_code, model, status, shape, known_mean, use_mbic)
    if (status /= 0) then
        result%status = cp_invalid_data
        return
    end if
    max_segments = min(max_segments, n / minseglen)
    allocate(dp_cost(0:max_segments, 0:n), prev(0:max_segments, 0:n))
    allocate(result%cost_by_ncpts(0:max_segments - 1))
    dp_cost = huge(1.0_dp)
    prev = -1
    dp_cost(0, 0) = 0.0_dp
    result%cost_by_ncpts = huge(1.0_dp)

    do m = 1, max_segments
        do t = m * minseglen, n
            best = huge(1.0_dp)
            best_s = -1
            do s = (m - 1) * minseglen, t - minseglen
                if (dp_cost(m - 1, s) >= huge(1.0_dp) / 4.0_dp) cycle
                value = dp_cost(m - 1, s) + cp_segment_cost(model, s, t)
                if (value < best) then
                    best = value
                    best_s = s
                end if
            end do
            dp_cost(m, t) = best
            prev(m, t) = best_s
        end do
        result%cost_by_ncpts(m - 1) = dp_cost(m, n)
    end do

    best_criterion = huge(1.0_dp)
    best_m = 1
    do m = 1, max_segments
        criterion = dp_cost(m, n) + real(m - 1, dp) * penalty
        if (criterion < best_criterion) then
            best_criterion = criterion
            best_m = m
        end if
    end do
    result%ncpts = best_m - 1
    allocate(result%cpts(result%ncpts))
    if (result%ncpts > 0) then
        allocate(rev(result%ncpts))
        last = n
        k = 0
        do m = best_m, 2, -1
            k = k + 1
            rev(k) = prev(m, last)
            last = prev(m, last)
        end do
        do i = 1, result%ncpts
            result%cpts(i) = rev(result%ncpts - i + 1)
        end do
    end if
end subroutine cp_segneigh

subroutine build_model_optional(data, code, model, status, shape, known_mean, mbic)
    real(dp), intent(in) :: data(:)
    integer, intent(in) :: code
    type(cp_cost_model), intent(out) :: model
    integer, intent(out) :: status
    real(dp), intent(in), optional :: shape, known_mean
    logical, intent(in) :: mbic

    if (present(shape) .and. present(known_mean)) then
        call cp_build_cost_model(data, code, model, shape=shape, known_mean=known_mean, mbic=mbic, status=status)
    else if (present(shape)) then
        call cp_build_cost_model(data, code, model, shape=shape, mbic=mbic, status=status)
    else if (present(known_mean)) then
        call cp_build_cost_model(data, code, model, known_mean=known_mean, mbic=mbic, status=status)
    else
        call cp_build_cost_model(data, code, model, mbic=mbic, status=status)
    end if
end subroutine build_model_optional

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

end module changepoint_multiple
