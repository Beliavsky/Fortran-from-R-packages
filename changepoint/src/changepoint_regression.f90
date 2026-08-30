! Computational translation of R package changepoint 2.3.
! Upstream license declaration: GPL (unversioned).
! See NOTICE.md and UPSTREAM.md for authorship and provenance.
module changepoint_regression
use r_kinds, only : dp
use r_linalg, only : least_squares_svd
use changepoint_types, only : changepoint_result, amoc_result
use changepoint_types, only : cp_invalid_argument, cp_linalg_failure
implicit none
private
public :: cp_regression_amoc
public :: cp_regression_pelt
public :: cp_regression_segment_fit

contains

subroutine cp_regression_segment_fit(y, x, beta, rss, status, rcond)
    real(dp), intent(in) :: y(:), x(:, :)
    real(dp), allocatable, intent(out) :: beta(:)
    real(dp), intent(out) :: rss
    integer, intent(out) :: status
    real(dp), intent(in), optional :: rcond
    integer :: p, rank, info
    real(dp), allocatable :: residual(:)

    status = 0
    p = size(x, 2)
    if (size(x, 1) /= size(y) .or. p < 1) then
        status = cp_invalid_argument
        allocate(beta(0))
        rss = huge(1.0_dp)
        return
    end if
    allocate(beta(p), residual(size(y)))
    if (present(rcond)) then
        call least_squares_svd(x, y, beta, rank, info, rcond)
    else
        call least_squares_svd(x, y, beta, rank, info)
    end if
    if (info /= 0) then
        status = cp_linalg_failure
        rss = huge(1.0_dp)
        return
    end if
    residual = y - matmul(x, beta)
    rss = dot_product(residual, residual)
end subroutine cp_regression_segment_fit

subroutine cp_regression_amoc(y, x, penalty, minseglen, shape, result, rcond)
    real(dp), intent(in) :: y(:), x(:, :)
    real(dp), intent(in) :: penalty, shape
    integer, intent(in) :: minseglen
    type(amoc_result), intent(out) :: result
    real(dp), intent(in), optional :: rcond
    integer :: n, tau, status, best_tau
    real(dp) :: null_cost, alt, best_alt

    result = amoc_result()
    n = size(y)
    if (size(x, 1) /= n .or. size(x, 2) < 1 .or. minseglen < 1 .or. n < 2 * minseglen .or. penalty < 0.0_dp) then
        result%status = cp_invalid_argument
        return
    end if
    null_cost = regression_cost(y, x, 1, n, shape, status, rcond)
    if (status /= 0) then
        result%status = status
        return
    end if
    best_alt = huge(1.0_dp)
    best_tau = n
    do tau = minseglen, n - minseglen
        alt = regression_cost(y, x, 1, tau, shape, status, rcond)
        if (status /= 0) cycle
        alt = alt + regression_cost(y, x, tau + 1, n, shape, status, rcond)
        if (status /= 0) cycle
        if (alt < best_alt) then
            best_alt = alt
            best_tau = tau
        end if
    end do
    result%cpt = best_tau
    result%null_cost = null_cost
    result%alt_cost = best_alt
    result%test_statistic = null_cost - best_alt
    result%changed = result%test_statistic >= penalty
    if (.not. result%changed) result%cpt = n
end subroutine cp_regression_amoc

subroutine cp_regression_pelt(y, x, penalty, minseglen, shape, result, rcond)
    real(dp), intent(in) :: y(:), x(:, :)
    real(dp), intent(in) :: penalty, shape
    integer, intent(in) :: minseglen
    type(changepoint_result), intent(out) :: result
    real(dp), intent(in), optional :: rcond
    real(dp), allocatable :: f(:), trial(:)
    integer, allocatable :: previous(:), nsegments(:), active(:), next_active(:), rev(:)
    integer :: n, t, i, s, nactive, nnext, best_i, new_s, last, k, status
    real(dp) :: best, val, segcost

    result = changepoint_result()
    n = size(y)
    if (size(x, 1) /= n .or. size(x, 2) < 1 .or. minseglen < 1 .or. n < minseglen .or. penalty < 0.0_dp) then
        result%status = cp_invalid_argument
        return
    end if
    allocate(f(0:n), previous(0:n), nsegments(0:n), active(0:n), next_active(0:n), trial(0:n))
    f = huge(1.0_dp)
    previous = -1
    nsegments = 0
    f(0) = -penalty
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
            segcost = regression_cost(y, x, s + 1, t, shape, status, rcond)
            if (status /= 0) cycle
            val = f(s) + segcost + penalty
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
        result%status = cp_linalg_failure
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
end subroutine cp_regression_pelt

real(dp) function regression_cost(y, x, first, last, shape, status, rcond) result(cost)
    real(dp), intent(in) :: y(:), x(:, :)
    integer, intent(in) :: first, last
    real(dp), intent(in) :: shape
    integer, intent(out) :: status
    real(dp), intent(in), optional :: rcond
    real(dp), allocatable :: beta(:)
    real(dp) :: rss, pi
    integer :: nn

    nn = last - first + 1
    if (present(rcond)) then
        call cp_regression_segment_fit(y(first:last), x(first:last, :), beta, rss, status, rcond)
    else
        call cp_regression_segment_fit(y(first:last), x(first:last, :), beta, rss, status)
    end if
    if (status /= 0) then
        cost = huge(1.0_dp)
        return
    end if
    pi = acos(-1.0_dp)
    if (shape < 0.0_dp) then
        cost = rss
    else if (abs(shape) <= epsilon(1.0_dp)) then
        rss = max(rss, tiny(1.0_dp))
        cost = real(nn, dp) + real(nn, dp) * log(2.0_dp * pi * rss) - real(nn, dp) * log(real(nn, dp))
    else
        cost = real(nn, dp) * log(2.0_dp * pi * shape) + rss / shape
    end if
end function regression_cost

end module changepoint_regression
