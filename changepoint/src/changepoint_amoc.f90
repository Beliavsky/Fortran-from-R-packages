! Computational translation of R package changepoint 2.3.
! Upstream license declaration: GPL (unversioned).
! See NOTICE.md and UPSTREAM.md for authorship and provenance.
module changepoint_amoc
use r_kinds, only : dp
use changepoint_costs, only : cp_cost_model, cp_build_cost_model, cp_segment_cost
use changepoint_types, only : amoc_result, cp_invalid_argument, cp_invalid_data
implicit none
private
public :: cp_amoc
public :: cp_amoc_cusum
public :: cp_amoc_css

contains

subroutine cp_amoc(data, cost_code, penalty, minseglen, result, shape, known_mean, mbic)
    real(dp), intent(in) :: data(:)
    integer, intent(in) :: cost_code
    real(dp), intent(in) :: penalty
    integer, intent(in) :: minseglen
    type(amoc_result), intent(out) :: result
    real(dp), intent(in), optional :: shape, known_mean
    logical, intent(in), optional :: mbic
    type(cp_cost_model) :: model
    integer :: n, tau, best_tau, status
    real(dp) :: null_cost, alt, best_alt
    logical :: use_mbic

    result = amoc_result()
    n = size(data)
    if (minseglen < 1 .or. n < 2 * minseglen .or. penalty < 0.0_dp) then
        result%status = cp_invalid_argument
        return
    end if
    use_mbic = .false.
    if (present(mbic)) use_mbic = mbic
    call build_model_optional(data, cost_code, model, status, shape, known_mean, .false.)
    if (status /= 0) then
        result%status = cp_invalid_data
        return
    end if
    null_cost = cp_segment_cost(model, 0, n)
    best_alt = huge(1.0_dp)
    best_tau = n
    do tau = minseglen, n - minseglen
        alt = cp_segment_cost(model, 0, tau) + cp_segment_cost(model, tau, n)
        if (use_mbic) alt = alt + log(real(tau, dp)) + log(real(n - tau + 1, dp))
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
end subroutine cp_amoc

subroutine cp_amoc_css(data, penalty, minseglen, result)
    real(dp), intent(in) :: data(:)
    real(dp), intent(in) :: penalty
    integer, intent(in) :: minseglen
    type(amoc_result), intent(out) :: result
    integer :: n, tau, best_tau
    real(dp), allocatable :: sx2(:)
    real(dp) :: d, best_d, total

    result = amoc_result()
    n = size(data)
    if (minseglen < 1 .or. n < 2 * minseglen .or. penalty < 0.0_dp) then
        result%status = cp_invalid_argument
        return
    end if
    allocate(sx2(0:n))
    sx2 = 0.0_dp
    do tau = 1, n
        sx2(tau) = sx2(tau - 1) + data(tau) * data(tau)
    end do
    total = sx2(n)
    if (total <= 0.0_dp) then
        result%status = cp_invalid_data
        return
    end if
    best_d = -1.0_dp
    best_tau = n
    do tau = minseglen, n - minseglen
        d = abs(sx2(tau) / total - real(tau, dp) / real(n, dp))
        if (d > best_d) then
            best_d = d
            best_tau = tau
        end if
    end do
    result%cpt = best_tau
    result%test_statistic = sqrt(real(n, dp) / 2.0_dp) * best_d
    result%null_cost = result%test_statistic
    result%changed = result%test_statistic >= penalty
    if (.not. result%changed) result%cpt = n
end subroutine cp_amoc_css

subroutine cp_amoc_cusum(data, penalty, minseglen, result)
    real(dp), intent(in) :: data(:)
    real(dp), intent(in) :: penalty
    integer, intent(in) :: minseglen
    type(amoc_result), intent(out) :: result
    integer :: n, i, best_tau
    real(dp) :: mu, acc, value, best_value

    result = amoc_result()
    n = size(data)
    if (minseglen < 1 .or. n < 2 * minseglen .or. penalty < 0.0_dp) then
        result%status = cp_invalid_argument
        return
    end if
    mu = sum(data) / real(n, dp)
    acc = 0.0_dp
    best_value = -1.0_dp
    best_tau = n
    do i = 1, n
        acc = acc + data(i) - mu
        if (i >= minseglen .and. i <= n - minseglen) then
            value = abs(acc / real(n, dp))
            if (value > best_value) then
                best_value = value
                best_tau = i
            end if
        end if
    end do
    result%cpt = best_tau
    result%test_statistic = best_value
    result%null_cost = best_value
    result%changed = best_value >= penalty
    if (.not. result%changed) result%cpt = n
end subroutine cp_amoc_cusum

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

end module changepoint_amoc
