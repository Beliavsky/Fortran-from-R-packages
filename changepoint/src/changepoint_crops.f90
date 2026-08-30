! Computational translation of R package changepoint 2.3.
! Upstream license declaration: GPL (unversioned).
! See NOTICE.md and UPSTREAM.md for authorship and provenance.
module changepoint_crops
use r_kinds, only : dp
use changepoint_multiple, only : cp_pelt
use changepoint_types, only : changepoint_result, crops_solution, cp_invalid_argument
implicit none
private
public :: cp_crops

contains

subroutine cp_crops(data, cost_code, beta_min, beta_max, minseglen, solutions, status, shape, known_mean, mbic)
    real(dp), intent(in) :: data(:)
    integer, intent(in) :: cost_code
    real(dp), intent(in) :: beta_min, beta_max
    integer, intent(in) :: minseglen
    type(crops_solution), allocatable, intent(out) :: solutions(:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: shape, known_mean
    logical, intent(in), optional :: mbic
    real(dp), allocatable :: betas(:), costs(:)
    integer, allocatable :: ks(:)
    type(changepoint_result), allocatable :: results(:)
    type(changepoint_result) :: r
    integer :: cap, nsol, i, j
    real(dp) :: beta_new, tol
    logical :: added, use_mbic

    status = 0
    if (beta_min < 0.0_dp .or. beta_max <= beta_min .or. minseglen < 1) then
        status = cp_invalid_argument
        allocate(solutions(0))
        return
    end if
    use_mbic = .false.
    if (present(mbic)) use_mbic = mbic
    cap = max(4, size(data) / max(1, minseglen) + 4)
    allocate(betas(cap), costs(cap), ks(cap), results(cap))
    nsol = 0
    call run_pelt_optional(data, cost_code, beta_min, minseglen, r, shape, known_mean, use_mbic)
    if (r%status /= 0) then
        status = r%status
        allocate(solutions(0))
        return
    end if
    call append_solution(beta_min, r)
    call run_pelt_optional(data, cost_code, beta_max, minseglen, r, shape, known_mean, use_mbic)
    if (r%status /= 0) then
        status = r%status
        allocate(solutions(0))
        return
    end if
    call append_solution(beta_max, r)

    tol = 100.0_dp * epsilon(1.0_dp)
    do
        call sort_by_beta()
        added = .false.
        do i = 1, nsol - 1
            if (abs(ks(i) - ks(i + 1)) <= 1) cycle
            beta_new = (costs(i + 1) - costs(i)) / real(ks(i) - ks(i + 1), dp)
            if (beta_new <= betas(i) + tol .or. beta_new >= betas(i + 1) - tol) cycle
            if (any(abs(betas(1:nsol) - beta_new) <= tol * max(1.0_dp, abs(beta_new)))) cycle
            call run_pelt_optional(data, cost_code, beta_new, minseglen, r, shape, known_mean, use_mbic)
            if (r%status /= 0) then
                status = r%status
                allocate(solutions(0))
                return
            end if
            call append_solution(beta_new, r)
            added = .true.
            exit
        end do
        if (.not. added) exit
        if (nsol >= cap) exit
    end do
    call sort_by_beta()

    do i = nsol, 2, -1
        if (ks(i) == ks(i - 1)) then
            do j = i, nsol - 1
                betas(j) = betas(j + 1)
                costs(j) = costs(j + 1)
                ks(j) = ks(j + 1)
                results(j) = results(j + 1)
            end do
            nsol = nsol - 1
        end if
    end do

    allocate(solutions(nsol))
    do i = 1, nsol
        solutions(i)%ncpts = ks(i)
        solutions(i)%unpenalized_cost = costs(i)
        solutions(i)%cpts = results(i)%cpts
        if (i == 1) then
            solutions(i)%beta_start = beta_min
        else
            solutions(i)%beta_start = solutions(i - 1)%beta_end
        end if
        if (i == nsol) then
            solutions(i)%beta_end = beta_max
        else
            if (ks(i) == ks(i + 1)) then
                solutions(i)%beta_end = betas(i + 1)
            else
                solutions(i)%beta_end = (costs(i + 1) - costs(i)) / real(ks(i) - ks(i + 1), dp)
            end if
        end if
    end do

contains

subroutine append_solution(beta, result)
    real(dp), intent(in) :: beta
    type(changepoint_result), intent(in) :: result
    if (nsol >= cap) return
    nsol = nsol + 1
    betas(nsol) = beta
    costs(nsol) = result%unpenalized_cost
    ks(nsol) = result%ncpts
    results(nsol) = result
end subroutine append_solution

subroutine sort_by_beta()
    real(dp) :: tb, tc
    integer :: tk
    type(changepoint_result) :: tr
    do i = 2, nsol
        tb = betas(i)
        tc = costs(i)
        tk = ks(i)
        tr = results(i)
        j = i - 1
        do while (j >= 1)
            if (betas(j) <= tb) exit
            betas(j + 1) = betas(j)
            costs(j + 1) = costs(j)
            ks(j + 1) = ks(j)
            results(j + 1) = results(j)
            j = j - 1
        end do
        betas(j + 1) = tb
        costs(j + 1) = tc
        ks(j + 1) = tk
        results(j + 1) = tr
    end do
end subroutine sort_by_beta

end subroutine cp_crops

subroutine run_pelt_optional(data, cost_code, penalty, minseglen, result, shape, known_mean, mbic)
    real(dp), intent(in) :: data(:)
    integer, intent(in) :: cost_code, minseglen
    real(dp), intent(in) :: penalty
    type(changepoint_result), intent(out) :: result
    real(dp), intent(in), optional :: shape, known_mean
    logical, intent(in) :: mbic

    if (present(shape) .and. present(known_mean)) then
        call cp_pelt(data, cost_code, penalty, minseglen, result, shape=shape, known_mean=known_mean, mbic=mbic)
    else if (present(shape)) then
        call cp_pelt(data, cost_code, penalty, minseglen, result, shape=shape, mbic=mbic)
    else if (present(known_mean)) then
        call cp_pelt(data, cost_code, penalty, minseglen, result, known_mean=known_mean, mbic=mbic)
    else
        call cp_pelt(data, cost_code, penalty, minseglen, result, mbic=mbic)
    end if
end subroutine run_pelt_optional

end module changepoint_crops
