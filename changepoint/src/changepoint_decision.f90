! Computational translation of R package changepoint 2.3.
! Upstream license declaration: GPL (unversioned).
! See NOTICE.md and UPSTREAM.md for authorship and provenance.
module changepoint_decision
use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
use r_kinds, only : dp
use changepoint_types, only : cp_ok, cp_invalid_argument
implicit none
private
public :: cp_decision

contains

subroutine cp_decision(tau, null_cost, penalty, n, cpts, status, alt_cost)
    integer, intent(in) :: tau(:)
    real(dp), intent(in) :: null_cost(:)
    real(dp), intent(in) :: penalty
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: cpts(:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: alt_cost(:)
    integer :: i
    real(dp) :: statistic

    status = cp_ok
    if (size(tau) /= size(null_cost) .or. penalty < 0.0_dp .or. n < 0) then
        status = cp_invalid_argument
        allocate(cpts(0))
        return
    end if
    if (present(alt_cost)) then
        if (size(alt_cost) /= size(tau)) then
            status = cp_invalid_argument
            allocate(cpts(0))
            return
        end if
    end if

    allocate(cpts(size(tau)))
    do i = 1, size(tau)
        if (present(alt_cost)) then
            if (ieee_is_finite(alt_cost(i))) then
                statistic = null_cost(i) - alt_cost(i)
            else
                statistic = null_cost(i)
            end if
        else
            statistic = null_cost(i)
        end if
        if (statistic >= penalty) then
            cpts(i) = tau(i)
        else
            cpts(i) = n
        end if
    end do
end subroutine cp_decision

end module changepoint_decision
