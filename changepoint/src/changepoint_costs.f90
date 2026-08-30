! Computational translation of R package changepoint 2.3.
! Upstream license declaration: GPL (unversioned).
! See NOTICE.md and UPSTREAM.md for authorship and provenance.
module changepoint_costs
use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
use r_kinds, only : dp
implicit none
private

integer, parameter, public :: cp_cost_mean_normal = 1
integer, parameter, public :: cp_cost_var_normal = 2
integer, parameter, public :: cp_cost_meanvar_normal = 3
integer, parameter, public :: cp_cost_exponential = 4
integer, parameter, public :: cp_cost_gamma = 5
integer, parameter, public :: cp_cost_poisson = 6

type, public :: cp_cost_model
    integer :: code = cp_cost_mean_normal
    real(dp) :: shape = 1.0_dp
    real(dp) :: reference_mean = 0.0_dp
    logical :: mbic = .false.
    real(dp), allocatable :: sx(:)
    real(dp), allocatable :: sx2(:)
    real(dp), allocatable :: scentered2(:)
end type cp_cost_model

public :: cp_build_cost_model
public :: cp_segment_cost
public :: cp_valid_cost_code

contains

subroutine cp_build_cost_model(data, code, model, shape, known_mean, mbic, status)
    real(dp), intent(in) :: data(:)
    integer, intent(in) :: code
    type(cp_cost_model), intent(out) :: model
    real(dp), intent(in), optional :: shape
    real(dp), intent(in), optional :: known_mean
    logical, intent(in), optional :: mbic
    integer, intent(out), optional :: status
    integer :: i, n
    real(dp) :: mu

    if (present(status)) status = 0
    if (.not. cp_valid_cost_code(code)) then
        if (present(status)) status = 1
        return
    end if
    n = size(data)
    model%code = code
    model%shape = 1.0_dp
    if (present(shape)) model%shape = shape
    model%mbic = .false.
    if (present(mbic)) model%mbic = mbic
    if (code == cp_cost_gamma .and. model%shape <= 0.0_dp) then
        if (present(status)) status = 2
        return
    end if
    if (.not. all(ieee_is_finite(data))) then
        if (present(status)) status = 3
        return
    end if
    if (code == cp_cost_exponential .and. any(data < 0.0_dp)) then
        if (present(status)) status = 4
        return
    end if
    if (code == cp_cost_gamma .and. any(data <= 0.0_dp)) then
        if (present(status)) status = 5
        return
    end if
    if (code == cp_cost_poisson) then
        if (any(data < 0.0_dp)) then
            if (present(status)) status = 6
            return
        end if
        if (any(abs(data - anint(data)) > 16.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(data)))) then
            if (present(status)) status = 7
            return
        end if
    end if

    allocate(model%sx(0:n), model%sx2(0:n), model%scentered2(0:n))
    model%sx = 0.0_dp
    model%sx2 = 0.0_dp
    model%scentered2 = 0.0_dp
    if (present(known_mean)) then
        mu = known_mean
    else if (n > 0) then
        mu = sum(data) / real(n, dp)
    else
        mu = 0.0_dp
    end if
    model%reference_mean = mu
    do i = 1, n
        model%sx(i) = model%sx(i - 1) + data(i)
        model%sx2(i) = model%sx2(i - 1) + data(i) * data(i)
        model%scentered2(i) = model%scentered2(i - 1) + (data(i) - mu) * (data(i) - mu)
    end do
end subroutine cp_build_cost_model

pure logical function cp_valid_cost_code(code)
    integer, intent(in) :: code
    cp_valid_cost_code = code >= cp_cost_mean_normal .and. code <= cp_cost_poisson
end function cp_valid_cost_code

pure real(dp) function cp_segment_cost(model, left, right) result(cost)
    type(cp_cost_model), intent(in) :: model
    integer, intent(in) :: left, right
    integer :: nseg
    real(dp) :: sumx, sumx2, ss, sigsq, pi
    real(dp), parameter :: tiny_variance = 1.0e-11_dp

    if (left < 0 .or. right <= left .or. right > ubound(model%sx, 1)) then
        cost = huge(1.0_dp)
        return
    end if
    nseg = right - left
    sumx = model%sx(right) - model%sx(left)
    sumx2 = model%sx2(right) - model%sx2(left)
    pi = acos(-1.0_dp)

    select case (model%code)
    case (cp_cost_mean_normal)
        cost = sumx2 - sumx * sumx / real(nseg, dp)
    case (cp_cost_var_normal)
        ss = model%scentered2(right) - model%scentered2(left)
        ss = max(ss, tiny_variance)
        cost = real(nseg, dp) * (log(2.0_dp * pi) + log(ss / real(nseg, dp)) + 1.0_dp)
    case (cp_cost_meanvar_normal)
        sigsq = (sumx2 - sumx * sumx / real(nseg, dp)) / real(nseg, dp)
        sigsq = max(sigsq, tiny_variance)
        cost = real(nseg, dp) * (log(2.0_dp * pi) + log(sigsq) + 1.0_dp)
    case (cp_cost_exponential)
        sumx = max(sumx, tiny(1.0_dp))
        cost = 2.0_dp * real(nseg, dp) * (log(sumx) - log(real(nseg, dp)))
    case (cp_cost_gamma)
        if (sumx <= 0.0_dp) then
            cost = huge(1.0_dp)
        else
            cost = 2.0_dp * real(nseg, dp) * model%shape * &
                (log(sumx) - log(real(nseg, dp) * model%shape))
        end if
    case (cp_cost_poisson)
        if (sumx <= 0.0_dp) then
            cost = 0.0_dp
        else
            cost = 2.0_dp * sumx * (log(real(nseg, dp)) - log(sumx))
        end if
    case default
        cost = huge(1.0_dp)
    end select

    if (model%mbic) then
        if (.not. (model%code == cp_cost_poisson .and. sumx <= 0.0_dp)) then
            cost = cost + log(real(nseg, dp))
        end if
    end if
end function cp_segment_cost

end module changepoint_costs
