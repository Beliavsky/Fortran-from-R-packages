! Computational translation of R package changepoint 2.3.
! Upstream license declaration: GPL (unversioned).
! See NOTICE.md and UPSTREAM.md for authorship and provenance.
module changepoint_penalties
use r_kinds, only : dp
implicit none
private
public :: cp_penalty_value

contains

subroutine cp_penalty_value(name, n, diffparam, value, status, manual_value, alpha, asymcheck)
    character(len=*), intent(in) :: name
    integer, intent(in) :: n, diffparam
    real(dp), intent(out) :: value
    integer, intent(out) :: status
    real(dp), intent(in), optional :: manual_value, alpha
    character(len=*), intent(in), optional :: asymcheck
    real(dp) :: a, alogn, blogn, pi, top, bottom
    character(len=:), allocatable :: check

    status = 0
    value = 0.0_dp
    if (n <= 1 .or. diffparam < 0) then
        status = 1
        return
    end if
    pi = acos(-1.0_dp)
    select case (trim(name))
    case ('SIC0', 'BIC0')
        value = real(diffparam, dp) * log(real(n, dp))
    case ('SIC', 'BIC')
        value = real(diffparam + 1, dp) * log(real(n, dp))
    case ('MBIC')
        value = real(diffparam + 2, dp) * log(real(n, dp))
    case ('AIC0')
        value = 2.0_dp * real(diffparam, dp)
    case ('AIC')
        value = 2.0_dp * real(diffparam + 1, dp)
    case ('Hannan-Quinn0')
        value = 2.0_dp * real(diffparam, dp) * log(log(real(n, dp)))
    case ('Hannan-Quinn')
        value = 2.0_dp * real(diffparam + 1, dp) * log(log(real(n, dp)))
    case ('None')
        value = 0.0_dp
    case ('Manual')
        if (.not. present(manual_value)) then
            status = 2
            return
        end if
        value = manual_value
    case ('Asymptotic')
        if (.not. present(alpha) .or. .not. present(asymcheck)) then
            status = 3
            return
        end if
        a = alpha
        if (a <= 0.0_dp .or. a > 1.0_dp) then
            status = 4
            return
        end if
        check = trim(asymcheck)
        select case (check)
        case ('mean.norm')
            alogn = (2.0_dp * log(log(real(n, dp))))**(-0.5_dp)
            blogn = 1.0_dp / alogn + 0.5_dp * alogn * log(log(log(real(n, dp))))
            value = (-alogn * log(log((1.0_dp - a + exp(-2.0_dp * sqrt(pi) * &
                exp(blogn / alogn)))**(-1.0_dp / (2.0_dp * sqrt(pi))))) + blogn)**2
        case ('var.norm')
            alogn = sqrt(2.0_dp * log(log(real(n, dp))))
            blogn = 2.0_dp * log(log(real(n, dp))) + 0.5_dp * log(log(log(real(n, dp)))) - log(gamma(0.5_dp))
            value = (-log(log((1.0_dp - a + exp(-2.0_dp * exp(blogn)))**(-0.5_dp))) / alogn + blogn / alogn)**2
        case ('meanvar.norm')
            alogn = sqrt(2.0_dp * log(log(real(n, dp))))
            blogn = 2.0_dp * log(log(real(n, dp))) + log(log(log(real(n, dp))))
            value = (-log(log((1.0_dp - a + exp(-2.0_dp * exp(blogn)))**(-0.5_dp))) / alogn + blogn / alogn)**2
        case ('reg.norm')
            top = -log(log((1.0_dp - a + exp(-2.0_dp * exp(2.0_dp * log(log(real(n, dp))) + &
                0.5_dp * real(diffparam, dp) * log(log(log(real(n, dp)))) - log(gamma(0.5_dp * real(diffparam, dp)))))) &
                **(-0.5_dp))) + 2.0_dp * log(log(real(n, dp))) + &
                0.5_dp * real(diffparam, dp) * log(log(log(real(n, dp)))) - log(gamma(0.5_dp * real(diffparam, dp)))
            bottom = sqrt(2.0_dp * log(log(real(n, dp))))
            value = (top / bottom)**2
        case ('var.css')
            call css_critical_value(a, value, status)
        case ('meanvar.exp')
            alogn = sqrt(2.0_dp * log(log(real(n, dp))))
            blogn = 2.0_dp * log(log(real(n, dp))) + 0.5_dp * log(log(log(real(n, dp)))) - 0.5_dp * log(pi)
            if (abs(a - 1.0_dp) <= epsilon(1.0_dp)) then
                value = 1.42417_dp
            else
                value = -log(-0.5_dp * log(1.0_dp - a)) / alogn + blogn
            end if
        case default
            status = 5
            return
        end select
    case default
        status = 6
        return
    end select
    if (value < 0.0_dp) status = 7
end subroutine cp_penalty_value

subroutine css_critical_value(alpha, value, status)
    real(dp), intent(in) :: alpha
    real(dp), intent(out) :: value
    integer, intent(out) :: status
    real(dp), parameter :: tol = 1.0e-12_dp

    status = 0
    if (abs(alpha - 0.01_dp) < tol) then
        value = 1.628_dp
    else if (abs(alpha - 0.05_dp) < tol) then
        value = 1.358_dp
    else if (abs(alpha - 0.10_dp) < tol) then
        value = 1.224_dp
    else if (abs(alpha - 0.25_dp) < tol) then
        value = 1.019_dp
    else if (abs(alpha - 0.50_dp) < tol) then
        value = 0.828_dp
    else if (abs(alpha - 0.75_dp) < tol) then
        value = 0.677_dp
    else if (abs(alpha - 0.90_dp) < tol) then
        value = 0.571_dp
    else if (abs(alpha - 0.95_dp) < tol) then
        value = 0.520_dp
    else
        value = 0.0_dp
        status = 8
    end if
end subroutine css_critical_value

end module changepoint_penalties
