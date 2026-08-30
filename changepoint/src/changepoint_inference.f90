! Computational translation of R package changepoint 2.3.
! Upstream license declaration: GPL (unversioned).
! See NOTICE.md and UPSTREAM.md for authorship and provenance.
module changepoint_inference
use r_kinds, only : dp
use changepoint_costs, only : cp_cost_mean_normal, cp_cost_var_normal, cp_cost_meanvar_normal
use changepoint_costs, only : cp_cost_exponential
implicit none
private
public :: cp_amoc_asymptotic_value

contains

subroutine cp_amoc_asymptotic_value(cost_code, n, null_cost, alt_cost, value, status)
    integer, intent(in) :: cost_code, n
    real(dp), intent(in) :: null_cost, alt_cost
    real(dp), intent(out) :: value
    integer, intent(out) :: status
    real(dp) :: a, b, difference, pi

    value = 0.0_dp
    status = 0
    if (n < 3) then
        status = 1
        return
    end if
    difference = sqrt(abs(null_cost - alt_cost))
    pi = acos(-1.0_dp)
    select case (cost_code)
    case (cp_cost_mean_normal)
        a = (2.0_dp * log(log(real(n, dp))))**(-0.5_dp)
        b = 1.0_dp / a + 0.5_dp * a * log(log(log(real(n, dp))))
        value = exp(-2.0_dp * sqrt(pi) * exp(-a * difference + b / a)) - &
            exp(-2.0_dp * sqrt(pi) * exp(b / a))
    case (cp_cost_var_normal)
        a = sqrt(2.0_dp * log(log(real(n, dp))))
        b = 2.0_dp * log(log(real(n, dp))) + 0.5_dp * log(log(log(real(n, dp)))) - log(gamma(0.5_dp))
        value = exp(-2.0_dp * exp(-a * difference + b)) - exp(-2.0_dp * exp(b))
    case (cp_cost_meanvar_normal)
        a = sqrt(2.0_dp * log(log(real(n, dp))))
        b = 2.0_dp * log(log(real(n, dp))) + log(log(log(real(n, dp))))
        value = exp(-2.0_dp * exp(-a * difference + b)) - exp(-2.0_dp * exp(b))
    case (cp_cost_exponential)
        a = sqrt(2.0_dp * log(log(real(n, dp))))
        b = 2.0_dp * log(log(real(n, dp))) + 0.5_dp * log(log(log(real(n, dp)))) - 0.5_dp * log(pi)
        value = exp(-2.0_dp * exp(-a * difference + a * b)) - exp(-2.0_dp * exp(a * b))
    case default
        status = 2
        return
    end select
end subroutine cp_amoc_asymptotic_value

end module changepoint_inference
