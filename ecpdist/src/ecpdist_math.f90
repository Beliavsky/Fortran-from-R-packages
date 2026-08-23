! SPDX-License-Identifier: GPL-3.0-only
module ecpdist_math
    use ecpdist_kinds, only: dp
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
    implicit none
    private
    public :: expm1_stable, log1p_stable, log_abs_expm1, logsumexp2
    public :: nan_dp, inf_dp

contains

    pure real(dp) function nan_dp() result(x)
        x = ieee_value(0.0_dp, ieee_quiet_nan)
    end function nan_dp

    pure real(dp) function inf_dp() result(x)
        x = ieee_value(0.0_dp, ieee_positive_inf)
    end function inf_dp

    pure real(dp) function expm1_stable(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: x2

        if (abs(x) > 1.0e-5_dp) then
            y = exp(x) - 1.0_dp
        else
            x2 = x*x
            y = x + 0.5_dp*x2 + x*x2/6.0_dp + x2*x2/24.0_dp + &
                x2*x2*x/120.0_dp + x2*x2*x2/720.0_dp
        end if
    end function expm1_stable

    pure real(dp) function log1p_stable(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: term, sumv
        integer :: k

        if (abs(x) > 1.0e-4_dp) then
            y = log(1.0_dp + x)
            return
        end if

        term = x
        sumv = 0.0_dp
        do k = 1, 12
            if (mod(k, 2) == 1) then
                sumv = sumv + term/real(k, dp)
            else
                sumv = sumv - term/real(k, dp)
            end if
            term = term*x
        end do
        y = sumv
    end function log1p_stable

    pure real(dp) function log_abs_expm1(x) result(y)
        real(dp), intent(in) :: x

        if (x > 50.0_dp) then
            y = x + log1p_stable(-exp(-x))
        else
            y = log(abs(expm1_stable(x)))
        end if
    end function log_abs_expm1

    pure real(dp) function logsumexp2(a, b) result(y)
        real(dp), intent(in) :: a, b
        real(dp) :: m

        m = max(a, b)
        y = m + log(exp(a - m) + exp(b - m))
    end function logsumexp2

end module ecpdist_math
