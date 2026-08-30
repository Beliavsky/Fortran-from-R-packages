module e1071_special
    use e1071_kinds, only: dp
    implicit none
    private

    public :: normal_cdf, chi_square_sf

contains

    elemental function normal_cdf(x) result(value)
        real(dp), intent(in) :: x !! Standard-normal quantile at which the lower-tail probability is evaluated.
        real(dp) :: value

        value = 0.5_dp * erfc(-x / sqrt(2.0_dp))
    end function normal_cdf

    function chi_square_sf(x, df) result(value)
        real(dp), intent(in) :: x !! Nonnegative chi-square statistic.
        real(dp), intent(in) :: df !! Positive chi-square degrees of freedom.
        real(dp) :: value

        if (df <= 0.0_dp) error stop "chi_square_sf: df must be positive"
        if (x <= 0.0_dp) then
            value = 1.0_dp
        else
            value = regularized_gamma_q(0.5_dp * df, 0.5_dp * x)
        end if
    end function chi_square_sf

    function regularized_gamma_q(a, x) result(value)
        real(dp), intent(in) :: a !! Positive shape parameter of the incomplete gamma ratio.
        real(dp), intent(in) :: x !! Nonnegative evaluation point of the incomplete gamma ratio.
        real(dp) :: value
        real(dp), parameter :: eps = 1.0e-14_dp
        real(dp), parameter :: fpmin = 1.0e-300_dp
        integer, parameter :: maxit = 10000
        real(dp) :: sum_value
        real(dp) :: del
        real(dp) :: ap
        real(dp) :: b
        real(dp) :: c
        real(dp) :: d
        real(dp) :: h
        real(dp) :: an
        integer :: n

        if (a <= 0.0_dp .or. x < 0.0_dp) error stop "regularized_gamma_q: invalid argument"
        if (x <= 0.0_dp) then
            value = 1.0_dp
            return
        end if
        if (x < a + 1.0_dp) then
            ap = a
            sum_value = 1.0_dp / a
            del = sum_value
            do n = 1, maxit
                ap = ap + 1.0_dp
                del = del * x / ap
                sum_value = sum_value + del
                if (abs(del) < abs(sum_value) * eps) exit
            end do
            value = 1.0_dp - sum_value * exp(-x + a * log(x) - log_gamma(a))
        else
            b = x + 1.0_dp - a
            c = 1.0_dp / fpmin
            d = 1.0_dp / b
            h = d
            do n = 1, maxit
                an = -real(n, dp) * (real(n, dp) - a)
                b = b + 2.0_dp
                d = an * d + b
                if (abs(d) < fpmin) d = fpmin
                c = b + an / c
                if (abs(c) < fpmin) c = fpmin
                d = 1.0_dp / d
                del = d * c
                h = h * del
                if (abs(del - 1.0_dp) < eps) exit
            end do
            value = h * exp(-x + a * log(x) - log_gamma(a))
        end if
        value = max(0.0_dp, min(1.0_dp, value))
    end function regularized_gamma_q

end module e1071_special
