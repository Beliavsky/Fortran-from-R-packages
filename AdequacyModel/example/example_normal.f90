! SPDX-License-Identifier: GPL-2.0-or-later
program example_normal
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf, ieee_negative_inf
    use adequacy_model, only: dp, goodness_result, goodness_from_mle
    implicit none
    real(dp) :: x(8), par(2), domain(2)
    type(goodness_result) :: fit

    x = [-1.1_dp, -0.4_dp, 0.0_dp, 0.2_dp, 0.6_dp, 0.9_dp, 1.2_dp, 1.8_dp]
    par = [0.4_dp, 0.9_dp]
    domain = [ieee_value(1.0_dp, ieee_negative_inf), ieee_value(1.0_dp, ieee_positive_inf)]
    call goodness_from_mle(pdf, cdf, par, x, fit, domain)
    print '(a,f10.5)', 'AIC = ', fit%aic
    print '(a,f10.5)', 'KS  = ', fit%ks
contains
    function pdf(p, z) result(v)
        real(dp), intent(in) :: p(:), z
        real(dp) :: v, q
        q = (z-p(1))/p(2)
        v = exp(-0.5_dp*q*q)/(sqrt(2.0_dp*acos(-1.0_dp))*p(2))
    end function pdf

    function cdf(p, z) result(v)
        real(dp), intent(in) :: p(:), z
        real(dp) :: v
        v = 0.5_dp*erfc(-(z-p(1))/(p(2)*sqrt(2.0_dp)))
    end function cdf
end program example_normal
