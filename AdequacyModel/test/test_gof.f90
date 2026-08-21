! SPDX-License-Identifier: GPL-2.0-or-later
program test_gof
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf, ieee_negative_inf
    use adequacy_model, only: dp, goodness_result, goodness_from_mle
    implicit none
    integer, parameter :: n = 80
    real(dp) :: data(n), par(2), domain(2), p
    type(goodness_result) :: g
    integer :: i

    par = [1.0_dp, 2.0_dp]
    do i = 1, n
        p = (real(i,dp)-0.5_dp)/real(n,dp)
        data(i) = par(1) + par(2)*invnorm(p)
    end do
    domain(1) = ieee_value(1.0_dp, ieee_negative_inf)
    domain(2) = ieee_value(1.0_dp, ieee_positive_inf)
    call goodness_from_mle(normal_pdf, normal_cdf_cb, par, data, g, domain)
    if (abs(g%pdf_integral - 1.0_dp) > 2.0e-5_dp) then
        print *, 'integral', g%pdf_integral
        error stop 1
    end if
    if (.not. g%cdf_endpoints_ok) error stop 'cdf endpoints'
    if (g%ks > 0.03_dp) then
        print *, 'ks', g%ks
        error stop 1
    end if
    if (g%aic <= 0.0_dp .or. g%bic <= 0.0_dp) error stop 'criteria'
    print '(a)', 'test_gof: PASS'
contains
    function normal_pdf(pars, x) result(v)
        real(dp), intent(in) :: pars(:), x
        real(dp) :: v, z
        if (pars(2) <= 0.0_dp) then
            v = 0.0_dp
            return
        end if
        z = (x-pars(1))/pars(2)
        v = exp(-0.5_dp*z*z)/(sqrt(2.0_dp*acos(-1.0_dp))*pars(2))
    end function normal_pdf

    function normal_cdf_cb(pars, x) result(v)
        real(dp), intent(in) :: pars(:), x
        real(dp) :: v
        v = 0.5_dp*erfc(-(x-pars(1))/(pars(2)*sqrt(2.0_dp)))
    end function normal_cdf_cb

    function invnorm(prob) result(x)
        real(dp), intent(in) :: prob
        real(dp) :: x, lo, hi, mid
        integer :: k
        lo = -10.0_dp
        hi = 10.0_dp
        do k = 1, 100
            mid = 0.5_dp*(lo+hi)
            if (0.5_dp*erfc(-mid/sqrt(2.0_dp)) < prob) then
                lo = mid
            else
                hi = mid
            end if
        end do
        x = 0.5_dp*(lo+hi)
    end function invnorm
end program test_gof
