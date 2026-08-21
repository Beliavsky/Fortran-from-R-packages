! SPDX-License-Identifier: GPL-2.0-or-later
program test_fit
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf
    use adequacy_model, only: dp, goodness_result, goodness_fit
    implicit none
    integer, parameter :: n = 100
    real(dp) :: data(n), starts(1), domain(2), expected
    type(goodness_result) :: g
    integer :: i

    do i = 1, n
        data(i) = -log(1.0_dp - (real(i,dp)-0.5_dp)/real(n,dp)) / 2.0_dp
    end do
    expected = real(n,dp)/sum(data)
    starts = [1.5_dp]
    domain = [0.0_dp, ieee_value(1.0_dp, ieee_positive_inf)]
    call goodness_fit(exp_pdf, exp_cdf, starts, data, g, method='BFGS', domain=domain, max_iter=500)
    if (abs(g%mle(1)-expected) > 2.0e-4_dp) then
        print *, 'mle', g%mle(1), expected, g%convergence
        error stop 1
    end if
    if (size(g%se) /= 1 .or. g%se(1) <= 0.0_dp) error stop 'se'
    if (abs(g%pdf_integral-1.0_dp) > 2.0e-5_dp) error stop 'integral'
    print '(a)', 'test_fit: PASS'
contains
    function exp_pdf(par, x) result(v)
        real(dp), intent(in) :: par(:), x
        real(dp) :: v
        if (par(1) <= 0.0_dp .or. x < 0.0_dp) then
            v = 0.0_dp
        else
            v = par(1)*exp(-par(1)*x)
        end if
    end function exp_pdf

    function exp_cdf(par, x) result(v)
        real(dp), intent(in) :: par(:), x
        real(dp) :: v
        if (x <= 0.0_dp) then
            v = 0.0_dp
        else
            v = 1.0_dp-exp(-par(1)*x)
        end if
    end function exp_cdf
end program test_fit
