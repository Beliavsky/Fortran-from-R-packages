! SPDX-License-Identifier: GPL-2.0-or-later
program block_maxima_example
    use evir, only : dp, gev, rlevel_gev, gev_fit_result
    implicit none
    type(gev_fit_result) :: fit
    real(dp) :: x(20)
    integer :: i

    do i = 1, size(x)
        x(i) = 10.0_dp + sin(real(i,dp)) + 0.03_dp*real(i,dp)
    end do
    fit = gev(x, block_size=5)
    print '(a,3f12.6)', 'xi, sigma, mu: ', fit%xi, fit%sigma, fit%mu
    print '(a,f12.6)', '20-block return level: ', rlevel_gev(fit,20.0_dp)
end program block_maxima_example
