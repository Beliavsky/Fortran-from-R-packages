! SPDX-License-Identifier: GPL-2.0-or-later
program cvar_demo
    use cvar, only : dp, normal_pdf, normal_cdf, normal_quantile, &
                     var_qf, var_cdf, es_qf, es_cdf, es_pdf
    implicit none

    real(dp), parameter :: levels(3) = [0.10_dp, 0.05_dp, 0.01_dp]
    real(dp), allocatable :: var_values(:), es_values(:)
    integer :: i

    var_values = var_qf(normal_quantile, levels)
    es_values = es_qf(normal_quantile, levels)

    print '(a)', "Standard-normal lower-tail risk"
    print '(a)', " p_loss             VaR              ES"
    do i = 1, size(levels)
        print '(f7.3,2(2x,f14.8))', levels(i), var_values(i), es_values(i)
    end do

    print '(/,a,f14.8)', "VaR from CDF at 5%: ", var_cdf(normal_cdf, 0.05_dp)
    print '(a,f14.8)', "ES from CDF at 5%:  ", es_cdf(normal_cdf, 0.05_dp)
    print '(a,f14.8)', "ES from PDF at 5%:  ", &
        es_pdf(normal_pdf, normal_quantile, 0.05_dp)
end program cvar_demo
