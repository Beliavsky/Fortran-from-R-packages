! SPDX-License-Identifier: LGPL-3.0-or-later
! Thin compatibility adapter from sadists' internal approximation names to
! the standalone PDQutils-fortran implementation.
module sadists_approximations
    use pdqutils, only : edgeworth_pdf => dapx_edgeworth, &
        edgeworth_cdf => papx_edgeworth, &
        cornish_fisher_quantile => qapx_cf, &
        as269
    implicit none
    private
    public :: edgeworth_pdf, edgeworth_cdf, cornish_fisher_quantile, as269
end module sadists_approximations
