! SPDX-License-Identifier: GPL-3.0-or-later
program basic
    use nbconv, only : dp, dnbconv_mu, qnbconv_mu, nbconv_params_mu, nbconv_summary
    implicit none
    real(dp), parameter :: mus(2) = [100.0_dp, 10.0_dp]
    real(dp), parameter :: phis(2) = [5.0_dp, 8.0_dp]
    integer, parameter :: counts(6) = [0, 25, 50, 100, 150, 200]
    real(dp), parameter :: probs(5) = [0.05_dp, 0.25_dp, 0.5_dp, 0.75_dp, 0.95_dp]
    real(dp), allocatable :: pmf(:)
    integer, allocatable :: q(:)
    type(nbconv_summary) :: s
    integer :: i

    pmf = dnbconv_mu(counts, mus, phis, method="exact")
    q = qnbconv_mu(probs, 500, mus, phis, method="exact")
    s = nbconv_params_mu(mus, phis)

    print '(a,f12.5)', "mean     = ", s%mean
    print '(a,f12.5)', "variance = ", s%variance
    print '(a)', "selected PMF values:"
    do i = 1, size(counts)
        print '(i5,2x,es14.6)', counts(i), pmf(i)
    end do
    print '(a,5(1x,i0))', "quantiles =", q
end program basic
