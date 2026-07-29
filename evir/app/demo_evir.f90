! SPDX-License-Identifier: GPL-2.0-or-later
program demo_evir
    use evir, only : dp, evir_rng, seed_rng, rgpd, gpd, riskmeasures, gpd_fit_result
    implicit none
    type(evir_rng) :: rng
    type(gpd_fit_result) :: fit
    real(dp), allocatable :: x(:)
    real(dp) :: p(3), q(3), es(3)

    allocate(x(2000))
    call seed_rng(rng, 20260728_8)
    call rgpd(size(x), x, rng, xi=0.2_dp, beta=1.5_dp)
    fit = gpd(x, threshold=1.0_dp, information='expected')
    p = [0.95_dp, 0.99_dp, 0.995_dp]
    call riskmeasures(fit, p, q, es)

    print '(a,f10.5)', 'xi estimate:   ', fit%xi
    print '(a,f10.5)', 'beta estimate: ', fit%beta
    print '(a)', 'probability       quantile       shortfall'
    print '(3f15.6)', p(1), q(1), es(1)
    print '(3f15.6)', p(2), q(2), es(2)
    print '(3f15.6)', p(3), q(3), es(3)
end program demo_evir
