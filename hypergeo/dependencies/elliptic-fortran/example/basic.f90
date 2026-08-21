! SPDX-License-Identifier: GPL-2.0-only
program basic
    use elliptic, only : dp, sn, elliptic_parameters, parameters_from_g, wp, &
        weierstrass_sigma, j_invariant
    implicit none
    type(elliptic_parameters) :: par
    complex(dp) :: g(2)

    print '(a,2f20.12)', 'sn(0.61802 | m=0.5) = ', &
        sn((0.61802_dp,0.0_dp),(0.5_dp,0.0_dp))

    g=[(10.0_dp,0.0_dp),(2.0_dp,0.0_dp)]
    par=parameters_from_g(g)
    print '(a,2es22.12)', 'P(0.07+0.1i; g2=10,g3=2) = ', &
        wp((0.07_dp,0.1_dp),params=par)

    print '(a,2f20.12)', 'J(i) = ',j_invariant((0.0_dp,1.0_dp))
    print '(a,2es22.12)', 'sigma(0.3+0.2i) = ', &
        weierstrass_sigma((0.3_dp,0.2_dp),params=par)
end program basic
