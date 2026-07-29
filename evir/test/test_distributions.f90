! SPDX-License-Identifier: GPL-2.0-or-later
program test_distributions
    use evir
    use test_support
    implicit none
    type(evir_rng)::rng
    real(dp)::p,x,d
    real(dp)::draws(1000)
    integer::i

    do i=1,9
        p=real(i,dp)/10.0_dp
        x=qgev(p,0.2_dp,1.0_dp,2.0_dp)
        call check_close(pgev(x,0.2_dp,1.0_dp,2.0_dp),p,1.0e-12_dp,'GEV inversion')
        d=dgev(x,0.2_dp,1.0_dp,2.0_dp)
        call check(d>0.0_dp,'GEV density positive')
        x=qgpd(p,-0.2_dp,1.0_dp,2.0_dp)
        call check_close(pgpd(x,-0.2_dp,1.0_dp,2.0_dp),p,1.0e-12_dp,'GPD inversion')
        call check(dgpd(x,-0.2_dp,1.0_dp,2.0_dp)>0.0_dp,'GPD density positive')
    end do
    call check_close(pgev(qgev(0.9_dp,0.0_dp),0.0_dp),0.9_dp,1.0e-12_dp,'Gumbel limit')
    call check_close(pgpd(qgpd(0.9_dp,0.0_dp),0.0_dp),0.9_dp,1.0e-12_dp,'exponential limit')
    call seed_rng(rng,123456_8)
    call rgpd(size(draws),draws,rng,0.2_dp)
    call check(minval(draws)>=0.0_dp,'random GPD support')
    call rgev(size(draws),draws,rng,0.0_dp)
    call check(all(draws<huge(1.0_dp)),'random GEV finite')
    call finish_tests()
end program test_distributions
