! SPDX-License-Identifier: GPL-2.0-or-later
program test_fitting
    use evir
    use test_support
    implicit none
    type(gev_fit_result)::gf,gg
    type(gpd_fit_result)::pf,pwm
    type(pot_fit_result)::pt
    type(profile_result)::pr
    real(dp),allocatable::xgev(:),xgum(:),xgpd(:),data_all(:),times(:)
    real(dp)::p,q,lo,se,hi,qs(2),es(2)
    integer::i,n,status

    n=240
    allocate(xgev(n),xgum(n))
    do i=1,n
        p=(real(i,dp)-0.5_dp)/real(n,dp)
        xgev(i)=qgev(p,0.15_dp,5.0_dp,2.0_dp)
        xgum(i)=qgev(p,0.0_dp,3.0_dp,1.25_dp)
    end do
    gf=gev(xgev)
    call check(gf%converged,'GEV converged')
    call check_close(gf%xi,0.1497925023923835_dp,2.0e-7_dp,'GEV SciPy xi reference')
    call check_close(gf%sigma,1.99375116_dp,2.0e-7_dp,'GEV SciPy sigma reference')
    call check_close(gf%mu,5.00003276_dp,2.0e-7_dp,'GEV SciPy mu reference')
    gg=gumbel(xgum)
    call check(gg%converged,'Gumbel converged')
    call check_close(gg%sigma,1.25_dp,3.0e-2_dp,'Gumbel sigma recovery')
    call check_close(gg%mu,3.0_dp,3.0e-2_dp,'Gumbel mu recovery')
    call check_close(rlevel_gev(gf,20.0_dp),qgev(0.95_dp,gf%xi,gf%mu,gf%sigma),1.0e-12_dp,'GEV return level')
    pr=rlevel_gev_profile(gf,20.0_dp)
    call check(pr%lower<pr%estimate.and.pr%upper>pr%estimate,'GEV profile interval')

    allocate(xgpd(160),data_all(240),times(240))
    do i=1,160
        p=(real(i,dp)-0.5_dp)/160.0_dp
        xgpd(i)=1.0_dp+qgpd(p,0.25_dp,beta=1.2_dp)
    end do
    do i=1,80
        data_all(i)=0.2_dp+0.79_dp*real(i-1,dp)/79.0_dp
    end do
    data_all(81:240)=xgpd
    times=[(real(i-1,dp),i=1,240)]
    pf=gpd(data_all,threshold=1.0_dp,information='expected')
    call check(pf%converged,'GPD ML converged')
    call check_close(pf%xi,0.23939195_dp,2.0e-7_dp,'GPD SciPy xi reference')
    call check_close(pf%beta,1.20947820_dp,2.0e-7_dp,'GPD SciPy beta reference')
    pwm=gpd(data_all,threshold=1.0_dp,method='pwm')
    call check(pwm%converged.and.pwm%beta>0.0_dp,'GPD PWM')
    call gpd_q_wald(pf,0.99_dp,lo,q,se,hi,status=status)
    call check(status==evir_ok.and.lo<q.and.hi>q,'GPD Wald quantile')
    call check(gpd_sfall(pf,0.99_dp)>gpd_q(pf,0.99_dp),'shortfall exceeds quantile')
    call riskmeasures(pf,[0.95_dp,0.99_dp],qs,es)
    call check(all(es>qs),'risk measure vectors')
    pr=gpd_q_profile(pf,0.99_dp,20.0_dp,n_grid=25)
    call check(pr%lower<pr%estimate.and.pr%upper>pr%estimate,'GPD quantile profile')

    pt=pot(data_all,times=times,threshold=1.0_dp)
    call check(pt%converged.and.pt%beta>0.0_dp,'POT point process fit')
    call check(pt%n_exceed==160,'POT exceedance count')
    call check_close(gev_negloglik_value([0.2_dp,1.3_dp,2.5_dp], &
        [1.2_dp,2.3_dp,3.1_dp,4.7_dp,2.9_dp]),8.261279140320084_dp,1.0e-13_dp,'GEV fixed likelihood')
    call check_close(gpd_negloglik_value([0.25_dp,1.1_dp], &
        [0.1_dp,0.4_dp,0.9_dp,1.7_dp,3.2_dp]),6.320623144514968_dp,1.0e-13_dp,'GPD fixed likelihood')
    call check_close(pot_negloglik_value([0.15_dp,0.9_dp,0.2_dp], &
        [1.2_dp,1.7_dp,2.4_dp,4.1_dp],1.0_dp,10.0_dp),13.046554714342406_dp,1.0e-13_dp,'POT fixed likelihood')
    call finish_tests()
end program test_fitting
