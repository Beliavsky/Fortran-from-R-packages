! SPDX-License-Identifier: MIT
! Copyright (c) 2014-2022 Robert L. McDonald
program test_simulation_asian
    use derivmkts, only: dp,simprice,simulation_result
    use derivmkts, only: arithasianmc,geomasianmc,arithavgpricecv,asian_mc_result
    implicit none
    type(simulation_result)::x,y,multi
    type(asian_mc_result)::a,g,cv
    real(dp)::cov(2,2),rates(2),divs(2),lam(2),alpha(2),jv(2)
    x=simprice(40.0_dp,0.30_dp,0.08_dp,1.0_dp,0.0_dp,50,12,.true.,2.0_dp,-0.2_dp,0.4_dp,123)
    y=simprice(40.0_dp,0.30_dp,0.08_dp,1.0_dp,0.0_dp,50,12,.true.,2.0_dp,-0.2_dp,0.4_dp,123)
    if(.not.x%valid .or. .not.y%valid)error stop 1
    if(maxval(abs(x%price-y%price))>epsilon(1.0_dp) .or. any(x%jumps/=y%jumps))error stop 1
    if(any(x%price<=0.0_dp) .or. any(x%jumps<0))error stop 1
    cov=reshape([0.04_dp,0.018_dp,0.018_dp,0.09_dp],[2,2])
    rates=[0.05_dp,0.04_dp]
    divs=[0.01_dp,0.02_dp]
    lam=[0.0_dp,0.0_dp]
    alpha=0.0_dp
    jv=0.0_dp
    multi=simprice(100.0_dp,cov,rates,0.5_dp,divs,20,6,.false.,lam,alpha,jv,17)
    if(.not.multi%valid)error stop 1
    if(any(shape(multi%price)/=[20,7,2]))error stop 1

    a=arithasianmc(40.0_dp,40.0_dp,0.30_dp,0.08_dp,0.25_dp,0.0_dp,3,30000,987)
    g=geomasianmc(40.0_dp,40.0_dp,0.30_dp,0.08_dp,0.25_dp,0.0_dp,3,30000,987)
    cv=arithavgpricecv(40.0_dp,40.0_dp,0.30_dp,0.08_dp,0.25_dp,0.0_dp,3,30000,987)
    if(a%avg_price_call<=0.0_dp .or. a%avg_price_put<=0.0_dp)error stop 1
    call assert_close(g%avg_price_call,g%exact_avg_price_call,0.025_dp)
    call assert_close(g%avg_price_put,g%exact_avg_price_put,0.025_dp)
    if(cv%sd_avg_price_call>=a%sd_avg_price_call)error stop 1
    if(abs(cv%beta)<=0.1_dp)error stop 1
    print '(a)', 'test_simulation_asian: PASS'
contains
    subroutine assert_close(actual,expected,tol)
        real(dp),intent(in)::actual,expected,tol
        if(abs(actual-expected)>tol*(1.0_dp+abs(expected)))then
            print '(a,3es24.16)','mismatch: ',actual,expected,abs(actual-expected)
            error stop 1
        end if
    end subroutine assert_close
end program test_simulation_asian
