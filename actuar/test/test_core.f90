program test_core
    use actuar
    implicit none
    real(dp), parameter :: tol = 5.0e-10_dp
    real(dp) :: prob(1), tmat(1,1), x

    call check_close(ppareto(3.0_dp, 2.0_dp, 5.0_dp), 1.0_dp-(1.0_dp+3.0_dp/5.0_dp)**(-2.0_dp), tol, 'pareto cdf')
    call check_close(qpareto(0.75_dp, 2.0_dp, 5.0_dp), 5.0_dp, tol, 'pareto quantile')
    call check_close(pburR_local(2.0_dp, 3.0_dp, 1.5_dp, 4.0_dp), pburr(2.0_dp,3.0_dp,1.5_dp,4.0_dp), tol, 'burr cdf')
    call check_close(pllogis(4.0_dp, 2.0_dp, 4.0_dp), 0.5_dp, tol, 'loglogistic median')
    call check_close(pgenpareto(2.0_dp, 3.0_dp, 2.0_dp, 4.0_dp), reg_beta(1.0_dp/3.0_dp,2.0_dp,3.0_dp), tol, 'genpareto cdf')
    call check_close(pgenbeta(2.0_dp,2.0_dp,3.0_dp,1.5_dp,5.0_dp), &
                     reg_beta((2.0_dp/5.0_dp)**1.5_dp,2.0_dp,3.0_dp), &
                     tol, 'genbeta cdf')
    call check_close(pgumbel(qgumbel(0.37_dp,1.2_dp,0.7_dp),1.2_dp,0.7_dp),0.37_dp,2.0e-12_dp,'gumbel roundtrip')
    call check_close(ptrbeta(qtrbeta(0.63_dp,2.0_dp,1.5_dp,3.0_dp,4.0_dp), &
                     2.0_dp,1.5_dp,3.0_dp,4.0_dp),0.63_dp,1.0e-10_dp, &
                     'transformed beta roundtrip')
    call check_close(pinvgamma(qinvgamma(0.6_dp,4.0_dp,3.0_dp),4.0_dp,3.0_dp),0.6_dp,1.0e-10_dp,'invgamma roundtrip')
    call check_close(pinvgauss(qinvgauss(0.42_dp,2.0_dp,0.3_dp),2.0_dp,0.3_dp),0.42_dp,2.0e-10_dp,'invgauss roundtrip')
    call check_close(minvgauss(1,2.0_dp,0.3_dp),2.0_dp,1.0e-12_dp,'invgauss mean')
    call check_close(minvgauss(2,2.0_dp,0.3_dp),6.4_dp,1.0e-12_dp,'invgauss second moment')
    call check_close(mgamma_act(2.0_dp,3.0_dp,2.0_dp),48.0_dp,1.0e-12_dp,'gamma moment')
    call check_close(mlnorm_act(1.0_dp,0.2_dp,0.5_dp),exp(0.325_dp),1.0e-12_dp,'lognormal mean')

    prob = [1.0_dp]
    tmat(1,1) = -2.0_dp
    call check_close(dphtype(1.25_dp,prob,tmat),2.0_dp*exp(-2.5_dp),1.0e-11_dp,'phase-type density')
    call check_close(pphtype(1.25_dp,prob,tmat),1.0_dp-exp(-2.5_dp),1.0e-11_dp,'phase-type cdf')
    call check_close(mphtype(1,prob,tmat),0.5_dp,1.0e-12_dp,'phase-type mean')
    call check_close(mgfphtype(0.5_dp,prob,tmat),2.0_dp/1.5_dp,1.0e-11_dp,'phase-type mgf')

    call check_close(dpoisinvgauss(0,2.0_dp,0.4_dp), exp((1.0_dp-sqrt(1.0_dp+2.0_dp*0.4_dp*4.0_dp))/(0.8_dp)), 1.0e-12_dp, 'PIG p0')
    call check_true(ppoisinvgauss(100,2.0_dp,0.4_dp) > 0.999999_dp, 'PIG normalization')

    print '(a)', 'test_core: PASS'
contains
    pure real(dp) function pburR_local(y,a,b,s) result(p)
        real(dp),intent(in)::y,a,b,s
        p=1.0_dp-(1.0_dp+(y/s)**b)**(-a)
    end function pburR_local
    subroutine check_close(got, expected, eps, name)
        real(dp),intent(in)::got,expected,eps
        character(*),intent(in)::name
        if(abs(got-expected)>eps*max(1.0_dp,abs(expected))) then
            print *, 'FAIL ',trim(name),got,expected
            error stop 1
        end if
    end subroutine check_close
    subroutine check_true(ok,name)
        logical,intent(in)::ok
        character(*),intent(in)::name
        if(.not.ok) then
            print *, 'FAIL ',trim(name)
            error stop 1
        end if
    end subroutine check_true
end program test_core
