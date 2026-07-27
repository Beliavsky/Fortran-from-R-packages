! SPDX-License-Identifier: MIT
! Copyright (c) 2014-2022 Robert L. McDonald
module test_greeks_support
    use derivmkts, only: dp,bscall
    implicit none
contains
    pure real(dp) function call_wrapper(s,k,v,r,tt,d) result(price)
        real(dp),intent(in)::s,k,v,r,tt,d
        price=bscall(s,k,v,r,tt,d)
    end function call_wrapper
end module test_greeks_support

program test_greeks_misc
    use derivmkts, only: dp,bs_call_greeks,numerical_greeks,greek_result
    use derivmkts, only: quincunx,quincunx_result,binormsdist_discrete
    use test_greeks_support, only: call_wrapper
    implicit none
    type(greek_result)::ga,gn
    type(quincunx_result)::q
    real(dp)::probs(0:10)
    ga=bs_call_greeks(40.0_dp,40.0_dp,0.30_dp,0.08_dp,0.25_dp,0.0_dp)
    gn=numerical_greeks(call_wrapper,40.0_dp,40.0_dp,0.30_dp,0.08_dp,0.25_dp,0.0_dp)
    call assert_close(gn%premium,ga%premium,1.0e-13_dp)
    call assert_close(gn%delta,ga%delta,2.0e-7_dp)
    call assert_close(gn%gamma,ga%gamma,2.0e-5_dp)
    call assert_close(gn%vega,ga%vega,2.0e-7_dp)
    call assert_close(gn%rho,ga%rho,2.0e-7_dp)
    call assert_close(gn%theta,ga%theta,3.0e-7_dp)
    call assert_close(gn%psi,ga%psi,2.0e-7_dp)
    q=quincunx(10,1000,0.6_dp,123)
    if(sum(q%counts)/=1000)error stop 1
    call assert_close(sum(q%expected),1000.0_dp,2.0e-12_dp)
    call binormsdist_discrete(10,0.6_dp,probs)
    call assert_close(sum(probs),1.0_dp,2.0e-14_dp)
    print '(a)', 'test_greeks_misc: PASS'
contains
    subroutine assert_close(actual,expected,tol)
        real(dp),intent(in)::actual,expected,tol
        if(abs(actual-expected)>tol*(1.0_dp+abs(expected)))then
            print '(a,3es24.16)','mismatch: ',actual,expected,abs(actual-expected)
            error stop 1
        end if
    end subroutine assert_close
end program test_greeks_misc
