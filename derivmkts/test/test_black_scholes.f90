! SPDX-License-Identifier: MIT
! Copyright (c) 2014-2022 Robert L. McDonald
program test_black_scholes
    use derivmkts, only: dp, bscall, bsput, assetcall, cashcall, assetput, cashput
    use derivmkts, only: bscallimpvol, bsputimpvol, bscallimps, bsputimps
    use derivmkts, only: bondpv, bondyield, duration, convexity
    use derivmkts, only: geomavgpricecall, geomavgpriceput, geomavgstrikecall, geomavgstrikeput
    implicit none
    real(dp) :: c,p,price,y
    logical :: ok

    c=bscall(40.0_dp,40.0_dp,0.30_dp,0.08_dp,0.25_dp,0.0_dp)
    p=bsput(40.0_dp,40.0_dp,0.30_dp,0.08_dp,0.25_dp,0.0_dp)
    call assert_close(c,2.784736657821661_dp,2.0e-13_dp)
    call assert_close(p,1.992683590091872_dp,2.0e-13_dp)
    call assert_close(c,assetcall(40.0_dp,40.0_dp,0.30_dp,0.08_dp,0.25_dp,0.0_dp)- &
        40.0_dp*cashcall(40.0_dp,40.0_dp,0.30_dp,0.08_dp,0.25_dp,0.0_dp),2.0e-13_dp)
    call assert_close(p,40.0_dp*cashput(40.0_dp,40.0_dp,0.30_dp,0.08_dp,0.25_dp,0.0_dp)- &
        assetput(40.0_dp,40.0_dp,0.30_dp,0.08_dp,0.25_dp,0.0_dp),2.0e-13_dp)

    call assert_close(bscallimpvol(40.0_dp,40.0_dp,0.08_dp,0.25_dp,0.0_dp,4.0_dp,ok=ok), &
        0.4554460631412654_dp,2.0e-11_dp)
    if(.not.ok)error stop 1
    call assert_close(bsputimpvol(40.0_dp,40.0_dp,0.08_dp,0.25_dp,0.0_dp,4.0_dp,ok=ok), &
        0.5568408186010122_dp,3.0e-11_dp)
    if(.not.ok)error stop 1
    call assert_close(bscallimps(40.0_dp,40.0_dp,0.30_dp,0.08_dp,0.25_dp,0.0_dp,c,ok=ok),40.0_dp,2.0e-10_dp)
    if(.not.ok)error stop 1
    call assert_close(bsputimps(40.0_dp,40.0_dp,0.30_dp,0.08_dp,0.25_dp,0.0_dp,p,ok=ok),40.0_dp,2.0e-10_dp)
    if(.not.ok)error stop 1

    price=bondpv(6.0_dp,20.0_dp,0.045_dp,100.0_dp,2)
    call assert_close(price,119.64514165431042_dp,2.0e-12_dp)
    y=bondyield(price,6.0_dp,20.0_dp,100.0_dp,2,ok)
    if(.not.ok)error stop 1
    call assert_close(y,0.045_dp,2.0e-10_dp)
    call assert_close(duration(price,6.0_dp,20.0_dp,100.0_dp,2,.false.),12.635375265069701_dp,3.0e-11_dp)
    call assert_close(duration(price,6.0_dp,20.0_dp,100.0_dp,2,.true.),12.357335222562055_dp,3.0e-11_dp)
    call assert_close(convexity(price,6.0_dp,20.0_dp,100.0_dp,2),205.97896747386878_dp,3.0e-10_dp)

    call assert_close(geomavgpricecall(40.0_dp,40.0_dp,0.30_dp,0.08_dp,0.25_dp,0.0_dp,3), &
        1.938526326081995_dp,3.0e-13_dp)
    call assert_close(geomavgpriceput(40.0_dp,40.0_dp,0.30_dp,0.08_dp,0.25_dp,0.0_dp,3), &
        1.478421552797165_dp,3.0e-13_dp)
    call assert_close(geomavgstrikecall(40.0_dp,40.0_dp,0.30_dp,0.08_dp,0.25_dp,0.0_dp,3), &
        1.200118417301642_dp,3.0e-13_dp)
    call assert_close(geomavgstrikeput(40.0_dp,40.0_dp,0.30_dp,0.08_dp,0.25_dp,0.0_dp,3), &
        0.8681701228566823_dp,3.0e-13_dp)
    print '(a)', 'test_black_scholes: PASS'
contains
    subroutine assert_close(actual,expected,tol)
        real(dp),intent(in)::actual,expected,tol
        if(abs(actual-expected)>tol*(1.0_dp+abs(expected)))then
            print '(a,3es24.16)', 'mismatch: ',actual,expected,abs(actual-expected)
            error stop 1
        end if
    end subroutine assert_close
end program test_black_scholes
