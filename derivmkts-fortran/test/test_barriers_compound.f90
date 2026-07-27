! SPDX-License-Identifier: MIT
! Copyright (c) 2014-2022 Robert L. McDonald
program test_barriers_compound
    use derivmkts, only: dp, pi, bscall, bsput, cashcall, cashput, assetcall, assetput
    use derivmkts, only: callupin,callupout,putupin,putupout,calldownin,calldownout,putdownin,putdownout
    use derivmkts, only: cashuicall,cashuocall,cashuiput,cashuoput,cashdicall,cashdocall,cashdiput,cashdoput
    use derivmkts, only: assetuicall,assetuocall,assetuiput,assetuoput,assetdicall,assetdocall,assetdiput,assetdoput
    use derivmkts, only: dr,ur,drdeferred,urdeferred,callperpetual,putperpetual
    use derivmkts, only: binormsdist,calloncall,putoncall,callonput,putonput,compound_result,perpetual_result
    implicit none
    real(dp),parameter::s=40.0_dp,k=40.0_dp,v=0.30_dp,r=0.08_dp,tt=2.0_dp,d=0.05_dp
    real(dp)::hu,hd
    type(compound_result)::cc,pc,cp,pp
    type(perpetual_result)::perp
    hu=45.0_dp
    hd=35.0_dp
    call assert_close(callupin(s,k,v,r,tt,d,hu)+callupout(s,k,v,r,tt,d,hu),bscall(s,k,v,r,tt,d),1e-11_dp)
    call assert_close(putupin(s,k,v,r,tt,d,hu)+putupout(s,k,v,r,tt,d,hu),bsput(s,k,v,r,tt,d),1e-11_dp)
    call assert_close(calldownin(s,k,v,r,tt,d,hd)+calldownout(s,k,v,r,tt,d,hd),bscall(s,k,v,r,tt,d),1e-11_dp)
    call assert_close(putdownin(s,k,v,r,tt,d,hd)+putdownout(s,k,v,r,tt,d,hd),bsput(s,k,v,r,tt,d),1e-11_dp)
    call assert_close(cashuicall(s,k,v,r,tt,d,hu)+cashuocall(s,k,v,r,tt,d,hu),cashcall(s,k,v,r,tt,d),1e-11_dp)
    call assert_close(cashuiput(s,k,v,r,tt,d,hu)+cashuoput(s,k,v,r,tt,d,hu),cashput(s,k,v,r,tt,d),1e-11_dp)
    call assert_close(cashdicall(s,k,v,r,tt,d,hd)+cashdocall(s,k,v,r,tt,d,hd),cashcall(s,k,v,r,tt,d),1e-11_dp)
    call assert_close(cashdiput(s,k,v,r,tt,d,hd)+cashdoput(s,k,v,r,tt,d,hd),cashput(s,k,v,r,tt,d),1e-11_dp)
    call assert_close(assetuicall(s,k,v,r,tt,d,hu)+assetuocall(s,k,v,r,tt,d,hu),assetcall(s,k,v,r,tt,d),1e-10_dp)
    call assert_close(assetuiput(s,k,v,r,tt,d,hu)+assetuoput(s,k,v,r,tt,d,hu),assetput(s,k,v,r,tt,d),1e-10_dp)
    call assert_close(assetdicall(s,k,v,r,tt,d,hd)+assetdocall(s,k,v,r,tt,d,hd),assetcall(s,k,v,r,tt,d),1e-10_dp)
    call assert_close(assetdiput(s,k,v,r,tt,d,hd)+assetdoput(s,k,v,r,tt,d,hd),assetput(s,k,v,r,tt,d),1e-10_dp)
    call assert_close(dr(30.0_dp,v,r,tt,d,hd),1.0_dp,0.0_dp)
    call assert_close(ur(50.0_dp,v,r,tt,d,hu),1.0_dp,0.0_dp)
    if(drdeferred(s,v,r,tt,d,hd)<=0.0_dp .or. urdeferred(s,v,r,tt,d,hu)<=0.0_dp)error stop 1

    perp=callperpetual(40.0_dp,40.0_dp,0.30_dp,0.08_dp,0.02_dp)
    if(perp%price<=0.0_dp .or. perp%barrier<=40.0_dp)error stop 1
    perp=putperpetual(40.0_dp,40.0_dp,0.30_dp,0.08_dp,0.02_dp)
    if(perp%price<=0.0_dp .or. perp%barrier>=40.0_dp)error stop 1

    call assert_close(binormsdist(0.0_dp,0.0_dp,0.5_dp),0.25_dp+asin(0.5_dp)/(2.0_dp*pi),2e-11_dp)
    cc=calloncall(40.0_dp,41.0_dp,1.75_dp,0.30_dp,0.08_dp,0.38_dp,0.75_dp,0.05_dp)
    pc=putoncall(40.0_dp,41.0_dp,1.75_dp,0.30_dp,0.08_dp,0.38_dp,0.75_dp,0.05_dp)
    cp=callonput(40.0_dp,41.0_dp,1.75_dp,0.30_dp,0.08_dp,0.38_dp,0.75_dp,0.05_dp)
    pp=putonput(40.0_dp,41.0_dp,1.75_dp,0.30_dp,0.08_dp,0.38_dp,0.75_dp,0.05_dp)
    call assert_close(cc%critical_spot,38.05831713551019_dp,2e-10_dp)
    call assert_close(cc%price,2.661926489803303_dp,5e-10_dp)
    call assert_close(pc%price,0.4151664393511885_dp,5e-10_dp)
    call assert_close(cp%critical_spot,43.56104477472397_dp,2e-10_dp)
    call assert_close(cp%price,2.618796392865772_dp,5e-10_dp)
    call assert_close(pp%price,0.2874671742923327_dp,5e-10_dp)
    print '(a)', 'test_barriers_compound: PASS'
contains
    subroutine assert_close(actual,expected,tol)
        real(dp),intent(in)::actual,expected,tol
        if(abs(actual-expected)>tol*(1.0_dp+abs(expected)))then
            print '(a,3es24.16)','mismatch: ',actual,expected,abs(actual-expected)
            error stop 1
        end if
    end subroutine assert_close
end program test_barriers_compound
