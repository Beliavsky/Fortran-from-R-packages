! SPDX-License-Identifier: MIT
! Copyright (c) 2014-2022 Robert L. McDonald
program test_binomial_jumps
    use derivmkts, only: dp,bscall,bsput,binomopt,binomial_result
    use derivmkts, only: mertonjump,assetjump,cashjump,option_pair
    implicit none
    type(binomial_result)::eurc,eurp,amp
    type(option_pair)::j,a,c
    real(dp),parameter::s=40.0_dp,k=40.0_dp,v=0.30_dp,r=0.08_dp,tt=2.0_dp,d=0.05_dp
    eurc=binomopt(s,k,v,r,tt,d,nstep=800,american=.false.,putopt=.false.)
    eurp=binomopt(s,k,v,r,tt,d,nstep=800,american=.false.,putopt=.true.)
    amp=binomopt(s,k,v,r,tt,d,nstep=800,american=.true.,putopt=.true.,returntrees=.true.)
    if(.not.eurc%valid .or. .not.eurp%valid .or. .not.amp%valid)error stop 1
    call assert_close(eurc%price,bscall(s,k,v,r,tt,d),4.0e-4_dp)
    call assert_close(eurp%price,bsput(s,k,v,r,tt,d),4.0e-4_dp)
    if(amp%price+1.0e-12_dp<eurp%price)error stop 1
    if(size(amp%stock_tree,1)/=801 .or. size(amp%option_tree,2)/=801)error stop 1
    call assert_close(sum(amp%probability_tree(:,801)),1.0_dp,2.0e-12_dp)

    j=mertonjump(s,k,v,r,tt,d,0.0_dp,-0.15_dp,0.20_dp)
    call assert_close(j%call,bscall(s,k,v,r,tt,d),1.0e-13_dp)
    call assert_close(j%put,bsput(s,k,v,r,tt,d),1.0e-13_dp)
    j=mertonjump(s,k,v,r,tt,d,0.75_dp,-0.05_dp,0.35_dp)
    a=assetjump(s,k,v,r,tt,d,0.75_dp,-0.05_dp,0.35_dp)
    c=cashjump(s,k,v,r,tt,d,0.75_dp,-0.05_dp,0.35_dp)
    call assert_close(j%call,a%call-k*c%call,2.0e-11_dp)
    call assert_close(j%put,k*c%put-a%put,2.0e-11_dp)
    call assert_close(j%call-j%put,s*exp(-d*tt)-k*exp(-r*tt),2.0e-10_dp)
    print '(a)', 'test_binomial_jumps: PASS'
contains
    subroutine assert_close(actual,expected,tol)
        real(dp),intent(in)::actual,expected,tol
        if(abs(actual-expected)>tol*(1.0_dp+abs(expected)))then
            print '(a,3es24.16)','mismatch: ',actual,expected,abs(actual-expected)
            error stop 1
        end if
    end subroutine assert_close
end program test_binomial_jumps
