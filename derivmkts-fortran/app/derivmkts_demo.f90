! SPDX-License-Identifier: MIT
! Copyright (c) 2014-2022 Robert L. McDonald
program derivmkts_demo
    use derivmkts, only: dp,bscall,bsput,bscallimpvol,binomopt,binomial_result
    use derivmkts, only: geomavgprice,option_pair,mertonjump
    implicit none
    type(binomial_result) :: tree
    type(option_pair) :: asian,jump
    real(dp) :: c,p,iv
    c=bscall(40.0_dp,40.0_dp,0.30_dp,0.08_dp,0.25_dp,0.0_dp)
    p=bsput(40.0_dp,40.0_dp,0.30_dp,0.08_dp,0.25_dp,0.0_dp)
    iv=bscallimpvol(40.0_dp,40.0_dp,0.08_dp,0.25_dp,0.0_dp,c)
    asian=geomavgprice(40.0_dp,40.0_dp,0.30_dp,0.08_dp,0.25_dp,0.0_dp,3)
    jump=mertonjump(40.0_dp,40.0_dp,0.30_dp,0.08_dp,2.0_dp,0.05_dp,0.75_dp,-0.05_dp,0.35_dp)
    tree=binomopt(40.0_dp,40.0_dp,0.30_dp,0.08_dp,0.25_dp,0.0_dp,200,.true.,.true.)
    print '(a,f12.6)', 'Black-Scholes call: ',c
    print '(a,f12.6)', 'Black-Scholes put:  ',p
    print '(a,f12.6)', 'Recovered vol:     ',iv
    print '(a,f12.6)', 'Geometric Asian:   ',asian%call
    print '(a,f12.6)', 'Merton jump call:  ',jump%call
    print '(a,f12.6)', 'American put tree: ',tree%price
end program derivmkts_demo
