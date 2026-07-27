! SPDX-License-Identifier: MIT
! Copyright (c) 2014-2022 Robert L. McDonald
program simulation_and_tree
    use derivmkts, only: dp,simprice,simulation_result,binomopt,binomial_result
    implicit none
    type(simulation_result) :: paths
    type(binomial_result) :: tree
    paths=simprice(100.0_dp,0.20_dp,0.05_dp,1.0_dp,0.01_dp,5,12,.true.,0.5_dp,-0.10_dp,0.25_dp,77)
    tree=binomopt(100.0_dp,100.0_dp,0.20_dp,0.05_dp,1.0_dp,0.01_dp,250,.true.,.true.)
    print '(a,5f11.3)', 'Terminal prices: ',paths%price(:,13,1)
    print '(a,f12.6)', 'American put:    ',tree%price
    print '(a,f12.6)', 'Initial delta:   ',tree%delta
end program simulation_and_tree
