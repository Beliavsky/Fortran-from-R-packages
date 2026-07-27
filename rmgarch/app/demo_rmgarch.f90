! SPDX-License-Identifier: GPL-3.0-only
!
! Experimental modern Fortran translation of computational methods from the
! R package rmgarch, copyright (C) 2008-2025 Alexios Galanos.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 3 only.
program demo_rmgarch
   use rmgarch
   implicit none

   integer, parameter :: n = 500, m = 3
   real(dp) :: qbar(m,m), z(n,m), q(m,m,n), r(m,m,n)
   real(dp) :: qforecast(m,m,5), rforecast(m,m,5), copula_ll
   type(dcc_spec) :: true_spec
   type(dcc_fit_result) :: fit
   real(dp), allocatable :: distances(:,:)
   logical :: ok

   call seed_rng(20260723)
   qbar = reshape([1.0_dp,0.45_dp,0.20_dp, &
                   0.45_dp,1.0_dp,0.35_dp, &
                   0.20_dp,0.35_dp,1.0_dp],[m,m])
   true_spec = make_dcc_spec([0.045_dp],[0.93_dp])
   call simulate_dcc(n,true_spec,qbar,z,q,r,burn=300)
   fit = fit_dcc11(z,max_iterations=500)
   call dcc_forecast(fit%spec,fit%qbar,fit%nbar,fit%q(:,:,n),z(n,:),5,qforecast,rforecast)
   copula_ll = dynamic_gaussian_copula_log_likelihood(z,fit%spec,ok)
   distances = correlation_distance_matrix(fit%r,stride=100,method='ma')

   print '(a,i0)', 'observations: ', n
   print '(a,i0)', 'assets: ', m
   print '(a,f10.6)', 'true DCC alpha: ', true_spec%alpha(1)
   print '(a,f10.6)', 'fitted DCC alpha: ', fit%spec%alpha(1)
   print '(a,f10.6)', 'true DCC beta:  ', true_spec%beta(1)
   print '(a,f10.6)', 'fitted DCC beta:  ', fit%spec%beta(1)
   print '(a,f14.4)', 'DCC log likelihood: ', fit%log_likelihood
   print '(a,f14.4)', 'Gaussian copula log likelihood: ', copula_ll
   print '(a,f10.6)', 'last correlation (1,2): ', fit%r(1,2,n)
   print '(a,f10.6)', 'one-step forecast correlation (1,2): ', rforecast(1,2,1)
   print '(a,i0,a,i0)', 'correlation-distance matrix dimensions: ', size(distances,1), ' x ', size(distances,2)
   print '(a,i0)', 'fit status: ', fit%status
end program demo_rmgarch
