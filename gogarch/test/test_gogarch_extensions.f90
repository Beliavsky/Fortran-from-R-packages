! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of gogarch, copyright (C) 2008-2026 Bernhard Pfaff.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
program test_gogarch_extensions
   use gogarch
   use gogarch_linalg, only : is_symmetric
   use test_helpers
   implicit none
   integer, parameter :: n = 180, m = 2
   real(dp) :: factors(n,m), factor_h(n,m), data(n,m), mixing(m,m)
   real(dp) :: mean_forecast(2,m), covariance_forecast(m,m,2)
   real(dp) :: means(m), omegas(m), arch21(m,2), lev0(m,0), garch21(m,1)
   real(dp) :: delta(m), shape(m), skew(m)
   real(dp) :: arch11(m,1), lev11(m,1), garch11(m,1)
   type(univariate_spec) :: high_spec, aparch_spec
   type(gogarch_fit) :: fits(4), aparch_fit
   integer :: i, t

   call seed_rng(778899)
   call simulate_garchpq(n,0.0_dp,0.035_dp,[0.07_dp,0.025_dp],[0.86_dp],'std',7.0_dp,1.0_dp, &
      factors(:,1),factor_h(:,1),burnin=500)
   call simulate_garchpq(n,0.0_dp,0.050_dp,[0.09_dp,0.020_dp],[0.82_dp],'std',9.0_dp,1.0_dp, &
      factors(:,2),factor_h(:,2),burnin=500)
   mixing = reshape([1.0_dp,-0.2_dp,0.4_dp,0.9_dp],[m,m])
   data = matmul(factors,transpose(mixing))

   high_spec%model = 'garch'
   high_spec%distribution = 'std'
   high_spec%p = 2
   high_spec%o = 0
   high_spec%q = 1
   high_spec%shape = 8.0_dp
   high_spec%fit_shape = .false.
   fits(1) = fit_gogarch_ica(data,max_ica_iterations=300,max_garch_iterations=180,factor_spec=high_spec)
   fits(2) = fit_gogarch_mm(data,lag_max=2,max_garch_iterations=180,factor_spec=high_spec)
   fits(3) = fit_gogarch_nls(data,max_outer_iterations=70,max_garch_iterations=150,factor_spec=high_spec)
   fits(4) = fit_gogarch_ml(data,max_outer_iterations=20,max_garch_iterations=95,factor_spec=high_spec)
   do i = 1, 4
      call assert_true(fits(i)%status <= 1,'higher-order Student GO-GARCH estimator')
      call assert_true(fits(i)%factor_models(1)%p == 2 .and. fits(i)%factor_models(1)%q == 1, &
         'GO-GARCH retains requested orders')
      call assert_true(trim(fits(i)%factor_models(1)%distribution) == 'std','GO-GARCH retains distribution')
      call forecast_gogarch(fits(i),2,mean_forecast,covariance_forecast)
      do t = 1, 2
         call assert_true(is_symmetric(covariance_forecast(:,:,t),1.0e-10_dp),'extended forecast covariance symmetry')
      end do
   end do
   call factor_coefficients_full(fits(1),means,omegas,arch21,lev0,garch21,delta,shape,skew)
   call assert_true(all(omegas > 0.0_dp) .and. all(arch21 >= 0.0_dp) .and. all(garch21 >= 0.0_dp), &
      'full higher-order coefficient extraction')
   call assert_true(all(abs(shape-8.0_dp) < 1.0e-12_dp),'fixed Student shape retained')

   call simulate_aparch(n,0.0_dp,0.04_dp,[0.08_dp],[0.18_dp],[0.84_dp],1.4_dp,'sstd',8.0_dp,1.25_dp, &
      factors(:,1),factor_h(:,1),burnin=600)
   call simulate_aparch(n,0.0_dp,0.05_dp,[0.10_dp],[-0.12_dp],[0.80_dp],1.6_dp,'sstd',10.0_dp,0.85_dp, &
      factors(:,2),factor_h(:,2),burnin=600)
   data = matmul(factors,transpose(mixing))
   aparch_spec%model = 'aparch'
   aparch_spec%distribution = 'sstd'
   aparch_spec%p = 1
   aparch_spec%o = 1
   aparch_spec%q = 1
   aparch_spec%delta = 1.5_dp
   aparch_spec%shape = 9.0_dp
   aparch_spec%skew = 1.1_dp
   aparch_spec%fit_delta = .true.
   aparch_spec%fit_shape = .false.
   aparch_spec%fit_skew = .false.
   aparch_fit = fit_gogarch_ica(data,max_ica_iterations=300,max_garch_iterations=300,factor_spec=aparch_spec)
   call assert_true(aparch_fit%status <= 1,'APARCH skew-Student GO-GARCH fit')
   call assert_true(aparch_fit%factor_models(1)%o == 1,'APARCH leverage order retained')
   call assert_true(aparch_fit%factor_models(1)%delta > 0.25_dp .and. aparch_fit%factor_models(1)%delta < 4.0_dp, &
      'APARCH delta estimated in GO-GARCH')
   call factor_coefficients_full(aparch_fit,means,omegas,arch11,lev11,garch11,delta,shape,skew)
   call assert_true(all(abs(lev11) < 1.0_dp),'APARCH leverage extraction')
   call assert_true(all(abs(shape-9.0_dp) < 1.0e-12_dp) .and. all(abs(skew-1.1_dp) < 1.0e-12_dp), &
      'fixed skew-Student parameters retained')

   write(*,'(a)') 'Extended GO-GARCH specification tests passed.'
end program test_gogarch_extensions
