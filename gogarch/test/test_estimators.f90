! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of gogarch, copyright (C) 2008-2026 Bernhard Pfaff.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
program test_estimators
   use gogarch
   use gogarch_linalg, only : is_orthogonal, is_symmetric, covariance_matrix, identity_matrix, outer_product
   use test_helpers
   implicit none
   integer, parameter :: n = 160, m = 2
   real(dp) :: factors(n,m), factor_h(n,m), data(n,m), mixing(m,m)
   real(dp) :: mean_forecast(3,m), covariance_forecast(m,m,3), factor_forecast(3,m)
   real(dp) :: residuals(n,m), simulated(25,m), simulated_factors(25,m), simulated_h(25,m)
   real(dp) :: variances(n,m), correlations(m,m,n), coefficients(m,4)
   real(dp) :: factor_cov(m,m), whitened(n,m), ssi(m,m,n), initial_nls(3)
   real(dp) :: initial_objective, final_objective, angle_negll
   type(gogarch_fit) :: fits(4), generic_fit, angle_fit
   character(len=3) :: names(4)
   logical :: ok
   integer :: i, t

   call seed_rng(13579)
   call simulate_garch11(n,0.0_dp,0.05_dp,0.08_dp,0.87_dp,factors(:,1),factor_h(:,1),burnin=400)
   call simulate_garch11(n,0.0_dp,0.08_dp,0.12_dp,0.80_dp,factors(:,2),factor_h(:,2),burnin=400)
   mixing = reshape([1.10_dp,-0.25_dp,0.45_dp,0.85_dp],[m,m])
   data = matmul(factors,transpose(mixing))

   fits(1) = fit_gogarch_ica(data,max_ica_iterations=350,max_garch_iterations=120)
   fits(2) = fit_gogarch_mm(data,lag_max=2,max_garch_iterations=120)
   fits(3) = fit_gogarch_nls(data,max_outer_iterations=90,max_garch_iterations=120)
   fits(4) = fit_gogarch_ml(data,max_outer_iterations=28,max_garch_iterations=70)
   names = ['ica','mm ','nls','ml ']

   do i = 1, 4
      call assert_true(fits(i)%status <= 1,trim(names(i))//' fit status')
      call assert_true(is_orthogonal(fits(i)%rotation,2.0e-8_dp),trim(names(i))//' orthogonal rotation')
      call assert_true(reconstruction_error(fits(i)) < 5.0e-9_dp,trim(names(i))//' exact data reconstruction')
      call assert_true(all(fits(i)%factor_variance > 0.0_dp),trim(names(i))//' positive factor variance')
      call assert_all_finite(fits(i)%covariance,trim(names(i))//' finite covariance path')
      do t = 1, n
         call assert_true(is_symmetric(fits(i)%covariance(:,:,t),1.0e-10_dp),trim(names(i))//' covariance symmetry')
      end do
      do t = 1, n
         call assert_true(fits(i)%covariance(1,1,t) > 0.0_dp .and. fits(i)%covariance(2,2,t) > 0.0_dp, &
            trim(names(i))//' positive covariance diagonal')
      end do
      call forecast_gogarch(fits(i),3,mean_forecast,covariance_forecast,factor_forecast)
      call assert_all_finite(mean_forecast,trim(names(i))//' finite mean forecast')
      call assert_true(all(factor_forecast > 0.0_dp),trim(names(i))//' positive factor forecasts')
      do t = 1, 3
         call assert_true(is_symmetric(covariance_forecast(:,:,t),1.0e-10_dp),trim(names(i))//' forecast symmetry')
      end do
      call standardized_residuals(fits(i),residuals,ok)
      call assert_true(ok,trim(names(i))//' standardized residual calculation')
      call assert_all_finite(residuals,trim(names(i))//' finite standardized residuals')
      call factor_coefficients(fits(i),coefficients)
      call assert_true(all(coefficients(:,2) > 0.0_dp),trim(names(i))//' positive fitted omega')
   end do

   call conditional_variances(fits(1)%covariance,variances)
   call conditional_correlations(fits(1)%covariance,correlations)
   call assert_true(all(variances > 0.0_dp),'conditional variances')
   do t = 1, n
      call assert_close(correlations(1,1,t),1.0_dp,1.0e-12_dp,'correlation diagonal')
      call assert_close(correlations(2,2,t),1.0_dp,1.0e-12_dp,'correlation diagonal')
   end do

   factor_cov = covariance_matrix(fits(1)%factors,center=.false.)
   call assert_true(maxval(abs(factor_cov-identity_matrix(m))) < 1.0e-8_dp,'ICA factors retain whitening covariance')
   call assert_true(allocated(fits(2)%mm_weights),'MM weights returned')
   call assert_close(sum(fits(2)%mm_weights),1.0_dp,1.0e-12_dp,'MM weights sum to one')

   whitened = matmul(data,fits(3)%covariance_invsqrt)
   do t = 1, n
      ssi(:,:,t) = outer_product(whitened(t,:),whitened(t,:))-identity_matrix(m)
   end do
   initial_nls = 0.1_dp
   initial_objective = gonls_objective(initial_nls,ssi)
   final_objective = gonls_objective(fits(3)%objective_parameters,ssi)
   call assert_true(final_objective <= initial_objective*(1.0_dp+1.0e-8_dp),'NLS objective must not worsen')
   call assert_true(all(fits(4)%objective_parameters > 0.0_dp) .and. &
      all(fits(4)%objective_parameters < 0.5_dp*acos(-1.0_dp)),'ML Euler-angle bounds')

   angle_fit = gogarch_from_angles(data,[0.4_dp],max_garch_iterations=90)
   angle_negll = gogarch_negloglik(data,[0.4_dp],max_garch_iterations=90)
   call assert_true(angle_fit%status <= 1,'direct angle model construction')
   call assert_close(angle_negll,-angle_fit%log_likelihood,1.0e-10_dp,'direct angle likelihood wrapper')

   call simulate_fitted_gogarch(fits(1),25,simulated,simulated_factors,simulated_h)
   call assert_all_finite(simulated,'fitted GO-GARCH simulation')
   call assert_true(all(simulated_h > 0.0_dp),'simulated factor variances')

   generic_fit = fit_gogarch(data,'mm',lag_max=1,max_garch_iterations=80)
   call assert_true(generic_fit%status <= 1 .and. trim(generic_fit%method) == 'mm','generic estimator dispatch')

   write(*,'(a)') 'GO-GARCH estimator and workflow tests passed.'
end program test_estimators
