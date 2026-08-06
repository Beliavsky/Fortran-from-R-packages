program test_fit
   use fgarch_kinds, only : dp
   use tvgarchkf
   use test_support
   implicit none
   type(tvgarch_spec) :: truth, initial
   type(tvgarch_simulation_result) :: simulation
   type(tvgarch_fit_result) :: fit
   type(tvgarch_filter_result) :: initial_filter

   truth = make_tvgarch_spec(make_tv_function([0.08_dp]), &
                             make_tv_function([0.12_dp]), &
                             make_tv_function([0.72_dp]))
   simulation = tvgarch_simulate(260,truth,seed=917,corrected_constraints=.true.)
   call assert_true(simulation%status == 0,'fit simulation')
   initial = make_tvgarch_spec(make_tv_function([0.16_dp]), &
                               make_tv_function([0.20_dp]), &
                               make_tv_function([0.55_dp]))
   initial_filter = tvgarch_kalman_filter(simulation%returns,initial,corrected_constraints=.true.)
   fit = tvgarch_kalman_fit(simulation%returns,initial,corrected_constraints=.true., &
                            max_iterations=1200,tolerance=2.0e-8_dp)
   call assert_true(fit%filter%status == 0,'fit produces valid filter')
   call assert_true(fit%criterion <= initial_filter%criterion+1.0e-7_dp,'fit improves criterion')
   call assert_true(fit%spec%omega%coefficients(1) > 0.0_dp,'fit omega positive')
   call assert_true(fit%spec%alpha%coefficients(1) >= 0.0_dp,'fit alpha nonnegative')
   call assert_true(fit%spec%beta%coefficients(1) >= 0.0_dp,'fit beta nonnegative')
   call assert_true(fit%spec%alpha%coefficients(1)+fit%spec%beta%coefficients(1) < 1.0_dp, &
                    'fit stationarity')
   call assert_true(size(fit%rounded_parameters) == 3,'rounded compatibility output')
   write(*,'(a,3f12.6)') 'test_fit: PASS estimates =',fit%parameters
end program test_fit
