program test_tv_parameter
   use fgarch_kinds, only : dp
   use tvgarchkf
   use test_support
   implicit none
   type(tvgarch_spec) :: spec
   type(tvgarch_simulation_result) :: simulation
   type(tv_parameter_result) :: local

   spec = make_tvgarch_spec(make_tv_function([0.05_dp]), &
                            make_tv_function([0.10_dp]), &
                            make_tv_function([0.82_dp]))
   simulation = tvgarch_simulate(180,spec,seed=4421,corrected_constraints=.true.)
   local = tv_parameter(simulation%returns,shift=50,window=80,max_iterations=700)
   call assert_true(local%status == 0,'tv_parameter status')
   call assert_true(size(local%midpoint) == 3,'tv_parameter block count')
   call assert_all_close(local%midpoint,[40.0_dp,90.0_dp,140.0_dp],1.0e-13_dp,'midpoints')
   call assert_true(all(local%omega > 0.0_dp),'local omega positive')
   call assert_true(all(local%alpha >= 0.0_dp),'local alpha nonnegative')
   call assert_true(all(local%beta >= 0.0_dp),'local beta nonnegative')
   call assert_true(all(local%alpha+local%beta < 1.0_dp),'local stationarity')
   write(*,'(a)') 'test_tv_parameter: PASS'
end program test_tv_parameter
