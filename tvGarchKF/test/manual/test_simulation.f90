program test_simulation
   use fgarch_kinds, only : dp
   use tvgarchkf
   use test_support
   implicit none
   type(tvgarch_spec) :: spec
   type(tvgarch_simulation_result) :: simulation, simulation2
   real(dp) :: z(4), expected_v(4), expected_y(4)

   spec = make_tvgarch_spec(make_tv_function([0.1_dp]), &
                            make_tv_function([0.2_dp]), &
                            make_tv_function([0.6_dp]))
   z = [0.0_dp,1.0_dp,-1.0_dp,0.5_dp]
   simulation = tvgarch_simulate(4,spec,innovations=z,corrected_constraints=.true.)
   call assert_true(simulation%status == 0,'simulation status')
   expected_v = [0.0_dp,0.1_dp,0.18_dp,0.244_dp]
   expected_y = [0.0_dp,sqrt(0.1_dp),-sqrt(0.18_dp),0.5_dp*sqrt(0.244_dp)]
   call assert_all_close(simulation%variance,expected_v,1.0e-13_dp,'simulation variance')
   call assert_all_close(simulation%returns,expected_y,1.0e-13_dp,'simulation returns')

   simulation = tvgarch_simulate(30,spec,seed=1234,corrected_constraints=.true.)
   simulation2 = tvgarch_simulate(30,spec,seed=1234,corrected_constraints=.true.)
   call assert_all_close(simulation%returns,simulation2%returns,0.0_dp,'seeded reproducibility')
   write(*,'(a)') 'test_simulation: PASS'
end program test_simulation
