program test_all
   use fgarch_kinds, only : dp
   use tvgarchkf
   implicit none
   type(tvgarch_spec) :: spec
   type(tvgarch_filter_result) :: filtered
   type(tvgarch_simulation_result) :: simulated
   real(dp), allocatable :: p(:)

   p = polynomial_values(2,[1.0_dp,2.0_dp])
   if (maxval(abs(p-[2.0_dp,3.0_dp])) > 1.0e-12_dp) error stop 'polynomial test failed'

   spec = make_tvgarch_spec(make_tv_function([0.1_dp]),make_tv_function([0.2_dp]),make_tv_function([0.6_dp]))
   filtered = tvgarch_kalman_filter([0.0_dp,1.0_dp,-0.5_dp],spec,corrected_constraints=.true.)
   if (filtered%status /= 0) error stop 'filter test failed'
   if (abs(filtered%criterion-0.831790649644333_dp) > 1.0e-10_dp) error stop 'criterion test failed'

   simulated = tvgarch_simulate(10,spec,seed=123,corrected_constraints=.true.)
   if (simulated%status /= 0 .or. any(simulated%variance < 0.0_dp)) error stop 'simulation test failed'
   write(*,'(a)') 'tvgarchkf fpm test: PASS'
end program test_all
