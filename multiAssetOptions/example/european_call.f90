program european_call
   use multi_asset_options
   implicit none

   type(pricing_config) :: config
   type(pricing_result) :: result
   type(status_type) :: status
   real(dp) :: value

   call initialize_config(config,1,status)
   if (.not. status%ok()) error stop status%message

   config%opt%pay_type = payoff_best_of
   config%opt%exercise_type = exercise_european
   config%opt%pc_flag = [option_call]
   config%opt%ttm = 1.0_dp
   config%opt%strike = [100.0_dp]
   config%opt%rf = 0.05_dp
   config%opt%q = [0.0_dp]
   config%opt%vol = [0.2_dp]
   config%fd%m = [120]
   config%fd%k_mult = [4.0_dp]
   config%fd%k_shift = [2]
   config%time%n_steps = 240

   call price_multi_asset(config,result,status)
   if (.not. status%ok()) error stop status%message
   call interpolate_value(result%grid,result%value(:,size(result%time)), &
      [100.0_dp],value,status)
   if (.not. status%ok()) error stop status%message

   print '(a,f12.6)', 'European call value at S=100: ', value
end program european_call
