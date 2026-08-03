program two_asset_digital
   use multi_asset_options
   implicit none

   type(pricing_config) :: config
   type(pricing_result) :: result
   type(status_type) :: status
   real(dp) :: value

   call initialize_config(config,2,status)
   if (.not. status%ok()) error stop status%message

   config%opt%pay_type = payoff_digital
   config%opt%exercise_type = exercise_european
   config%opt%pc_flag = [option_put,option_call]
   config%opt%ttm = 0.5_dp
   config%opt%strike = [110.0_dp,90.0_dp]
   config%opt%rf = 0.10_dp
   config%opt%q = [0.05_dp,0.04_dp]
   config%opt%vol = [0.20_dp,0.25_dp]
   config%opt%rho = reshape([1.0_dp,-0.5_dp,-0.5_dp,1.0_dp],[2,2])
   config%fd%m = [20,20]
   config%fd%k_mult = [0.0_dp,0.0_dp]
   config%fd%density = [5.0_dp,5.0_dp]
   config%fd%k_shift = [1,1]
   config%time%n_steps = 80

   call price_multi_asset(config,result,status)
   if (.not. status%ok()) error stop status%message
   call interpolate_value(result%grid,result%value(:,size(result%time)), &
      [100.0_dp,100.0_dp],value,status)
   if (.not. status%ok()) error stop status%message

   print '(a,f12.6)', 'Two-asset digital value: ', value
   print '(a,i0)', 'Spatial nodes: ', result%grid%n_nodes
end program two_asset_digital
