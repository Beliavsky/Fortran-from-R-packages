program demo_multi_asset_options
   use multi_asset_options
   implicit none

   type(pricing_config) :: config
   type(pricing_result) :: result
   type(status_type) :: status
   real(dp) :: value

   call initialize_config(config,2,status)
   if (.not. status%ok()) error stop status%message

   config%opt%pay_type = payoff_best_of
   config%opt%exercise_type = exercise_european
   config%opt%pc_flag = [option_call,option_call]
   config%opt%ttm = 1.0_dp
   config%opt%strike = [100.0_dp,100.0_dp]
   config%opt%rf = 0.04_dp
   config%opt%q = [0.01_dp,0.02_dp]
   config%opt%vol = [0.20_dp,0.30_dp]
   config%opt%rho = reshape([1.0_dp,0.35_dp,0.35_dp,1.0_dp],[2,2])
   config%fd%m = [28,28]
   config%fd%k_mult = [4.0_dp,4.0_dp]
   config%fd%density = [5.0_dp,5.0_dp]
   config%fd%k_shift = [1,1]
   config%fd%theta = 0.5_dp
   config%fd%max_smooth = 2
   config%time%n_steps = 112

   call price_multi_asset(config,result,status)
   if (.not. status%ok()) error stop status%message
   call interpolate_value(result%grid,result%value(:,size(result%time)), &
      [100.0_dp,100.0_dp],value,status)
   if (.not. status%ok()) error stop status%message

   print '(a)', 'multiAssetOptions-fortran demonstration'
   print '(a,i0)', 'Assets: ', config%opt%n_asset
   print '(a,i0)', 'Grid nodes: ', result%grid%n_nodes
   print '(a,i0)', 'Stored time levels: ', size(result%time)
   print '(a,f12.6)', 'Best-of call value at (100,100): ', value
end program demo_multi_asset_options
