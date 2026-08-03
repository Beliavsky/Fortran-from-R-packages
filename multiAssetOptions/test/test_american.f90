program test_american
   use multi_asset_options
   implicit none

   type(pricing_config) :: config
   type(pricing_result) :: european, american
   type(status_type) :: status
   real(dp), allocatable :: intrinsic(:)
   real(dp) :: european_value, american_value

   call initialize_config(config,1,status)
   call require(status%ok(),status_message(status))
   config%opt%pay_type = payoff_best_of
   config%opt%pc_flag = [option_put]
   config%opt%ttm = 1.0_dp
   config%opt%strike = [100.0_dp]
   config%opt%rf = 0.05_dp
   config%opt%q = [0.0_dp]
   config%opt%vol = [0.2_dp]
   config%fd%m = [120]
   config%fd%left_bound = [0.0_dp]
   config%fd%k_mult = [4.0_dp]
   config%fd%density = [5.0_dp]
   config%fd%k_shift = [2]
   config%fd%theta = 0.5_dp
   config%fd%max_smooth = 2
   config%fd%tol = 1.0e-7_dp
   config%fd%max_iter = 10
   config%time%n_steps = 240

   config%opt%exercise_type = exercise_european
   call price_multi_asset(config,european,status)
   call require(status%ok(),status_message(status))
   call interpolate_value(european%grid,european%value(:,size(european%time)), &
      [100.0_dp],european_value,status)
   call require(status%ok(),status_message(status))

   config%opt%exercise_type = exercise_american
   call price_multi_asset(config,american,status)
   call require(status%ok(),status_message(status))
   call interpolate_value(american%grid,american%value(:,size(american%time)), &
      [100.0_dp],american_value,status)
   call require(status%ok(),status_message(status))
   call payoff_values(config%opt%pay_type,config%opt%pc_flag,config%opt%strike, &
      american%grid,intrinsic,status)
   call require(status%ok(),status_message(status))

   call require(american_value >= european_value-1.0e-5_dp, &
      'American put is below the European put')
   call require(american_value > 6.0_dp .and. american_value < 6.2_dp, &
      'American put is outside its expected benchmark range')
   call require(minval(american%value(:,size(american%time))-intrinsic) > -2.0e-6_dp, &
      'American solution violates the exercise payoff floor')
   call require(maxval(american%penalty_iterations) >= 2, &
      'American penalty iteration was not exercised')

   print '(a)', 'test_american: PASS'

contains

   subroutine require(condition,message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         write(*,'(a)') 'FAIL: '//trim(message)
         error stop 1
      end if
   end subroutine require

   function status_message(status_value) result(message)
      type(status_type), intent(in) :: status_value
      character(len=:), allocatable :: message
      if (allocated(status_value%message)) then
         message = status_value%message
      else
         message = 'unexpected status failure'
      end if
   end function status_message

end program test_american
