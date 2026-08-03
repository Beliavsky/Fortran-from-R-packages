program test_multiasset
   use multi_asset_options
   implicit none

   type(pricing_config) :: config
   type(pricing_result) :: result
   type(status_type) :: status
   real(dp) :: value, exact, d21, d22

   call initialize_config(config,2,status)
   call require(status%ok(),status_message(status))
   config%opt%pay_type = payoff_digital
   config%opt%exercise_type = exercise_european
   config%opt%pc_flag = [option_call,option_call]
   config%opt%ttm = 0.5_dp
   config%opt%strike = [100.0_dp,100.0_dp]
   config%opt%rf = 0.05_dp
   config%opt%q = [0.0_dp,0.0_dp]
   config%opt%vol = [0.2_dp,0.25_dp]
   config%opt%rho = reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp],[2,2])
   config%fd%m = [30,30]
   config%fd%left_bound = [0.0_dp,0.0_dp]
   config%fd%k_mult = [4.0_dp,4.0_dp]
   config%fd%density = [5.0_dp,5.0_dp]
   config%fd%k_shift = [1,1]
   config%fd%theta = 0.5_dp
   config%fd%max_smooth = 2
   config%time%n_steps = 120

   call price_multi_asset(config,result,status)
   call require(status%ok(),status_message(status))
   call interpolate_value(result%grid,result%value(:,size(result%time)), &
      [100.0_dp,100.0_dp],value,status)
   call require(status%ok(),status_message(status))

   d21 = (log(1.0_dp)+(0.05_dp-0.5_dp*0.2_dp**2)*0.5_dp) / &
      (0.2_dp*sqrt(0.5_dp))
   d22 = (log(1.0_dp)+(0.05_dp-0.5_dp*0.25_dp**2)*0.5_dp) / &
      (0.25_dp*sqrt(0.5_dp))
   exact = exp(-0.05_dp*0.5_dp)*normal_cdf(d21)*normal_cdf(d22)

   call require(value >= 0.0_dp .and. value <= exp(-0.025_dp)+1.0e-5_dp, &
      'two-asset digital price violates cash-payoff bounds')
   call require(abs(value-exact) < 0.04_dp, &
      'independent two-asset digital price is outside tolerance')
   call require(result%grid%n_nodes == 31*31, &
      'two-asset grid dimension is incorrect')

   print '(a)', 'test_multiasset: PASS'

contains

   pure real(dp) function normal_cdf(x) result(p)
      real(dp), intent(in) :: x
      p = 0.5_dp*(1.0_dp+erf(x/sqrt(2.0_dp)))
   end function normal_cdf

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

end program test_multiasset
