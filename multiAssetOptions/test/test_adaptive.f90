program test_adaptive
   use multi_asset_options
   implicit none

   type(pricing_config) :: config
   type(pricing_result) :: result
   type(status_type) :: status
   real(dp) :: value, exact, d1, d2

   call initialize_config(config,1,status)
   call require(status%ok(),status_message(status))
   config%opt%pay_type = payoff_best_of
   config%opt%exercise_type = exercise_european
   config%opt%pc_flag = [option_call]
   config%opt%ttm = 1.0_dp
   config%opt%strike = [100.0_dp]
   config%opt%rf = 0.05_dp
   config%opt%q = [0.0_dp]
   config%opt%vol = [0.2_dp]
   config%fd%m = [100]
   config%fd%k_mult = [4.0_dp]
   config%fd%k_shift = [2]
   config%time%ts_type = timestep_adaptive
   config%time%dt_init = 0.01_dp
   config%time%d_norm = 0.08_dp
   config%time%scale_d = 0.05_dp

   call price_multi_asset(config,result,status)
   call require(status%ok(),status_message(status))
   call require(abs(result%time(size(result%time))) < 1.0e-12_dp, &
      'adaptive stepping did not reach valuation time zero')
   call require(all(result%time(2:) < result%time(:size(result%time)-1)), &
      'adaptive time history is not strictly decreasing')
   call require(size(result%time) > 2 .and. size(result%time) < 10000, &
      'adaptive timestep count is unreasonable')

   call interpolate_value(result%grid,result%value(:,size(result%time)), &
      [100.0_dp],value,status)
   call require(status%ok(),status_message(status))
   d1 = (0.05_dp+0.5_dp*0.2_dp**2)/0.2_dp
   d2 = d1-0.2_dp
   exact = 100.0_dp*normal_cdf(d1) - &
      100.0_dp*exp(-0.05_dp)*normal_cdf(d2)
   call require(abs(value-exact) < 0.08_dp, &
      'adaptive European call is outside tolerance')

   print '(a)', 'test_adaptive: PASS'

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

end program test_adaptive
