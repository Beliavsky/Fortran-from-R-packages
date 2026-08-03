program test_european
   use multi_asset_options
   implicit none

   type(pricing_config) :: config
   type(pricing_result) :: call_result, put_result
   type(status_type) :: status
   real(dp) :: call_value, put_value, call_exact, put_exact

   call initialize_config(config,1,status)
   call require(status%ok(),status_message(status))
   call set_common_inputs(config)

   config%opt%pc_flag = [option_call]
   call price_multi_asset(config,call_result,status)
   call require(status%ok(),status_message(status))
   call interpolate_value(call_result%grid, &
      call_result%value(:,size(call_result%time)),[100.0_dp],call_value,status)
   call require(status%ok(),status_message(status))

   config%opt%pc_flag = [option_put]
   call price_multi_asset(config,put_result,status)
   call require(status%ok(),status_message(status))
   call interpolate_value(put_result%grid, &
      put_result%value(:,size(put_result%time)),[100.0_dp],put_value,status)
   call require(status%ok(),status_message(status))

   call_exact = black_scholes(.true.,100.0_dp,100.0_dp,0.05_dp,0.0_dp,0.2_dp,1.0_dp)
   put_exact = black_scholes(.false.,100.0_dp,100.0_dp,0.05_dp,0.0_dp,0.2_dp,1.0_dp)
   call require(abs(call_value-call_exact) < 0.03_dp, &
      'European call is outside the finite-difference tolerance')
   call require(abs(put_value-put_exact) < 0.03_dp, &
      'European put is outside the finite-difference tolerance')
   call require(abs((call_value-put_value) - &
      (100.0_dp-100.0_dp*exp(-0.05_dp))) < 0.04_dp, &
      'put-call parity is not satisfied')

   print '(a)', 'test_european: PASS'

contains

   subroutine set_common_inputs(c)
      type(pricing_config), intent(inout) :: c
      c%opt%pay_type = payoff_best_of
      c%opt%exercise_type = exercise_european
      c%opt%ttm = 1.0_dp
      c%opt%strike = [100.0_dp]
      c%opt%rf = 0.05_dp
      c%opt%q = [0.0_dp]
      c%opt%vol = [0.2_dp]
      c%fd%m = [120]
      c%fd%left_bound = [0.0_dp]
      c%fd%k_mult = [4.0_dp]
      c%fd%density = [5.0_dp]
      c%fd%k_shift = [2]
      c%fd%theta = 0.5_dp
      c%fd%max_smooth = 2
      c%time%ts_type = timestep_constant
      c%time%n_steps = 240
   end subroutine set_common_inputs

   real(dp) function black_scholes(is_call,s,k,r,q,sigma,t) result(value)
      logical, intent(in) :: is_call
      real(dp), intent(in) :: s,k,r,q,sigma,t
      real(dp) :: d1,d2
      d1 = (log(s/k)+(r-q+0.5_dp*sigma**2)*t)/(sigma*sqrt(t))
      d2 = d1-sigma*sqrt(t)
      if (is_call) then
         value = s*exp(-q*t)*normal_cdf(d1)-k*exp(-r*t)*normal_cdf(d2)
      else
         value = k*exp(-r*t)*normal_cdf(-d2)-s*exp(-q*t)*normal_cdf(-d1)
      end if
   end function black_scholes

   pure real(dp) function normal_cdf(x) result(value)
      real(dp), intent(in) :: x
      value = 0.5_dp*(1.0_dp+erf(x/sqrt(2.0_dp)))
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

end program test_european
