program test_grid_payoff
   use multi_asset_options
   implicit none

   type(status_type) :: status
   type(grid_set) :: grid
   real(dp), allocatable :: x(:), values(:)
   real(dp), parameter :: tol = 1.0e-12_dp

   call node_spacer(50.0_dp,0.0_dp,100.0_dp,5,0.0_dp,0,x,status)
   call require(status%ok(),status_message(status))
   call require(maxval(abs(x-[0.0_dp,25.0_dp,50.0_dp,75.0_dp,100.0_dp])) < tol, &
      'uniform node spacing is incorrect')

   call node_spacer(50.0_dp,0.0_dp,150.0_dp,12,5.0_dp,2,x,status)
   call require(status%ok(),status_message(status))
   call require(all(x(2:) > x(:size(x)-1)),'nonuniform nodes are not increasing')
   call require(minval(abs(x-50.0_dp)) < 1.0e-10_dp, &
      'k_shift=2 did not place strike on a node')

   allocate(grid%asset(2),grid%dims(2),grid%strides(2))
   grid%dims = [2,2]
   grid%strides = [1,2]
   grid%n_nodes = 4
   grid%asset(1)%x = [80.0_dp,120.0_dp]
   grid%asset(2)%x = [70.0_dp,130.0_dp]

   call payoff_values(payoff_digital,[option_call,option_call], &
      [100.0_dp,100.0_dp],grid,values,status)
   call require(status%ok(),status_message(status))
   call require(maxval(abs(values-[0.0_dp,0.0_dp,0.0_dp,1.0_dp])) < tol, &
      'digital payoff is incorrect')

   call payoff_values(payoff_best_of,[option_call,option_call], &
      [100.0_dp,100.0_dp],grid,values,status)
   call require(status%ok(),status_message(status))
   call require(maxval(abs(values-[0.0_dp,20.0_dp,30.0_dp,30.0_dp])) < tol, &
      'best-of payoff is incorrect')

   call payoff_values(payoff_worst_of,[option_call,option_call], &
      [100.0_dp,100.0_dp],grid,values,status)
   call require(status%ok(),status_message(status))
   call require(maxval(abs(values-[0.0_dp,0.0_dp,0.0_dp,20.0_dp])) < tol, &
      'worst-of payoff is incorrect')

   print '(a)', 'test_grid_payoff: PASS'

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

end program test_grid_payoff
