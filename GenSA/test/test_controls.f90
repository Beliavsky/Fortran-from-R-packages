program test_controls
   use gensa
   implicit none

   type(gensa_control) :: control
   type(gensa_result) :: result
   real(dp) :: lower(2), upper(2), initial(2)
   integer :: i

   lower = -1.0_dp
   upper = 1.0_dp
   initial = 0.0_dp
   control%maxit = 1000
   control%max_call = 25
   control%local_search = .false.
   control%trace = .true.
   control%seed = -123_i8

   call gensa_minimize(flat_objective, lower, upper, result, control, initial)

   call assert_true(result%status == gensa_max_calls, 'max-call status')
   call assert_true(result%counts == control%max_call, 'exact maximum call count')
   call assert_true(result%trace%n >= 2, 'trace collected')
   do i = 2, result%trace%n
      call assert_true(result%trace%best_value(i) <= result%trace%best_value(i-1), 'monotone best trace')
      call assert_true(result%trace%step(i) > result%trace%step(i-1), 'increasing trace steps')
   end do

   print '(a)', 'PASS test_controls'

contains

   function flat_objective(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value
      value = 1.0_dp + 0.0_dp*sum(x)
   end function flat_objective

   subroutine assert_true(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*, '(a)') 'FAIL: ' // label
         error stop 1
      end if
   end subroutine assert_true

end program test_controls
