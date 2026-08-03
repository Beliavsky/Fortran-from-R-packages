program trace_example
   use gensa
   implicit none
   type(gensa_control) :: control
   type(gensa_result) :: result
   real(dp) :: lower(2), upper(2)
   integer :: i

   lower = -4.0_dp
   upper = 4.0_dp
   control%maxit = 40
   control%local_search = .false.
   control%seed = -55_i8

   call gensa_minimize(objective, lower, upper, result, control)
   write(*, '(a)') ' step    temperature       current          best'
   do i = 1, result%trace%n
      write(*, '(i5,3(1x,es15.7))') result%trace%step(i), result%trace%temperature(i), &
         result%trace%current_value(i), result%trace%best_value(i)
   end do
contains
   function objective(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value
      value = sin(3.0_dp*x(1))*sin(4.0_dp*x(2)) + 0.05_dp*sum(x*x)
   end function objective
end program trace_example
