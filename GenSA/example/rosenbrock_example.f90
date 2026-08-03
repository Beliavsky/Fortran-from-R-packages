program rosenbrock_example
   use gensa
   implicit none
   type(gensa_control) :: control
   type(gensa_result) :: result
   real(dp) :: lower(4), upper(4), initial(4)

   lower = -3.0_dp
   upper = 3.0_dp
   initial = [-1.2_dp, 1.0_dp, -1.2_dp, 1.0_dp]
   control%maxit = 800
   control%has_threshold = .true.
   control%threshold_stop = 1.0e-10_dp

   call gensa_minimize(rosenbrock, lower, upper, result, control, initial)
   write(*, '(a,es14.6)') 'minimum value: ', result%value
   write(*, '(a,*(f11.6,1x))') 'parameters:    ', result%par
contains
   function rosenbrock(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value
      integer :: i
      value = 0.0_dp
      do i = 1, size(x) - 1
         value = value + 100.0_dp*(x(i+1) - x(i)**2)**2 + (1.0_dp - x(i))**2
      end do
   end function rosenbrock
end program rosenbrock_example
