program nonsmooth_example
   use gensa
   implicit none
   type(gensa_control) :: control
   type(gensa_result) :: result
   real(dp) :: lower(2), upper(2)

   lower = -10.0_dp
   upper = 10.0_dp
   control%smooth = .false.
   control%maxit = 700
   control%seed = -2026_i8

   call gensa_minimize(objective, lower, upper, result, control)
   write(*, '(a,es14.6)') 'minimum value: ', result%value
   write(*, '(a,*(f11.6,1x))') 'parameters:    ', result%par
contains
   function objective(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value
      value = max(abs(x(1) - 2.0_dp), 3.0_dp*abs(x(2) + 1.0_dp))
   end function objective
end program nonsmooth_example
