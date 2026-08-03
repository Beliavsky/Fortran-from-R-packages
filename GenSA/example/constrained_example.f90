program constrained_example
   use gensa
   implicit none
   type(gensa_control) :: control
   type(gensa_result) :: result
   real(dp) :: lower(2), upper(2), initial(2)

   lower = -2.0_dp
   upper = 2.0_dp
   initial = [-1.0_dp, 1.0_dp]
   control%smooth = .false.
   control%maxit = 1000

   call gensa_minimize(objective, lower, upper, result, control, initial, feasible)
   write(*, '(a,es14.6)') 'minimum value: ', result%value
   write(*, '(a,*(f11.6,1x))') 'parameters:    ', result%par
   write(*, '(a,l1)') 'feasible:      ', feasible(result%par)
contains
   function objective(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value
      value = (x(1) - 1.0_dp)**2 + (x(2) - 2.0_dp)**2
   end function objective
   logical function feasible(x)
      real(dp), intent(in) :: x(:)
      feasible = sum(x*x) <= 1.0_dp
   end function feasible
end program constrained_example
