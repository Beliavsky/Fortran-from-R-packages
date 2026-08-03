program rastrigin_example
   use gensa
   implicit none
   type(gensa_control) :: control
   type(gensa_result) :: result
   real(dp) :: lower(5), upper(5)

   lower = -5.12_dp
   upper = 5.12_dp
   control%maxit = 2500
   control%max_call = 500000
   control%has_threshold = .true.
   control%threshold_stop = 1.0e-8_dp
   control%seed = -1234_i8

   call gensa_minimize(rastrigin, lower, upper, result, control)
   write(*, '(a,es14.6)') 'minimum value: ', result%value
   write(*, '(a,*(f11.6,1x))') 'parameters:    ', result%par
   write(*, '(a,i0)') 'function calls: ', result%counts
   write(*, '(a,a)') 'termination: ', result%message
contains
   function rastrigin(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value
      value = sum(x*x - 10.0_dp*cos(2.0_dp*acos(-1.0_dp)*x)) + 10.0_dp*real(size(x), dp)
   end function rastrigin
end program rastrigin_example
