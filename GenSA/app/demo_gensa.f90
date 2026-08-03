program demo_gensa
   use gensa
   implicit none
   type(gensa_control) :: control
   type(gensa_result) :: result
   real(dp) :: lower(3), upper(3)

   lower = -5.12_dp
   upper = 5.12_dp
   control%maxit = 1800
   control%max_call = 400000
   control%has_threshold = .true.
   control%threshold_stop = 1.0e-8_dp
   control%seed = -100377_i8

   call gensa_minimize(rastrigin, lower, upper, result, control)

   write(*, '(a)') 'GenSA modern Fortran demonstration'
   write(*, '(a,es14.6)') 'best value:     ', result%value
   write(*, '(a,*(f11.6,1x))') 'best parameters:', result%par
   write(*, '(a,i0)') 'function calls: ', result%counts
   write(*, '(a,i0)') 'iterations:     ', result%iterations
   write(*, '(a,a)') 'termination:    ', result%message
contains
   function rastrigin(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value
      value = sum(x*x - 10.0_dp*cos(2.0_dp*acos(-1.0_dp)*x)) + 10.0_dp*real(size(x), dp)
   end function rastrigin
end program demo_gensa
