program test_quadratic
   use gensa
   implicit none

   type(gensa_control) :: control
   type(gensa_result) :: result
   real(dp) :: lower(3), upper(3), initial(3)

   lower = -10.0_dp
   upper = 10.0_dp
   initial = [8.0_dp, -7.0_dp, 5.0_dp]
   control%maxit = 200
   control%max_call = 20000
   control%has_threshold = .true.
   control%threshold_stop = 1.0e-12_dp
   control%seed = -4321_i8

   call gensa_minimize(quadratic, lower, upper, result, control, initial)

   call assert_true(result%value < 1.0e-10_dp, 'quadratic objective')
   call assert_true(maxval(abs(result%par - [1.0_dp, -2.0_dp, 0.5_dp])) < 2.0e-5_dp, 'quadratic minimizer')
   call assert_true(result%status == gensa_threshold_reached, 'threshold status')
   call assert_true(result%counts > 0, 'positive function count')

   print '(a)', 'PASS test_quadratic'

contains

   function quadratic(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value
      value = (x(1) - 1.0_dp)**2 + 2.0_dp * (x(2) + 2.0_dp)**2 + 0.5_dp * (x(3) - 0.5_dp)**2
   end function quadratic

   subroutine assert_true(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*, '(a)') 'FAIL: ' // label
         error stop 1
      end if
   end subroutine assert_true

end program test_quadratic
