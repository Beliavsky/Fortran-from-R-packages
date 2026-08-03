program test_nonsmooth
   use gensa
   implicit none

   type(gensa_control) :: control
   type(gensa_result) :: result
   real(dp) :: lower(2), upper(2), initial(2)

   lower = -5.0_dp
   upper = 5.0_dp
   initial = [4.0_dp, -4.0_dp]
   control%smooth = .false.
   control%maxit = 500
   control%max_call = 100000
   control%has_threshold = .true.
   control%threshold_stop = 1.0e-6_dp
   control%seed = -77_i8

   call gensa_minimize(nonsmooth, lower, upper, result, control, initial)

   call assert_true(result%value < 2.0e-5_dp, 'nonsmooth objective')
   call assert_true(maxval(abs(result%par - [0.25_dp, -1.5_dp])) < 2.0e-4_dp, 'nonsmooth minimizer')

   print '(a)', 'PASS test_nonsmooth'

contains

   function nonsmooth(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value
      value = abs(x(1) - 0.25_dp) + 2.0_dp*abs(x(2) + 1.5_dp)
   end function nonsmooth

   subroutine assert_true(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*, '(a)') 'FAIL: ' // label
         error stop 1
      end if
   end subroutine assert_true

end program test_nonsmooth
