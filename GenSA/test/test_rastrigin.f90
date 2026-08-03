program test_rastrigin
   use gensa
   implicit none

   type(gensa_control) :: control
   type(gensa_result) :: result
   real(dp) :: lower(2), upper(2), initial(2)

   lower = -5.12_dp
   upper = 5.12_dp
   initial = [3.7_dp, -4.2_dp]
   control%maxit = 1200
   control%max_call = 250000
   control%has_threshold = .true.
   control%threshold_stop = 1.0e-8_dp
   control%seed = -12345_i8

   call gensa_minimize(rastrigin, lower, upper, result, control, initial)

   call assert_true(result%value < 1.0e-8_dp, 'Rastrigin global value')
   call assert_true(maxval(abs(result%par)) < 2.0e-4_dp, 'Rastrigin global minimizer')
   call assert_true(result%iterations < control%maxit, 'early threshold stopping')

   print '(a)', 'PASS test_rastrigin'

contains

   function rastrigin(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value
      value = sum(x*x - 10.0_dp*cos(2.0_dp*acos(-1.0_dp)*x)) + 10.0_dp*real(size(x), dp)
   end function rastrigin

   subroutine assert_true(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*, '(a)') 'FAIL: ' // label
         error stop 1
      end if
   end subroutine assert_true

end program test_rastrigin
