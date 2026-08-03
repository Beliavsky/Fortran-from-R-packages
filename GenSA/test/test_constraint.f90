program test_constraint
   use gensa
   implicit none

   type(gensa_control) :: control
   type(gensa_result) :: result
   real(dp) :: lower(2), upper(2), initial(2)

   lower = -2.0_dp
   upper = 2.0_dp
   initial = [-1.0_dp, 1.0_dp]
   control%smooth = .false.
   control%maxit = 1200
   control%max_call = 200000
   control%seed = -909_i8
   control%simple_function = .true.

   call gensa_minimize(objective, lower, upper, result, control, initial, feasible)

   call assert_true(feasible(result%par), 'returned point is feasible')
   call assert_true(result%value < 2.02_dp, 'constrained objective')
   call assert_true(abs(sum(result%par) - 1.0_dp) < 2.0e-2_dp, 'active nonlinear boundary')
   call assert_true(maxval(abs(result%par - [0.0_dp, 1.0_dp])) < 4.0e-2_dp, 'constrained minimizer')

   print '(a)', 'PASS test_constraint'

contains

   function objective(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value
      value = (x(1) - 1.0_dp)**2 + (x(2) - 2.0_dp)**2
   end function objective

   logical function feasible(x)
      real(dp), intent(in) :: x(:)
      feasible = x(1) + x(2) <= 1.0_dp + 1.0e-12_dp
   end function feasible

   subroutine assert_true(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*, '(a)') 'FAIL: ' // label
         error stop 1
      end if
   end subroutine assert_true

end program test_constraint
