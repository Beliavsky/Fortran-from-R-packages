module test_stosoo_det_callbacks
   use oor, only : dp, guirland
   implicit none
contains
   function objective(x) result(y)
      real(dp), intent(in) :: x(:)
      real(dp) :: y
      y = -guirland(x(1))
   end function objective
end module test_stosoo_det_callbacks

program test_stosoo_deterministic
   use oor, only : dp, stosoo, stosoo_options, stosoo_result
   use test_stosoo_det_callbacks, only : objective
   implicit none
   type(stosoo_options) :: options
   type(stosoo_result) :: result
   real(dp) :: lower(1), upper(1)

   lower = 0.0_dp
   upper = 1.0_dp
   options%stochastic = .false.
   call stosoo(objective, lower, upper, 500, result, options)
   if (abs(result%par(1) - acos(-1.0_dp)/6.0_dp) > 2.0e-3_dp) then
      error stop "deterministic StoSOO parameter failure"
   end if
   if (result%value > -0.99_dp) error stop "deterministic StoSOO value failure"
   if (result%evaluations < 500) error stop "deterministic StoSOO evaluation accounting failure"
   print *, "test_stosoo_deterministic: PASS"
end program test_stosoo_deterministic
