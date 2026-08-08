module test_stosoo_sto_callbacks
   use oor, only : dp
   implicit none
contains
   function quadratic(x) result(y)
      real(dp), intent(in) :: x(:)
      real(dp) :: y
      y = (x(1) - 0.2_dp)**2 + (x(2) - 0.7_dp)**2
   end function quadratic
end module test_stosoo_sto_callbacks

program test_stosoo_stochastic
   use oor, only : dp, stosoo, stosoo_options, stosoo_result
   use test_stosoo_sto_callbacks, only : quadratic
   implicit none
   type(stosoo_options) :: options
   type(stosoo_result) :: result
   real(dp) :: lower(2), upper(2)

   lower = 0.0_dp
   upper = 1.0_dp
   options%stochastic = .true.
   options%k_max = 3
   options%h_max = 12
   call stosoo(quadratic, lower, upper, 500, result, options)
   if (result%value > 0.01_dp) error stop "stochastic StoSOO value failure"
   if (result%evaluations < 499 .or. result%evaluations > 501) then
      error stop "stochastic StoSOO evaluation accounting failure"
   end if
   print *, "test_stosoo_stochastic: PASS"
end program test_stosoo_stochastic
