! SPDX-License-Identifier: MIT
program strategies_example
   use etrm
   implicit none

   real(dp) :: futures(6)
   type(strategy_result) :: result
   integer :: status
   character(len=:), allocatable :: message

   futures = [100.0_dp, 95.0_dp, 90.0_dp, 92.0_dp, 88.0_dp, 85.0_dp]
   call dppi(10.0_dp, futures, 0.10_dp, 0.12_dp, result, status, message)
   if (status /= etrm_ok) error stop message

   print '(a,6f10.3)', "Dynamic target: ", result%target
   print '(a,6f10.3)', "Hedge position: ", result%position
   print '(a,6f10.3)', "Portfolio price: ", result%portfolio
end program strategies_example
