! SPDX-License-Identifier: MIT
program msfc_example
   use etrm
   implicit none

   logical :: include(3)
   character(len=5) :: names(3)
   integer :: start_day(3), end_day(3), status, i
   real(dp) :: price(3)
   type(msfc_result) :: result
   character(len=:), allocatable :: message

   include = .true.
   names = [character(len=5) :: "JUL", "AUG", "SEP"]
   start_day = [14, 45, 76]
   end_day = [43, 75, 105]
   price = [32.55_dp, 32.50_dp, 32.50_dp]

   call msfc(include, names, start_day, end_day, price, result, status, message)
   if (status /= etrm_ok) error stop message

   print '(a,es12.4)', "Maximum contract repricing error: ", &
      maxval(abs(result%computed_price - result%market_price))
   print '(a)', "First ten daily curve values:"
   do i = 1, min(10, result%n_days)
      print '(i5,f14.6)', result%day(i), result%curve(i)
   end do
end program msfc_example
