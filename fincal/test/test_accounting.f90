! SPDX-License-Identifier: GPL-2.0-or-later
program test_accounting
   use fincal
   implicit none
   real(dp), parameter :: tol = 2.0e-12_dp
   real(dp), allocatable :: depreciation(:)
   type(inventory_result) :: result_value
   integer :: status

   call assert_close(eps(10000.0_dp, 1000.0_dp, 11000.0_dp), 9.0_dp / 11.0_dp, tol, 'eps')
   call assert_close(diluted_eps(115600.0_dp, 10000.0_dp, 200000.0_dp, &
      convertible_debt_interest = 42000.0_dp, tax_rate = 0.4_dp, convertible_debt_shares = 60000.0_dp), &
      0.503076923076923_dp, tol, 'diluted_eps')
   call assert_close(weighted_average_shares([10000.0_dp, 2000.0_dp], [12.0_dp, 6.0_dp]), &
      11000.0_dp, tol, 'weighted_average_shares')
   call assert_close(issuable_shares(20.0_dp, 15.0_dp, 10000.0_dp), 2500.0_dp, tol, 'issuable_shares')
   call assert_close(straight_line_depreciation(1200.0_dp, 200.0_dp, 5), 200.0_dp, tol, 'slde')
   call assert_close(slde(1200.0_dp, 200.0_dp, 5), 200.0_dp, tol, 'slde alias')
   call assert_close(was([10000.0_dp, 2000.0_dp], [12.0_dp, 6.0_dp]), 11000.0_dp, tol, 'was alias')
   call assert_close(iss(20.0_dp, 15.0_dp, 10000.0_dp), 2500.0_dp, tol, 'iss alias')

   depreciation = double_declining_balance(1200.0_dp, 200.0_dp, 5, status)
   call assert_equal_int(status, fincal_ok, 'ddb status')
   call assert_array_close(depreciation, [480.0_dp, 288.0_dp, 172.8_dp, 59.2_dp, 0.0_dp], 1.0e-10_dp, 'ddb')

   result_value = cogs(2.0_dp, 2.0_dp, [3.0_dp, 5.0_dp], [3.0_dp, 5.0_dp], 7.0_dp, 'FIFO')
   call assert_equal_int(result_value%status, fincal_ok, 'fifo status')
   call assert_close(result_value%cost_of_goods, 23.0_dp, tol, 'fifo cogs')
   call assert_close(result_value%ending_inventory, 15.0_dp, tol, 'fifo ending')

   result_value = cogs(2.0_dp, 2.0_dp, [3.0_dp, 5.0_dp], [3.0_dp, 5.0_dp], 7.0_dp, 'LIFO')
   call assert_close(result_value%cost_of_goods, 31.0_dp, tol, 'lifo cogs')
   call assert_close(result_value%ending_inventory, 7.0_dp, tol, 'lifo ending')

   result_value = cogs(2.0_dp, 2.0_dp, [3.0_dp, 5.0_dp], [3.0_dp, 5.0_dp], 7.0_dp, 'WAC')
   call assert_close(result_value%cost_of_goods, 26.6_dp, tol, 'wac cogs')
   call assert_close(result_value%ending_inventory, 11.4_dp, tol, 'wac ending')

   result_value = cogs(2.0_dp, 2.0_dp, [3.0_dp], [3.0_dp], 6.0_dp, 'FIFO')
   call assert_equal_int(result_value%status, fincal_insufficient_inventory, 'insufficient inventory')
   print '(a)', 'test_accounting: PASS'
contains
   subroutine assert_close(actual, expected, tolerance, message)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: message
      if (abs(actual - expected) > tolerance) then
         print *, trim(message), actual, expected, abs(actual - expected)
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_array_close(actual, expected, tolerance, message)
      real(dp), intent(in) :: actual(:), expected(:), tolerance
      character(len=*), intent(in) :: message
      if (size(actual) /= size(expected) .or. any(abs(actual - expected) > tolerance)) then
         print *, trim(message), actual, expected
         error stop 1
      end if
   end subroutine assert_array_close

   subroutine assert_equal_int(actual, expected, message)
      integer, intent(in) :: actual, expected
      character(len=*), intent(in) :: message
      if (actual /= expected) then
         print *, trim(message), actual, expected
         error stop 1
      end if
   end subroutine assert_equal_int
end program test_accounting
