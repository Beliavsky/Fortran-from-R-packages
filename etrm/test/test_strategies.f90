! SPDX-License-Identifier: MIT
program test_strategies
   use etrm
   implicit none

   real(dp), parameter :: tol = 2.0e-10_dp
   real(dp) :: f(8), falling(6), seller_prices(8)
   type(strategy_result) :: result
   type(strategy_summary) :: summary
   integer :: status
   character(len=:), allocatable :: message

   f = [100.0_dp, 102.0_dp, 105.0_dp, 103.0_dp, 108.0_dp, 112.0_dp, 109.0_dp, 115.0_dp]

   call cppi(10.0_dp, f, 0.10_dp, 0.12_dp, result, status, message, &
      transaction_cost=0.05_dp, integer_trades=.true.)
   call assert_status(status, message, "CPPI")
   call assert_vector(result%position, [2.0_dp, 2.0_dp, 3.0_dp, 5.0_dp, &
      4.0_dp, 6.0_dp, 8.0_dp, 7.0_dp], tol, "CPPI positions")
   call assert_vector(result%portfolio, [100.01_dp, 101.61_dp, 104.015_dp, 102.625_dp, &
      105.13_dp, 107.54_dp, 106.35_dp, 107.555_dp], tol, "CPPI portfolio")
   call summarize_strategy(result, summary)
   call assert_close(summary%churn, 1.1_dp, tol, "CPPI churn")

   falling = [100.0_dp, 95.0_dp, 90.0_dp, 92.0_dp, 88.0_dp, 85.0_dp]
   call dppi(10.0_dp, falling, 0.10_dp, 0.12_dp, result, status, message, &
      integer_trades=.true.)
   call assert_status(status, message, "DPPI")
   call assert_vector(result%position, [2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 4.0_dp, 2.0_dp], &
      tol, "DPPI positions")
   call assert_vector(result%target, [110.0_dp, 105.6_dp, 101.2_dp, 101.2_dp, &
      99.44_dp, 97.46_dp], tol, "DPPI target")
   call assert_vector(result%portfolio, [100.0_dp, 96.0_dp, 92.0_dp, 93.6_dp, &
      90.4_dp, 88.6_dp], tol, "DPPI portfolio")

   call obpi(10.0_dp, f, 0.20_dp, 8, result, status, message, &
      transaction_cost=0.05_dp, integer_trades=.true.)
   call assert_status(status, message, "OBPI")
   call assert_vector(result%position, [5.0_dp, 7.0_dp, 9.0_dp, 9.0_dp, &
      10.0_dp, 10.0_dp, 10.0_dp, 10.0_dp], tol, "OBPI positions")
   call assert_vector(result%portfolio, [100.025_dp, 101.035_dp, 101.945_dp, 101.745_dp, &
      102.25_dp, 102.25_dp, 102.25_dp, 102.25_dp], tol, "OBPI portfolio")
   call assert_close(result%target(1), 101.4272231739467_dp, 2.0e-12_dp, "OBPI target")

   call slpi(10.0_dp, f, 0.10_dp, result, status, message, transaction_cost=0.05_dp)
   call assert_status(status, message, "SLPI")
   call assert_vector(result%position, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 10.0_dp, 10.0_dp, 10.0_dp], tol, "SLPI positions")
   call assert_vector(result%portfolio, [100.0_dp, 102.0_dp, 105.0_dp, 103.0_dp, &
      108.0_dp, 112.05_dp, 112.05_dp, 112.05_dp], tol, "SLPI portfolio")

   call shpi(10.0_dp, f, 8, 0.10_dp, result, status, message, &
      transaction_cost=0.05_dp, integer_trades=.true.)
   call assert_status(status, message, "SHPI")
   call assert_vector(result%position, [1.0_dp, 2.0_dp, 4.0_dp, 5.0_dp, &
      6.0_dp, 8.0_dp, 9.0_dp, 10.0_dp], tol, "SHPI positions")
   call assert_vector(result%portfolio, [100.005_dp, 101.81_dp, 104.22_dp, 103.025_dp, &
      105.53_dp, 107.14_dp, 106.545_dp, 107.15_dp], tol, "SHPI portfolio")

   seller_prices = [100.0_dp, 98.0_dp, 95.0_dp, 97.0_dp, 92.0_dp, &
      88.0_dp, 91.0_dp, 85.0_dp]
   call obpi(-10.0_dp, seller_prices, 0.20_dp, 8, result, status, message, &
      transaction_cost=0.05_dp, integer_trades=.true.)
   call assert_status(status, message, "seller OBPI")
   call assert_vector(result%position, [-5.0_dp, -7.0_dp, -9.0_dp, -9.0_dp, &
      -10.0_dp, -10.0_dp, -10.0_dp, -10.0_dp], tol, "seller OBPI positions")
   call assert_close(result%target(1), 98.5727768260533_dp, 2.0e-12_dp, &
      "seller OBPI target")

   call slpi(-10.0_dp, seller_prices, -0.10_dp, result, status, message, &
      transaction_cost=0.05_dp)
   call assert_status(status, message, "seller SLPI")
   call assert_vector(result%position, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, -10.0_dp, -10.0_dp, -10.0_dp], tol, "seller SLPI positions")

   call cppi(0.0_dp, f, 0.10_dp, 0.12_dp, result, status, message)
   if (status /= etrm_err_argument) error stop "Expected invalid-volume status"

   call obpi(10.0_dp, f, 0.20_dp, 7, result, status, message)
   if (status /= etrm_err_argument) error stop "Expected days_left validation status"

   print '(a)', "strategy tests passed"

contains

   subroutine assert_status(code, text, label)
      integer, intent(in) :: code
      character(len=*), intent(in) :: text, label
      if (code /= etrm_ok) then
         print '(a,2a)', trim(label), ": ", trim(text)
         error stop "Unexpected status"
      end if
   end subroutine assert_status

   subroutine assert_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: label
      if (abs(actual - expected) > tolerance) then
         print '(a,2es24.15)', trim(label)//" actual/expected: ", actual, expected
         error stop "Scalar assertion failed"
      end if
   end subroutine assert_close

   subroutine assert_vector(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual(:), expected(:), tolerance
      character(len=*), intent(in) :: label
      if (size(actual) /= size(expected)) error stop "Vector size assertion failed"
      if (maxval(abs(actual - expected)) > tolerance) then
         print '(a,es24.15)', trim(label)//" maximum error: ", maxval(abs(actual - expected))
         error stop "Vector assertion failed"
      end if
   end subroutine assert_vector

end program test_strategies
