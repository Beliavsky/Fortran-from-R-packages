! SPDX-License-Identifier: MIT
program test_msfc
   use etrm
   implicit none

   logical :: include(10)
   character(len=6) :: names(10)
   integer :: start_day(10), end_day(10), status, i
   real(dp) :: price(10)
   real(dp), allocatable :: prior(:)
   type(msfc_result) :: curve, curve_with_prior
   character(len=:), allocatable :: message

   include = .true.
   names = [character(len=6) :: "JUL-21", "AUG-21", "SEP-21", "OCT-21", "NOV-21", &
      "DEC-21", "Q1-22 ", "Q2-22 ", "Q3-22 ", "Q4-22 "]
   start_day = [ &
      day_offset(2021, 7, 1, 2021, 6, 17), day_offset(2021, 8, 1, 2021, 6, 17), &
      day_offset(2021, 9, 1, 2021, 6, 17), day_offset(2021, 10, 1, 2021, 6, 17), &
      day_offset(2021, 11, 1, 2021, 6, 17), day_offset(2021, 12, 1, 2021, 6, 17), &
      day_offset(2022, 1, 1, 2021, 6, 17), day_offset(2022, 4, 1, 2021, 6, 17), &
      day_offset(2022, 7, 1, 2021, 6, 17), day_offset(2022, 10, 1, 2021, 6, 17)]
   end_day = [ &
      day_offset(2021, 7, 30, 2021, 6, 17), day_offset(2021, 8, 31, 2021, 6, 17), &
      day_offset(2021, 9, 30, 2021, 6, 17), day_offset(2021, 10, 31, 2021, 6, 17), &
      day_offset(2021, 11, 30, 2021, 6, 17), day_offset(2021, 12, 31, 2021, 6, 17), &
      day_offset(2022, 3, 31, 2021, 6, 17), day_offset(2022, 6, 30, 2021, 6, 17), &
      day_offset(2022, 9, 30, 2021, 6, 17), day_offset(2022, 12, 31, 2021, 6, 17)]
   price = [32.55_dp, 32.50_dp, 32.50_dp, 32.08_dp, 36.88_dp, &
      39.80_dp, 39.40_dp, 25.20_dp, 21.15_dp, 29.50_dp]

   call msfc(include, names, start_day, end_day, price, curve, status, message)
   call assert_status(status, message)
   if (curve%n_days /= 563) error stop "Unexpected MSFC daily length"
   if (curve%n_polynomials /= 20) error stop "Unexpected polynomial count"
   if (maxval(abs(curve%computed_price - price)) > 1.0e-9_dp) then
      error stop "MSFC failed contract repricing"
   end if
   call assert_vector(curve%curve(1:5), [32.7850385852_dp, 32.7766459646_dp, &
      32.7682457028_dp, 32.7598443497_dp, 32.7514468574_dp], 2.0e-9_dp, &
      "MSFC fixed reference")

   allocate(prior(curve%n_days))
   do i = 1, size(prior)
      prior(i) = 5.0_dp + 2.0_dp * sin(2.0_dp * acos(-1.0_dp) * real(i - 1, dp) / 365.0_dp)
   end do
   call msfc(include, names, start_day, end_day, price, prior, curve_with_prior, status, message)
   call assert_status(status, message)
   if (maxval(abs(curve_with_prior%computed_price - price)) > 1.0e-9_dp) then
      error stop "Prior-adjusted MSFC failed contract repricing"
   end if
   if (maxval(abs(curve_with_prior%curve - curve%curve)) < 1.0e-5_dp) then
      error stop "Nonconstant prior did not alter daily curve"
   end if

   include(10) = .false.
   call msfc(include, names, start_day, end_day, price, curve, status, message)
   call assert_status(status, message)
   if (curve%n_contracts /= 9) error stop "Contract inclusion filter failed"

   print '(a)', "MSFC tests passed"

contains

   subroutine assert_status(code, text)
      integer, intent(in) :: code
      character(len=*), intent(in) :: text
      if (code /= etrm_ok) then
         print '(a)', trim(text)
         error stop "Unexpected MSFC status"
      end if
   end subroutine assert_status

   subroutine assert_vector(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual(:), expected(:), tolerance
      character(len=*), intent(in) :: label
      if (size(actual) /= size(expected)) error stop "Vector size assertion failed"
      if (maxval(abs(actual - expected)) > tolerance) then
         print '(a,es24.15)', trim(label)//" maximum error: ", maxval(abs(actual - expected))
         error stop "Vector assertion failed"
      end if
   end subroutine assert_vector

end program test_msfc
