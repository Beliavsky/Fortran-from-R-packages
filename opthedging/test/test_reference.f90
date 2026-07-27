! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2013 Bruno Remillard
! Modern Fortran translation copyright (C) 2026 OpenAI
program test_reference
   use opthedging, only : dp, hedging_iid, hedging_result
   implicit none

   real(dp) :: factors(9)
   real(dp) :: expected_a(5)
   real(dp) :: expected_c(5)
   type(hedging_result) :: result

   factors = [0.78_dp, 0.86_dp, 0.94_dp, 0.99_dp, 1.01_dp, 1.05_dp, &
      1.11_dp, 1.19_dp, 1.28_dp]
   expected_c = [8.440248535662011_dp, 7.101498478371031_dp, &
      5.625407483587852_dp, 3.9582748208191756_dp, 2.0997637529399498_dp]
   expected_a = [-23.149186951809092_dp, -23.242665805554214_dp, &
      -22.706493093052355_dp, -20.831440666589580_dp, -15.444704669654602_dp]

   result = hedging_iid(log(factors), 0.75_dp, 105.0_dp, 0.04_dp, .true., &
      5, 501, 30.0_dp, 200.0_dp)
   if (.not. result%ok) error stop 1

   call check_close(result%rho, 1.045296167247387_dp, 2.0e-14_dp, "rho")
   call check_vector(result%c(:, 251), expected_c, 2.0e-11_dp, "C reference")
   call check_vector(result%a(:, 251), expected_a, 2.0e-11_dp, "a reference")
   call check_close(result%phi1(251), -0.2780151860587127_dp, &
      2.0e-12_dp, "phi reference")

   print '(a)', 'test_reference: PASS'

contains

   subroutine check_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual
      real(dp), intent(in) :: expected
      real(dp), intent(in) :: tolerance
      character(len=*), intent(in) :: label

      if (abs(actual - expected) > tolerance) then
         print '(a)', 'failed: ' // label
         print '(a,es24.16)', 'actual:   ', actual
         print '(a,es24.16)', 'expected: ', expected
         error stop 1
      end if
   end subroutine check_close

   subroutine check_vector(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual(:)
      real(dp), intent(in) :: expected(:)
      real(dp), intent(in) :: tolerance
      character(len=*), intent(in) :: label

      if (size(actual) /= size(expected)) error stop 1
      if (maxval(abs(actual - expected)) > tolerance) then
         print '(a)', 'failed: ' // label
         print '(a,es24.16)', 'maximum error: ', maxval(abs(actual - expected))
         error stop 1
      end if
   end subroutine check_vector

end program test_reference
