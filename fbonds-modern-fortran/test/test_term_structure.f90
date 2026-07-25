! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 Modern Fortran translation contributors
! Based on fBonds, copyright its original authors.
! This file is free software: you may redistribute it and/or modify it
! under the terms of the GNU General Public License, version 2 or later.
program test_term_structure
   use fbonds_kinds, only : dp
   use fbonds_term_structure, only : term_structure_fit, nelson_siegel_curve, &
      svensson_curve, fit_nelson_siegel, fit_svensson
   implicit none

   call test_curve_values()
   call test_nelson_siegel_fit()
   call test_svensson_fit()
   call test_reference_data()
   print '(a)', 'Term-structure tests passed.'
contains
   subroutine test_curve_values()
      real(dp) :: maturity(4), ns_par(4), sv_par(6), ns_rate(4), sv_rate(4)
      maturity = [0.0_dp, 0.5_dp, 2.0_dp, 10.0_dp]
      ns_par = [0.05_dp, -0.02_dp, 0.03_dp, 2.0_dp]
      sv_par = [0.04_dp, -0.01_dp, 0.02_dp, -0.005_dp, 1.5_dp, 6.0_dp]
      ns_rate = nelson_siegel_curve(maturity, ns_par)
      sv_rate = svensson_curve(maturity, sv_par)
      call assert_close(ns_rate(1), 0.03_dp, 1.0e-13_dp, 'NS zero-maturity limit')
      call assert_close(sv_rate(1), 0.03_dp, 1.0e-13_dp, 'Svensson zero-maturity limit')
      call assert_true(all(ns_rate < 0.08_dp) .and. all(ns_rate > -0.02_dp), &
         'NS curve finite range')
      call assert_true(all(sv_rate < 0.08_dp) .and. all(sv_rate > -0.02_dp), &
         'Svensson curve finite range')
   end subroutine test_curve_values

   subroutine test_nelson_siegel_fit()
      integer, parameter :: n = 40
      real(dp) :: maturity(n), truth(4), rate(n)
      type(term_structure_fit) :: fit = term_structure_fit()
      integer :: i
      truth = [0.047_dp, -0.026_dp, 0.038_dp, 2.4_dp]
      do i = 1, n
         maturity(i) = 0.25_dp + real(i - 1, dp) * 0.5_dp
      end do
      rate = nelson_siegel_curve(maturity, truth)
      call fit_nelson_siegel(rate, maturity, fit, max_iterations=3000, tolerance=1.0e-11_dp)
      call assert_true(fit%converged, 'NS fit convergence')
      call assert_true(fit%rmse < 1.0e-8_dp, 'NS exact-data RMSE')
      call assert_true(maxval(abs(fit%parameters - truth)) < 2.0e-4_dp, &
         'NS parameter recovery')
   end subroutine test_nelson_siegel_fit

   subroutine test_svensson_fit()
      integer, parameter :: n = 48
      real(dp) :: maturity(n), truth(6), rate(n)
      type(term_structure_fit) :: fit_sse = term_structure_fit(), fit_l1 = term_structure_fit()
      integer :: i
      truth = [0.041_dp, -0.017_dp, 0.024_dp, -0.011_dp, 1.8_dp, 7.0_dp]
      do i = 1, n
         maturity(i) = 0.25_dp + real(i - 1, dp) * 0.5_dp
      end do
      rate = svensson_curve(maturity, truth)
      call fit_svensson(rate, maturity, fit_sse, objective='sse', &
         max_iterations=5000, tolerance=1.0e-10_dp)
      call assert_true(fit_sse%rmse < 2.0e-5_dp, 'Svensson SSE exact-data RMSE')
      call fit_svensson(rate, maturity, fit_l1, max_iterations=5000, tolerance=1.0e-10_dp)
      call assert_true(fit_l1%mae < 2.0e-5_dp, 'Svensson source L1 exact-data MAE')
      call assert_true(all(fit_l1%parameters(5:6) > 0.0_dp), 'Svensson positive taus')
   end subroutine test_svensson_fit

   subroutine test_reference_data()
      real(dp), parameter :: rate(48) = [ &
         0.04984_dp, 0.05283_dp, 0.05549_dp, 0.05777_dp, 0.05961_dp, 0.06102_dp, &
         0.06216_dp, 0.06314_dp, 0.06403_dp, 0.06488_dp, 0.06568_dp, 0.06644_dp, &
         0.06717_dp, 0.06786_dp, 0.06852_dp, 0.06913_dp, 0.06969_dp, 0.07020_dp, &
         0.07134_dp, 0.07205_dp, 0.07339_dp, 0.07500_dp, 0.07710_dp, 0.07860_dp, &
         0.08011_dp, 0.08114_dp, 0.08194_dp, 0.08274_dp, 0.08355_dp, 0.08434_dp, &
         0.08512_dp, 0.08588_dp, 0.08662_dp, 0.08731_dp, 0.08794_dp, 0.08851_dp, &
         0.08900_dp, 0.08939_dp, 0.08967_dp, 0.08980_dp, 0.08976_dp, 0.08954_dp, &
         0.08910_dp, 0.08843_dp, 0.08748_dp, 0.08626_dp, 0.08474_dp, 0.08291_dp]
      real(dp), parameter :: maturity(48) = [ &
         0.083_dp, 0.167_dp, 0.250_dp, 0.333_dp, 0.417_dp, 0.500_dp, 0.583_dp, &
         0.667_dp, 0.750_dp, 0.833_dp, 0.917_dp, 1.000_dp, 1.083_dp, 1.167_dp, &
         1.250_dp, 1.333_dp, 1.417_dp, 1.500_dp, 1.750_dp, 2.000_dp, 2.500_dp, &
         3.000_dp, 4.000_dp, 5.000_dp, 6.000_dp, 7.000_dp, 8.000_dp, 9.000_dp, &
         10.000_dp, 11.000_dp, 12.000_dp, 13.000_dp, 14.000_dp, 15.000_dp, &
         16.000_dp, 17.000_dp, 18.000_dp, 19.000_dp, 20.000_dp, 21.000_dp, &
         22.000_dp, 23.000_dp, 24.000_dp, 25.000_dp, 26.000_dp, 27.000_dp, &
         28.000_dp, 29.000_dp]
      type(term_structure_fit) :: ns_fit = term_structure_fit(), sv_fit = term_structure_fit()
      real(dp) :: constant_sse
      constant_sse = sum((rate - sum(rate) / real(size(rate), dp))**2)
      call fit_nelson_siegel(rate, maturity, ns_fit, max_iterations=4000)
      call fit_svensson(rate, maturity, sv_fit, objective='sse', max_iterations=5000)
      call assert_true(ns_fit%sse < 0.15_dp * constant_sse, 'NS reference-data improvement')
      call assert_true(sv_fit%sse < ns_fit%sse, 'Svensson reference-data improvement')
      call assert_true(all(ns_fit%parameters(4:4) > 0.0_dp), 'NS positive tau')
      call assert_true(all(sv_fit%parameters(5:6) > 0.0_dp), 'Svensson positive taus')
   end subroutine test_reference_data

   subroutine assert_close(actual, expected, tolerance, message)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: message
      if (abs(actual - expected) > tolerance) then
         print '(a,2(1x,es24.15))', trim(message), actual, expected
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         print '(a)', trim(message)
         error stop 1
      end if
   end subroutine assert_true
end program test_term_structure
