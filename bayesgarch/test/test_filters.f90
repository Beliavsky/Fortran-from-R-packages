! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2008-2021 David Ardia
! Modern Fortran translation of computational routines from bayesGARCH.
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
program test_filters
   use bayesgarch_garch, only : filter_alpha, filter_alpha_asymmetric, filter_w, &
      filter_w_asymmetric, garch11_filter, garch11_student_loglik, garch_filter, &
      quasi_difference, simulate_garch, simulate_garch11_student, simulate_threshold_garch, &
      threshold_garch_filter
   use bayesgarch_kinds, only : dp
   use bayesgarch_rng, only : seed_rng
   use test_support, only : assert_close, assert_true
   implicit none

   real(dp), parameter :: y(5) = [0.2_dp, -0.1_dp, 0.3_dp, -0.4_dp, 0.25_dp]
   real(dp), parameter :: theta(3) = [0.05_dp, 0.12_dp, 0.8_dp]
   real(dp) :: h(6)
   real(dp) :: x(5, 2)
   real(dp) :: xa(5, 3)
   real(dp) :: w(5)
   real(dp) :: wa(5)
   real(dp) :: q(5)
   real(dp) :: ysim(200)
   real(dp) :: hsim(201)
   real(dp) :: hcheck(201)
   real(dp) :: hgeneral(6)
   real(dp) :: hthreshold(6)
   real(dp) :: ll
   real(dp), parameter :: innovations(5) = [0.5_dp, -1.0_dp, 0.25_dp, -0.75_dp, 1.2_dp]
   real(dp) :: ygeneral(5)
   real(dp) :: ythreshold(5)
   real(dp) :: hsim_general(6)
   real(dp) :: hsim_threshold(6)

   call garch11_filter(y, theta, h)
   call assert_close(h, [0.05_dp, 0.0948_dp, 0.12704_dp, 0.162432_dp, &
      0.1991456_dp, 0.21681648_dp], 1.0e-13_dp, "GARCH(1,1) recursion")

   call filter_alpha(y, theta(3), x)
   call assert_close(x(:, 1), [1.0_dp, 1.8_dp, 2.44_dp, 2.952_dp, 3.3616_dp], &
      1.0e-13_dp, "symmetric alpha lstar")
   call assert_close(x(:, 2), [0.0_dp, 0.04_dp, 0.042_dp, 0.1236_dp, 0.25888_dp], &
      1.0e-13_dp, "symmetric alpha vstar")

   call filter_alpha_asymmetric(y, theta(3), xa)
   call assert_close(xa(:, 1), [0.0_dp, 1.0_dp, 1.8_dp, 2.44_dp, 2.952_dp], &
      1.0e-13_dp, "asymmetric alpha lstar")
   call assert_close(xa(:, 2), [0.0_dp, 0.04_dp, 0.032_dp, 0.1156_dp, 0.09248_dp], &
      1.0e-13_dp, "asymmetric alpha positive")
   call assert_close(xa(:, 3), [0.0_dp, 0.04_dp, 0.042_dp, 0.0336_dp, 0.18688_dp], &
      1.0e-13_dp, "asymmetric alpha negative")

   call filter_w(y, theta(1:2), theta(3), w)
   call assert_close(w, [-0.01_dp, -0.0848_dp, -0.03704_dp, -0.002432_dp, &
      -0.1366456_dp], 1.0e-13_dp, "symmetric W filter")
   call quasi_difference(y**2 - w, -theta(3), q)
   call assert_close(q, [0.0_dp, 0.05_dp, 0.1348_dp, 0.23488_dp, 0.350336_dp], &
      1.0e-13_dp, "quasi difference")

   call filter_w_asymmetric(y, [0.05_dp, 0.12_dp, 0.18_dp], 0.8_dp, wa)
   call assert_close(wa, [-0.01_dp, -0.0848_dp, -0.03764_dp, -0.002912_dp, &
      -0.1466296_dp], 1.0e-13_dp, "asymmetric W filter")

   call garch_filter(y, theta(1), [theta(2)], [theta(3)], hgeneral)
   call assert_close(hgeneral, h, 1.0e-13_dp, "general GARCH recursion")
   call threshold_garch_filter(y, 0.05_dp, [0.12_dp], [0.8_dp], [0.18_dp], hthreshold)
   call assert_close(hthreshold, [0.05_dp, 0.0948_dp, 0.12764_dp, 0.162912_dp, &
      0.2091296_dp, 0.22480368_dp], 1.0e-13_dp, "threshold GARCH recursion")

   call simulate_garch(innovations, 0.05_dp, [0.12_dp], [0.8_dp], ygeneral, hsim_general)
   call garch_filter(ygeneral, 0.05_dp, [0.12_dp], [0.8_dp], hgeneral)
   call assert_close(hsim_general, hgeneral, 1.0e-13_dp, "general GARCH simulation/filter identity")
   call simulate_threshold_garch(innovations, 0.05_dp, [0.12_dp], [0.8_dp], [0.18_dp], &
      ythreshold, hsim_threshold)
   call threshold_garch_filter(ythreshold, 0.05_dp, [0.12_dp], [0.8_dp], [0.18_dp], hthreshold)
   call assert_close(hsim_threshold, hthreshold, 1.0e-13_dp, &
      "threshold GARCH simulation/filter identity")

   ll = garch11_student_loglik(y, theta, 8.0_dp)
   call assert_close(ll, -0.665609339519487_dp, 1.0e-12_dp, "Student-t log likelihood")

   call seed_rng(24680)
   call simulate_garch11_student(size(ysim), theta, 8.0_dp, ysim, hsim)
   call garch11_filter(ysim, theta, hcheck)
   call assert_close(hsim, hcheck, 1.0e-12_dp, "simulation/filter identity")
   call assert_true(all(hsim > 0.0_dp), "simulated variances are positive")

   write(*, '(a)') "Filter, likelihood, and simulation tests passed."
end program test_filters
