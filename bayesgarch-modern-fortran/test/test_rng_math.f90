! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2008-2021 David Ardia
! Modern Fortran translation of computational routines from bayesGARCH.
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
program test_rng_math
   use bayesgarch_kinds, only : dp
   use bayesgarch_math, only : digamma
   use bayesgarch_rng, only : random_exponential, random_gamma, random_normal, &
      random_standardized_student, seed_rng
   use test_support, only : assert_close, assert_true
   implicit none

   integer, parameter :: n = 100000
   real(dp), parameter :: euler_gamma = 0.5772156649015328606_dp
   real(dp) :: sum_x
   real(dp) :: sum_x2
   real(dp) :: x
   real(dp) :: mean
   real(dp) :: variance
   integer :: i

   call assert_close(digamma(1.0_dp), -euler_gamma, 2.0e-12_dp, "digamma at one")
   call assert_close(digamma(0.5_dp), -euler_gamma - 2.0_dp * log(2.0_dp), &
      2.0e-12_dp, "digamma at one half")

   call seed_rng(10101)
   sum_x = 0.0_dp
   sum_x2 = 0.0_dp
   do i = 1, n
      x = random_normal()
      sum_x = sum_x + x
      sum_x2 = sum_x2 + x * x
   end do
   mean = sum_x / real(n, dp)
   variance = sum_x2 / real(n, dp) - mean * mean
   call assert_true(abs(mean) < 0.015_dp, "normal RNG mean")
   call assert_true(abs(variance - 1.0_dp) < 0.025_dp, "normal RNG variance")

   call seed_rng(20202)
   sum_x = 0.0_dp
   sum_x2 = 0.0_dp
   do i = 1, n
      x = random_gamma(2.0_dp, 3.0_dp)
      sum_x = sum_x + x
      sum_x2 = sum_x2 + x * x
   end do
   mean = sum_x / real(n, dp)
   variance = sum_x2 / real(n, dp) - mean * mean
   call assert_true(abs(mean - 2.0_dp / 3.0_dp) < 0.012_dp, "gamma RNG mean")
   call assert_true(abs(variance - 2.0_dp / 9.0_dp) < 0.012_dp, "gamma RNG variance")

   call seed_rng(30303)
   sum_x = 0.0_dp
   do i = 1, n
      sum_x = sum_x + random_exponential(2.0_dp)
   end do
   mean = sum_x / real(n, dp)
   call assert_true(abs(mean - 0.5_dp) < 0.008_dp, "exponential RNG mean")

   call seed_rng(40404)
   sum_x = 0.0_dp
   sum_x2 = 0.0_dp
   do i = 1, n
      x = random_standardized_student(8.0_dp)
      sum_x = sum_x + x
      sum_x2 = sum_x2 + x * x
   end do
   mean = sum_x / real(n, dp)
   variance = sum_x2 / real(n, dp) - mean * mean
   call assert_true(abs(mean) < 0.02_dp, "standardized Student-t RNG mean")
   call assert_true(abs(variance - 1.0_dp) < 0.04_dp, "standardized Student-t RNG variance")

   write(*, '(a)') "Random-number and special-function tests passed."
end program test_rng_math
