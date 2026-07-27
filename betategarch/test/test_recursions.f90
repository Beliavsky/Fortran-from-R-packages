! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2011-2025 Genaro Sucarrat
! Copyright (C) 2026 contributors to the Modern Fortran translation
!
! This file is part of betategarch-modern-fortran.
! It is free software: you can redistribute it and/or modify it under the
! terms of the GNU General Public License version 2 only.

program test_recursions
  use betategarch, only : dp, tegarch_parameters, tegarch_filter_result, &
    tegarch_filter, tegarch_loglik, tegarch_simulate
  implicit none

  real(dp), parameter :: y(5) = [0.2_dp, -0.5_dp, 0.1_dp, -1.2_dp, 0.7_dp]
  real(dp), parameter :: lambda1_expected(5) = [0.03_dp, -0.05113463828788703_dp, &
    -0.08369688552525649_dp, -0.153226225376797_dp, -0.042354109062008635_dp]
  real(dp), parameter :: dagger1_expected(5) = [0.0_dp, -0.08113463828788703_dp, &
    -0.11369688552525649_dp, -0.183226225376797_dp, -0.07235410906200863_dp]
  real(dp), parameter :: score1_expected(4) = [-1.0283659571971755_dp, -0.6722975922179846_dp, &
    -1.0224757101016546_dp, 0.43791244814257224_dp]
  real(dp), parameter :: lambda2_expected(5) = [0.02_dp, -0.06554301929241725_dp, &
    -0.10854302355217747_dp, -0.17026859598916128_dp, -0.00221176383308613_dp]
  real(dp), parameter :: dagger21_expected(5) = [0.0_dp, -0.028089865411035968_dp, &
    -0.04895094769359137_dp, -0.07566800427144288_dp, -0.0419137164868196_dp]
  real(dp), parameter :: dagger22_expected(5) = [0.0_dp, -0.05745315388138129_dp, &
    -0.0795920758585861_dp, -0.11460059171771839_dp, 0.019701952653733472_dp]
  real(dp), parameter :: score2_expected(4) = [-0.9363288470345323_dp, -0.7421858517702401_dp, &
    -0.9721534654177029_dp, 0.9990295857017042_dp]

  type(tegarch_parameters) :: p1, p2
  type(tegarch_filter_result) :: f1, f2, sim1, sim2, refiltered

  p1%components = 1
  p1%asym = .true.
  p1%skewed = .true.
  p1%omega = 0.03_dp
  p1%phi1 = 0.9_dp
  p1%kappa1 = 0.08_dp
  p1%kappastar = 0.04_dp
  p1%df = 8.0_dp
  p1%skew = 0.8_dp
  call tegarch_filter(y, p1, f1)
  call assert_vector(f1%lambda, lambda1_expected, 5.0e-14_dp, "one-component lambda")
  call assert_vector(f1%lambda1_dagger, dagger1_expected, 5.0e-14_dp, "one-component dagger")
  call assert_vector(f1%score(1:4), score1_expected, 5.0e-14_dp, "one-component score")
  call assert_close(tegarch_loglik(y, p1), -6.045430005781165_dp, 7.0e-14_dp, "one-component loglik")

  p2%components = 2
  p2%asym = .true.
  p2%skewed = .true.
  p2%omega = 0.02_dp
  p2%phi1 = 0.95_dp
  p2%phi2 = 0.7_dp
  p2%kappa1 = 0.03_dp
  p2%kappa2 = 0.06_dp
  p2%kappastar = 0.02_dp
  p2%df = 10.0_dp
  p2%skew = 1.1_dp
  call tegarch_filter(y, p2, f2)
  call assert_vector(f2%lambda, lambda2_expected, 6.0e-14_dp, "two-component lambda")
  call assert_vector(f2%lambda1_dagger, dagger21_expected, 6.0e-14_dp, "two-component dagger 1")
  call assert_vector(f2%lambda2_dagger, dagger22_expected, 6.0e-14_dp, "two-component dagger 2")
  call assert_vector(f2%score(1:4), score2_expected, 6.0e-14_dp, "two-component score")
  call assert_close(tegarch_loglik(y, p2), -5.893354473366062_dp, 8.0e-14_dp, "two-component loglik")

  call tegarch_simulate(300, p1, sim1, seed=1001)
  call tegarch_filter(sim1%y, p1, refiltered)
  call assert_vector(refiltered%lambda, sim1%lambda, 2.0e-13_dp, "simulate/filter one-component")

  call tegarch_simulate(300, p2, sim2, seed=1002)
  call tegarch_filter(sim2%y, p2, refiltered)
  call assert_vector(refiltered%lambda, sim2%lambda, 3.0e-13_dp, "simulate/filter two-component")

  print '(a)', "Recursion, simulation, and likelihood tests passed."

contains

  subroutine assert_close(actual, expected, tolerance, label)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: label

    if (abs(actual - expected) > tolerance) then
      print '(a,1x,es24.16)', trim(label)//" actual:", actual
      print '(a,1x,es24.16)', trim(label)//" expected:", expected
      error stop "assert_close failed"
    end if
  end subroutine assert_close

  subroutine assert_vector(actual, expected, tolerance, label)
    real(dp), intent(in) :: actual(:), expected(:), tolerance
    character(len=*), intent(in) :: label

    if (size(actual) /= size(expected)) error stop "assert_vector size mismatch"
    if (maxval(abs(actual - expected)) > tolerance) then
      print '(a,1x,es24.16)', trim(label)//" maximum error:", maxval(abs(actual - expected))
      error stop "assert_vector failed"
    end if
  end subroutine assert_vector

end program test_recursions
