! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2011-2025 Genaro Sucarrat
! Copyright (C) 2026 contributors to the Modern Fortran translation
!
! This file is part of betategarch-modern-fortran.
! It is free software: you can redistribute it and/or modify it under the
! terms of the GNU General Public License version 2 only.

program test_skew_t
  use betategarch, only : dp, set_random_seed, skew_t_random, skew_t_pdf, skew_t_logpdf, &
    skew_t_mean, skew_t_variance, skew_t_skewness, skew_t_kurtosis
  implicit none

  integer, parameter :: n = 200000
  real(dp), allocatable :: x(:), x2(:)
  real(dp) :: sample_mean, sample_variance

  call assert_close(skew_t_logpdf(-2.0_dp, 8.0_dp, 1.2_dp, 0.8_dp), &
    -2.060140911114238_dp, 2.0e-14_dp, "log density -2")
  call assert_close(skew_t_pdf(0.5_dp, 8.0_dp, 1.2_dp, 0.8_dp), &
    0.27058183397714936_dp, 2.0e-14_dp, "density 0.5")
  call assert_close(skew_t_mean(8.0_dp, 0.8_dp), -0.3977475644174328_dp, &
    2.0e-14_dp, "mean")
  call assert_close(skew_t_variance(8.0_dp, 0.8_dp), 1.4451302083333317_dp, &
    3.0e-14_dp, "variance")
  call assert_close(skew_t_skewness(8.0_dp, 0.8_dp), -0.5848399454966141_dp, &
    5.0e-14_dp, "skewness")
  call assert_close(skew_t_kurtosis(8.0_dp, 0.8_dp), 4.872060165576732_dp, &
    8.0e-14_dp, "kurtosis")

  allocate(x(n), x2(20))
  call set_random_seed(12345)
  call skew_t_random(x, 8.0_dp, 0.8_dp)
  sample_mean = sum(x)/real(n, dp)
  sample_variance = sum((x - sample_mean)**2)/real(n - 1, dp)
  call assert_close(sample_mean, skew_t_mean(8.0_dp, 0.8_dp), 0.02_dp, "random mean")
  call assert_close(sample_variance, skew_t_variance(8.0_dp, 0.8_dp), 0.05_dp, "random variance")

  call set_random_seed(2468)
  call skew_t_random(x2, 7.0_dp, 1.3_dp)
  x(1:20) = x2
  call set_random_seed(2468)
  call skew_t_random(x2, 7.0_dp, 1.3_dp)
  if (maxval(abs(x(1:20) - x2)) > 0.0_dp) error stop "seed reproducibility failed"

  print '(a)', "Skew-t tests passed."

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

end program test_skew_t
