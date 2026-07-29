! RiskPortfolios Fortran, derived from RiskPortfolios 2.1.7.
! Original code Copyright (C) 2013-2021 David Ardia.
! Original authors: David Ardia, Kris Boudt, Jean-Philippe Gagnon-Fleury.
! SPDX-License-Identifier: GPL-2.0-or-later
program test_reference_values
   use riskportfolios
   implicit none

   real(dp) :: rets(8, 3), expected(3, 3), tolerance
   real(dp), allocatable :: mu(:), semidev(:), sigma(:, :)
   integer :: info, failures

   failures = 0
   tolerance = 5.0e-13_dp
   rets = reshape([ &
      0.010_dp, 0.005_dp, -0.008_dp, 0.013_dp, -0.004_dp, 0.009_dp, -0.006_dp, 0.011_dp, &
     -0.020_dp, 0.012_dp,  0.007_dp,-0.003_dp,  0.018_dp, 0.002_dp, -0.009_dp, 0.015_dp, &
      0.015_dp,-0.004_dp,  0.020_dp, 0.006_dp, -0.011_dp, 0.014_dp,  0.003_dp,-0.007_dp], [8, 3])

   call mean_estimation(rets, mu, MEAN_NAIVE, info=info)
   call check(maxval(abs(mu - [0.00375_dp, 0.00275_dp, 0.00450_dp])) < tolerance, &
      'naive mean reference')

   call mean_estimation(rets, mu, MEAN_EWMA, lambda=0.9_dp, info=info)
   call check(maxval(abs(mu - [0.0037257826717931371_dp, 0.0036202452891254950_dp, &
      0.0034793588442905985_dp])) < 5.0e-12_dp, 'EWMA mean reference')

   call mean_estimation(rets, mu, MEAN_BAYES_STEIN, info=info)
   call check(maxval(abs(mu - [0.0036929668839812_dp, 0.0036827319744295_dp, &
      0.0037006430661449_dp])) < tolerance, 'Bayes-Stein mean reference')

   call semideviation_estimation(rets, semidev, SEMIDEV_NAIVE, info=info)
   call check(maxval(abs(semidev - [0.0098858063235462319_dp, 0.013126785592825076_dp, &
      0.010571187255932988_dp])) < 5.0e-13_dp, 'naive semideviation reference')

   expected = reshape([ &
      7.13571428571428511e-05_dp, -2.17857142857142814e-05_dp, -1.71428571428571512e-06_dp, &
     -2.17857142857142814e-05_dp,  1.67928571428571449e-04_dp, -8.95714285714285628e-05_dp, &
     -1.71428571428571512e-06_dp, -8.95714285714285628e-05_dp,  1.27142857142857108e-04_dp], [3, 3])
   call covariance_estimation(rets, sigma, COV_NAIVE, info=info)
   call check(maxval(abs(sigma - expected)) < tolerance, 'naive covariance reference')

   expected = reshape([ &
      9.41734617691698659e-05_dp, -1.35351763672025348e-05_dp, -3.36124407519759246e-06_dp, &
     -1.35351763672025348e-05_dp,  2.07431925442429286e-04_dp, -1.13338012268602439e-04_dp, &
     -3.36124407519759246e-06_dp, -1.13338012268602439e-04_dp,  1.63063494847789746e-04_dp], [3, 3])
   call covariance_estimation(rets, sigma, COV_EWMA, lambda=0.9_dp, info=info)
   call check(maxval(abs(sigma - expected)) < tolerance, 'EWMA covariance reference')

   expected = reshape([ &
      7.13571428571428511e-05_dp, -3.02863042491570752e-05_dp, -2.63529942609987423e-05_dp, &
     -3.02863042491570752e-05_dp,  1.67928571428571421e-04_dp, -4.04271883537212488e-05_dp, &
     -2.63529942609987423e-05_dp, -4.04271883537212488e-05_dp,  1.27142857142857108e-04_dp], [3, 3])
   call covariance_estimation(rets, sigma, COV_CONSTANT, info=info)
   call check(maxval(abs(sigma - expected)) < tolerance, 'constant-correlation covariance reference')

   expected = reshape([ &
      6.24374999999999930e-05_dp, -8.83650225314606028e-06_dp,  1.97133182725117197e-06_dp, &
     -8.83650225314606028e-06_dp,  1.46937499999999987e-04_dp, -5.24814984198689161e-05_dp, &
      1.97133182725117197e-06_dp, -5.24814984198689161e-05_dp,  1.11249999999999967e-04_dp], [3, 3])
   call covariance_estimation(rets, sigma, COV_LEDOIT_WOLF, info=info)
   call check(maxval(abs(sigma - expected)) < tolerance, 'Ledoit-Wolf covariance reference')

   expected = reshape([ &
      6.24374999999999930e-05_dp, -2.65005162180124353e-05_dp, -2.30588699783738953e-05_dp, &
     -2.65005162180124353e-05_dp,  1.46937500000000014e-04_dp, -3.53737898095060808e-05_dp, &
     -2.30588699783738953e-05_dp, -3.53737898095060808e-05_dp,  1.11249999999999980e-04_dp], [3, 3])
   call covariance_estimation(rets, sigma, COV_COR_SHRINKAGE, info=info)
   call check(maxval(abs(sigma - expected)) < tolerance, 'correlation shrinkage reference')

   expected = reshape([ &
      6.24374999999999930e-05_dp, -9.25937674244965627e-06_dp, -7.28606694487842309e-07_dp, &
     -9.25937674244965627e-06_dp,  1.46937500000000014e-04_dp, -3.80696997869897391e-05_dp, &
     -7.28606694487842309e-07_dp, -3.80696997869897391e-05_dp,  1.11249999999999980e-04_dp], [3, 3])
   call covariance_estimation(rets, sigma, COV_DIAGONAL, info=info)
   call check(maxval(abs(sigma - expected)) < tolerance, 'diagonal shrinkage reference')

   expected = reshape([ &
      9.17021798892733015e-05_dp, -6.50873436536095555e-06_dp, -5.12162704159551007e-07_dp, &
     -6.50873436536095555e-06_dp,  1.20554012223594680e-04_dp, -2.67605012923365235e-05_dp, &
     -5.12162704159551007e-07_dp, -2.67605012923365235e-05_dp,  1.08368807887132020e-04_dp], [3, 3])
   call covariance_estimation(rets, sigma, COV_ONE_PARAMETER, info=info)
   call check(maxval(abs(sigma - expected)) < tolerance, 'one-parameter shrinkage reference')

   expected = reshape([ &
      6.24374999999999930e-05_dp, -8.79084227261255538e-06_dp,  1.98683162957670878e-06_dp, &
     -8.79084227261255538e-06_dp,  1.46937500000000014e-04_dp, -5.23658816548978152e-05_dp, &
      1.98683162957670878e-06_dp, -5.23658816548978152e-05_dp,  1.11249999999999980e-04_dp], [3, 3])
   call covariance_estimation(rets, sigma, COV_LARGE, info=info)
   call check(maxval(abs(sigma - expected)) < tolerance, 'large-dimensional covariance reference')

   expected = reshape([ &
      7.33606924842079042e-05_dp, -1.99013282484294325e-05_dp,  1.95778979999318009e-07_dp, &
     -1.99013282484294325e-05_dp,  1.70055671035224576e-04_dp, -8.77737651159699446e-05_dp, &
      1.95778979999318009e-07_dp, -8.77737651159699446e-05_dp,  1.29217776987420902e-04_dp], [3, 3])
   call covariance_estimation(rets, sigma, COV_BAYES_STEIN, info=info)
   call check(maxval(abs(sigma - expected)) < tolerance, 'Bayes-Stein covariance reference')

   if (failures > 0) then
      write(*, '(a,1x,i0)') 'FAILED reference tests:', failures
      error stop 1
   end if
   write(*, '(a)') 'All independent reference-value tests passed.'

contains

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         failures = failures + 1
         write(*, '(a)') 'FAIL: ' // label
      end if
   end subroutine check

end program test_reference_values
