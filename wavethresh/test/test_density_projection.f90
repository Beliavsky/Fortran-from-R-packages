! SPDX-License-Identifier: GPL-2.0-or-later
program test_density_projection
   use wavethresh, only : dp, density_grid_t, density_projection_t, density_wavelet_t, scaling_function_t
   use wavethresh, only : chires5, chires6, cwavde, dencvwd, denproj, denwd, evaluate_density
   implicit none

   type(density_projection_t) :: projection
   type(density_projection_t) :: covariance_projection
   type(density_projection_t) :: invalid_projection
   type(density_grid_t) :: grid
   type(density_wavelet_t) :: density_wavelet
   type(scaling_function_t) :: scaling_function
   type(scaling_function_t) :: wavelet_function
   real(dp), parameter :: root_two = sqrt(2.0_dp)
   real(dp) :: observations(4)
   real(dp) :: locations(2)

   observations = [0.125_dp, 0.375_dp, 0.625_dp, 0.875_dp]
   locations = [0.25_dp, 0.75_dp]

   projection = chires5(observations, resolution_level=1, filter_number=1.0_dp, family="DaubExPhase")
   call assert_true(projection%ok, "Chires5 Haar status")
   call assert_true(projection%k_min == 0 .and. projection%k_max == 1, "Chires5 Haar translation bounds")
   call assert_close(projection%coefficients, [0.5_dp * root_two, 0.5_dp * root_two], 1.0e-13_dp, &
      "Chires5 Haar coefficients")
   call assert_true(.not. projection%has_covariance, "Chires5 covariance flag")

   grid = evaluate_density(projection, locations)
   call assert_true(grid%ok, "density evaluation status")
   call assert_close(grid%x, locations, 0.0_dp, "density evaluation grid")
   call assert_close(grid%y, [1.0_dp, 1.0_dp], 1.0e-13_dp, "Haar uniform density")

   covariance_projection = chires6(observations, resolution_level=1, filter_number=1.0_dp, family="DaubExPhase")
   call assert_true(covariance_projection%ok, "Chires6 Haar status")
   call assert_true(covariance_projection%has_covariance, "Chires6 covariance flag")
   call assert_true(all(shape(covariance_projection%covariance) == [2, 1]), "Chires6 covariance shape")
   call assert_close(covariance_projection%covariance(:, 1), [0.25_dp, 0.25_dp], 1.0e-13_dp, &
      "Chires6 Haar covariance diagonal")

   projection = denproj(observations, 1, filter_number=1.0_dp, family="DaubExPhase", covariance=.true.)
   call assert_true(projection%ok .and. projection%has_covariance, "denproj covariance dispatch")
   call assert_close(projection%coefficients, covariance_projection%coefficients, 0.0_dp, "denproj coefficients")

   projection = chires5([0.13_dp, 0.37_dp, 0.61_dp, 0.88_dp], resolution_level=2, &
      filter_number=2.0_dp, family="DaubExPhase", n_iterations=24)
   call assert_true(projection%ok, "Chires5 Daubechies-2 status")
   call assert_true(projection%k_min == -2 .and. projection%k_max == 3, "Chires5 R translation bounds")
   call assert_close(projection%coefficients, &
      [0.017099165665497962_dp, 0.052626208504778976_dp, 0.500328788835714855_dp, &
       0.521128257887379487_dp, 0.448705732017985082_dp, 0.460111847119365924_dp], &
      2.0e-14_dp, "Chires5 upstream R coefficients")
   call assert_scalar_close(sum(projection%coefficients) / sqrt(projection%primary_resolution), 1.0_dp, &
      2.0e-7_dp, "Daubechies-2 coefficient integral")
   density_wavelet = denwd(projection)
   call assert_true(density_wavelet%ok .and. density_wavelet%nlevels == 2, "denwd status and levels")
   call assert_true(all(density_wavelet%scaling_first == [-2, -2, -2]), "denwd scaling first indices")
   call assert_true(all(density_wavelet%scaling_last == [0, 1, 3]), "denwd scaling last indices")
   call assert_true(all(density_wavelet%detail_first == [-1, -1]), "denwd detail first indices")
   call assert_true(all(density_wavelet%detail_last == [1, 2]), "denwd detail last indices")
   call assert_close(density_wavelet%scaling(0)%values, &
      [-0.013218545076352896_dp, 0.162911696891345770_dp, 0.850306848201648036_dp], &
      3.0e-14_dp, "denwd upstream R coarse scaling coefficients")
   call assert_close(density_wavelet%detail(0)%values, &
      [-0.049332281827228583_dp, 0.289220240438053555_dp, -0.227839033293969040_dp], &
      3.0e-14_dp, "denwd upstream R coarse detail coefficients")
   call assert_close(density_wavelet%scaling(1)%values, &
      [-0.0029776593837171502_dp, 0.0969876152616868520_dp, 0.7186043173791327732_dp, &
       0.6015992891386214847_dp], 3.0e-14_dp, "denwd upstream R fine scaling coefficients")
   call assert_close(density_wavelet%detail(1)%values, &
      [-0.011112776107738156_dp, 0.152838930739656431_dp, -0.028422310794173483_dp, &
       -0.161198043691628434_dp], 3.0e-14_dp, "denwd upstream R fine detail coefficients")
   covariance_projection = chires6([0.13_dp, 0.37_dp, 0.61_dp, 0.88_dp], resolution_level=2, &
      filter_number=2.0_dp, family="DaubExPhase", n_iterations=24)
   density_wavelet = dencvwd(covariance_projection)
   call assert_true(density_wavelet%ok .and. density_wavelet%nlevels == 2, "dencvwd status and levels")
   call assert_close(density_wavelet%detail(0)%values, &
      [0.0040027637886049135_dp, 0.1850699665599081500_dp, 0.0145985959236172337_dp], &
      5.0e-14_dp, "dencvwd upstream R coarse variances")
   call assert_close(density_wavelet%detail(1)%values, &
      [0.00021852939385404906_dp, 0.18714353440085357283_dp, 0.17401344697698684083_dp, &
       0.01429100194550118018_dp], 5.0e-14_dp, "dencvwd upstream R fine variances")

   invalid_projection = chires5([real(dp) ::], resolution_level=1)
   call assert_true(.not. invalid_projection%ok, "empty density sample rejection")
   invalid_projection = chires5(observations, tau=0.0_dp, resolution_level=1)
   call assert_true(.not. invalid_projection%ok, "nonpositive tau rejection")

   scaling_function%x = [0.0_dp, 0.5_dp, 1.0_dp, 1.5_dp]
   scaling_function%y = [0.0_dp, 1.0_dp, 0.5_dp, 0.0_dp]
   scaling_function%ok = .true.
   wavelet_function%x = [-0.5_dp, 0.0_dp, 0.5_dp, 1.0_dp]
   wavelet_function%y = [0.0_dp, 1.0_dp, -1.0_dp, 0.0_dp]
   wavelet_function%ok = .true.
   grid = cwavde([0.1_dp, 0.35_dp, 0.6_dp, 0.9_dp], 2, scaling_function, wavelet_function, &
      threshold=0.05_dp, n_output=9, filter_number=1.0_dp, family="DaubExPhase")
   call assert_true(grid%ok, "CWavDE status")
   call assert_close(grid%x, [-1.0_dp, -0.625_dp, -0.25_dp, 0.125_dp, 0.5_dp, 0.875_dp, &
      1.25_dp, 1.625_dp, 2.0_dp], 2.0e-15_dp, "CWavDE upstream R grid")
   call assert_close(grid%y, [0.0_dp, 0.17083333333333334_dp, 0.34166666666666667_dp, &
      0.91875000000000018_dp, -0.16805555555555590_dp, 0.40486111111111101_dp, &
      0.38333333333333330_dp, 0.0_dp, 0.0_dp], 3.0e-14_dp, "CWavDE upstream R density")
   call assert_true(all(grid%scaling_indices == [-1, 0, 1]), "CWavDE scaling translations")
   call assert_true(all(grid%wavelet_k_min == [-1, -1]), "CWavDE wavelet lower translations")
   call assert_true(all(grid%wavelet_k_max == [1, 1]), "CWavDE wavelet upper translations")

   print *, "test_density_projection: PASS"

contains

   subroutine assert_true(condition, label)
      logical, intent(in) :: condition !! Condition that must be true.
      character(len=*), intent(in) :: label !! Diagnostic label printed on failure.

      if (.not. condition) then
         print *, "FAIL: ", trim(label)
         error stop 1
      end if
   end subroutine assert_true

   subroutine assert_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual(:) !! Values produced by the translated procedure.
      real(dp), intent(in) :: expected(:) !! Reference values required by the test.
      real(dp), intent(in) :: tolerance !! Maximum permitted absolute error.
      character(len=*), intent(in) :: label !! Diagnostic label printed on failure.

      if (size(actual) /= size(expected)) then
         print *, "FAIL: ", trim(label), " shape"
         error stop 1
      end if
      if (any(abs(actual - expected) > tolerance)) then
         print *, "FAIL: ", trim(label), maxval(abs(actual - expected))
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_scalar_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual !! Scalar value produced by the translated procedure.
      real(dp), intent(in) :: expected !! Scalar reference value required by the test.
      real(dp), intent(in) :: tolerance !! Maximum permitted absolute error.
      character(len=*), intent(in) :: label !! Diagnostic label printed on failure.

      if (abs(actual - expected) > tolerance) then
         print *, "FAIL: ", trim(label), actual, expected
         error stop 1
      end if
   end subroutine assert_scalar_close

end program test_density_projection
