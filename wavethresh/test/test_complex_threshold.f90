program test_complex_threshold
   use wavethresh, only : complex_threshold_parameters_t, complex_threshold_result_t, cthresh, cwd_t, dp
   use wavethresh, only : complex_wd, complex_wr, find_parameters, lina_mayrand_31, wt_filter_t
   implicit none

   integer, parameter :: reference_indices(6) = [1, 2, 8, 16, 24, 32]
   complex(dp), parameter :: reference_values(6) = [ &
      cmplx(0.47079558459643051_dp, -0.0045714117757979893_dp, dp), &
      cmplx(0.48840775455349694_dp, -0.0054877082511772951_dp, dp), &
      cmplx(1.2937346742596543_dp, -0.0042785067423359857_dp, dp), &
      cmplx(0.025117672498552141_dp, -0.010368665268024062_dp, dp), &
      cmplx(-0.79147933851431929_dp, -0.016273768711619098_dp, dp), &
      cmplx(0.49009761880525882_dp, 0.0043778227174188274_dp, dp) ]
   complex(dp), parameter :: level2_reference(4) = [ &
      cmplx(0.7675033319219261_dp, -0.1279742776698391_dp, dp), &
      cmplx(-0.33517498016037561_dp, -0.34445708783125734_dp, dp), &
      cmplx(0.046234358935228065_dp, -0.40141745272401919_dp, dp), &
      cmplx(0.78263151325592628_dp, -0.26208082722433873_dp, dp) ]
   type(wt_filter_t) :: filter
   type(cwd_t) :: transform
   type(complex_threshold_result_t) :: result
   type(complex_threshold_parameters_t) :: parameters
   real(dp) :: data(32)
   real(dp) :: sigma(0:4, 2, 2)
   complex(dp), allocatable :: reconstructed(:)
   integer :: i

   do i = 1, size(data)
      data(i) = sin(real(i, dp) * 0.21_dp) + 0.3_dp * cos(real(i, dp) * 0.73_dp) + &
         real(i, dp) / 200.0_dp
   end do
   filter = lina_mayrand_31()
   call check(filter%ok, "Lina-Mayrand 3.1 filter status")
   transform = complex_wd(data, filter)
   call check(transform%ok, "complex DWT status")
   call check_close(maxval(abs(transform%detail(2)%values - level2_reference)), 0.0_dp, 3.0e-14_dp, &
      "complex DWT upstream coefficients")
   reconstructed = complex_wr(transform)
   call check_close(maxval(abs(reconstructed - cmplx(data, 0.0_dp, dp))), 0.0_dp, 6.0e-7_dp, &
      "complex DWT reconstruction")

   result = cthresh(data, j0=2, rule="hard", policy="mws")
   call check(result%ok, "cthresh mws status")
   call check_close(maxval(abs(result%estimate(reference_indices) - reference_values)), 0.0_dp, 3.0e-10_dp, &
      "cthresh upstream mws estimate")

   sigma = 0.0_dp
   sigma(:, 1, 1) = 0.4_dp
   sigma(:, 2, 2) = 0.3_dp
   sigma(:, 1, 2) = 0.05_dp
   sigma(:, 2, 1) = 0.05_dp
   parameters = find_parameters(transform, sigma, 2, 0.02_dp)
   call check(parameters%ok, "find.parameters status")
   call check(all(parameters%nonzero_probability(2:) >= 0.02_dp), "find.parameters probability lower bound")
   call check(all(parameters%nonzero_probability(2:) <= 0.98_dp), "find.parameters probability upper bound")
   call check(all(parameters%signal_covariance(2:, 1, 1) > 0.0_dp), "find.parameters real variances")
   call check(all(parameters%signal_covariance(2:, 2, 2) > 0.0_dp), "find.parameters imaginary variances")

   result = cthresh(data, j0=2, rule="mean", policy="ebayes", tolerance=0.02_dp)
   call check(result%ok, "cthresh ebayes status")
   call check(all(abs(result%estimate) < huge(1.0_dp)), "cthresh ebayes finite estimate")
   print *, "test_complex_threshold: PASS"

contains

   subroutine check(condition, label)
      !! Stops the test program when a logical assertion fails.
      logical, intent(in) :: condition !! Assertion value.
      character(len=*), intent(in) :: label !! Human-readable assertion label.
      if (.not. condition) then
         print *, "FAIL: ", trim(label)
         error stop 1
      end if
   end subroutine check

   subroutine check_close(value, target, tolerance, label)
      !! Stops the test program when a scalar differs from its target beyond tolerance.
      real(dp), intent(in) :: value !! Computed value.
      real(dp), intent(in) :: target !! Expected value.
      real(dp), intent(in) :: tolerance !! Maximum permitted absolute error.
      character(len=*), intent(in) :: label !! Human-readable assertion label.
      if (abs(value - target) > tolerance) then
         print *, "FAIL: ", trim(label), value, target
         error stop 1
      end if
   end subroutine check_close

end program test_complex_threshold
