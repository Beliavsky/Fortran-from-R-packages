program test_numerical
   use, intrinsic :: ieee_arithmetic, only: ieee_negative_inf, ieee_positive_inf, ieee_value
   use r_stats, only: dp, integrate, integrate_result_t, optimize, optimize_result_t, &
                      root_extend_both, uniroot, uniroot_result_t
   use stats_test_assertions, only: assert_close, assert_true
   implicit none

   type(integrate_result_t) :: quadrature
   type(optimize_result_t) :: optimum
   type(uniroot_result_t) :: root

   root = uniroot(root_objective, 0.0_dp, 1.0_dp, tolerance=1.0e-12_dp)
   call assert_true(root%converged, "uniroot convergence")
   call assert_close(root%root, 0.7390851332151559_dp, "R uniroot root", 2.0e-12_dp)
   call assert_close(root%function_value, 0.0_dp, "R uniroot function value", 2.0e-12_dp)

   root = uniroot(distant_root_objective, 0.0_dp, 1.0_dp, tolerance=1.0e-12_dp, &
                  extend_interval=root_extend_both)
   call assert_true(root%converged, "extended uniroot convergence")
   call assert_close(root%root, 10.0_dp, "R extended uniroot root", 2.0e-12_dp)

   optimum = optimize(quadratic_objective, -4.0_dp, 6.0_dp, tolerance=1.0e-10_dp)
   call assert_true(optimum%converged, "optimize minimum convergence")
   call assert_close(optimum%location, 2.25_dp, "R optimize minimum", 8.0e-10_dp)
   call assert_close(optimum%objective, 1.5_dp, "R optimize minimum objective", 1.0e-12_dp)

   optimum = optimize(sine_function, 0.0_dp, 4.0_dp, maximum=.true., tolerance=1.0e-10_dp)
   call assert_true(optimum%converged, "optimize maximum convergence")
   call assert_close(optimum%location, 1.5707963267948966_dp, "R optimize maximum", 8.0e-9_dp)
   call assert_close(optimum%objective, 1.0_dp, "R optimize maximum objective", 1.0e-12_dp)

   quadrature = integrate(gaussian_function, -2.0_dp, 2.0_dp, &
                          abs_tolerance=1.0e-11_dp, rel_tolerance=1.0e-11_dp)
   call assert_true(quadrature%converged, "Gaussian integration convergence")
   call assert_close(quadrature%value, 1.7641627815248431_dp, &
                     "R Gaussian integral", 2.0e-11_dp)

   quadrature = integrate(sine_function, 0.0_dp, acos(-1.0_dp), &
                          abs_tolerance=1.0e-12_dp, rel_tolerance=1.0e-12_dp)
   call assert_true(quadrature%converged, "sine integration convergence")
   call assert_close(quadrature%value, 2.0_dp, "R sine integral", 2.0e-12_dp)

   quadrature = integrate(right_exponential, 0.0_dp, &
                          ieee_value(0.0_dp, ieee_positive_inf), &
                          abs_tolerance=1.0e-10_dp, rel_tolerance=1.0e-10_dp)
   call assert_true(quadrature%converged, "right-tail integration convergence")
   call assert_close(quadrature%value, 1.0_dp, "R right-tail integral", 2.0e-10_dp)

   quadrature = integrate(left_exponential, ieee_value(0.0_dp, ieee_negative_inf), &
                          0.0_dp, abs_tolerance=1.0e-10_dp, rel_tolerance=1.0e-10_dp)
   call assert_true(quadrature%converged, "left-tail integration convergence")
   call assert_close(quadrature%value, 1.0_dp, "R left-tail integral", 2.0e-10_dp)

   quadrature = integrate(gaussian_function, ieee_value(0.0_dp, ieee_negative_inf), &
                          ieee_value(0.0_dp, ieee_positive_inf), &
                          abs_tolerance=1.0e-10_dp, rel_tolerance=1.0e-10_dp)
   call assert_true(quadrature%converged, "whole-line integration convergence")
   call assert_close(quadrature%value, sqrt(acos(-1.0_dp)), &
                     "R whole-line Gaussian integral", 2.0e-10_dp)

   print *, "test_numerical: PASS"

contains

   pure function root_objective(x) result(value)
      !! Evaluates the root-test objective.
      real(dp), intent(in) :: x !! Scalar argument.
      real(dp) :: value !! Difference between cosine and the identity.

      value = cos(x) - x
   end function root_objective

   pure function quadratic_objective(x) result(value)
      !! Evaluates the bounded-minimization fixture.
      real(dp), intent(in) :: x !! Scalar argument.
      real(dp) :: value !! Shifted quadratic objective.

      value = (x - 2.25_dp)**2 + 1.5_dp
   end function quadratic_objective

   pure function distant_root_objective(x) result(value)
      !! Evaluates a root outside the initial test bracket.
      real(dp), intent(in) :: x !! Scalar argument.
      real(dp) :: value !! Difference from the target root ten.

      value = x - 10.0_dp
   end function distant_root_objective

   pure function sine_function(x) result(value)
      !! Evaluates sine for optimization and quadrature fixtures.
      real(dp), intent(in) :: x !! Scalar argument in radians.
      real(dp) :: value !! Sine of `x`.

      value = sin(x)
   end function sine_function

   pure function gaussian_function(x) result(value)
      !! Evaluates an unnormalized standard Gaussian kernel.
      real(dp), intent(in) :: x !! Scalar argument.
      real(dp) :: value !! Exponential of negative squared argument.

      value = exp(-x*x)
   end function gaussian_function

   pure function right_exponential(x) result(value)
      !! Evaluates a decaying exponential on the positive half-line.
      real(dp), intent(in) :: x !! Nonnegative scalar argument.
      real(dp) :: value !! Exponential of negative `x`.

      value = exp(-x)
   end function right_exponential

   pure function left_exponential(x) result(value)
      !! Evaluates a decaying exponential on the negative half-line.
      real(dp), intent(in) :: x !! Nonpositive scalar argument.
      real(dp) :: value !! Exponential of `x`.

      value = exp(x)
   end function left_exponential

end program test_numerical
