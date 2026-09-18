program test_spline_interpolation
   use r_stats, only: cubic_spline_t, dp, predict_splinefun, smooth_xy_t, spline, &
                      spline_method_natural, splinefun
   use stats_test_assertions, only: assert_vector_close
   implicit none

   real(dp), parameter :: x(5) = [0.0_dp, 1.0_dp, 2.0_dp, 4.0_dp, 5.0_dp]
   real(dp), parameter :: xout(6) = [-1.0_dp, 0.5_dp, 1.5_dp, 3.0_dp, 5.0_dp, 6.0_dp]
   real(dp), parameter :: y(5) = [0.0_dp, 1.0_dp, 0.0_dp, 2.0_dp, 1.0_dp]
   type(cubic_spline_t) :: interpolator
   type(smooth_xy_t) :: result

   interpolator = splinefun(x, y)
   call assert_vector_close(predict_splinefun(interpolator, xout), &
      [-6.04098360655737743_dp, 0.97387295081967207_dp, 0.45850409836065564_dp, &
       0.78483606557377028_dp, 1.0_dp, -4.34426229508196649_dp], "R FMM spline", 2.0e-11_dp)
   call assert_vector_close(predict_splinefun(interpolator, xout, derivative=1), &
      [9.39480874316939918_dp, 0.89583333333333304_dp, -1.23941256830601110_dp, &
       1.46243169398907114_dp, -2.83879781420765021_dp, -8.18306010928961669_dp], &
      "R FMM first derivative", 2.0e-11_dp)
   call assert_vector_close(predict_splinefun(interpolator, xout, derivative=2), &
      [-7.54098360655737743_dp, -3.79098360655737743_dp, 0.33196721311475397_dp, &
       0.43032786885245944_dp, -4.34426229508196649_dp, -6.34426229508196649_dp], &
      "R FMM second derivative", 2.0e-11_dp)

   interpolator = splinefun(x, y, spline_method_natural)
   call assert_vector_close(predict_splinefun(interpolator, xout), &
      [-1.65573770491803285_dp, 0.74590163934426235_dp, 0.51229508196721318_dp, &
       0.87704918032786860_dp, 1.0_dp, -0.54098360655737698_dp], &
      "R natural spline", 2.0e-11_dp)
   call assert_vector_close(predict_splinefun(interpolator, xout, derivative=1), &
      [1.6557377049180328_dp, 1.1639344262295082_dp, -1.3196721311475410_dp, &
       1.5819672131147540_dp, -1.5409836065573770_dp, -1.5409836065573770_dp], &
      "R natural first derivative", 2.0e-11_dp)
   call assert_vector_close(predict_splinefun(interpolator, xout, derivative=2), &
      [0.0_dp, -1.967213114754098546_dp, -0.098360655737705027_dp, &
       0.245901639344262346_dp, 0.0_dp, 0.0_dp], "R natural second derivative", 2.0e-11_dp)

   result = spline([2.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 3.0_dp], &
                   [4.0_dp, 0.0_dp, 1.0_dp, 3.0_dp, 9.0_dp], &
                   [0.5_dp, 1.0_dp, 1.5_dp], spline_method_natural)
   call assert_vector_close(result%y, [1.075_dp, 2.0_dp, 2.775_dp], &
                            "R natural spline with averaged ties", 2.0e-12_dp)

   print *, "test_spline_interpolation: PASS"

end program test_spline_interpolation
