program test_interpolation
   use r_stats, only: approx, approxfun, dp, interpolation_constant, interpolation_rule_constant, &
                      interpolation_t, isoreg, isoreg_t, predict_approxfun, smooth_xy_t
   use stats_test_assertions, only: assert_integer_vector, assert_vector_close
   implicit none

   real(dp), parameter :: x(5) = [4.0_dp, 1.0_dp, 2.0_dp, 2.0_dp, 6.0_dp]
   real(dp), parameter :: xout(7) = [0.0_dp, 1.0_dp, 1.5_dp, 2.0_dp, 3.0_dp, 6.0_dp, 7.0_dp]
   real(dp), parameter :: y(5) = [8.0_dp, 1.0_dp, 4.0_dp, 6.0_dp, 12.0_dp]
   type(interpolation_t) :: interpolator
   type(isoreg_t) :: isotonic
   type(smooth_xy_t) :: result

   result = approx(x, y, xout, left_rule=interpolation_rule_constant, &
                   right_rule=interpolation_rule_constant)
   call assert_vector_close(result%y, [1.0_dp, 1.0_dp, 3.0_dp, 5.0_dp, 6.5_dp, 12.0_dp, 12.0_dp], &
                            "R linear approx")

   interpolator = approxfun(x, y, method=interpolation_constant, &
                            left_rule=interpolation_rule_constant, &
                            right_rule=interpolation_rule_constant)
   call assert_vector_close(interpolator%x, [1.0_dp, 2.0_dp, 4.0_dp, 6.0_dp], &
                            "approxfun sorted distinct knots")
   call assert_vector_close(interpolator%y, [1.0_dp, 5.0_dp, 8.0_dp, 12.0_dp], &
                            "approxfun averaged ties")
   call assert_vector_close(predict_approxfun(interpolator, xout), &
                            [1.0_dp, 1.0_dp, 1.0_dp, 5.0_dp, 5.0_dp, 12.0_dp, 12.0_dp], &
                            "R constant approx")

   interpolator = approxfun(x, y, method=interpolation_constant, &
                            left_rule=interpolation_rule_constant, &
                            right_rule=interpolation_rule_constant, constant_mix=0.25_dp)
   call assert_vector_close(predict_approxfun(interpolator, xout), &
                            [1.0_dp, 1.0_dp, 2.0_dp, 5.0_dp, 5.75_dp, 12.0_dp, 12.0_dp], &
                            "R mixed constant approx")

   isotonic = isoreg([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp], &
                     [3.0_dp, 2.0_dp, 4.0_dp, 5.0_dp, 1.0_dp, 6.0_dp, 8.0_dp, 7.0_dp])
   call assert_vector_close(isotonic%fitted, &
                            [2.5_dp, 2.5_dp, 10.0_dp/3.0_dp, 10.0_dp/3.0_dp, &
                             10.0_dp/3.0_dp, 6.0_dp, 7.5_dp, 7.5_dp], "R isotonic fit")
   call assert_integer_vector(isotonic%knots, [2, 5, 6, 8], "isotonic block endpoints")

   print *, "test_interpolation: PASS"

end program test_interpolation
