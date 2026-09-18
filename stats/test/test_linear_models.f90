program test_linear_models
   use stats_test_assertions, only: assert_close, assert_equal, assert_matrix_close, &
                                      assert_true, assert_vector_close
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
   use r_stats, only: dp, lm_confint, lm_fit, lm_fit_t, predict_lm
   implicit none

   real(dp) :: predictors(5, 2), response(5), design(5, 3), identity(3, 3)
   real(dp) :: new_data(2, 2), no_intercept_x(4, 1)
   real(dp) :: rank_x(6, 2), weighted_x(6, 2), weighted_y(6), weights(6)
   real(dp), allocatable :: actual(:), intervals(:, :)
   type(lm_fit_t) :: fit

   predictors(:, 1) = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
   predictors(:, 2) = [1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp]
   response = 2.0_dp + 3.0_dp*predictors(:, 1) - 4.0_dp*predictors(:, 2)

   fit = lm_fit(response, predictors)
   call assert_true(fit%has_intercept, "linear-model intercept flag")
   call assert_equal(fit%df, 2, "linear-model residual degrees of freedom")
   call assert_vector_close(fit%coef, [2.0_dp, 3.0_dp, -4.0_dp], "linear-model coefficients")
   call assert_vector_close(fit%fitted, response, "linear-model fitted values")
   call assert_vector_close(fit%resid, spread(0.0_dp, 1, 5), "linear-model residuals")
   call assert_close(fit%r_squared, 1.0_dp, "linear-model R-squared")

   design(:, 1) = 1.0_dp
   design(:, 2:) = predictors
   identity = 0.0_dp
   identity(1, 1) = 1.0_dp
   identity(2, 2) = 1.0_dp
   identity(3, 3) = 1.0_dp
   call assert_matrix_close(matmul(matmul(transpose(design), design), fit%cov_unscaled), identity, &
                            "linear-model unscaled covariance")

   new_data(1, :) = [5.0_dp, 0.0_dp]
   new_data(2, :) = [5.0_dp, 1.0_dp]
   actual = predict_lm(fit, new_data)
   call assert_vector_close(actual, [17.0_dp, 13.0_dp], "linear-model predictions")

   no_intercept_x(:, 1) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
   fit = lm_fit(2.0_dp*no_intercept_x(:, 1), no_intercept_x, intercept=.false.)
   call assert_true(.not. fit%has_intercept, "no-intercept flag")
   call assert_vector_close(fit%coef, [2.0_dp], "no-intercept coefficient")

   intervals = lm_confint(fit, 0.90_dp)
   call assert_true(all(shape(intervals) == [1, 2]), "confidence-interval shape")
   call assert_true(intervals(1, 1) <= fit%coef(1), "confidence-interval lower bound")
   call assert_true(intervals(1, 2) >= fit%coef(1), "confidence-interval upper bound")

   weighted_x(:, 1) = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
   weighted_x(:, 2) = [1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp]
   weighted_y = [-1.8_dp, 4.9_dp, 4.4_dp, 10.2_dp, 10.1_dp, 15.3_dp]
   weights = [1.0_dp, 2.0_dp, 1.0_dp, 3.0_dp, 2.0_dp, 1.0_dp]
   fit = lm_fit(weighted_y, weighted_x, weights=weights)
   call assert_true(fit%weighted, "weighted linear-model flag")
   call assert_equal(fit%rank, 3, "weighted linear-model rank")
   call assert_equal(fit%df, 3, "weighted linear-model residual degrees of freedom")
   call assert_vector_close(fit%coef, &
      [1.86044776119403_dp, 2.78358208955224_dp, -3.11940298507463_dp], &
      "weighted linear-model coefficients")
   call assert_close(sum(weights*fit%resid**2), 0.762313432835819_dp, &
                     "weighted linear-model RSS")
   call assert_matrix_close(fit%cov_unscaled, reshape( &
      [0.485074626865671_dp, -0.119402985074627_dp, -0.186567164179104_dp, &
       -0.119402985074627_dp, 0.0447761194029851_dp, 0.00746268656716411_dp, &
       -0.186567164179104_dp, 0.00746268656716413_dp, 0.417910447761194_dp], &
      [3, 3]), "weighted linear-model covariance")

   rank_x(:, 1) = weighted_x(:, 1)
   rank_x(:, 2) = 2.0_dp*rank_x(:, 1)
   fit = lm_fit(weighted_y, rank_x)
   call assert_equal(fit%rank, 2, "rank-deficient linear-model rank")
   call assert_equal(fit%df, 4, "rank-deficient linear-model degrees of freedom")
   call assert_true(fit%aliased(3), "rank-deficient aliased coefficient")
   call assert_vector_close(fit%coef, &
      [-0.452380952380950_dp, 3.05428571428571_dp, 0.0_dp], &
      "rank-deficient linear-model coefficients")
   call assert_close(sum(fit%resid**2), 13.8967619047619_dp, &
                     "rank-deficient linear-model RSS")
   intervals = lm_confint(fit)
   call assert_true(all(ieee_is_nan(intervals(3, :))), "aliased confidence interval")

   print *, "test_linear_models: PASS"

end program test_linear_models
