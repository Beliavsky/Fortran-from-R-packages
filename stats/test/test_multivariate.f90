program test_multivariate
   use stats_test_assertions, only: assert_integer, assert_matrix_close, assert_vector_close
   use r_stats, only: dp, mahalanobis, prcomp, prcomp_fit_t, scale, scale_columns, scale_result_t
   implicit none

   real(dp), parameter :: tolerance = 8.0e-10_dp
   real(dp) :: covariance(3, 3), data(6, 3), expected_rotation(3, 3)
   real(dp) :: expected_scaled(6, 3), expected_scores(6, 3), inverse_covariance(3, 3)
   real(dp), allocatable :: actual(:), scaled_values(:, :)
   type(prcomp_fit_t) :: fit
   type(scale_result_t) :: scaled

   data(1, :) = [2.0_dp, 0.0_dp, 1.0_dp]
   data(2, :) = [3.0_dp, 1.0_dp, 4.0_dp]
   data(3, :) = [4.0_dp, 1.0_dp, 2.0_dp]
   data(4, :) = [5.0_dp, 3.0_dp, 5.0_dp]
   data(5, :) = [7.0_dp, 4.0_dp, 3.0_dp]
   data(6, :) = [8.0_dp, 6.0_dp, 7.0_dp]

   expected_scaled(1, :) = [-1.2230532340215274_dp, -1.1070186069251191_dp, -1.2344267996967351_dp]
   expected_scaled(2, :) = [-0.79138738671981179_dp, -0.66421116415507153_dp, 0.15430334996209197_dp]
   expected_scaled(3, :) = [-0.3597215394180962_dp, -0.66421116415507153_dp, -0.77151674981045937_dp]
   expected_scaled(4, :) = [0.071944307883619385_dp, 0.22140372138502384_dp, 0.61721339984836765_dp]
   expected_scaled(5, :) = [0.93527600248705056_dp, 0.66421116415507153_dp, -0.30860669992418371_dp]
   expected_scaled(6, :) = [1.3669418497887662_dp, 1.5498260496951668_dp, 1.543033499620919_dp]
   scaled = scale_columns(data)
   call assert_integer(scaled%status, 0, "scaling status")
   call assert_vector_close(scaled%center, &
      [4.833333333333333_dp, 2.5_dp, 3.6666666666666665_dp], "scaling centers", &
      allowed_error=tolerance)
   call assert_vector_close(scaled%scale, &
      [2.3166067138525408_dp, 2.2583179581272428_dp, 2.1602468994692869_dp], &
      "scaling divisors", allowed_error=tolerance)
   call assert_matrix_close(scaled%values, expected_scaled, "scaled values", &
      allowed_error=tolerance)
   scaled_values = scale(data)
   call assert_matrix_close(scaled_values, expected_scaled, "scale generic", &
      allowed_error=tolerance)

   expected_rotation(:, 1) = [0.60314758019344572_dp, 0.60902394462605891_dp, 0.5150755589026782_dp]
   expected_rotation(:, 2) = [-0.5240273905359355_dp, -0.18428296866483651_dp, 0.83152575512017268_dp]
   expected_rotation(:, 3) = [-0.60133874852273517_dp, 0.77144644812986873_dp, -0.20799540185523502_dp]
   expected_scores(:, 1) = [-4.6050128291870518_dp, -1.8476146276595122_dp, &
                             -2.2746181652714235_dp, 1.0918039808821749_dp, &
                              1.8769719680897692_dp, 5.7584696731460454_dp]
   expected_scores(:, 2) = [-0.27195031880655196_dp, 1.5143165873531941_dp, &
                             -0.6727623134230869_dp, 0.92922162407182263_dp, &
                             -1.9661676359052302_dp, 0.46734205670985229_dp]
   expected_scores(:, 3) = [0.32983140543703759_dp, -0.12404710052153389_dp, &
                            -0.30939504533379891_dp, 0.008172896837498278_dp, &
                            -0.0070673483676334525_dp, 0.10250519194842911_dp]
   fit = prcomp(data)
   call assert_integer(fit%status, 0, "PCA status")
   call assert_vector_close(fit%sdev, &
      [3.6788243623557007_dp, 1.2463865716612978_dp, 0.21472111739161795_dp], &
      "PCA standard deviations", allowed_error=tolerance)
   call assert_vector_close(fit%center, &
      [4.833333333333333_dp, 2.5_dp, 3.6666666666666665_dp], "PCA centers", allowed_error=tolerance)
   call assert_vector_close(fit%scale, [1.0_dp, 1.0_dp, 1.0_dp], "PCA scales", &
      allowed_error=tolerance)
   call assert_matrix_close(fit%rotation, expected_rotation, "PCA rotations", &
      allowed_error=tolerance)
   call assert_matrix_close(fit%x, expected_scores, "PCA scores", allowed_error=tolerance)

   fit = prcomp(data, rank=2)
   call assert_integer(size(fit%sdev), 2, "rank-limited PCA component count")
   call assert_matrix_close(fit%rotation, expected_rotation(:, :2), "rank-limited PCA rotations", &
      allowed_error=tolerance)

   covariance(1, :) = [5.3666666666666663_dp, 5.0999999999999996_dp, 3.5333333333333332_dp]
   covariance(2, :) = [5.0999999999999996_dp, 5.0999999999999996_dp, 4.0_dp]
   covariance(3, :) = [3.5333333333333332_dp, 4.0_dp, 4.666666666666667_dp]
   inverse_covariance(1, :) = [8.0467675378267174_dp, -9.9724896836313963_dp, 2.4552957359009717_dp]
   inverse_covariance(2, :) = [-9.9724896836314016_dp, 12.957359009628657_dp, -3.5557083906465041_dp]
   inverse_covariance(3, :) = [2.4552957359009731_dp, -3.555708390646505_dp, 1.4030261348005533_dp]
   actual = mahalanobis(data, [4.833333333333333_dp, 2.5_dp, 3.6666666666666665_dp], covariance)
   call assert_vector_close(actual, &
      [3.9740944520862307_dp, 2.0621274644658478_dp, 2.7498853736818045_dp, &
       0.64534617148097151_dp, 2.7498853736818094_dp, 2.8186611646034119_dp], &
      "Mahalanobis distances", allowed_error=tolerance)
   actual = mahalanobis(data, [4.833333333333333_dp, 2.5_dp, 3.6666666666666665_dp], &
                        inverse_covariance, inverted=.true.)
   call assert_vector_close(actual, &
      [3.9740944520862307_dp, 2.0621274644658478_dp, 2.7498853736818045_dp, &
       0.64534617148097151_dp, 2.7498853736818094_dp, 2.8186611646034119_dp], &
      "inverted Mahalanobis distances", allowed_error=tolerance)

   print *, "test_multivariate: PASS"

end program test_multivariate
