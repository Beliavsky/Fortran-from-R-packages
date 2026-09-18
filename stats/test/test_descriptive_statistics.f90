program test_descriptive_statistics
   use r_stats, only: cor, cov, cov2cor, dp, iqr, mad, median, quantile, sd, var, weighted_mean
   use stats_test_assertions, only: assert_close, assert_matrix_close, assert_vector_close
   implicit none

   real(dp), parameter :: x(6) = [1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp, 16.0_dp, 32.0_dp]
   real(dp), parameter :: y(6) = [3.0_dp, 1.0_dp, 4.0_dp, 1.0_dp, 5.0_dp, 9.0_dp]
   real(dp), parameter :: z(6) = [2.0_dp, 7.0_dp, 1.0_dp, 8.0_dp, 2.0_dp, 8.0_dp]
   real(dp) :: data(6, 3), covariance(3, 3)

   data(:, 1) = x
   data(:, 2) = y
   data(:, 3) = z

   call assert_close(median(x), 6.0_dp, "R median")
   call assert_vector_close(quantile(x, [0.1_dp, 0.25_dp, 0.5_dp, 0.75_dp, 0.9_dp]), &
                            [1.5_dp, 2.5_dp, 6.0_dp, 14.0_dp, 24.0_dp], "R type-7 quantiles")
   call assert_close(sd(x), 11.861703081766969_dp, "R standard deviation")
   call assert_close(var(x), 140.69999999999999_dp, "R variance")
   call assert_close(cov(x, y), 31.100000000000001_dp, "R vector covariance")
   call assert_close(cor(x, y), 0.87558402174921179_dp, "R vector correlation")

   covariance = reshape([140.69999999999999_dp, 31.100000000000001_dp, 15.599999999999999_dp, &
                         31.100000000000001_dp, 8.9666666666666668_dp, -0.066666666666666666_dp, &
                         15.599999999999999_dp, -0.066666666666666666_dp, &
                         11.066666666666666_dp], [3, 3])
   call assert_matrix_close(var(data), covariance, "R matrix variance")
   call assert_matrix_close(cov(data), covariance, "R covariance matrix")
   call assert_matrix_close(cor(data), &
      reshape([1.0_dp, 0.87558402174921179_dp, 0.39533852277393483_dp, &
               0.87558402174921179_dp, 1.0_dp, -0.0066924450472417548_dp, &
               0.39533852277393483_dp, -0.0066924450472417548_dp, 1.0_dp], [3, 3]), &
      "R correlation matrix", 2.0e-12_dp)
   call assert_matrix_close(cov2cor(covariance), cor(data), "R cov2cor", 2.0e-12_dp)

   call assert_close(iqr(x), 11.5_dp, "R interquartile range")
   call assert_close(mad(x), 6.6717_dp, "R median absolute deviation")
   call assert_close(mad(x, center=5.0_dp, constant=1.0_dp), 3.5_dp, "R custom MAD")
   call assert_close(weighted_mean(x, [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]), &
                     15.285714285714286_dp, "R weighted mean")

   print *, "test_descriptive_statistics: PASS"

end program test_descriptive_statistics
