program test_empirical_distributions
   use r_stats, only: bw_nrd, bw_nrd0, density, density_result_t, dp, ecdf, ecdf_t, &
                      histogram, histogram_t, predict_ecdf
   use stats_test_assertions, only: assert_close, assert_integer_vector, assert_true, &
                                    assert_vector_close
   implicit none

   real(dp), parameter :: sample(7) = [1.0_dp, 1.0_dp, 2.0_dp, 4.0_dp, 4.0_dp, 4.0_dp, 7.0_dp]
   real(dp), parameter :: bandwidth_sample(6) = [1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp, 16.0_dp, 32.0_dp]
   type(density_result_t) :: estimate
   type(ecdf_t) :: empirical
   type(histogram_t) :: bins

   empirical = ecdf(sample)
   call assert_integer_vector([empirical%n, size(empirical%values)], [7, 4], "ECDF dimensions")
   call assert_vector_close(empirical%values, [1.0_dp, 2.0_dp, 4.0_dp, 7.0_dp], "ECDF knots")
   call assert_vector_close(empirical%probabilities, &
                            [2.0_dp/7.0_dp, 3.0_dp/7.0_dp, 6.0_dp/7.0_dp, 1.0_dp], &
                            "ECDF cumulative probabilities")
   call assert_vector_close(predict_ecdf(empirical, [0.0_dp, 1.0_dp, 1.5_dp, 4.0_dp, 10.0_dp]), &
                            [0.0_dp, 2.0_dp/7.0_dp, 2.0_dp/7.0_dp, 6.0_dp/7.0_dp, 1.0_dp], &
                            "R ECDF evaluation")

   bins = histogram(sample, [0.0_dp, 2.0_dp, 5.0_dp, 8.0_dp])
   call assert_integer_vector(bins%counts, [3, 3, 1], "R right-closed histogram counts")
   call assert_vector_close(bins%density, [3.0_dp/14.0_dp, 1.0_dp/7.0_dp, 1.0_dp/21.0_dp], &
                            "R histogram densities")
   call assert_true(.not. bins%equidistant, "unequal histogram widths")

   bins = histogram(sample, [0.0_dp, 2.0_dp, 5.0_dp, 8.0_dp], right=.false.)
   call assert_integer_vector(bins%counts, [2, 4, 1], "R left-closed histogram counts")

   call assert_close(bw_nrd0(bandwidth_sample), 5.397657223347645_dp, "R bw.nrd0")
   call assert_close(bw_nrd(bandwidth_sample), 6.3572407297205604_dp, "R bw.nrd")

   estimate = density(sample, bandwidth=1.25_dp, number_points=6, from=-1.0_dp, to=9.0_dp)
   call assert_vector_close(estimate%x, [-1.0_dp, 1.0_dp, 3.0_dp, 5.0_dp, 7.0_dp, 9.0_dp], &
                            "density grid")
   call assert_vector_close(estimate%y, &
      [0.027958593238841450_dp, 0.131972990371196308_dp, 0.158056224983577909_dp, &
       0.115103794247048724_dp, 0.053287728469703810_dp, 0.012722558651495231_dp], &
      "direct Gaussian density", 2.0e-12_dp)
   call assert_vector_close(estimate%y, &
      [0.027967386336783025_dp, 0.131966232707846326_dp, 0.158052950389171304_dp, &
       0.115102459565100040_dp, 0.053288814402006182_dp, 0.012725758166249612_dp], &
      "R FFT Gaussian density", 1.0e-5_dp)
   call assert_close(estimate%bandwidth, 1.25_dp, "stored density bandwidth")

   print *, "test_empirical_distributions: PASS"

end program test_empirical_distributions
