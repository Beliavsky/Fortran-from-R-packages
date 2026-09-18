! SPDX-License-Identifier: MIT
! SPDX-FileComment: Deterministic multivariate linear-model and canonical-correlation tests.
program test_manova
   use stats_test_assertions, only: assert_equal, assert_matrix_close, assert_vector_close
   use r_stats, only: cancor, cancor_result_t, dp, manova_lm, manova_result_t
   use r_stats, only: mlm_fit, mlm_fit_t
   implicit none

   real(dp) :: cx(10, 2), cy(10, 2), full_x(10, 2), reduced_x(10, 1), responses(10, 3)
   real(dp) :: noise1(10), noise2(10), noise3(10), z(10)
   integer :: i
   type(cancor_result_t) :: canonical
   type(manova_result_t) :: tests
   type(mlm_fit_t) :: full, reduced

   z = [0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp]
   noise1 = [0.1_dp, -0.1_dp, 0.2_dp, 0.0_dp, -0.2_dp, 0.1_dp, 0.0_dp, -0.1_dp, 0.1_dp, -0.1_dp]
   noise2 = [0.0_dp, 0.2_dp, -0.1_dp, 0.1_dp, 0.0_dp, -0.2_dp, 0.1_dp, 0.0_dp, -0.1_dp, 0.1_dp]
   noise3 = [0.2_dp, 0.0_dp, -0.2_dp, 0.1_dp, -0.1_dp, 0.2_dp, -0.1_dp, 0.0_dp, 0.1_dp, -0.2_dp]
   do i = 1, 10
      reduced_x(i, 1) = real(i - 1, dp)
   end do
   full_x(:, 1) = reduced_x(:, 1)
   full_x(:, 2) = z
   responses(:, 1) = 1.0_dp + 0.5_dp*reduced_x(:, 1) + 1.2_dp*z + noise1
   responses(:, 2) = 2.0_dp - 0.3_dp*reduced_x(:, 1) + 0.8_dp*z + noise2
   responses(:, 3) = -1.0_dp + 0.2_dp*reduced_x(:, 1) - 0.5_dp*z + noise3

   full = mlm_fit(responses, full_x)
   reduced = mlm_fit(responses, reduced_x)
   call assert_equal(full%status, 0, "multivariate linear-model status")
   call assert_matrix_close(full%coefficients, reshape( &
      [1.07_dp, 0.4925_dp, 1.1275_dp, 2.01_dp, -0.3075_dp, 0.8675_dp, &
       -0.96_dp, 0.185_dp, -0.445_dp], [3, 3]), "multivariate coefficients")
   call assert_matrix_close(full%residual_sscp, reshape( &
      [0.1195_dp, -0.0725_dp, 0.049_dp, -0.0725_dp, 0.1155_dp, -0.065_dp, &
       0.049_dp, -0.065_dp, 0.178_dp], [3, 3]), "residual SSCP")

   tests = manova_lm(reduced, full)
   call assert_equal(tests%status, 0, "MANOVA status")
   call assert_vector_close(tests%statistic, &
      [0.990803895356027_dp, 0.00919610464397838_dp, &
       107.741694305838_dp, 107.741694305838_dp], "MANOVA criteria")
   call assert_vector_close(tests%approximate_f, spread(179.56949050973_dp, 1, 4), &
                            "MANOVA approximate F")
   call assert_vector_close(tests%p_value, spread(1.646677620432e-5_dp, 1, 4), &
                            "MANOVA p values")

   cx(:, 1) = reduced_x(:, 1)
   cx(:, 2) = [2.0_dp, 1.0_dp, 4.0_dp, 3.0_dp, 7.0_dp, &
               5.0_dp, 8.0_dp, 6.0_dp, 10.0_dp, 9.0_dp]
   cy(:, 1) = [1.1_dp, 1.8_dp, 2.9_dp, 3.7_dp, 5.2_dp, &
               5.8_dp, 7.1_dp, 7.9_dp, 9.2_dp, 9.8_dp]
   cy(:, 2) = [4.0_dp, 2.0_dp, 5.0_dp, 3.0_dp, 8.0_dp, &
               6.0_dp, 9.0_dp, 7.0_dp, 11.0_dp, 10.0_dp]
   canonical = cancor(cx, cy)
   call assert_equal(canonical%status, 0, "canonical-correlation status")
   call assert_vector_close(canonical%correlation, &
      [0.999610333714639_dp, 0.957678339966711_dp], "canonical correlations")
   call assert_matrix_close(abs(canonical%x_coefficients), abs(reshape( &
      [-0.108164530578339_dp, -0.00213506184944830_dp, &
       -0.232347041514549_dp, 0.256281397855356_dp], [2, 2])), &
      "canonical x coefficients")
   call assert_matrix_close(abs(canonical%y_coefficients), abs(reshape( &
      [-0.118963065786905_dp, 0.0110334213078566_dp, &
       -0.210866907274288_dp, 0.243621615630978_dp], [2, 2])), &
      "canonical y coefficients")
   call assert_vector_close(canonical%x_center, [4.5_dp, 5.5_dp], "canonical x centers")
   call assert_vector_close(canonical%y_center, [5.45_dp, 6.5_dp], "canonical y centers")

   print *, "test_manova: PASS"

end program test_manova
