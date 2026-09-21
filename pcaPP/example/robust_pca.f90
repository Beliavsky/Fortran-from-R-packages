program robust_pca_example
  use pcapp_api, only: dp, pca_result, pca_grid
  implicit none

  real(dp) :: x(8, 3)
  type(pca_result) :: fit

  x(:, 1) = [-4.0_dp, -3.0_dp, -2.0_dp, -1.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
  x(:, 2) = [0.2_dp, -0.1_dp, 0.1_dp, -0.2_dp, 0.2_dp, -0.1_dp, 0.1_dp, -0.2_dp]
  x(:, 3) = [0.0_dp, 0.1_dp, -0.1_dp, 0.0_dp, 0.0_dp, -0.1_dp, 0.1_dp, 0.0_dp]

  fit = pca_grid(x, k=2, method=0, maxiter=6, splitcircle=31)
  write (*, '(a,2f12.6)') 'robust component scales:', fit%sdev
  write (*, '(a,3f12.6)') 'first loading:', fit%loadings(:, 1)
end program robust_pca_example
