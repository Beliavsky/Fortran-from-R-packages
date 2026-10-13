program covariance_demo
   use burstfin
   implicit none
   real(dp) :: returns(8, 3), intensity
   real(dp), allocatable :: covariance(:, :), augmented(:, :), relative(:, :), weights(:)
   type(statistical_factor_model) :: model
   integer :: i, j, status
   do j = 1, 3
      do i = 1, 8
         returns(i, j) = 0.02_dp * sin(real(i * (j + 1), dp)) + 0.001_dp * real(i - j, dp)
      end do
   end do
   call slide_weight(8, weights, status)
   if (status /= 0) error stop 'time weights failed'
   call var_shrink_eqcor(returns, covariance, intensity, status, weights = weights)
   if (status /= 0) error stop 'shrinkage failed'
   print '(a,f9.6)', 'Estimated shrinkage intensity: ', intensity
   call var_add_benchmark(covariance, [0.5_dp, 0.3_dp, 0.2_dp], augmented, status)
   if (status /= 0) error stop 'benchmark augmentation failed'
   call var_relative_benchmark(augmented, 4, relative, status)
   if (status /= 0) error stop 'relative covariance failed'
   print '(a,3es14.5)', 'Benchmark-relative asset variances: ', (relative(i, i), i = 1, 3)
   call factor_model_stat(returns, model, status)
   if (status /= 0) error stop 'factor model failed'
   call factor_covariance(model, covariance, status)
   if (status /= 0) error stop 'factor reconstruction failed'
   print '(a,i0)', 'Selected factors: ', size(model%loadings, 2)
   print '(a,3es14.5)', 'Factor-model asset variances: ', (covariance(i, i), i = 1, 3)
end program covariance_demo
