program example_sparse_index_tracking
   use sparse_index_tracking, only : dp, sparse_index_fit, fit_sparse_index_tracking, sit_success
   implicit none

   integer, parameter :: observations = 180, assets = 8
   real(dp) :: x(observations, assets), benchmark(observations), time
   type(sparse_index_fit) :: fit
   integer :: i

   do i = 1, observations
      time = real(i, dp)
      x(i, 1) = 0.010_dp * sin(0.11_dp * time) + 0.002_dp * cos(0.03_dp * time)
      x(i, 2) = 0.008_dp * cos(0.07_dp * time) - 0.001_dp * sin(0.17_dp * time)
      x(i, 3) = 0.006_dp * sin(0.19_dp * time + 0.3_dp)
      x(i, 4) = 0.007_dp * cos(0.13_dp * time + 0.5_dp)
      x(i, 5) = 0.005_dp * sin(0.23_dp * time)
      x(i, 6) = 0.004_dp * cos(0.29_dp * time)
      x(i, 7) = 0.003_dp * sin(0.31_dp * time + 0.2_dp)
      x(i, 8) = 0.004_dp * cos(0.37_dp * time - 0.1_dp)
   end do
   benchmark = 0.60_dp * x(:, 1) + 0.25_dp * x(:, 2) + 0.15_dp * x(:, 4)
   benchmark(75) = benchmark(75) + 0.12_dp

   call fit_sparse_index_tracking(x, benchmark, 5.0e-7_dp, fit, upper_bound=0.7_dp, &
                                  measure='hete', huber_parameter=0.02_dp)
   if (fit%info /= sit_success) then
      print '(a)', trim(fit%message)
      error stop 1
   end if

   print '(a,i0)', 'MM iterations: ', fit%iterations
   print '(a,i0)', 'Selected assets: ', fit%cardinality
   print '(a)', 'Weights:'
   do i = 1, assets
      print '(i3,2x,f12.8)', i, fit%weights(i)
   end do
   print '(a,es14.6)', 'Tracking RMSE: ', &
      norm2(matmul(x, fit%weights) - benchmark) / sqrt(real(observations, dp))
end program example_sparse_index_tracking
