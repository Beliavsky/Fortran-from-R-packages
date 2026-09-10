program regression_example
   use grf, only : dp, grf_options, grf_forest, regression_forest, predict_regression_forest
   implicit none

   integer, parameter :: n = 60
   real(dp) :: x(n,2)
   real(dp) :: y(n)
   real(dp), allocatable :: prediction(:)
   type(grf_options) :: options
   type(grf_forest) :: forest
   integer :: i
   integer :: info

   do i = 1, n
      x(i,1) = -1.0_dp + 2.0_dp * real(i - 1,dp) / real(n - 1,dp)
      x(i,2) = sin(0.31_dp * real(i,dp))
      y(i) = 1.5_dp * x(i,1) - 0.4_dp * x(i,2)
   end do
   options%num_trees = 40
   options%min_node_size = 4
   options%sample_fraction = 0.8_dp
   options%seed = 2026
   call regression_forest(x, y, forest, info, options)
   if (info /= 0) error stop 'regression_forest failed'
   call predict_regression_forest(forest, x(1:5,:), prediction)
   print '(a)', 'First five regression-forest predictions:'
   do i = 1, size(prediction)
      print '(i3,2f12.6)', i, y(i), prediction(i)
   end do
end program regression_example
