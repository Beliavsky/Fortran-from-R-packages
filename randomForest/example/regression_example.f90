program regression_example
   use randomforest, only : dp, rf_options, rf_regression_forest
   use randomforest, only : fit_regression, predict_regression
   implicit none

   integer, parameter :: n = 30
   real(dp) :: x(n, 2), y(n), prediction(n), mse
   integer :: i, status
   type(rf_options) :: options
   type(rf_regression_forest) :: forest
   character(len=256) :: message

   do i = 1, n
      x(i, 1) = real(i - 15, dp) / 5.0_dp
      x(i, 2) = real(mod(i, 5), dp)
      y(i) = 1.0_dp + 2.0_dp * x(i, 1) - 0.25_dp * x(i, 2)
   end do

   options%ntree = 101
   options%mtry = 2
   options%seed = 24680
   call fit_regression(x, y, forest, options=options, status=status, message=message)
   if (status /= 0) error stop trim(message)
   call predict_regression(forest, x, prediction)
   mse = sum((prediction - y) ** 2) / real(n, dp)

   print '(a,f12.6)', 'training MSE: ', mse
   print '(a,f12.6)', 'final OOB MSE: ', forest%mse_curve(size(forest%mse_curve))
end program regression_example
