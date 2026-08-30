program classification_example
   use randomforest, only : dp, rf_options, rf_classification_forest
   use randomforest, only : fit_classification, predict_classification
   implicit none

   integer, parameter :: n = 20
   real(dp) :: x(n, 2), probability(n, 2)
   integer :: y(n), prediction(n), i, status
   type(rf_options) :: options
   type(rf_classification_forest) :: forest
   character(len=256) :: message

   do i = 1, n
      x(i, 1) = real(i - 10, dp)
      x(i, 2) = real(mod(i, 4), dp)
      y(i) = merge(1, 2, i <= 10)
   end do

   options%ntree = 101
   options%mtry = 2
   options%seed = 12345
   call fit_classification(x, y, forest, options=options, status=status, message=message)
   if (status /= 0) error stop trim(message)
   call predict_classification(forest, x, prediction, probabilities=probability)

   print '(a,i0)', 'training errors: ', count(prediction /= y)
   print '(a,2f10.4)', 'first-case class probabilities: ', probability(1, :)
end program classification_example
