! SPDX-License-Identifier: GPL-3.0-only
! See NOTICE.md for ranger upstream provenance.
program classification_example
   use r_kinds, only : dp
   use ranger, only : ranger_options, ranger_classification_forest
   use ranger, only : fit_ranger_classification, predict_ranger_classification
   implicit none
   real(dp) :: x(12, 2), votes(12, 2)
   integer :: y(12), prediction(12), i, status
   type(ranger_options) :: options
   type(ranger_classification_forest) :: forest

   do i = 1, 6
      x(i, :) = [-real(7 - i, dp), real(mod(i, 2), dp)]
      y(i) = 1
   end do
   do i = 7, 12
      x(i, :) = [real(i - 6, dp), real(mod(i, 2), dp)]
      y(i) = 2
   end do

   options%num_trees = 25
   options%mtry = 2
   options%seed = 2026
   call fit_ranger_classification(x, y, forest, options=options, status=status)
   if (status /= 0) error stop 'classification fit failed'
   call predict_ranger_classification(forest, x, prediction, votes)

   print '(a,i0)', 'training errors: ', count(prediction /= y)
   print '(a,2f10.5)', 'first case vote fractions: ', votes(1, :)
end program classification_example
