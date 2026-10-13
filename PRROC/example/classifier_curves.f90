! SPDX-License-Identifier: GPL-3.0-only
program classifier_curves
   use prroc_api, only : dp, prroc_result, pr_curve, roc_curve, prroc_ok
   implicit none
   type(prroc_result) :: pr, roc
   real(dp), parameter :: scores(4) = [4.0_dp, 4.0_dp, 2.0_dp, 1.0_dp]
   real(dp), parameter :: positive_probability(4) = [0.2_dp, 0.9_dp, 0.7_dp, 0.1_dp]
   integer :: status
   call pr_curve(scores, pr, weights_class0 = positive_probability, curve = .true., &
      max_compute = .true., min_compute = .true., rand_compute = .true., status = status)
   if (status /= prroc_ok) error stop 'PR calculation failed'
   call roc_curve(scores, roc, weights_class0 = positive_probability, curve = .true., status = status)
   if (status /= prroc_ok) error stop 'ROC calculation failed'
   print '(a,f12.9)', 'Nonlinear PR AUC: ', pr%auc
   print '(a,f12.9)', 'Maximum PR AUC:   ', pr%maximum%auc
   print '(a,f12.9)', 'Minimum PR AUC:   ', pr%minimum%auc
   print '(a,f12.9)', 'Random PR AUC:    ', pr%random%auc
   print '(a,f12.9)', 'ROC AUC:          ', roc%auc
   print '(a,i0)', 'PR curve rows:   ', size(pr%curve, 1)
end program classifier_curves
