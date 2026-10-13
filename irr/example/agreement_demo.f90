! SPDX-License-Identifier: GPL-2.0-or-later
program agreement_demo
   use irr
   implicit none
   real(dp) :: table(2,2), ratings(5,3), n
   integer :: status
   type(irr_result) :: result
   table = reshape([20.0_dp,10.0_dp,5.0_dp,15.0_dp], [2,2])
   result = kappa2_table(table)
   if (result%status /= irr_ok) error stop 'Kappa failed'
   print '(a,f8.4)', 'Cohen kappa: ', result%value
   ratings(:,1) = [1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp]
   ratings(:,2) = [1.0_dp,3.0_dp,3.0_dp,4.0_dp,4.0_dp]
   ratings(:,3) = [2.0_dp,2.0_dp,3.0_dp,5.0_dp,5.0_dp]
   result = icc(ratings, model = 'twoway', kind = 'agreement')
   if (result%status /= irr_ok) error stop 'ICC failed'
   print '(a,f8.4,a,2f8.4)', 'ICC: ', result%value, '; 95% interval: ', result%lbound, result%ubound
   call n2_cohen_kappa([0.2_dp,0.3_dp,0.5_dp], 0.6_dp, 0.2_dp, n, status)
   if (status /= irr_ok) error stop 'Sample size failed'
   print '(a,f8.0)', 'Three-category sample size: ', n
end program
