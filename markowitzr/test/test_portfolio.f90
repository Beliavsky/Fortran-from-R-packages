! SPDX-License-Identifier: LGPL-3.0-or-later
program test_portfolio
   use markowitzr, only: dp, markowitz_result, mp_vcov
   use markowitzr, only: covariance_empirical, weights_upstream, weights_all_columns
   implicit none
   real(dp) :: x(6,2), feat(6,1), weights(6), jmat(1,2), gmat(1,2)
   real(dp) :: sigma_n(2,2), hedge_covariance
   type(markowitz_result) :: plain, constrained, hedged, conditional
   type(markowitz_result) :: upstream_weighted, all_weighted

   x = reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp, &
                2.0_dp,0.0_dp,1.0_dp,-1.0_dp,3.0_dp,2.0_dp],[6,2])
   feat(:,1) = [-1.0_dp,0.0_dp,1.0_dp,2.0_dp,-2.0_dp,3.0_dp]
   weights = [1.0_dp,1.1_dp,0.9_dp,1.2_dp,0.8_dp,1.05_dp]

   plain = mp_vcov(x,covariance_method=covariance_empirical)
   if (plain%status /= 0) error stop 1
   if (any(shape(plain%w) /= [2,1])) error stop 1
   if (maxval(abs(plain%w(:,1)-[1.144736842105263_dp, &
      0.276315789473684_dp])) > 2.0e-13_dp) error stop 1
   if (maxval(abs(plain%w_covariance-transpose(plain%w_covariance))) > 0.0_dp) error stop 1

   jmat = reshape([1.0_dp,0.0_dp],[1,2])
   constrained = mp_vcov(x,jmat=jmat)
   if (constrained%status /= 0) error stop 1
   if (maxval(abs(constrained%w(:,1)-[1.2_dp,0.0_dp])) > 2.0e-13_dp) error stop 1

   gmat = jmat
   hedged = mp_vcov(x,gmat=gmat)
   if (hedged%status /= 0) error stop 1
   sigma_n = reshape([2.916666666666667_dp,0.583333333333333_dp, &
                      0.583333333333333_dp,1.805555555555556_dp],[2,2])
   hedge_covariance = dot_product(gmat(1,:),matmul(sigma_n,hedged%w(:,1)))
   if (abs(hedge_covariance) > 2.0e-13_dp) error stop 1

   conditional = mp_vcov(x,feat=feat)
   if (conditional%status /= 0) error stop 1
   if (any(shape(conditional%w) /= [2,2])) error stop 1
   if (any(shape(conditional%w_covariance) /= [4,4])) error stop 1

   upstream_weighted = mp_vcov(x,feat=feat,weights=weights, &
      weight_mode=weights_upstream)
   all_weighted = mp_vcov(x,feat=feat,weights=weights, &
      weight_mode=weights_all_columns)
   if (upstream_weighted%status /= 0 .or. all_weighted%status /= 0) error stop 1
   if (maxval(abs(upstream_weighted%w-all_weighted%w)) < 1.0e-8_dp) error stop 1

   print '(a)', 'test_portfolio: PASS'
end program test_portfolio
