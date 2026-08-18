program test_model_utils
   use mcmcpack
   implicit none
   real(dp)::lm(3),pr(3),p(3),bf(2,3),tr(3,3),a(3,3)
   integer::br(2),votes(3,4),gam(6,3)
   type(model_prob_result)::tm
   type(bayes_factor_result)::bfr
   lm=[0.0_dp,-1.0_dp,-2.0_dp];pr=[1.0_dp,1.0_dp,2.0_dp];p=post_prob_mod(lm,pr)
   if(abs(sum(p)-1.0_dp)>1.0e-12_dp.or.any(p<=0.0_dp)) error stop 'post_prob_mod'
   bf(1,:)=[0.0_dp,-2.0_dp,-4.0_dp];bf(2,:)=[0.0_dp,-0.1_dp,-3.0_dp];br=make_breaklist(bf,3.0_dp)
   if(any(br/=[0,0])) error stop 'make_breaklist'
   tr=transition_prior(2,30)
   if(abs(tr(1,2)-0.1_dp)>1e-12_dp.or.tr(3,3)/=1.0_dp) error stop 'transition_prior'
   votes=reshape([1,1,0,1, 1,1,0,0, 0,1,0,1],[3,4],order=[2,1])
   votes(1,:)=[1,1,0,1];votes(2,:)=[1,1,0,0];votes(3,:)=[0,1,0,1]
   a=agreement_matrix(votes)
   if(abs(a(1,2)-0.75_dp)>1e-12_dp) error stop 'agreement_matrix'
   gam(1,:)=[1,0,0];gam(2,:)=[1,0,0];gam(3,:)=[1,1,0];gam(4,:)=[1,0,0];gam(5,:)=[0,0,0];gam(6,:)=[1,1,0]
   tm=top_models(gam,2)
   if(size(tm%probability)/=2.or.abs(tm%probability(1)-0.5_dp)>1e-12_dp) error stop 'top_models'
   bfr=bayes_factor([0.0_dp,log(2.0_dp)])
   if(abs(bfr%factor(2,1)-2.0_dp)>1.0e-12_dp) error stop 'bayes_factor'
   print '(a)','test_model_utils: PASS'
end program
