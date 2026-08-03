program test_ffp
  use ffp_mod
  implicit none
  integer::fails,info
  real(dp)::p(5),x(5,1),mu(1),s(1,1),a(2,5),b(2),post(5),stats(6,1)
  real(dp)::kprob(5),ksig(1,1),ddm(1),dds(1,1)
  fails=0
  call exp_decay_probabilities(5,0.1_dp,p)
  call check(abs(sum(p)-1.0_dp)<1.0e-12_dp,'decay normalization')
  call check(abs(half_life(0.1_dp)-7.0_dp)<1.0e-12_dp,'half life')
  x(:,1)=[-2.0_dp,-1.0_dp,0.0_dp,1.0_dp,2.0_dp]; p=0.2_dp
  call weighted_moments(x,p,mu,s)
  call check(abs(mu(1))<1.0e-12_dp,'weighted mean')
  call check(abs(s(1,1)-2.0_dp)<1.0e-12_dp,'weighted variance')
  a(1,:)=1.0_dp; b(1)=1.0_dp; a(2,:)=x(:,1); b(2)=0.5_dp
  call entropy_pool_equalities(p,a,b,post,info)
  call check(info==0,'entropy solver convergence')
  call weighted_moments(x,post,mu,s)
  call check(abs(mu(1)-0.5_dp)<1.0e-8_dp,'entropy mean view')
  call check(relative_entropy(p,post)>=0.0_dp,'relative entropy')
  call check(abs(effective_scenarios(p)-5.0_dp)<1.0e-12_dp,'effective scenarios')
  ksig(1,1)=1.0_dp; call normal_kernel_probabilities(x,[0.0_dp],ksig,kprob,info)
  call check(info==0 .and. maxloc(kprob,dim=1)==3,'normal kernel')
  call weighted_empirical_stats(x,p,0.2_dp,stats)
  call check(abs(stats(1,1))<1.0e-12_dp,'empirical stats')
  call double_decay_covariance(x,0.05_dp,0.20_dp,ddm,dds)
  call check(dds(1,1)>0.0_dp,'double decay covariance')
  if(fails/=0) then; print '(a,i0)','FAILED: ',fails; error stop 1; end if
  print '(a)','All FFP tests passed.'
contains
  subroutine check(ok,name)
    logical,intent(in)::ok; character(*),intent(in)::name
    if(.not.ok) then; fails=fails+1; print '(a,a)','FAIL: ',name; end if
  end subroutine
end program test_ffp
