program test_workflow
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use portvine
   use rugarch, only : garch_spec, dist_norm, seed_rng
   use rvinecopulib, only : bicop_gaussian, bicop_indep
   implicit none
   integer,parameter::d=3,n=70
   real(dp)::r(d,n),weights(d,4),u1,u2
   type(marginal_settings_type)::ms
   type(vine_settings_type)::vs
   type(portvine_roll_result)::ans
   type(garch_spec)::spec
   integer::i,status
   call seed_rng(77123)
   do i=1,n
      call random_number(u1);call random_number(u2)
      r(1,i)=0.012_dp*(u1-0.5_dp)
      r(2,i)=0.65_dp*r(1,i)+0.008_dp*(u2-0.5_dp)
      call random_number(u1);r(3,i)=0.010_dp*(u1-0.5_dp)
   end do
   spec=make_portvine_spec(0,0,1,1,cond_dist=dist_norm)
   ms=make_marginal_settings(50,10,d,spec);ms%max_iterations=140
   vs=make_vine_settings(30,5,vine_dvine,[bicop_indep,bicop_gaussian])
   weights=0.0_dp;weights(2,:)=0.6_dp;weights(3,:)=0.4_dp
   call estimate_risk_roll(r,ms,vs,[0.05_dp],[risk_var,risk_es_mean],40,ans, &
      weights=weights,cond_indices=[1],cond_u=[0.05_dp,0.5_dp], &
      n_mc_samples=50,prior_residual_strategy=.true.,status=status)
   if(status/=0 .or. any(shape(ans%overall)/=[2,1,20]))error stop 1
   if(any(shape(ans%conditional)/=[2,1,3,20]))error stop 2
   if(maxval(abs(ans%realized-(0.6_dp*r(2,51:70)+0.4_dp*r(3,51:70))))>1.0e-12_dp)error stop 3
   if(.not.all(ieee_is_finite(ans%overall)) .or. &
      .not.all(ieee_is_finite(ans%conditional)))error stop 4
   print '(a)', 'test_workflow: PASS'
end program test_workflow
