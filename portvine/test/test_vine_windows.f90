program test_vine_windows
   use portvine
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use rugarch, only : garch_spec, dist_norm, seed_rng
   use rvinecopulib, only : dvine_model, cvine_model, bicop_gaussian, bicop_indep
   implicit none
   integer,parameter::d=3,n=70
   real(dp)::r(d,n),u1,u2
   type(marginal_settings_type)::ms
   type(vine_settings_type)::vs
   type(asset_marginal_result),allocatable::marg(:)
   type(dvine_model),allocatable::dv(:)
   type(cvine_model),allocatable::cv(:)
   type(garch_spec)::spec
   integer::i,status
   call seed_rng(1881)
   do i=1,n
      call random_number(u1);call random_number(u2)
      r(1,i)=0.01_dp*(u1-0.5_dp)
      r(2,i)=0.7_dp*r(1,i)+0.006_dp*(u2-0.5_dp)
      call random_number(u1);r(3,i)=0.008_dp*(u1-0.5_dp)
   end do
   spec=make_portvine_spec(0,0,1,1,cond_dist=dist_norm)
   ms=make_marginal_settings(50,10,d,spec);ms%max_iterations=140
   vs=make_vine_settings(30,5,vine_dvine,[bicop_indep,bicop_gaussian])
   call fit_rolling_marginals(r,ms,vs%train_size,marg,status)
   if(status/=0)error stop 1
   call fit_vine_windows(marg,ms,vs,n,[1],dv,cv,status)
   if(status/=0 .or. size(dv)/=4 .or. dv(1)%order(1)/=1)error stop 2
   if(.not.ieee_is_finite(dv(1)%loglik))error stop 3
   print '(a)', 'test_vine_windows: PASS'
end program test_vine_windows
