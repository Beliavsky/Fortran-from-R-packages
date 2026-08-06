program test_fit
   use tvgarch, only : dp, tvgarch_spec, tvgarch_simulation, tvgarch_fit, &
                       make_tvgarch_spec, make_tv_spec, tvgarch_simulate, fit_tvgarch, &
                       tvgarch_forecast, tvgarch_quantile_path, set_random_seed
   implicit none
   integer,parameter::n=320
   type(tvgarch_spec)::spec
   type(tvgarch_simulation)::sim
   type(tvgarch_fit)::fit
   real(dp)::ph(3)
   real(dp),allocatable::x(:),fc(:),q(:,:)
   integer::i,st
   allocate(x(n));do i=1,n;x(i)=real(i,dp)/real(n,dp);end do
   call make_tvgarch_spec(spec,order_g=[1],order_h=[1,1,0],speed_option=2,status=st)
   call make_tv_spec(spec%tv,[1],1.1_dp,[1.8_dp],[log(16.0_dp)],[0.52_dp],2)
   ph=[0.10_dp,0.10_dp,0.82_dp]
   call set_random_seed(1201)
   call tvgarch_simulate(n,spec,ph,sim,xtv=x)
   call check(sim%status==0,'source simulation')
   call fit_tvgarch(sim%y,spec,fit,xtv=x,initial_g=[1.0_dp,1.5_dp,log(12.0_dp),0.5_dp], &
                    initial_h=[0.12_dp,0.08_dp,0.80_dp],max_outer=8,rel_tol=2e-3_dp,turbo=.true.)
   call check(allocated(fit%sigma2) .and. all(fit%sigma2>0.0_dp),'fit variance')
   call check(fit%par_g(2)>0.3_dp,'positive transition size recovery')
   call check(abs(fit%par_g(4)-0.52_dp)<0.25_dp,'location recovery')
   call tvgarch_forecast(fit,5,fc,st,n_sim=100)
   call check(st==0 .and. size(fc)==5 .and. all(fc>0.0_dp),'forecast')
   call tvgarch_quantile_path(fit,[0.05_dp,0.95_dp],q,st)
   call check(st==0 .and. all(q(:,1)<q(:,2)),'quantile paths')
   print '(a)', 'test_fit: PASS'
contains
   subroutine check(ok,msg)
      logical,intent(in)::ok;character(len=*),intent(in)::msg
      if(.not.ok)then;write(*,'(a)')'FAIL: '//msg;error stop 1;end if
   end subroutine check
end program test_fit
