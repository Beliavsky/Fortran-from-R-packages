program test_multivariate
   use tvgarch, only : dp, tvgarch_spec, mtvgarch_simulation, mtvgarch_fit, dcc_fit, &
                       make_tvgarch_spec, make_tv_spec, mtvgarch_simulate, fit_mtvgarch, fit_mtvgarch_spillover, &
                       fit_dcc, lower_triangle_series, set_random_seed
   implicit none
   integer,parameter::n=260,m=2
   type(tvgarch_spec)::specs(m),spill_specs(m)
   type(mtvgarch_simulation)::sim
   type(mtvgarch_fit)::fit,spill_fit
   type(dcc_fit)::dfit
   real(dp)::ph(3,m),r(m,m),dccp(2)
   real(dp),allocatable::lt(:,:)
   integer::j,st
   logical::order_x(m,m)
   do j=1,m
      call make_tvgarch_spec(specs(j),order_g=[1],order_h=[1,1,0],speed_option=2,status=st)
      call make_tv_spec(specs(j)%tv,[1],1.0_dp,[0.8_dp],[log(10.0_dp)],[0.5_dp],2)
   end do
   ph(:,1)=[0.10_dp,0.08_dp,0.85_dp];ph(:,2)=[0.12_dp,0.10_dp,0.82_dp]
   r=reshape([1.0_dp,0.45_dp,0.45_dp,1.0_dp],[m,m]);dccp=[0.05_dp,0.90_dp]
   call set_random_seed(711)
   call mtvgarch_simulate(n,specs,ph,r,sim,dcc_par=dccp)
   call check(sim%status==0 .and. all(sim%sigma2>0.0_dp),'multivariate simulation')
   call fit_dcc(sim%innovations,dfit,sim%sigma2,turbo=.true.)
   call check(dfit%status==0 .and. dfit%alpha>=0.0_dp .and. dfit%beta>=0.0_dp,'dcc fit')
   call lower_triangle_series(dfit%correlations,lt,st)
   call check(st==0 .and. size(lt,2)==1 .and. all(abs(lt)<1.0_dp),'lower triangle')
   call fit_mtvgarch(sim%y,specs,fit,dynamic_correlation=.false.,turbo=.true.)
   call check(fit%status==0 .and. abs(fit%correlation(1,2))>0.1_dp,'mtvgarch fit')
   do j=1,m
      call make_tvgarch_spec(spill_specs(j),order_h=[1,1,0],status=st)
   end do
   order_x=.false.;order_x(1,2)=.true.;order_x(2,1)=.true.
   call fit_mtvgarch_spillover(sim%y,spill_specs,order_x,spill_fit,max_outer=2,turbo=.true.)
   call check(spill_fit%status==0 .and. size(spill_fit%sigma2,2)==m,'spillover fit')
   print '(a)', 'test_multivariate: PASS'
contains
   subroutine check(ok,msg)
      logical,intent(in)::ok;character(len=*),intent(in)::msg
      if(.not.ok)then;write(*,'(a)')'FAIL: '//msg;error stop 1;end if
   end subroutine check
end program test_multivariate
