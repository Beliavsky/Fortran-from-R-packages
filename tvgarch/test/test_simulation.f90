program test_simulation
   use tvgarch, only : dp, tvgarch_spec, tvgarch_simulation, make_tvgarch_spec, &
                       make_tv_spec, tvgarch_simulate, set_random_seed
   implicit none
   type(tvgarch_spec) :: spec
   type(tvgarch_simulation) :: sim
   integer :: st,n,i
   real(dp),allocatable::z(:),x(:)
   real(dp)::ph(3)
   n=180;allocate(z(n),x(n));z=0.0_dp
   do i=1,n;x(i)=real(i,dp)/real(n,dp);z(i)=merge(0.6_dp,-0.4_dp,mod(i,2)==0);end do
   call make_tvgarch_spec(spec,order_g=[1],order_h=[1,1,0],speed_option=2,status=st)
   call make_tv_spec(spec%tv,[1],1.0_dp,[2.0_dp],[log(18.0_dp)],[0.55_dp],2)
   ph=[0.12_dp,0.08_dp,0.86_dp]
   call tvgarch_simulate(n,spec,ph,sim,xtv=x,innovations=z)
   call check(sim%status==0,'simulation status')
   call check(all(sim%sigma2>0.0_dp),'positive variance')
   call check(maxval(abs(sim%y-sqrt(sim%sigma2)*z))<1e-12_dp,'simulation identity')
   call check(sim%g(n)>sim%g(1)+1.0_dp,'tv transition')
   print '(a)', 'test_simulation: PASS'
contains
   subroutine check(ok,msg)
      logical,intent(in)::ok;character(len=*),intent(in)::msg
      if(.not.ok)then;write(*,'(a)')'FAIL: '//msg;error stop 1;end if
   end subroutine check
end program test_simulation
