program test_lm
   use tvgarch, only : dp, tvgarch_spec, tvgarch_simulation, tvgarch_test_result, &
                       make_tvgarch_spec, tvgarch_simulate, tvgarch_test, set_random_seed
   implicit none
   integer,parameter::n=220
   type(tvgarch_spec)::spec
   type(tvgarch_simulation)::sim
   type(tvgarch_test_result)::test
   real(dp)::ph(3)
   integer::st
   call make_tvgarch_spec(spec,order_h=[1,1,0],status=st)
   ph=[0.1_dp,0.08_dp,0.86_dp]
   call set_random_seed(999)
   call tvgarch_simulate(n,spec,ph,sim)
   call tvgarch_test(sim%y,test,turbo=.true.)
   call check(test%status==0,'test status')
   call check(all(test%nonrobust(:,2)>=0.0_dp) .and. all(test%nonrobust(:,2)<=1.0_dp),'nonrobust pvalues')
   call check(all(test%robust(:,2)>=0.0_dp) .and. all(test%robust(:,2)<=1.0_dp),'robust pvalues')
   call check(test%selected_order>=0 .and. test%selected_order<=3,'selected order')
   print '(a)', 'test_lm: PASS'
contains
   subroutine check(ok,msg)
      logical,intent(in)::ok;character(len=*),intent(in)::msg
      if(.not.ok)then;write(*,'(a)')'FAIL: '//msg;error stop 1;end if
   end subroutine check
end program test_lm
