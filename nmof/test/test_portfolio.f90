! SPDX-License-Identifier: GPL-3.0-only
program test_portfolio
   use nmof
   implicit none
   integer :: failures,status
   real(dp) :: cov5(5,5),w5(5),groups(2,5),ret(240,5),target_ret(240,6),mu(5)
   real(dp) :: w3(3),cov3(3,3),m3(3)

   failures=0
   cov5=reshape([ &
      0.000988087100677907_dp,-0.0000179669410403153_dp,0.000368923882626859_dp,0.000208303611101873_dp,0.000262742052359594_dp, &
      -0.0000179669410403153_dp,0.00171852167358765_dp,0.0000857467457561209_dp,0.0000215059246610556_dp,0.0000283532159921211_dp, &
      0.000368923882626859_dp,0.0000857467457561209_dp,0.00075871953281751_dp,0.000194002299424151_dp,0.000188824454515841_dp, &
      0.000208303611101873_dp,0.0000215059246610556_dp,0.000194002299424151_dp,0.000265780633005374_dp,0.000132611196599808_dp, &
      0.000262742052359594_dp,0.0000283532159921211_dp,0.000188824454515841_dp,0.000132611196599808_dp,0.00025948420130626_dp],[5,5])

   call minimum_variance(cov5,w5,lower=[0.0_dp],upper=[0.25_dp],status=status)
   call check_true('bounded minvar status',status==nmof_ok)
   call check_close('bounded minvar budget',sum(w5),1.0_dp,5.0e-9_dp)
   call check_true('bounded minvar limits',all(w5>=-1.0e-10_dp).and.all(w5<=0.25_dp+1.0e-9_dp))

   groups=0.0_dp; groups(1,1)=1.0_dp; groups(2,4:5)=1.0_dp
   call minimum_variance(cov5,w5,lower=[0.0_dp],upper=[0.4_dp],groups=groups, &
      group_lower=[0.29_dp,0.10_dp],group_upper=[0.30_dp,0.20_dp],status=status)
   call check_true('group minvar status',status==nmof_ok)
   call check_true('group1 lower',w5(1)>=0.29_dp-1.0e-8_dp)
   call check_true('group1 upper',w5(1)<=0.30_dp+1.0e-8_dp)
   call check_true('group2 lower',sum(w5(4:5))>=0.10_dp-1.0e-8_dp)
   call check_true('group2 upper',sum(w5(4:5))<=0.20_dp+1.0e-8_dp)

   cov3=reshape([0.04_dp,0.006_dp,0.004_dp,0.006_dp,0.09_dp,0.012_dp,0.004_dp,0.012_dp,0.16_dp],[3,3])
   m3=[0.06_dp,0.09_dp,0.12_dp]
   call mean_variance_portfolio(m3,cov3,0.08_dp,w3,status=status)
   call check_true('MV status',status==nmof_ok)
   call check_close('MV budget',sum(w3),1.0_dp,5.0e-9_dp)
   call check_true('MV return',dot_product(m3,w3)>=0.08_dp-1.0e-8_dp)

   call maximum_sharpe(m3,cov3,w3,status=status)
   call check_true('max Sharpe status',status==nmof_ok)
   call check_close('max Sharpe budget',sum(w3),1.0_dp,5.0e-10_dp)

   ret=random_returns(5,240,[0.015_dp],mean=[0.001_dp,0.0015_dp,0.0008_dp,0.0012_dp,0.001_dp],rho=0.25_dp,exact=.true.,seed=901_i8,status=status)
   mu=column_means(ret)
   call minimum_cvar(ret,w5,q=0.10_dp,lower=[0.05_dp],upper=[0.5_dp],min_return=0.0008_dp,mean_returns=mu,status=status,maxiter=4000,tol=1.0e-7_dp)
   call check_true('CVaR status',status==nmof_ok.or.status==nmof_numerical_failure)
   call check_close('CVaR budget',sum(w5),1.0_dp,2.0e-7_dp)
   call check_true('CVaR limits',all(w5>=0.05_dp-1.0e-7_dp).and.all(w5<=0.5_dp+1.0e-7_dp))

   call minimum_mad(ret,w5,lower=[0.05_dp],upper=[0.5_dp],min_return=0.0008_dp,mean_returns=mu,status=status,maxiter=4000,tol=1.0e-7_dp)
   call check_true('MAD status',status==nmof_ok.or.status==nmof_numerical_failure)
   call check_close('MAD budget',sum(w5),1.0_dp,2.0e-7_dp)
   call check_true('MAD limits',all(w5>=0.05_dp-1.0e-7_dp).and.all(w5<=0.5_dp+1.0e-7_dp))

   target_ret(:,2:6)=ret
   target_ret(:,1)=0.2_dp*ret(:,1)+0.3_dp*ret(:,2)+0.5_dp*ret(:,3)
   call tracking_portfolio(w=w5,lower=[0.0_dp],upper=[1.0_dp],returns=target_ret,objective='sum.of.squares',status=status)
   call check_true('tracking status',status==nmof_ok)
   call check_close('tracking budget',sum(w5),1.0_dp,5.0e-9_dp)
   call check_close('tracking w1',w5(1),0.2_dp,2.0e-5_dp)
   call check_close('tracking w2',w5(2),0.3_dp,2.0e-5_dp)
   call check_close('tracking w3',w5(3),0.5_dp,2.0e-5_dp)

   if(failures>0) then
      write(*,'(a,i0)') 'test_portfolio failures: ',failures
      error stop 1
   end if
   write(*,'(a)') 'test_portfolio: PASS'
contains
   subroutine check_close(name,actual,expected,tol)
      character(len=*),intent(in)::name
      real(dp),intent(in)::actual,expected,tol
      if(abs(actual-expected)>tol) then
         failures=failures+1; write(*,'(a,2(1x,es16.8))') trim(name)//' failed:',actual,expected
      end if
   end subroutine check_close
   subroutine check_true(name,condition)
      character(len=*),intent(in)::name
      logical,intent(in)::condition
      if(.not.condition) then; failures=failures+1; write(*,'(a)') trim(name)//' failed'; end if
   end subroutine check_true
end program test_portfolio
