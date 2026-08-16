program test_batch
   use rfast
   implicit none
   real(dp) :: x(6,2), pars(2,2), tt(3,2), av(4,2), y(6), cr(2)
   integer :: g(6)
   type(test_result) :: a

   x(:,1)=[1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp]
   x(:,2)=2.0_dp*x(:,1)
   pars=colnormal_mle(x)
   call assert_close(pars(1,1),3.5_dp,1e-12_dp,'col normal mean')
   call assert_close(pars(1,2),7.0_dp,1e-12_dp,'col normal mean 2')
   tt=column_ttests(x)
   call assert_true(all(tt(3,:)<0.05_dp),'column t tests')
   g=[1,1,1,2,2,2]
   y=[1.0_dp,1.2_dp,0.8_dp,5.0_dp,5.2_dp,4.8_dp]
   a=one_way_anova(y,g)
   call assert_true(a%statistic>100.0_dp .and. a%pvalue<1e-3_dp,'anova')
   av=one_way_anovas(x,g)
   call assert_true(all(av(1,:)>0.0_dp),'anovas')
   cr=column_correlations(x,x(:,1))
   call assert_close(cr(2),1.0_dp,1e-12_dp,'column cor')
   print *, 'test_batch: PASS'
contains
   subroutine assert_close(got,want,tol,msg)
      real(dp),intent(in)::got,want,tol
      character(*),intent(in)::msg
      if(abs(got-want)>tol)then
         print *, 'FAIL ',trim(msg),got,want
         error stop 1
      end if
   end subroutine
   subroutine assert_true(ok,msg)
      logical,intent(in)::ok
      character(*),intent(in)::msg
      if(.not.ok)then
         print *, 'FAIL ',trim(msg)
         error stop 1
      end if
   end subroutine
end program test_batch
