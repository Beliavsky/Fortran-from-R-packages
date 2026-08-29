program test_reference
   use evd
   implicit none
   real(dp) :: dat(12,2), v, marg(2,3), x(2)
   integer :: fails
   fails=0
   dat(:,1)=[0.1_dp,0.4_dp,0.2_dp,0.9_dp,1.1_dp,0.5_dp,1.7_dp,0.3_dp,1.5_dp,0.7_dp,2.0_dp,1.3_dp]
   dat(:,2)=[0.2_dp,0.1_dp,0.6_dp,0.8_dp,0.4_dp,1.2_dp,1.5_dp,0.5_dp,1.7_dp,0.9_dp,1.8_dp,1.4_dp]
   v=bvpot_censored_nll(dat,[0.5_dp,0.5_dp],'log',[0.7_dp],[1.0_dp,1.2_dp],[0.1_dp,-0.05_dp])
   call check(abs(v-25.56903897339505_dp)<1.0e-7_dp,'upstream censored logistic likelihood',fails)
   v=bvpot_poisson_nll(dat,[0.5_dp,0.5_dp],'log',[0.7_dp],[1.0_dp,1.2_dp],[0.1_dp,-0.05_dp])
   call check(abs(v-23.524363780270296_dp)<1.0e-10_dp,'upstream poisson logistic likelihood',fails)
   marg=0.0_dp
   marg(:,2)=1.0_dp
   x=[0.2_dp,-0.1_dp]
   call check(abs(dmvlog(x,0.7_dp,marg,.true.)-dbvlog(x(1),x(2),0.7_dp,marg(1,:),marg(2,:),.true.))<1e-12_dp, &
      'multivariate density agrees with bivariate density in d=2',fails)
   if(fails>0) then
   write(*,'(a,i0)') 'test_reference: FAIL ',fails
   error stop 1
   end if
   write(*,'(a)') 'test_reference: PASS'
contains
   subroutine check(ok,msg,fails)
      logical,intent(in)::ok
      character(len=*),intent(in)::msg
      integer,intent(inout)::fails
      if(.not.ok) then
      write(*,'(a)') 'FAIL: '//msg
      fails=fails+1
      end if
   end subroutine
end program
