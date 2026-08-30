program test_survival
   use rpart
   implicit none
   type(rpart_model) :: model
   type(rpart_control) :: ctl
   real(dp) :: x(6,1), time(6), status(6), newtime(6), pred(6)
   real(dp) :: start(3), stop(3), status2(3), newtime2(3)
   real(dp), parameter :: tol=2.0e-12_dp
   integer :: i,stat

   call rpart_exp_transform_right([1.0_dp,2.0_dp,3.0_dp],[1.0_dp,0.0_dp,1.0_dp],newtime(1:3),stat=stat)
   call check(stat==0,'right transform status')
   call check(all(abs(newtime(1:3)-[1.0_dp/3.0_dp,2.0_dp/3.0_dp,1.0_dp])<tol),'right transform values')

   start=[0.5_dp,1.0_dp,0.5_dp];stop=[1.0_dp,2.0_dp,3.0_dp];status2=[1.0_dp,0.0_dp,1.0_dp]
   call rpart_exp_transform_startstop(start,stop,status2,newtime2,stat=stat)
   call check(stat==0,'start-stop transform status')
   call check(all(abs(newtime2-[0.5_dp,1.0_dp/3.0_dp,7.0_dp/6.0_dp])<tol),'start-stop transform values')
   start(1)=0.0_dp
   call rpart_exp_transform_startstop(start,stop,status2,newtime2,stat=stat)
   call check(stat/=0,'zero start rejected like rpart.exp')

   do i=1,6;x(i,1)=real(i,dp);end do
   time=[1.0_dp,1.2_dp,2.0_dp,2.2_dp,3.0_dp,3.2_dp]
   status=[1.0_dp,1.0_dp,0.0_dp,1.0_dp,0.0_dp,1.0_dp]
   ctl%cp=0.0_dp;ctl%xval=0;ctl%minsplit=2;ctl%minbucket=1;ctl%maxdepth=2
   call rpart_fit_survival(x,time,status,model,ctl,stat=stat)
   call check(stat==0,'survival fit status')
   call check(model%method==RPART_EXP,'survival model method')
   call rpart_predict_values(model,x,pred)
   call check(all(pred>0.0_dp),'survival rates positive')

   start=[0.2_dp,0.4_dp,0.6_dp];stop=[1.0_dp,2.0_dp,3.0_dp]
   status2=[1.0_dp,0.0_dp,1.0_dp]
   call rpart_fit_survival_startstop(reshape([1.0_dp,2.0_dp,3.0_dp],[3,1]),start,stop,status2,model,ctl,stat=stat)
   call check(stat==0,'start-stop fit status')

   print '(a)','test_survival: PASS'
contains
   subroutine check(ok,msg)
      logical,intent(in)::ok
      character(*),intent(in)::msg
      if(.not.ok)then
         write(*,'(a)')'FAIL: '//trim(msg)
         error stop 1
      end if
   end subroutine check
end program test_survival
