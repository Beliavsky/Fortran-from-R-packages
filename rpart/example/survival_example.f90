program survival_example
   use rpart
   implicit none
   integer,parameter::n=8
   real(dp)::x(n,1),time(n),status(n),rate(n)
   integer::i,stat
   type(rpart_control)::control
   type(rpart_model)::model

   do i=1,n;x(i,1)=real(i,dp);end do
   time=[1.0_dp,1.3_dp,1.8_dp,2.0_dp,2.4_dp,2.8_dp,3.2_dp,3.7_dp]
   status=[1.0_dp,1.0_dp,0.0_dp,1.0_dp,0.0_dp,1.0_dp,0.0_dp,1.0_dp]
   control=rpart_make_control(minsplit=3,minbucket=1,cp=0.0_dp,xval=0,stat=stat)
   if(stat/=0)error stop 'invalid control values'
   call rpart_fit_survival(x,time,status,model,control=control,stat=stat)
   if(stat/=0)error stop 'survival fit failed'
   call rpart_predict_values(model,x,rate)
   write(*,'(a,*(f7.3,1x))')'relative event rates: ',rate
end program survival_example
