program basic_example
   use rpart
   implicit none
   integer,parameter::n=20
   real(dp)::x(n,2),y(n),pred(n)
   integer::i,stat
   type(rpart_control)::control
   type(rpart_model)::model

   do i=1,n
      x(i,1)=real(i,dp)
      x(i,2)=sin(real(i,dp))
      if(i<=10)then
         y(i)=1.0_dp
      else
         y(i)=5.0_dp
      end if
   end do

   control=rpart_make_control(minsplit=4,minbucket=2,cp=0.0_dp,xval=0,stat=stat)
   if(stat/=0)error stop 'invalid control values'
   call rpart_fit_regression(x,y,model,control=control,stat=stat)
   if(stat/=0)error stop 'fit failed'
   call rpart_predict_values(model,x,pred)

   write(*,'(a,i0)')'nodes: ',count_nodes(model%root)
   write(*,'(a,f8.3)')'first prediction: ',pred(1)
   write(*,'(a,f8.3)')'last prediction:  ',pred(n)
end program basic_example
