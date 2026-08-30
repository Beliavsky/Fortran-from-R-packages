program test_priors
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use rpart
   implicit none
   type(rpart_model) :: model
   type(rpart_control) :: ctl
   real(dp) :: x(15,3), dummy(15), prior(3), loss(3,3), nanv
   real(dp), parameter :: tol = 2.0e-12_dp
   integer :: y(15), yorig(15), i, stat

   dummy=[3,1,4,1,5,9,2,6,5,3,5,8,9,7,9]/5.0_dp
   do i=1,15
      y(i)=1+mod(i-1,3)
      yorig(i)=y(i)
      x(i,1)=real(i,dp)
   end do
   y(15)=1
   x(:,2)=[(real(1+mod(i-1,6),dp),i=1,12),1.0_dp,2.0_dp,3.0_dp]
   x(:,3)=(real(yorig,dp)+dummy)*10.0_dp
   nanv=ieee_value(0.0_dp,ieee_quiet_nan)
   x(1,3)=nanv;x(5,3)=nanv;x(10,3)=nanv

   prior=[0.2_dp,0.3_dp,0.5_dp]
   loss=0.0_dp
   loss(1,:)=[0.0_dp,2.0_dp,2.0_dp]
   loss(2,:)=[2.0_dp,0.0_dp,6.0_dp]
   loss(3,:)=[1.0_dp,1.0_dp,0.0_dp]
   ctl%cp=0.0_dp;ctl%xval=0;ctl%minsplit=5;ctl%minbucket=2;ctl%maxdepth=1

   call rpart_fit_classification(x,y,model,ctl,prior=prior,loss=loss,stat=stat)
   call check(stat==0,'classification priors fit status')
   call check(abs(model%root%risk-13.5_dp)<tol,'root loss')
   call check(model%root%primary(1)%var==3,'best variable')
   call check(abs(model%root%primary(1)%improve-3.9876304669118356_dp)<tol,'x3 split improvement')
   call check(abs(model%root%primary(2)%improve-1.9121161978304828_dp)<tol,'x2 split improvement')
   call check(abs(model%root%primary(3)%improve-0.29049457620886265_dp)<tol,'x1 split improvement')
   call check(abs(model%root%left%risk-4.0_dp)<tol,'left loss')
   call check(abs(model%root%right%risk-7.4_dp)<tol,'right loss')
   call check(all(abs(model%root%left%response(2:4)-[4.0_dp,4.0_dp,0.0_dp])<tol),'left counts')
   call check(all(abs(model%root%right%response(2:4)-[2.0_dp,1.0_dp,4.0_dp])<tol),'right counts')

   ctl%maxdepth=2
   call rpart_fit_classification(x,y,model,ctl,prior=prior,loss=loss,split_rule=RPART_INFORMATION,stat=stat)
   call check(stat==0,'information split fit status')
   call check(allocated(model%root%left),'information split produces tree')
   call check(model%root%primary(1)%improve>0.0_dp,'information improvement positive')

   print '(a)','test_priors: PASS'
contains
   subroutine check(ok,msg)
      logical,intent(in)::ok
      character(*),intent(in)::msg
      if(.not.ok)then
         write(*,'(a)')'FAIL: '//trim(msg)
         error stop 1
      end if
   end subroutine check
end program test_priors
