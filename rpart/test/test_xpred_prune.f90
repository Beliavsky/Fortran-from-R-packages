program test_xpred_prune
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use rpart
   implicit none
   type(rpart_model) :: model
   type(rpart_control) :: ctl
   real(dp) :: x(24,2),y(24),pred(24,3),cp(3),p0(24),p1(24),cost(2),expect
   integer :: groups(24),i,stat,n_before,n_after

   do i=1,24
      x(i,1)=real(i,dp)
      x(i,2)=x(i,1)
      if(i<=8)then;y(i)=1.0_dp
      else if(i<=16)then;y(i)=5.0_dp
      else;y(i)=10.0_dp
      end if
      groups(i)=1+mod(i-1,4)
   end do
   ctl%cp=0.0_dp;ctl%xval=0;ctl%minsplit=4;ctl%minbucket=2;ctl%maxdepth=4
   cost=[1.0_dp,2.0_dp]
   call rpart_fit_regression(x,y,model,ctl,cost=cost,stat=stat)
   call check(stat==0,'regression fit')
   call check(model%root%primary(1)%var==1,'cost-adjusted primary split')

   cp=[0.0_dp,0.05_dp,1.0_dp]
   call rpart_xpred_regression(model,x,y,groups,pred,cp=cp,stat=stat)
   call check(stat==0,'xpred status')
   call check(all(ieee_is_finite(pred)),'xpred finite')
   do i=1,24
      expect=sum(pack(y,groups/=groups(i)))/real(count(groups/=groups(i)),dp)
      call check(abs(pred(i,3)-expect)<1.0e-12_dp,'root-only xpred equals fold mean')
   end do

   call rpart_predict_values(model,x,p0)
   n_before=count_nodes(model%root)
   call prune_model(model,0.20_dp)
   n_after=count_nodes(model%root)
   call check(n_after<=n_before,'pruning reduces tree')
   call rpart_predict_values(model,x,p1)
   call check(all(ieee_is_finite(p1)),'pruned predictions finite')

   print '(a)','test_xpred_prune: PASS'
contains
   subroutine check(ok,msg)
      logical,intent(in)::ok
      character(*),intent(in)::msg
      if(.not.ok)then
         write(*,'(a)')'FAIL: '//trim(msg)
         error stop 1
      end if
   end subroutine check
end program test_xpred_prune
