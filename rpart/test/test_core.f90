program test_core
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
   use rpart_kinds, only : dp, i8
   use rpart_types
   use rpart_fit
   use rpart_predict
   use rpart_cp, only : prune_model
   use rpart_control_api, only : rpart_make_control
   implicit none
   integer, parameter :: n=20
   real(dp) :: x(n,2), y(n), pred(n), prob(n,2), time(n), event(n)
   real(dp) :: nan
   integer :: cls(n), pc(n), i, stat
   type(rpart_model) :: m
   type(rpart_control) :: ctl
   integer :: ncat2(2)

   nan=ieee_value(0.0_dp,ieee_quiet_nan)
   ctl%minsplit=2;ctl%minbucket=2;ctl%cp=0.0_dp;ctl%xval=0;ctl%maxdepth=5
   do i=1,n
      x(i,1)=real(i,dp);x(i,2)=real(i,dp)+0.1_dp
      if(i<=10)then;y(i)=0.0_dp;cls(i)=1;event(i)=1.0_dp
      else;y(i)=10.0_dp;cls(i)=2;event(i)=5.0_dp
      end if
      time(i)=1.0_dp
   end do

   call rpart_fit_regression(x,y,m,control=ctl,stat=stat)
   call check(stat==0,'regression fit status')
   call check(allocated(m%root%left),'regression root split')
   call check(abs(m%root%primary(1)%spoint-10.5_dp)<1.0e-12_dp,'regression split point')
   call rpart_predict_values(m,x,pred)
   call check(maxval(abs(pred-y))<1.0e-12_dp,'regression predictions')
   call check(size(m%cptable)>=2,'cp table')

   call rpart_fit_classification(x,cls,m,control=ctl,stat=stat)
   call check(stat==0,'classification fit status')
   call rpart_predict_class(m,x,pc)
   call check(all(pc==cls),'classification predictions')
   call rpart_predict_proba(m,x,prob)
   call check(minval(prob([(i,i=1,10)],1))>0.999999_dp,'classification probabilities left')
   call check(minval(prob([(i,i=11,20)],2))>0.999999_dp,'classification probabilities right')

   call rpart_fit_poisson(x,time,event,m,control=ctl,shrink=0.0_dp,stat=stat)
   call check(stat==0,'poisson fit status')
   call rpart_predict_values(m,x,pred)
   call check(maxval(abs(pred(1:10)-1.0_dp))<1.0e-12_dp,'poisson left rate')
   call check(maxval(abs(pred(11:20)-5.0_dp))<1.0e-12_dp,'poisson right rate')

   ! Categorical predictor: levels 1 and 3 share response 0, level 2 response 10.
   do i=1,n
      x(i,1)=real(1+mod(i-1,3),dp)
      if(nint(x(i,1))==2)then;y(i)=10.0_dp;else;y(i)=0.0_dp;end if
      x(i,2)=real(i,dp)
   end do
   ncat2=[3,0]
   call rpart_fit_regression(x,y,m,control=ctl,ncat=ncat2,cost=[1.0_dp,100.0_dp],stat=stat)
   call check(stat==0,'categorical fit status')
   call check(m%root%primary(1)%ncat==3,'categorical split selected')
   call rpart_predict_values(m,x,pred)
   call check(maxval(abs(pred-y))<1.0e-12_dp,'categorical predictions')

   ! Missing primary values should be recovered by the second predictor as surrogate.
   do i=1,n
      x(i,1)=real(i,dp);x(i,2)=real(i,dp);y(i)=merge(0.0_dp,10.0_dp,i<=10)
   end do
   x(3,1)=nan;x(18,1)=nan
   ctl%maxsurrogate=2;ctl%usesurrogate=2
   call rpart_fit_regression(x,y,m,control=ctl,cost=[1.0_dp,2.0_dp],stat=stat)
   call check(stat==0,'surrogate fit status')
   call check(allocated(m%root%surrogate),'surrogate allocated')
   call rpart_predict_values(m,x,pred)
   call check(abs(pred(3)-0.0_dp)<1.0e-12_dp.and.abs(pred(18)-10.0_dp)<1.0e-12_dp,'surrogate routing')

   ctl%xval=4;ctl%cp=0.001_dp
   call rpart_fit_regression(x,y,m,control=ctl,seed=123_i8,stat=stat)
   call check(stat==0,'xval fit status')
   call check(all([(ieee_is_finite(m%cptable(i)%xerror),i=1,size(m%cptable))]),'xval errors finite')

   ! Zero-risk trees must stay valid under pruning/prediction CP comparisons.
   x(:,1)=[(real(i,dp),i=1,n)];x(:,2)=x(:,1);y=3.0_dp
   ctl%minsplit=2;ctl%minbucket=1;ctl%cp=0.0_dp;ctl%xval=0
   call rpart_fit_regression(x,y,m,control=ctl,stat=stat)
   call check(stat==0,'zero-risk fit status')
   call rpart_predict_values(m,x,pred,cp=1.0_dp)
   call check(maxval(abs(pred-3.0_dp))<1.0e-12_dp,'zero-risk prediction')

   ctl=rpart_make_control(minsplit=9,xval=0,stat=stat)
   call check(stat==0,'control constructor status')
   call check(ctl%minbucket==3,'minbucket follows minsplit')
   ctl=rpart_make_control(minbucket=4,xval=0,stat=stat)
   call check(ctl%minsplit==12,'minsplit follows explicit minbucket')

   write(*,'(a)') 'test_core: PASS'
contains
   subroutine check(ok,msg)
      logical,intent(in)::ok
      character(*),intent(in)::msg
      if(.not.ok)then
         write(*,'(a)') 'FAIL: '//trim(msg)
         error stop 1
      end if
   end subroutine check
end program test_core
