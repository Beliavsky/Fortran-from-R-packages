program test_new
   use grf
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   implicit none
   integer, parameter :: n = 80
   real(dp) :: x(n,2), y(n), w(n), s_time(n), priorities(n)
   real(dp) :: wm(n,2)
   integer :: event(n), clusters(n)
   integer :: i, j, info
   type(grf_options) :: opt
   type(grf_forest) :: f, fm, fs, fc
   type(grf_ate_result) :: ate
   type(grf_rate_result) :: rate
   real(dp), allocatable :: scores(:,:), pred(:,:)
   real(dp) :: grid(4), eval(4)

   do i=1,n
      x(i,1)=real(i,dp)/real(n,dp)
      x(i,2)=sin(0.17_dp*real(i,dp))
      y(i)=1.0_dp+2.0_dp*x(i,1)+0.1_dp*x(i,2)
      w(i)=merge(1.0_dp,0.0_dp,mod(i,2)==0)
      s_time(i)=0.6_dp+0.012_dp*real(i,dp)
      event(i)=merge(1,0,mod(i,4)/=0)
      clusters(i)=1+(i-1)/10
      priorities(i)=x(i,1)
      wm(i,:)=0.0_dp
      select case(mod(i,3))
      case(1)
         wm(i,1)=1.0_dp
      case(2)
         wm(i,2)=1.0_dp
      end select
   end do
   opt%num_trees=12
   opt%honesty=.false.
   opt%sample_fraction=0.5_dp
   opt%clusters=clusters
   opt%equalize_cluster_weights=.true.
   call regression_forest(x,y,f,info,opt)
   if(info/=0) error stop 'cluster regression'
   if(f%n_clusters/=8 .or. f%samples_per_cluster/=10) error stop 'cluster metadata'
   do i=1,size(f%trees)
      if(any([(any(f%trees(i)%inbag((j-1)*10+1:j*10) .neqv. f%trees(i)%inbag((j-1)*10+1)), j=1,8)])) then
         error stop 'cluster inbag'
      end if
   end do

   opt%equalize_cluster_weights=.false.
   opt%clusters=clusters
   opt%num_trees=20
   opt%honesty=.true.
   call causal_forest(x,y,w,fc,info,opt)
   if(info/=0) error stop 'cluster causal'
   call average_treatment_effect(fc,ate,info)
   if(info/=0 .or. any(ieee_is_nan(ate%std_err))) error stop 'cluster ate'
   call rank_average_treatment_effect(fc,priorities,rate,info,bootstrap_reps=10,seed=9)
   if(info/=0 .or. ieee_is_nan(rate%std_err)) error stop 'cluster rate'

   opt=grf_options()
   opt%num_trees=20
   opt%honesty=.false.
   call multi_arm_causal_forest(x,y,wm,fm,info,opt)
   if(info/=0) error stop 'multi arm'
   if(.not.fm%multi_arm_categorical) error stop 'categorical flag'
   if(.not.allocated(fm%arm_propensity_hat)) error stop 'propensity'
   if(maxval(abs(sum(fm%arm_propensity_hat,dim=2)-1.0_dp))>1.0e-8_dp) error stop 'prop sum'
   call get_scores_multi_arm_causal_forest(fm,scores)
   if(size(scores,2)/=2) error stop 'multi scores'

   grid=[0.7_dp,1.0_dp,1.3_dp,1.6_dp]
   eval=[0.5_dp,0.9_dp,1.2_dp,1.5_dp]
   opt=grf_options()
   opt%num_trees=20
   opt%honesty=.false.
   opt%fast_logrank=.true.
   call survival_forest(x,s_time,event,fs,info,opt,failure_times=grid)
   if(info/=0) error stop 'survival'
   if(minval(fs%survival_time_index)<0 .or. maxval(fs%survival_time_index)>4) error stop 'surv index'
   call predict_survival_forest(fs,x(1:4,:),pred,failure_times=eval,individual_times=.true.)
   if(any(shape(pred)/=[4,1])) error stop 'surv shape'
   if(any(pred<0.0_dp) .or. any(pred>1.0_dp)) error stop 'surv range'
   print *, 'new feature smoke tests passed'
end program test_new
