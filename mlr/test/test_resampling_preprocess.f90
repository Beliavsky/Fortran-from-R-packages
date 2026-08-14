program test_resampling_preprocess
  use mlr_kinds, only : dp, i8
  use mlr_rng, only : rng_state, rng_seed
  use mlr_types, only : resample_plan, scaler_model
  use mlr_resampling
  use mlr_preprocess
  use mlr_smote
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  implicit none
  type(rng_state)::rng
  type(resample_plan)::plan
  type(scaler_model)::sc
  real(dp)::x(10,2),xm(4,2),minor(4,2)
  real(dp),allocatable::z(:,:),syn(:,:)
  integer,allocatable::nn(:,:)
  logical::isnum(2)
  integer::i
  do i=1,10;x(i,1)=real(i,dp);x(i,2)=2.0_dp*real(i,dp);end do
  call standardize(x,z,sc)
  if(maxval(abs(sum(z,dim=1)))>1.0e-12_dp)error stop 'standardize center'
  xm=reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp,10.0_dp,20.0_dp,30.0_dp,40.0_dp],[4,2])
  xm(2,1)=ieee_value(0.0_dp,ieee_quiet_nan);call impute_mean(xm)
  if(abs(xm(2,1)-(1.0_dp+3.0_dp+4.0_dp)/3.0_dp)>1.0e-12_dp)error stop 'impute'
  call rng_seed(rng,77_i8);call make_kfold(10,5,rng,plan)
  if(size(plan%test)/=5)error stop 'kfold size'
  do i=1,5
    if(size(plan%test(i)%idx)/=2.or.size(plan%train(i)%idx)/=8)error stop 'kfold fold'
  end do
  minor=reshape([0.0_dp,0.1_dp,0.2_dp,0.3_dp,1.0_dp,1.0_dp,2.0_dp,2.0_dp],[4,2])
  isnum=[.true.,.false.];call nearest_neighbors(minor,2,nn);call smote_generate(minor,isnum,nn,20,rng,syn)
  if(any(syn(:,1)<0.0_dp).or.any(syn(:,1)>0.3_dp))error stop 'smote numeric'
  if(any((abs(syn(:,2)-1.0_dp)>1.0e-12_dp).and.(abs(syn(:,2)-2.0_dp)>1.0e-12_dp)))error stop 'smote factor'
  print *, 'test_resampling_preprocess: PASS'
end program test_resampling_preprocess
