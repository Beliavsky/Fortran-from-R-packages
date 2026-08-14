program test_learners
  use mlr_kinds, only : dp, i8
  use mlr_rng, only : rng_state, rng_seed
  use mlr_types, only : linear_model, logistic_model, kmeans_model, knn_model
  use mlr_learners
  implicit none
  real(dp) :: x(8,2), y(8), lx(10,1)
  integer :: ly(10), i
  type(linear_model) :: lm
  type(logistic_model) :: logm
  type(kmeans_model) :: km
  type(knn_model) :: knn
  type(rng_state) :: rng
  real(dp), allocatable :: pred(:), prob(:)
  integer, allocatable :: cls(:)
  do i=1,8
    x(i,1)=real(i-1,dp);x(i,2)=real(mod(i,3),dp);y(i)=2.0_dp+3.0_dp*x(i,1)-0.5_dp*x(i,2)
  end do
  call fit_linear_regression(x,y,lm)
  call predict_linear_regression(lm,x,pred)
  if(maxval(abs(pred-y))>1.0e-9_dp)error stop 'linear regression'
  do i=1,10;lx(i,1)=real(i-5,dp);ly(i)=merge(1,0,lx(i,1)>0.0_dp);end do
  call fit_logistic_regression(lx,ly,logm,ridge=1.0e-4_dp,maxiter=200)
  call predict_logistic_probability(logm,lx,prob)
  if(prob(1)>=0.5_dp.or.prob(10)<=0.5_dp)error stop 'logistic'
  x(1:4,1)=0.0_dp;x(1:4,2)=[0.0_dp,0.1_dp,-0.1_dp,0.05_dp]
  x(5:8,1)=5.0_dp;x(5:8,2)=[5.0_dp,5.1_dp,4.9_dp,5.05_dp]
  call rng_seed(rng,123_i8);call fit_kmeans(x,2,km,rng,maxiter=100)
  call predict_kmeans(km,x,cls)
  if(count(cls==cls(1))/=4)error stop 'kmeans group1'
  if(count(cls==cls(5))/=4.or.cls(1)==cls(5))error stop 'kmeans group2'
  y=[1.0_dp,1.1_dp,0.9_dp,1.05_dp,5.0_dp,5.1_dp,4.9_dp,5.05_dp]
  call fit_knn_regression(x,y,2,knn);call predict_knn_regression(knn,x,pred)
  if(maxval(abs(pred-y))>0.2_dp)error stop 'knn regression'
  print *, 'test_learners: PASS'
end program test_learners
