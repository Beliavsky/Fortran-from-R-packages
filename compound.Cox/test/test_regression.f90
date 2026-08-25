program test_regression
  use compound_cox
  implicit none
  integer, parameter :: n=18, p=2
  real(dp) :: t(n), x(n,p)
  integer :: d(n), i
  type(compound_result) :: cr
  type(depend_cox_result) :: dr
  type(depend_cv_result) :: dcv
  type(selection_result) :: sr

  do i=1,n
    x(i,1)=real(i-9,dp)/5.0_dp
    x(i,2)=sin(0.7_dp*real(i,dp))
    t(i)=20.0_dp-exp(0.35_dp*x(i,1))+0.12_dp*real(mod(i,3),dp)
    d(i)=1
    if(mod(i,5)==0)d(i)=0
  end do

  call compound_reg(t,d,x,cr,kfold=3,delta_a=0.2_dp,a0=0.0_dp,with_variance=.true.)
  if(size(cr%beta)/=p) error stop 'compound beta size'
  if(cr%a<0.0_dp .or. cr%a>1.0_dp) error stop 'compound a'
  if(any(abs(cr%beta)>100.0_dp)) error stop 'compound beta finite'
  if(.not.allocated(cr%se)) error stop 'compound variance path'
  if(any(cr%se < 0.0_dp)) error stop 'compound se'

  call depend_cox_reg(t,d,x(:,1),2.0_dp,dr,with_variance=.true.,baseline=.true.)
  if(.not.allocated(dr%baseline)) error stop 'depend baseline'
  if(any(dr%baseline(2:)<dr%baseline(:n-1))) error stop 'depend baseline monotone'
  if(abs(dr%beta)>100.0_dp) error stop 'depend beta finite'
  if(dr%se <= 0.0_dp) error stop 'depend se'

  call depend_cox_reg_cv(t,d,x,dcv,kfold=3,ngrid=3)
  if(dcv%alpha<=0.0_dp) error stop 'depend cv alpha'
  if(dcv%c_index<0.0_dp .or. dcv%c_index>1.0_dp) error stop 'depend cv cindex'

  call uni_selection(t,d,x,sr,p_value=1.0_dp,kfold=3,use_score=.true.,permutation=.false.)
  if(sr%n_selected/=p) error stop 'selection all variables'
  if(sr%c_index_full<0.0_dp .or. sr%c_index_full>1.0_dp) error stop 'selection cindex'

  print '(a)', 'test_regression: PASS'
end program test_regression
