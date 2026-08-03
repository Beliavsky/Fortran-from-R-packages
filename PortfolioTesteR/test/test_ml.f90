! SPDX-License-Identifier: MIT
! PortfolioTesteR modern Fortran translation
program test_ml
  use portfolio_tester
  implicit none
  real(dp)::x(100,2),y(100)
  real(dp),allocatable::pred(:),scores(:,:),labels(:,:),ic(:),prices(:,:),mom(:,:),lagmom(:,:)
  real(dp),allocatable::features(:,:,:)
  type(linear_model)::model
  integer::i,st
  do i=1,100
    x(i,1)=real(i,dp)/100.0_dp;x(i,2)=sin(0.1_dp*real(i,dp))
    y(i)=1.5_dp+2.0_dp*x(i,1)-0.75_dp*x(i,2)
  end do
  call fit_linear_model(x,y,model,lambda=0.0_dp,status=st)
  call assert_true(st==0.and.model%fitted,'linear fit')
  call predict_linear_model(model,x,pred)
  call assert_true(maxval(abs(pred-y))<1.0e-10_dp,'linear exact prediction')
  call generate_sample_prices(100,5,prices,999_i8)
  call calc_momentum(prices,4,mom);call panel_lag(mom,1,lagmom)
  call make_labels(prices,2,1,labels)
  allocate(features(100,5,1));features(:,:,1)=lagmom
  call rolling_fit_predict(features,labels,40,5,5,scores,lambda=0.1_dp)
  call ic_series(scores,labels,ic,.true.)
  call assert_true(count(is_finite(scores))>0,'rolling predictions')
  call assert_true(count(is_finite(ic))>0,'information coefficients')
  print '(a)','test_ml: PASS'
contains
  subroutine assert_true(ok,msg)
    logical,intent(in)::ok;character(len=*),intent(in)::msg
    if(.not.ok)then;write(*,'(a)')'FAIL: '//msg;error stop 1;end if
  end subroutine
end program test_ml
