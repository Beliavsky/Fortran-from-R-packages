! SPDX-License-Identifier: GPL-3.0-only
program test_ordinal_bandwidth
  use mass
  use test_support
  implicit none
  integer, parameter :: n=72
  real(dp) :: x(n,1), sample(n), h_ucv, h_bcv, h_sj
  integer :: y(n), predicted_correct, i, status
  integer, allocatable :: predicted(:)
  real(dp), allocatable :: probabilities(:,:)
  type(ordinal_model) :: model
  real(dp) :: parameters(3), exact_y(n)

  do i=1,n
    x(i,1)=-3.0_dp+6.0_dp*real(i-1,dp)/real(n-1,dp)
    sample(i)=sin(0.31_dp*real(i,dp))+0.15_dp*cos(0.91_dp*real(i,dp))
    if(x(i,1)+0.35_dp*sin(real(i,dp)) < -0.8_dp)then
      y(i)=1
    else if(x(i,1)+0.35_dp*sin(real(i,dp)) < 0.8_dp)then
      y(i)=2
    else
      y(i)=3
    end if
  end do
  call polr_fit(x,y,model,maxit=600)
  call assert_true(model%status==mass_success,'polr fit status')
  call polr_predict(model,x,predicted,probabilities,status=status)
  predicted_correct=count(predicted==y)
  call assert_true(status==mass_success .and. predicted_correct>=60,'polr prediction')
  call assert_true(maxval(abs(sum(probabilities,dim=2)-1.0_dp))<1.0e-10_dp, &
    'polr probabilities sum one')

  call ucv_bandwidth(sample,h_ucv,status)
  call assert_true(status==mass_success .and. h_ucv>0.0_dp,'ucv bandwidth')
  call bcv_bandwidth(sample,h_bcv,status)
  call assert_true(status==mass_success .and. h_bcv>0.0_dp,'bcv bandwidth')
  call width_sj(sample,h_sj,status,method='dpi')
  call assert_true(status==mass_success .and. h_sj>0.0_dp,'SJ bandwidth')

  exact_y=negative_exponential(x(:,1)+4.0_dp,2.0_dp,5.0_dp,1.8_dp)
  call negexp_initial(x(:,1)+4.0_dp,exact_y,parameters,status)
  call assert_true(status==mass_success,'negative exponential initialization')
  call assert_close(parameters(1),2.0_dp,0.2_dp,'negative exponential baseline')
  call assert_close(parameters(3),1.8_dp,0.6_dp,'negative exponential theta')
  write(*,'(a)') 'test_ordinal_bandwidth: PASS'
end program test_ordinal_bandwidth
