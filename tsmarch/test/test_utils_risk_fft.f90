program test_utils_risk_fft
  use ghyp_kinds, only : dp, i8
  use tsmarch
  use test_support
  implicit none
  real(dp) :: x(100,2), r(2,2), sigma(2), w(2), probs(2), dens(3,2)
  real(dp), allocatable :: cov3(:,:,:), lw(:,:), sim(:,:,:), var(:,:), es(:,:), draws(:,:), comb_real(:), cov2(:,:)
  integer, allocatable :: comb(:,:)
  type(escc_result) :: test
  type(fft_distribution) :: fd
  integer :: i

  do i=1,100
    x(i,1)=sin(0.15_dp*real(i,dp))
    x(i,2)=cos(0.19_dp*real(i,dp))
  end do
  cov3=ewma_covariance(x,0.94_dp)
  call assert_true(all([(cov3(i,i,100)>0.0_dp,i=1,2)]),'EWMA positive diagonal')
  lw=lw_covariance(x)
  call assert_true(lw(1,1)>0.0_dp .and. lw(2,2)>0.0_dp,'LW covariance')
  r=reshape([1.0_dp,0.4_dp,0.4_dp,1.0_dp],[2,2])
  sigma=[0.2_dp,0.3_dp]
  cov2=cor2cov(r,sigma)
  call assert_close(cov2(1,2),0.024_dp,1.0e-12_dp,'correlation to covariance')
  draws=rmvnorm(200,[0.0_dp,0.0_dp],r,seed=7_i8)
  call assert_true(size(draws,1)==200,'multivariate normal draws')

  comb=combn_fast(5,3)
  call assert_true(size(comb,2)==10,'combination count')
  allocate(comb_real(size(comb,2)))
  comb_real=real(comb(1,:),dp)
  call assert_close(minval(comb_real),1.0_dp,0.0_dp,'combination values')

  test=escc_test(x,lags=1)
  call assert_true(test%status==tsm_success,'ESCC test')
  call assert_true(test%p_value>=0.0_dp .and. test%p_value<=1.0_dp,'ESCC p-value')

  allocate(sim(2,2,1000))
  do i=1,1000
    sim(1,1,i)=-1.0_dp+2.0_dp*real(i-1,dp)/999.0_dp
    sim(1,2,i)=0.5_dp*sim(1,1,i)
    sim(2,:,i)=sim(1,:,i)+0.1_dp
  end do
  probs=[0.05_dp,0.10_dp]
  w=[0.6_dp,0.4_dp]
  var=value_at_risk(sim,probs,w)
  es=expected_shortfall(sim,probs,w)
  call assert_true(all(es<=var+1.0e-12_dp),'expected shortfall below VaR')

  dens(:,1)=[0.25_dp,0.5_dp,0.25_dp]
  dens(:,2)=[0.25_dp,0.5_dp,0.25_dp]
  fd=fft_convolution(dens,-1.0_dp,1.0_dp)
  call assert_true(fd%status==tsm_success,'FFT convolution')
  call assert_close(sum(fd%density),1.0_dp,1.0e-10_dp,'FFT density normalized')
  call assert_true(pfft(fd,0.0_dp)>0.5_dp-0.2_dp .and. pfft(fd,0.0_dp)<0.8_dp,'FFT CDF')
  call assert_true(abs(qfft(fd,0.5_dp))<1.0_dp,'FFT median')

  call finish_test('test_utils_risk_fft')
end program test_utils_risk_fft
