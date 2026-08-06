program test_ica_gogarch
  use ghyp_kinds, only : dp, i8
  use tsmarch
  use test_support
  implicit none
  type(ica_result) :: ica, rad
  real(dp), allocatable :: source(:, :), mixed(:, :), reconstructed(:, :), cov(:, :)
  real(dp) :: mixing(2,2), fv(2), fs(2), fk(2), w(2)
  real(dp), allocatable :: c3(:,:,:), c4(:,:,:,:), gcov(:,:)
  integer :: i,n

  n=300
  allocate(source(n,2),mixed(n,2))
  do i=1,n
    source(i,1)=sin(0.23_dp*real(i,dp))+0.15_dp*sin(0.71_dp*real(i,dp))
    source(i,2)=sign(1.0_dp,sin(0.13_dp*real(i,dp)))+0.2_dp*cos(0.31_dp*real(i,dp))
  end do
  source(:,1)=source(:,1)-sum(source(:,1))/real(n,dp)
  source(:,2)=source(:,2)-sum(source(:,2))/real(n,dp)
  mixing=reshape([1.0_dp,0.45_dp,0.35_dp,1.2_dp],[2,2])
  mixed=matmul(source,transpose(mixing))
  ica=fastica(mixed,max_iterations=800,tolerance=1.0e-8_dp,seed=41_i8)
  call assert_true(ica%status==tsm_success .or. ica%status==tsm_no_convergence,'FastICA returns a solution')
  reconstructed=spread(ica%mean,1,n)+matmul(ica%components,transpose(ica%mixing))
  call assert_true(maxval(abs(reconstructed-mixed))<1.0e-7_dp,'FastICA reconstructs data')
  cov=sample_covariance(ica%components)
  call assert_true(abs(cov(1,2))<1.0e-6_dp,'FastICA components decorrelated')

  rad=radical(mixed,sweeps=2,angle_grid=21,seed=19_i8)
  call assert_true(rad%status==tsm_success,'RADICAL returns a solution')
  reconstructed=spread(rad%mean,1,n)+matmul(rad%components,transpose(rad%mixing))
  call assert_true(maxval(abs(reconstructed-mixed))<1.0e-7_dp,'RADICAL reconstructs data')

  fv=[0.5_dp,1.2_dp]
  fs=[0.3_dp,-0.2_dp]
  fk=[1.5_dp,0.8_dp]
  gcov=gogarch_covariance(fv,mixing)
  call assert_true(gcov(1,1)>0.0_dp .and. gcov(2,2)>0.0_dp,'GO-GARCH covariance positive diagonal')
  c3=gogarch_coskewness(mixing,fs,fv)
  c4=gogarch_cokurtosis(mixing,fk,fv)
  w=[0.6_dp,0.4_dp]
  call assert_true(portfolio_variance(gcov,w)>0.0_dp,'portfolio variance positive')
  call assert_true(abs(portfolio_skewness(c3,gcov,w))<10.0_dp,'portfolio skewness finite')
  call assert_true(portfolio_kurtosis(c4,gcov,w)>0.0_dp,'portfolio kurtosis positive')

  call finish_test('test_ica_gogarch')
end program test_ica_gogarch
