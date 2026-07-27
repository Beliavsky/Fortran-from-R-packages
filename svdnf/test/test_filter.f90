! SPDX-License-Identifier: GPL-3.0-only
program test_filter
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use svdnf
  use test_support
  implicit none
  type(svm_dynamics) :: dynamics, capm
  type(simulation_result) :: simulated
  type(filter_result) :: filtered, factor_filtered
  type(percentile_result) :: median, prediction
  real(dp), allocatable :: factors(:,:), returns(:)
  integer :: i

  dynamics=dynamics_svm('Taylor',phi=0.93_dp,theta=-1.2_dp,sigma=0.22_dp)
  simulated=model_simulate(dynamics,60,initial_volatility=-1.2_dp,seed=2026)
  call assert_true(simulated%ok,'Taylor simulation')
  filtered=dnf_filter(dynamics,simulated%returns,n=24)
  call assert_true(filtered%ok,'Taylor filtering')
  call assert_true(ieee_is_finite(filtered%log_likelihood),'finite log likelihood')
  do i=1,size(filtered%filter_grid,2)
    call assert_close(sum(filtered%filter_grid(:,i)),1.0_dp,2.0e-12_dp,'filter normalization')
  end do
  call assert_true(all(filtered%likelihoods>0.0_dp),'positive likelihoods')

  median=extract_vol_percentile(filtered,0.5_dp)
  prediction=extract_vol_percentile(filtered,0.5_dp,.true.)
  call assert_true(median%ok .and. prediction%ok,'percentile extraction')
  call assert_true(size(median%values)==61,'filter percentile length')
  call assert_true(size(prediction%values)==60,'prediction percentile length')
  call assert_true(all(median%values>=filtered%grids%var_mid_points(1)),'percentile lower bound')
  call assert_true(all(median%values<=filtered%grids%var_mid_points(24)),'percentile upper bound')

  allocate(factors(60,2),returns(60))
  do i=1,60
    factors(i,1)=1.0_dp
    factors(i,2)=sin(0.1_dp*real(i,dp))
  end do
  capm=dynamics_svm('CAPM_SV',phi=0.9_dp,theta=-1.0_dp,sigma=0.2_dp,coefs=[0.001_dp,0.2_dp])
  simulated=model_simulate(capm,60,initial_volatility=-1.0_dp,seed=77)
  returns=simulated%returns+matmul(factors,capm%coefs)
  factor_filtered=dnf_filter(capm,returns,factors=factors,n=20)
  call assert_true(factor_filtered%ok,'factor-adjusted filtering')

  write(*,'(a)') 'test_filter: PASS'
end program test_filter
