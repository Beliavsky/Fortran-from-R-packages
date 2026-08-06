program test_dcc
  use ghyp_kinds, only : dp, i8
  use tsmarch
  use test_support
  implicit none
  type(dcc_spec) :: spec
  type(dcc_parameters) :: p
  type(dcc_simulation) :: sim
  type(dcc_filter_result) :: filtered, constant_fit
  real(dp) :: qbar(2,2)
  real(dp), allocatable :: sigma(:, :), z(:, :)
  integer :: t

  qbar = reshape([1.0_dp, 0.45_dp, 0.45_dp, 1.0_dp], [2,2])
  spec%distribution = 'mvn'
  spec%alpha_order = 1
  spec%gamma_order = 1
  spec%beta_order = 1
  allocate(p%alpha(1), p%gamma(1), p%beta(1))
  p%alpha = 0.04_dp
  p%gamma = 0.02_dp
  p%beta = 0.90_dp
  p%shape = 8.0_dp

  sim = simulate_dcc_innovations(spec, p, qbar, 180, paths=1, burn=100, seed=17_i8)
  call assert_true(sim%status == tsm_success, 'DCC simulation succeeds')
  z = sim%innovations(:,:,1)
  allocate(sigma(size(z,1),size(z,2)))
  sigma = 1.0_dp
  filtered = dcc_filter(z, sigma, spec, p)
  call assert_true(filtered%status == tsm_success, 'dynamic DCC filter succeeds')
  call assert_true(filtered%log_likelihood < 0.0_dp, 'dynamic DCC log likelihood is finite')
  do t = 1, size(z,1)
    call assert_close(filtered%correlation(1,1,t), 1.0_dp, 1.0e-10_dp, 'DCC diagonal one')
    call assert_close(filtered%correlation(2,2,t), 1.0_dp, 1.0e-10_dp, 'DCC diagonal one')
    call assert_true(abs(filtered%correlation(1,2,t)) < 1.0_dp, 'DCC correlation bounded')
  end do
  call assert_true(dcc_stationarity(p, filtered%qbar, filtered%nbar) < 1.0_dp, 'DCC stationarity')

  spec%constant_correlation = .true.
  constant_fit = dcc_filter(z, sigma, spec, p)
  call assert_true(constant_fit%status == tsm_success, 'constant correlation filter succeeds')
  call assert_close(constant_fit%correlation(1,2,1), constant_fit%correlation(1,2,size(z,1)), &
    1.0e-12_dp, 'constant correlation remains constant')

  call finish_test('test_dcc')
end program test_dcc
