! SPDX-License-Identifier: Apache-2.0
program test_fixed_acceleration
  use intraday_model
  use test_support, only : check, make_parameters
  implicit none
  type(volume_parameters) :: truth
  type(volume_model_spec) :: spec
  type(volume_model) :: model, one_step
  type(volume_fit_control) :: control
  real(dp), allocatable :: volume(:, :), states(:, :)

  truth = make_parameters(8)
  call simulate_intraday_volume(truth, 70, volume, states, seed=77)
  call initialize_volume_spec(spec, 8)
  spec%fixed%a_mu = truth%a_mu
  spec%fixed%var_mu = truth%var_mu
  spec%is_fixed%a_mu = .true.
  spec%is_fixed%var_mu = .true.
  spec%initial%a_eta = 0.5_dp
  spec%has_initial%a_eta = .true.

  control%acceleration = .true.
  control%maxit = 250
  control%abstol = 2.0e-5_dp
  control%save_history = .true.
  call fit_volume(volume, model, spec=spec, control=control)

  call check(model%status == intraday_ok, 'accelerated EM status')
  call check(model%converged, 'accelerated EM convergence')
  call check(abs(model%par%a_mu - truth%a_mu) < 1.0e-15_dp, 'fixed a_mu unchanged')
  call check(abs(model%par%var_mu - truth%var_mu) < 1.0e-15_dp, 'fixed var_mu unchanged')
  call check(abs(model%par%a_eta - truth%a_eta) < 0.16_dp, 'accelerated daily persistence')

  control%maxit = 1
  control%abstol = 1.0e-14_dp
  control%save_history = .false.
  call fit_volume(volume, one_step, spec=spec, control=control)
  call check(one_step%status == intraday_not_converged, 'maxit nonconvergence status')
  call check(.not. one_step%converged, 'maxit nonconvergence flag')
  call check(one_step%iterations == 1, 'one-iteration result is safe')

  write(*, '(a)') 'test_fixed_acceleration: PASS'
end program test_fixed_acceleration
