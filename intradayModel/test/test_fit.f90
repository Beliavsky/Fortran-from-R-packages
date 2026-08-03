! SPDX-License-Identifier: Apache-2.0
program test_fit
  use intraday_model
  use test_support, only : check, make_parameters
  implicit none
  type(volume_parameters) :: truth
  type(volume_model) :: model
  type(volume_fit_control) :: control
  real(dp), allocatable :: volume(:, :), states(:, :)

  truth = make_parameters(8)
  call simulate_intraday_volume(truth, 80, volume, states, seed=123)
  control%acceleration = .false.
  control%maxit = 600
  control%abstol = 2.0e-4_dp
  control%save_history = .true.
  call fit_volume(volume, model, control=control)

  call check(model%status == intraday_ok, 'ordinary EM status')
  call check(model%converged, 'ordinary EM convergence')
  call check(model%iterations > 5 .and. model%iterations <= control%maxit, 'iteration count')
  call check(allocated(model%history), 'parameter history allocated')
  call check(size(model%history) == model%iterations + 1, 'parameter history length')
  call check(abs(model%par%a_eta - truth%a_eta) < 0.12_dp, 'daily persistence estimate')
  call check(abs(model%par%a_mu - truth%a_mu) < 0.18_dp, 'intraday persistence estimate')
  call check(abs(model%par%var_eta - truth%var_eta) < 0.012_dp, 'daily variance estimate')
  call check(abs(model%par%var_mu - truth%var_mu) < 0.012_dp, 'dynamic variance estimate')
  call check(sqrt(sum((model%par%phi - truth%phi)**2) / real(size(truth%phi), dp)) < 0.04_dp, &
             'seasonal profile estimate')

  write(*, '(a)') 'test_fit: PASS'
end program test_fit
