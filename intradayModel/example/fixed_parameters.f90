! SPDX-License-Identifier: Apache-2.0
program fixed_parameters
  use intraday_model
  implicit none
  type(volume_parameters) :: truth
  type(volume_model_spec) :: spec
  type(volume_model) :: model
  type(volume_fit_control) :: control
  real(dp), allocatable :: volume(:, :)
  integer :: i

  allocate(truth%phi(6))
  do i = 1, 6
    truth%phi(i) = 0.12_dp * sin(2.0_dp * acos(-1.0_dp) * real(i - 1, dp) / 6.0_dp)
  end do
  truth%x0 = [2.8_dp, 0.0_dp]
  truth%a_eta = 0.82_dp
  truth%a_mu = 0.35_dp
  truth%var_eta = 0.018_dp
  truth%var_mu = 0.009_dp
  truth%r = 0.006_dp

  call simulate_intraday_volume(truth, 60, volume, seed=909)
  call initialize_volume_spec(spec, 6)
  spec%fixed%a_mu = truth%a_mu
  spec%fixed%var_mu = truth%var_mu
  spec%is_fixed%a_mu = .true.
  spec%is_fixed%var_mu = .true.
  spec%initial%a_eta = 0.6_dp
  spec%has_initial%a_eta = .true.

  control%acceleration = .true.
  control%maxit = 200
  control%abstol = 3.0e-5_dp
  control%save_history = .false.
  call fit_volume(volume, model, spec=spec, control=control)

  write(*, '(a,l1)') 'converged: ', model%converged
  write(*, '(a,f10.5)') 'estimated a_eta: ', model%par%a_eta
  write(*, '(a,f10.5)') 'fixed a_mu:      ', model%par%a_mu
  write(*, '(a,f10.5)') 'fixed var_mu:    ', model%par%var_mu
end program fixed_parameters
