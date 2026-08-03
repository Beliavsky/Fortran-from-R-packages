! SPDX-License-Identifier: Apache-2.0
program basic_fit
  use intraday_model
  implicit none
  type(volume_parameters) :: truth
  type(volume_model) :: model
  type(volume_fit_control) :: control
  real(dp), allocatable :: volume(:, :)
  integer :: i

  allocate(truth%phi(10))
  do i = 1, 10
    truth%phi(i) = -0.20_dp * cos(2.0_dp * acos(-1.0_dp) * real(i - 1, dp) / 10.0_dp)
  end do
  truth%phi = truth%phi - sum(truth%phi) / 10.0_dp
  truth%x0 = [3.5_dp, 0.0_dp]
  truth%a_eta = 0.88_dp
  truth%a_mu = 0.40_dp
  truth%var_eta = 0.012_dp
  truth%var_mu = 0.008_dp
  truth%r = 0.004_dp

  call simulate_intraday_volume(truth, 70, volume, seed=1001)
  control%acceleration = .true.
  control%maxit = 250
  control%abstol = 2.0e-5_dp
  control%save_history = .false.
  call fit_volume(volume, model, control=control)

  write(*, '(a,l1)') 'converged: ', model%converged
  write(*, '(a,i0)') 'iterations: ', model%iterations
  write(*, '(a,2f12.6)') 'a_eta, a_mu: ', model%par%a_eta, model%par%a_mu
  write(*, '(a,3f12.6)') 'var_eta, var_mu, r: ', &
    model%par%var_eta, model%par%var_mu, model%par%r
end program basic_fit
