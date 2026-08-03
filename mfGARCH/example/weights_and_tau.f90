program weights_and_tau
  use mfgarch
  implicit none
  type(mfgarch_model) :: model
  real(dp), allocatable :: weights(:), tau(:)
  real(dp) :: covariate(8)
  integer :: period(16), status

  model%k = 4
  model%asymmetric = .false.
  model%gamma = 0.0_dp
  model%m = -0.1_dp
  model%theta = 0.5_dp
  model%w1 = 1.0_dp
  model%w2 = 3.0_dp
  covariate = [0.1_dp,0.2_dp,-0.1_dp,0.4_dp,0.3_dp,0.0_dp,-0.2_dp,0.1_dp]
  period = [1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8]

  call beta_weights(model%k, model%w1, model%w2, weights, status)
  call build_tau(model, period, covariate, tau, status)
  write(*,'(a)') 'MIDAS weights:'
  write(*,'(*(f10.6,1x))') weights
  write(*,'(a)') 'Expanded tau values (leading values are NaN because lags are unavailable):'
  write(*,'(*(es12.4,1x))') tau
end program weights_and_tau
