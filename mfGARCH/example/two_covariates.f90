program two_covariates
  use mfgarch
  implicit none
  type(mfgarch_model) :: model
  real(dp), allocatable :: tau(:)
  real(dp) :: macro_one(6), macro_two(4)
  integer :: period_one(12), period_two(12), status

  model%k = 2
  model%k_two = 1
  model%has_second = .true.
  model%asymmetric = .true.
  model%m = 0.0_dp
  model%theta = 0.3_dp
  model%w1 = 1.0_dp
  model%w2 = 3.0_dp
  model%theta_two = -0.15_dp
  model%w1_two = 1.0_dp
  model%w2_two = 1.0_dp
  macro_one = [0.2_dp,0.1_dp,-0.1_dp,0.3_dp,0.4_dp,0.2_dp]
  macro_two = [1.0_dp,0.8_dp,1.2_dp,1.1_dp]
  period_one = [1,1,2,2,3,3,4,4,5,5,6,6]
  period_two = [1,1,1,2,2,2,3,3,3,4,4,4]

  call build_tau(model, period_one, macro_one, tau, status, period_two, macro_two)
  if (status /= mfgarch_success) error stop 'tau construction failed'
  write(*,'(a)') 'Two-covariate long-run component:'
  write(*,'(*(es12.4,1x))') tau
end program two_covariates
