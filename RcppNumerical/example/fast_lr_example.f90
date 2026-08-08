program fast_lr_example
  use rcppnumerical
  implicit none

  type(logistic_fit_t) :: fit
  real(dp) :: x(10,2), y(10)

  x(:,1) = 1.0_dp
  x(:,2) = [-2.0_dp, -1.5_dp, -1.0_dp, -0.5_dp, 0.0_dp, &
            0.5_dp, 1.0_dp, 1.5_dp, 2.0_dp, 2.5_dp]
  y = [0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
       1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp]

  call fast_lr(x, y, fit)
  write(*,'(a,*(1x,f12.7))') 'coefficients:', fit%coefficients
  write(*,'(a,f14.8)') 'log likelihood: ', fit%log_likelihood
  write(*,'(a,l1)') 'converged: ', fit%converged
end program fast_lr_example
