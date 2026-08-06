! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of the RPEGLMEN computational code.

program periodogram_glm_demo
  use rpeglmen, only : dp, enet_options, fit_result, glmnet_exp
  implicit none

  integer, parameter :: n = 120, degree = 5
  real(dp) :: frequency(n), spectrum(n), design(n, degree + 1)
  real(dp) :: estimated_standard_error
  type(enet_options) :: options
  type(fit_result) :: fit
  integer :: i, j

  do i = 1, n
    frequency(i) = 0.5_dp * real(i, dp) / real(n + 1, dp)
    spectrum(i) = exp(-2.2_dp + 1.1_dp * frequency(i) - 0.8_dp * frequency(i)**2) &
      * (1.0_dp + 0.08_dp * sin(0.8_dp * real(i, dp)))
    design(i, 1) = 1.0_dp
    do j = 1, degree
      design(i, j + 1) = frequency(i)**j
    end do
  end do

  options%num_lambda = 20
  options%k_fold = 5
  options%k_fold_iter = 2
  options%max_iter = 1000
  options%cv_metric = 'nll'
  call glmnet_exp(design, spectrum, fit, options)

  if (.not. allocated(fit%coefficients)) error stop 'fit failed'
  estimated_standard_error = sqrt(exp(fit%coefficients(1)) / real(n, dp))

  print '(a,es13.5)', 'selected lambda: ', fit%selected_lambda
  print '(a,es13.5)', 'intercept:       ', fit%coefficients(1)
  print '(a,es13.5)', 'estimated SE:    ', estimated_standard_error
end program periodogram_glm_demo
