! SPDX-License-Identifier: GPL-3.0-only
program demo_mass
  use mass
  implicit none
  integer, parameter :: n = 50
  real(dp) :: design(n,2), predictors(n,3), response(n), contaminated(n)
  integer :: i
  type(regression_result) :: ordinary, robust
  type(model_selection_result) :: selection
  real(dp) :: bandwidth
  integer :: status

  do i = 1, n
    design(i,1) = 1.0_dp
    design(i,2) = -1.0_dp + 2.0_dp * real(i-1,dp) / real(n-1,dp)
    predictors(i,1) = design(i,2)
    predictors(i,2) = sin(0.37_dp * real(i,dp))
    predictors(i,3) = cos(0.23_dp * real(i,dp))
    response(i) = 1.0_dp + 2.0_dp * design(i,2) + 0.03_dp * sin(real(i,dp))
  end do
  contaminated = response
  contaminated(n) = contaminated(n) + 40.0_dp

  call linear_model_fit(design, contaminated, ordinary)
  call rlm_fit(design, contaminated, robust, psi='bisquare')
  call step_aic_linear(predictors, response, selection)
  call ucv_bandwidth(response, bandwidth, status)

  write(*,'(a,2f12.6)') 'OLS coefficients:    ', ordinary%coefficients
  write(*,'(a,2f12.6)') 'Robust coefficients: ', robust%coefficients
  write(*,'(a,3l3)') 'AIC-selected columns:', selection%selected
  write(*,'(a,f12.6)') 'UCV bandwidth:       ', bandwidth
end program demo_mass
