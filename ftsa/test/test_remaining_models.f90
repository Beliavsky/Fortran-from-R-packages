program test_remaining_models
   use ftsa, only : dp, ftsa_gaev_forecast, ftsa_mfdm, ftsa_ok, ftsa_skew_t_fun
   implicit none

   real(dp) :: average(5, 12), female(5, 12), first_fraction(3), male(5, 12)
   real(dp) :: gaev_data(10, 8), gaev_q(2), gaev_shape
   real(dp) :: grid(61), skew_data(12, 31), skew_forecast_parameters(4)
   real(dp), allocatable :: female_forecast(:, :), gaev_density(:, :), gaev_location(:), gaev_scale(:)
   real(dp), allocatable :: male_forecast(:, :), skew_density(:), skew_parameters(:, :)
   real(dp) :: common_score, female_score, male_score, z
   integer :: i, j, k_loc, k_scale, lag, n_common, n_female, n_male, status, t

   do t = 1, 12
      common_score = 0.15_dp * real(t, dp)
      female_score = 0.08_dp * sin(0.5_dp * real(t, dp))
      male_score = -0.06_dp * cos(0.4_dp * real(t, dp))
      do i = 1, 5
         average(i, t) = 1.0_dp + 0.2_dp * real(i, dp) + common_score * real(i - 3, dp)
         female(i, t) = average(i, t) + 0.3_dp + female_score * real(i - 1, dp)
         male(i, t) = average(i, t) - 0.2_dp + male_score * real(5 - i, dp)
      end do
   end do
   call ftsa_mfdm(female, male, average, 0.95_dp, 0.95_dp, 3, female_forecast, male_forecast, &
      n_common, n_female, n_male, first_fraction, status)
   if (status /= ftsa_ok) error stop "MFDM status"
   if (any(shape(female_forecast) /= [5, 3]) .or. any(shape(male_forecast) /= [5, 3])) &
      error stop "MFDM forecast shape"
   if (n_common < 1 .or. n_female < 1 .or. n_male < 1) error stop "MFDM components"
   if (any(first_fraction <= 0.0_dp) .or. any(first_fraction > 1.0_dp)) error stop "MFDM fractions"

   do t = 1, 12
      do j = 1, 31
         z = real(j - 16, dp) / 5.5_dp
         skew_data(t, j) = 0.04_dp * real(t, dp) + z + 0.18_dp * z * z + &
            0.03_dp * sin(0.7_dp * real(t * j, dp))
      end do
   end do
   do j = 1, 61
      grid(j) = -4.0_dp + 0.15_dp * real(j - 1, dp)
   end do
   call ftsa_skew_t_fun(skew_data, grid, skew_density, skew_forecast_parameters, skew_parameters, lag, status, &
      max_var_lag=1)
   if (status /= ftsa_ok) error stop "skew_t_fun status"
   if (size(skew_density) /= 61 .or. any(skew_density < 0.0_dp)) error stop "skew_t_fun density"
   if (any(shape(skew_parameters) /= [4, 12])) error stop "skew_t_fun parameter shape"
   if (skew_forecast_parameters(2) <= 0.0_dp .or. skew_forecast_parameters(3) <= 2.0_dp .or. &
      skew_forecast_parameters(4) <= 0.0_dp) error stop "skew_t_fun parameters"
   if (lag /= 1) error stop "skew_t_fun VAR lag"

   do t = 1, 10
      do i = 1, 8
         gaev_data(t, i) = 5.0_dp + 0.25_dp * real(i, dp) + 0.04_dp * real(t, dp) + &
            0.18_dp * sin(0.7_dp * real(i + t, dp)) + 0.05_dp * cos(0.3_dp * real(i * t, dp))
      end do
   end do
   gaev_q = [6.0_dp, 7.0_dp]
   call ftsa_gaev_forecast(gaev_data, gaev_q, 3, 3, gaev_density, gaev_location, gaev_scale, gaev_shape, &
      k_loc, k_scale, status, max_var_lag=0)
   if (status /= ftsa_ok) error stop "GAEVforecast status"
   if (k_loc /= 3 .or. k_scale /= 3) error stop "GAEVforecast selected basis"
   if (any(shape(gaev_density) /= [2, 8])) error stop "GAEVforecast density shape"
   if (size(gaev_location) /= 8 .or. size(gaev_scale) /= 8) error stop "GAEVforecast parameter shape"
   if (any(gaev_scale <= 0.0_dp) .or. any(gaev_density < 0.0_dp)) error stop "GAEVforecast positivity"
   if (gaev_shape < -0.45_dp .or. gaev_shape > 0.45_dp) error stop "GAEVforecast shape range"

   print '(a)', 'test_remaining_models: PASS'
end program test_remaining_models
