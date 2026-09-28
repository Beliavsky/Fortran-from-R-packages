program test_remaining
   use ftsa, only : dp, ftsa_coda_bayes_nw, ftsa_coda_fpca, ftsa_farforecast, ftsa_horta_ziegelmann_fpca
   use ftsa, only : ftsa_lqdt_fpca, ftsa_ok, ftsa_t_stationary
   implicit none

   real(dp), allocatable :: basis(:, :), coda_forecast(:), eigenvalues(:), forecasts(:, :)
   real(dp), allocatable :: horta_forecast(:), horta_scores(:, :)
   real(dp), allocatable :: lqd_diff(:), lqd_forecast(:), lqd_mass(:), lqd_reconstruction(:, :)
   real(dp) :: cv, p_value, statistic
   real(dp) :: densities(5, 12), density_grid(5), sample(5, 20), y(3, 12)
   real(dp) :: score(12), pi
   integer :: i, j, k_selected, ncomp, psel, dsel, status, var_lag

   pi = acos(-1.0_dp)
   score(1) = -1.0_dp
   do j = 2, 12
      score(j) = 0.8_dp * score(j - 1) + 0.15_dp
   end do
   do j = 1, 12
      y(1, j) = 2.0_dp + score(j)
      y(2, j) = 1.0_dp + 2.0_dp * score(j)
      y(3, j) = -0.5_dp - score(j)
   end do
   call ftsa_farforecast(y, 3, 1, 1, forecasts, psel, dsel, status)
   if (status /= ftsa_ok) error stop "farforecast status"
   if (dsel /= 1) error stop "farforecast dimension"
   if (psel < 0 .or. psel > 1) error stop "farforecast lag"
   if (any(shape(forecasts) /= [3, 3])) error stop "farforecast shape"
   if (maxval(abs(forecasts)) > 100.0_dp) error stop "farforecast finite scale"

   do j = 1, 12
      densities(:, j) = [0.12_dp, 0.18_dp, 0.30_dp, 0.24_dp, 0.16_dp] + &
         0.02_dp * sin(0.4_dp * real(j, dp)) * [-1.0_dp, 0.5_dp, 1.0_dp, 0.5_dp, -1.0_dp]
      densities(:, j) = densities(:, j) / sum(densities(:, j))
   end do
   call ftsa_coda_fpca(densities, 0.90_dp, coda_forecast, ncomp, status, drift=.true.)
   if (status /= ftsa_ok) error stop "CoDa_FPCA status"
   if (ncomp < 1) error stop "CoDa_FPCA components"
   if (any(coda_forecast <= 0.0_dp)) error stop "CoDa_FPCA positivity"
   if (abs(sum(coda_forecast) - 1.0_dp) > 1.0e-10_dp) error stop "CoDa_FPCA closure"

   call ftsa_coda_bayes_nw(densities, coda_forecast, k_selected, cv, status)
   if (status /= ftsa_ok) error stop "CoDa_BayesNW status"
   if (k_selected < 1) error stop "CoDa_BayesNW k"
   if (cv < 0.0_dp) error stop "CoDa_BayesNW cv"
   if (any(coda_forecast <= 0.0_dp)) error stop "CoDa_BayesNW positivity"
   if (abs(sum(coda_forecast) - 1.0_dp) > 1.0e-10_dp) error stop "CoDa_BayesNW closure"

   call ftsa_horta_ziegelmann_fpca(densities, 0.25_dp, 2, 1, horta_forecast, basis, horta_scores, &
      eigenvalues, var_lag, status, max_var_lag=2)
   if (status /= ftsa_ok) error stop "Horta_Ziegelmann_FPCA status"
   if (size(basis, 2) /= 1 .or. size(horta_scores, 1) /= 1) error stop "Horta dimensions"
   if (var_lag < 0 .or. var_lag > 2) error stop "Horta VAR lag"
   if (any(horta_forecast < 0.0_dp)) error stop "Horta positivity"
   if (abs(sum(horta_forecast) * 0.25_dp - 1.0_dp) > 1.0e-10_dp) error stop "Horta normalization"

   density_grid = [-2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp]
   call ftsa_lqdt_fpca(densities, density_grid, 0.85_dp, lqd_forecast, lqd_reconstruction, &
      lqd_diff, lqd_mass, ncomp, status, lqd_size=11, drift=.true.)
   if (status /= ftsa_ok) error stop "LQDT_FPCA status"
   if (ncomp < 1) error stop "LQDT_FPCA components"
   if (size(lqd_reconstruction, 2) /= 12) error stop "LQDT reconstruction shape"
   if (any(lqd_diff < 0.0_dp) .or. any(lqd_mass <= 0.0_dp)) error stop "LQDT diagnostics"
   if (any(lqd_forecast < 0.0_dp)) error stop "LQDT forecast positivity"
   if (sum(lqd_forecast) <= 0.0_dp) error stop "LQDT forecast mass"

   do j = 1, 20
      do i = 1, 5
         sample(i, j) = sin(2.0_dp * pi * real(i, dp) / 5.0_dp) * &
            (1.0_dp + 0.15_dp * sin(0.7_dp * real(j, dp))) + &
            0.05_dp * cos(0.4_dp * real(i * j, dp))
      end do
   end do
   call set_seed(314159)
   call ftsa_t_stationary(sample, statistic, p_value, ncomp, status, l_basis=5, j_terms=20, mc_rep=80)
   if (status /= ftsa_ok) error stop "T_stationary status"
   if (statistic < 0.0_dp) error stop "T_stationary statistic"
   if (p_value < 0.0_dp .or. p_value > 1.0_dp) error stop "T_stationary p value"
   if (ncomp < 1 .or. ncomp > 5) error stop "T_stationary components"

   call set_seed(271828)
   call ftsa_t_stationary(sample, statistic, p_value, ncomp, status, l_basis=5, j_terms=20, mc_rep=80, pivotal=.true.)
   if (status /= ftsa_ok) error stop "T_stationary pivotal status"
   if (statistic < 0.0_dp) error stop "T_stationary pivotal statistic"
   if (p_value < 0.0_dp .or. p_value > 1.0_dp) error stop "T_stationary pivotal p value"

   print '(a)', 'test_remaining: PASS'

contains

   subroutine set_seed(seed_value)
      integer, intent(in) :: seed_value !! Deterministic scalar expanded across the compiler RNG seed vector.
      integer, allocatable :: seed(:)
      integer :: k, n

      call random_seed(size=n)
      allocate(seed(n))
      do k = 1, n
         seed(k) = modulo(seed_value + 104729 * k, huge(1) - 1)
      end do
      call random_seed(put=seed)
   end subroutine set_seed

end program test_remaining
