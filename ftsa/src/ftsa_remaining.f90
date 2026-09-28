! SPDX-License-Identifier: GPL-3.0-only
! Derived from R package ftsa 6.7 by Rob Hyndman and Han Lin Shang.
module ftsa_remaining
   use r_kinds, only : dp
   use r_linalg, only : solve_spd, symmetric_eigen, thin_svd
   use ftsa_models, only : ftsa_forecast_ftsm, ftsa_ftsm
   use ftsa_statistics, only : ftsa_invalid_input, ftsa_ok
   use ftsa_types, only : ftsm_fit
   implicit none
   private

   public :: ftsa_coda_bayes_nw, ftsa_coda_fpca, ftsa_farforecast, ftsa_horta_ziegelmann_fpca
   public :: ftsa_lqdt_fpca, ftsa_mfdm, ftsa_skew_t_fun, ftsa_t_stationary

contains



   subroutine ftsa_mfdm(mort_female, mort_male, mort_ave, percent_common, percent_specific, horizon, &
      female_forecast, male_forecast, ncomp_common, ncomp_female, ncomp_male, explained_first, status)
      real(dp), intent(in) :: mort_female(:, :) !! Grid-by-time female functional series used in the multilevel decomposition.
      real(dp), intent(in) :: mort_male(:, :) !! Grid-by-time male functional series with the same dimensions as `mort_female`.
      real(dp), intent(in) :: mort_ave(:, :) !! Grid-by-time common/average functional series with matching dimensions.
      real(dp), intent(in) :: percent_common !! Target cumulative variance fraction for common components in `(0,1]`.
      real(dp), intent(in) :: percent_specific !! Target cumulative variance fraction for sex-specific components in `(0,1]`.
      integer, intent(in) :: horizon !! Positive number of future functional observations to forecast.
      real(dp), allocatable, intent(out) :: female_forecast(:, :) !! Grid-by-horizon deterministic female point forecasts.
      real(dp), allocatable, intent(out) :: male_forecast(:, :) !! Grid-by-horizon deterministic male point forecasts.
      integer, intent(out) :: ncomp_common !! Number of retained common components.
      integer, intent(out) :: ncomp_female !! Number of retained female-specific components.
      integer, intent(out) :: ncomp_male !! Number of retained male-specific components.
      real(dp), intent(out) :: explained_first(3) !! First-component variance fractions for common, female, and male decompositions.
      integer, intent(out) :: status !! Zero on success; nonzero for dimension errors or failed SVDs.
      real(dp), allocatable :: ave_center(:, :), basis_common(:, :), basis_female(:, :), basis_male(:, :)
      real(dp), allocatable :: female_center(:, :), female_resid(:, :), female_score(:, :), female_score_fore(:, :)
      real(dp), allocatable :: male_center(:, :), male_resid(:, :), male_score(:, :), male_score_fore(:, :)
      real(dp), allocatable :: mean_ave(:), mean_female(:), mean_male(:), score_common(:, :), score_common_fore(:, :)
      integer :: g, local_status, n

      allocate(female_forecast(0, 0), male_forecast(0, 0))
      ncomp_common = 0
      ncomp_female = 0
      ncomp_male = 0
      explained_first = 0.0_dp
      g = size(mort_female, 1)
      n = size(mort_female, 2)
      if (g < 2 .or. n < 3 .or. horizon < 1 .or. size(mort_male, 1) /= g .or. &
         size(mort_male, 2) /= n .or. size(mort_ave, 1) /= g .or. size(mort_ave, 2) /= n .or. &
         percent_common <= 0.0_dp .or. percent_common > 1.0_dp .or. percent_specific <= 0.0_dp .or. &
         percent_specific > 1.0_dp) then
         status = ftsa_invalid_input
         return
      end if
      allocate(mean_ave(g), mean_female(g), mean_male(g), ave_center(g, n), female_center(g, n), male_center(g, n))
      mean_ave = sum(mort_ave, dim=2) / real(n, dp)
      mean_female = sum(mort_female, dim=2) / real(n, dp)
      mean_male = sum(mort_male, dim=2) / real(n, dp)
      ave_center = mort_ave - spread(mean_ave, 2, n)
      female_center = mort_female - spread(mean_female, 2, n)
      male_center = mort_male - spread(mean_male, 2, n)

      call pca_by_target(ave_center, percent_common, basis_common, score_common, ncomp_common, explained_first(1), local_status)
      if (local_status /= ftsa_ok) then
         status = local_status
         return
      end if
      allocate(female_resid(g, n), male_resid(g, n))
      female_resid = female_center - matmul(basis_common, transpose(score_common))
      male_resid = male_center - matmul(basis_common, transpose(score_common))
      call pca_by_target(female_resid, percent_specific, basis_female, female_score, ncomp_female, &
         explained_first(2), local_status)
      if (local_status /= ftsa_ok) then
         status = local_status
         return
      end if
      call pca_by_target(male_resid, percent_specific, basis_male, male_score, ncomp_male, explained_first(3), local_status)
      if (local_status /= ftsa_ok) then
         status = local_status
         return
      end if

      call drift_score_forecast(score_common, horizon, score_common_fore)
      call drift_score_forecast(female_score, horizon, female_score_fore)
      call drift_score_forecast(male_score, horizon, male_score_fore)
      deallocate(female_forecast, male_forecast)
      allocate(female_forecast(g, horizon), male_forecast(g, horizon))
      female_forecast = spread(mean_female, 2, horizon) + matmul(basis_common, transpose(score_common_fore)) + &
         matmul(basis_female, transpose(female_score_fore))
      male_forecast = spread(mean_male, 2, horizon) + matmul(basis_common, transpose(score_common_fore)) + &
         matmul(basis_male, transpose(male_score_fore))
      status = ftsa_ok
   end subroutine ftsa_mfdm

   subroutine pca_by_target(centered, target, basis, scores, ncomp, first_fraction, status)
      real(dp), intent(in) :: centered(:, :) !! Grid-by-time centered functional observations to decompose by SVD.
      real(dp), intent(in) :: target !! Target cumulative squared-singular-value fraction in `(0,1]`.
      real(dp), allocatable, intent(out) :: basis(:, :) !! Retained grid-by-component right singular vectors.
      real(dp), allocatable, intent(out) :: scores(:, :) !! Time-by-component principal-component scores.
      integer, intent(out) :: ncomp !! Component count whose cumulative fraction is nearest to `target`.
      real(dp), intent(out) :: first_fraction !! Fraction of total squared singular value in the first component.
      integer, intent(out) :: status !! Zero on success; nonzero for degenerate data or failed SVD.
      real(dp), allocatable :: singular(:), u(:, :), vt(:, :)
      real(dp) :: cumulative, distance, best_distance, total
      integer :: i, info, r

      allocate(basis(0, 0), scores(0, 0))
      ncomp = 0
      first_fraction = 0.0_dp
      call thin_svd(transpose(centered), u, singular, vt, info)
      if (info /= 0 .or. size(singular) < 1) then
         status = ftsa_invalid_input
         return
      end if
      total = sum(singular**2)
      if (total <= tiny(1.0_dp)) then
         status = ftsa_invalid_input
         return
      end if
      first_fraction = singular(1)**2 / total
      cumulative = 0.0_dp
      best_distance = huge(1.0_dp)
      r = size(singular)
      do i = 1, r
         cumulative = cumulative + singular(i)**2 / total
         distance = abs(cumulative - target)
         if (distance < best_distance) then
            best_distance = distance
            ncomp = i
         end if
      end do
      deallocate(basis, scores)
      allocate(basis(size(centered, 1), ncomp), scores(size(centered, 2), ncomp))
      basis = transpose(vt(1:ncomp, :))
      scores = matmul(transpose(centered), basis)
      status = ftsa_ok
   end subroutine pca_by_target

   pure subroutine drift_score_forecast(scores, horizon, forecasts)
      real(dp), intent(in) :: scores(:, :) !! Historical time-by-component scores to extrapolate by a random walk with drift.
      integer, intent(in) :: horizon !! Positive forecast horizon.
      real(dp), allocatable, intent(out) :: forecasts(:, :) !! Horizon-by-component deterministic score forecasts.
      real(dp) :: drift(size(scores, 2))
      integer :: h, n

      n = size(scores, 1)
      allocate(forecasts(horizon, size(scores, 2)))
      if (n > 1) then
         drift = (scores(n, :) - scores(1, :)) / real(n - 1, dp)
      else
         drift = 0.0_dp
      end if
      do h = 1, horizon
         forecasts(h, :) = scores(n, :) + real(h, dp) * drift
      end do
   end subroutine drift_score_forecast

   subroutine ftsa_skew_t_fun(data, gridpoints, density_forecast, forecast_parameters, &
      estimated_parameters, selected_lag, status, max_var_lag)
      real(dp), intent(in) :: data(:, :) !! Time-by-observation samples fitted independently by standardized skew-Student laws.
      real(dp), intent(in) :: gridpoints(:) !! Support locations where the one-step forecast density is evaluated.
      real(dp), allocatable, intent(out) :: density_forecast(:) !! One-step standardized skew-Student density forecast.
      real(dp), intent(out) :: forecast_parameters(4) !! Forecast mean, standard deviation, degrees of freedom, and skew parameter.
      real(dp), allocatable, intent(out) :: estimated_parameters(:, :) !! Four-by-time fitted skew-Student parameters.
      integer, intent(out) :: selected_lag !! VAR lag selected by the translated AIC criterion for transformed parameters.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid dimensions, failed fits, or VAR failure.
      integer, intent(in), optional :: max_var_lag !! Maximum parameter-VAR lag; default is the largest feasible lag up to 10.
      real(dp), allocatable :: coefficients(:, :), residuals(:, :), series(:, :), one_step(:, :)
      real(dp), allocatable :: transformed(:, :)
      real(dp) :: candidate, df_fore, lower_df, upper_df
      integer :: i, local_status, max_lag, n_time

      allocate(density_forecast(0), estimated_parameters(0, 0))
      forecast_parameters = 0.0_dp
      selected_lag = 0
      n_time = size(data, 1)
      if (n_time < 5 .or. size(data, 2) < 5 .or. size(gridpoints) < 2) then
         status = ftsa_invalid_input
         return
      end if
      deallocate(estimated_parameters)
      allocate(estimated_parameters(4, n_time))
      do i = 1, n_time
         call fit_sstd(data(i, :), estimated_parameters(:, i), local_status)
         if (local_status /= ftsa_ok) then
            status = local_status
            return
         end if
      end do

      allocate(transformed(4, n_time))
      transformed(1, :) = estimated_parameters(1, :)
      transformed(2, :) = log(max(estimated_parameters(2, :)**2, tiny(1.0_dp)))
      do i = 1, n_time
         candidate = student_t_cdf(-2.0_dp, estimated_parameters(3, i))
         transformed(3, i) = log(max(candidate, tiny(1.0_dp)))
      end do
      transformed(4, :) = estimated_parameters(4, :)
      allocate(series(n_time - 1, 4))
      series(:, 1) = transformed(1, 2:n_time)
      series(:, 2) = transformed(2, 2:n_time) - transformed(2, 1:n_time - 1)
      series(:, 3) = transformed(3, 2:n_time) - transformed(3, 1:n_time - 1)
      series(:, 4) = transformed(4, 2:n_time)

      max_lag = min(10, max(1, (size(series, 1) - 2) / 5))
      if (present(max_var_lag)) max_lag = min(max_lag, max_var_lag)
      if (max_lag < 1) then
         status = ftsa_invalid_input
         return
      end if
      call select_var_aic(series, max_lag, selected_lag, coefficients, residuals, local_status)
      if (local_status /= ftsa_ok) then
         status = local_status
         return
      end if
      call forecast_score_var(series, coefficients, selected_lag, 1, .true., one_step, local_status)
      if (local_status /= ftsa_ok) then
         status = local_status
         return
      end if

      forecast_parameters(1) = one_step(1, 1)
      forecast_parameters(2) = sqrt(exp(one_step(1, 2) + transformed(2, n_time)))
      candidate = exp(one_step(1, 3) + transformed(3, n_time))
      lower_df = max(2.000001_dp, minval(estimated_parameters(3, :)) - 1.0_dp)
      upper_df = max(lower_df + 0.25_dp, maxval(estimated_parameters(3, :)) + 1.0_dp)
      call invert_student_tail(candidate, lower_df, upper_df, df_fore, local_status)
      if (local_status /= ftsa_ok) df_fore = estimated_parameters(3, n_time)
      forecast_parameters(3) = df_fore
      forecast_parameters(4) = max(one_step(1, 4), 0.05_dp)

      deallocate(density_forecast)
      allocate(density_forecast(size(gridpoints)))
      do i = 1, size(gridpoints)
         density_forecast(i) = sstd_density(gridpoints(i), forecast_parameters(1), forecast_parameters(2), &
            forecast_parameters(3), forecast_parameters(4))
      end do
      if (any(density_forecast < 0.0_dp)) then
         status = ftsa_invalid_input
      else
         status = ftsa_ok
      end if
   end subroutine ftsa_skew_t_fun

   subroutine fit_sstd(x, parameters, status)
      real(dp), intent(in) :: x(:) !! Observations fitted by maximum likelihood under the standardized skew-Student law.
      real(dp), intent(out) :: parameters(4) !! Fitted mean, standard deviation, degrees of freedom, and positive skew parameter.
      integer, intent(out) :: status !! Zero when the derivative-free likelihood optimization returns finite parameters.
      integer, parameter :: ndim = 4
      real(dp) :: centroid(ndim), contracted(ndim), expanded(ndim), reflected(ndim)
      real(dp) :: f(ndim + 1), simplex(ndim, ndim + 1), step(ndim), theta(ndim)
      real(dp) :: fc, fe, fr, mean_x, scale_x
      integer :: iter, j

      if (size(x) < 5) then
         parameters = 0.0_dp
         status = ftsa_invalid_input
         return
      end if
      mean_x = sum(x) / real(size(x), dp)
      scale_x = sqrt(sum((x - mean_x)**2) / real(max(1, size(x) - 1), dp))
      if (scale_x <= sqrt(tiny(1.0_dp))) then
         parameters = 0.0_dp
         status = ftsa_invalid_input
         return
      end if
      theta = [mean_x, log(scale_x), log(6.0_dp), 0.0_dp]
      step = [max(0.20_dp * scale_x, 1.0e-3_dp), 0.18_dp, 0.25_dp, 0.18_dp]
      simplex(:, 1) = theta
      do j = 1, ndim
         simplex(:, j + 1) = theta
         simplex(j, j + 1) = simplex(j, j + 1) + step(j)
      end do
      do j = 1, ndim + 1
         f(j) = sstd_negloglik(x, simplex(:, j))
      end do

      do iter = 1, 350
         call sort_simplex(simplex, f)
         if (maxval(abs(f(2:ndim + 1) - f(1))) < 1.0e-9_dp .and. &
            maxval(abs(simplex(:, 2:ndim + 1) - spread(simplex(:, 1), 2, ndim))) < 1.0e-6_dp) exit
         centroid = sum(simplex(:, 1:ndim), dim=2) / real(ndim, dp)
         reflected = centroid + (centroid - simplex(:, ndim + 1))
         fr = sstd_negloglik(x, reflected)
         if (fr < f(1)) then
            expanded = centroid + 2.0_dp * (reflected - centroid)
            fe = sstd_negloglik(x, expanded)
            if (fe < fr) then
               simplex(:, ndim + 1) = expanded
               f(ndim + 1) = fe
            else
               simplex(:, ndim + 1) = reflected
               f(ndim + 1) = fr
            end if
         else if (fr < f(ndim)) then
            simplex(:, ndim + 1) = reflected
            f(ndim + 1) = fr
         else
            if (fr < f(ndim + 1)) then
               contracted = centroid + 0.5_dp * (reflected - centroid)
               fc = sstd_negloglik(x, contracted)
               if (fc <= fr) then
                  simplex(:, ndim + 1) = contracted
                  f(ndim + 1) = fc
                  cycle
               end if
            else
               contracted = centroid + 0.5_dp * (simplex(:, ndim + 1) - centroid)
               fc = sstd_negloglik(x, contracted)
               if (fc < f(ndim + 1)) then
                  simplex(:, ndim + 1) = contracted
                  f(ndim + 1) = fc
                  cycle
               end if
            end if
            do j = 2, ndim + 1
               simplex(:, j) = simplex(:, 1) + 0.5_dp * (simplex(:, j) - simplex(:, 1))
               f(j) = sstd_negloglik(x, simplex(:, j))
            end do
         end if
      end do
      call sort_simplex(simplex, f)
      parameters(1) = simplex(1, 1)
      parameters(2) = exp(simplex(2, 1))
      parameters(3) = 2.0_dp + exp(simplex(3, 1))
      parameters(4) = exp(simplex(4, 1))
      if (f(1) >= huge(1.0_dp) / 10.0_dp .or. parameters(2) <= 0.0_dp .or. &
         parameters(3) <= 2.0_dp .or. parameters(4) <= 0.0_dp) then
         status = ftsa_invalid_input
      else
         status = ftsa_ok
      end if
   end subroutine fit_sstd

   pure real(dp) function sstd_negloglik(x, theta) result(value)
      real(dp), intent(in) :: x(:) !! Observations entering the skew-Student likelihood.
      real(dp), intent(in) :: theta(4) !! Unconstrained mean, log-scale, log(df-2), and log-skew parameters.
      real(dp) :: log_density, mu, nu, sigma, xi
      integer :: i

      mu = theta(1)
      if (theta(2) > 30.0_dp .or. theta(2) < -30.0_dp .or. theta(3) > 12.0_dp .or. &
         theta(3) < -12.0_dp .or. abs(theta(4)) > 8.0_dp) then
         value = huge(1.0_dp)
         return
      end if
      sigma = exp(theta(2))
      nu = 2.0_dp + exp(theta(3))
      xi = exp(theta(4))
      value = 0.0_dp
      do i = 1, size(x)
         log_density = sstd_log_density(x(i), mu, sigma, nu, xi)
         if (log_density <= -0.5_dp * huge(1.0_dp)) then
            value = huge(1.0_dp)
            return
         end if
         value = value - log_density
      end do
   end function sstd_negloglik

   pure real(dp) function sstd_log_density(x, mu, sigma, nu, xi) result(log_pdf)
      real(dp), intent(in) :: x !! Evaluation point for the standardized Fernandez-Steel skew-Student log density.
      real(dp), intent(in) :: mu !! Location parameter of the standardized skew-Student distribution.
      real(dp), intent(in) :: sigma !! Positive scale parameter of the standardized skew-Student distribution.
      real(dp), intent(in) :: nu !! Student degrees of freedom, which must exceed two.
      real(dp), intent(in) :: xi !! Positive Fernandez-Steel skew parameter; one is symmetric.
      real(dp) :: beta_value, g, m1, skew_mu, skew_sigma, t_scale, xxi, z, zz

      log_pdf = -huge(1.0_dp)
      if (sigma <= 0.0_dp .or. nu <= 2.0_dp .or. xi <= 0.0_dp) return
      beta_value = exp(log_gamma(0.5_dp) + log_gamma(0.5_dp * nu) - log_gamma(0.5_dp * (nu + 1.0_dp)))
      m1 = 2.0_dp * sqrt(nu - 2.0_dp) / (nu - 1.0_dp) / beta_value
      skew_mu = m1 * (xi - 1.0_dp / xi)
      skew_sigma = sqrt(max(tiny(1.0_dp), (1.0_dp - m1 * m1) * (xi * xi + 1.0_dp / (xi * xi)) + &
         2.0_dp * m1 * m1 - 1.0_dp))
      z = (x - mu) / sigma
      zz = z * skew_sigma + skew_mu
      if (zz < 0.0_dp) then
         xxi = 1.0_dp / xi
      else if (zz > 0.0_dp) then
         xxi = xi
      else
         xxi = 1.0_dp
      end if
      g = 2.0_dp / (xi + 1.0_dp / xi)
      t_scale = sqrt(nu / (nu - 2.0_dp))
      log_pdf = log(g) + log(t_scale) + student_t_log_pdf((zz / xxi) * t_scale, nu) + &
         log(skew_sigma) - log(sigma)
   end function sstd_log_density

   pure real(dp) function sstd_density(x, mu, sigma, nu, xi) result(pdf)
      real(dp), intent(in) :: x !! Evaluation point for the standardized Fernandez-Steel skew-Student density.
      real(dp), intent(in) :: mu !! Location parameter of the standardized skew-Student distribution.
      real(dp), intent(in) :: sigma !! Positive scale parameter of the standardized skew-Student distribution.
      real(dp), intent(in) :: nu !! Student degrees of freedom, which must exceed two.
      real(dp), intent(in) :: xi !! Positive Fernandez-Steel skew parameter; one is symmetric.
      real(dp) :: log_pdf

      pdf = 0.0_dp
      log_pdf = sstd_log_density(x, mu, sigma, nu, xi)
      if (log_pdf < log(tiny(1.0_dp))) return
      pdf = exp(log_pdf)
   end function sstd_density

   pure real(dp) function student_t_log_pdf(x, nu) result(log_pdf)
      real(dp), intent(in) :: x !! Ordinary Student-t evaluation point.
      real(dp), intent(in) :: nu !! Positive Student degrees of freedom.

      if (nu <= 0.0_dp) then
         log_pdf = -huge(1.0_dp)
         return
      end if
      log_pdf = log_gamma(0.5_dp * (nu + 1.0_dp)) - log_gamma(0.5_dp * nu) - &
         0.5_dp * log(acos(-1.0_dp) * nu) - 0.5_dp * (nu + 1.0_dp) * log(1.0_dp + x * x / nu)
   end function student_t_log_pdf

   pure real(dp) function student_t_cdf(x, nu) result(probability)
      real(dp), intent(in) :: x !! Ordinary Student-t evaluation point.
      real(dp), intent(in) :: nu !! Positive Student degrees of freedom.
      real(dp) :: ibeta, z

      if (nu <= 0.0_dp) then
         probability = 0.0_dp
         return
      end if
      if (abs(x) <= epsilon(1.0_dp)) then
         probability = 0.5_dp
         return
      end if
      z = nu / (nu + x * x)
      ibeta = regularized_beta(z, 0.5_dp * nu, 0.5_dp)
      if (x < 0.0_dp) then
         probability = 0.5_dp * ibeta
      else
         probability = 1.0_dp - 0.5_dp * ibeta
      end if
   end function student_t_cdf

   pure real(dp) function regularized_beta(x, a, b) result(value)
      real(dp), intent(in) :: x !! Argument in `[0,1]` of the regularized incomplete beta function.
      real(dp), intent(in) :: a !! Positive first beta-shape parameter.
      real(dp), intent(in) :: b !! Positive second beta-shape parameter.
      real(dp) :: bt

      if (x <= 0.0_dp) then
         value = 0.0_dp
         return
      end if
      if (x >= 1.0_dp) then
         value = 1.0_dp
         return
      end if
      bt = exp(log_gamma(a + b) - log_gamma(a) - log_gamma(b) + a * log(x) + b * log(1.0_dp - x))
      if (x < (a + 1.0_dp) / (a + b + 2.0_dp)) then
         value = bt * beta_continued_fraction(a, b, x) / a
      else
         value = 1.0_dp - bt * beta_continued_fraction(b, a, 1.0_dp - x) / b
      end if
      value = min(1.0_dp, max(0.0_dp, value))
   end function regularized_beta

   pure real(dp) function beta_continued_fraction(a, b, x) result(value)
      real(dp), intent(in) :: a !! Positive first beta-shape parameter for the continued fraction.
      real(dp), intent(in) :: b !! Positive second beta-shape parameter for the continued fraction.
      real(dp), intent(in) :: x !! Argument in `(0,1)` for the beta continued fraction.
      real(dp), parameter :: eps_cf = 3.0e-14_dp
      real(dp), parameter :: fpmin = 1.0e-300_dp
      real(dp) :: aa, c, d, del, h, qab, qam, qap
      integer :: m, m2

      qab = a + b
      qap = a + 1.0_dp
      qam = a - 1.0_dp
      c = 1.0_dp
      d = 1.0_dp - qab * x / qap
      if (abs(d) < fpmin) d = fpmin
      d = 1.0_dp / d
      h = d
      do m = 1, 200
         m2 = 2 * m
         aa = real(m, dp) * (b - real(m, dp)) * x / &
            ((qam + real(m2, dp)) * (a + real(m2, dp)))
         d = 1.0_dp + aa * d
         if (abs(d) < fpmin) d = fpmin
         c = 1.0_dp + aa / c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp / d
         h = h * d * c
         aa = -(a + real(m, dp)) * (qab + real(m, dp)) * x / &
            ((a + real(m2, dp)) * (qap + real(m2, dp)))
         d = 1.0_dp + aa * d
         if (abs(d) < fpmin) d = fpmin
         c = 1.0_dp + aa / c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp / d
         del = d * c
         h = h * del
         if (abs(del - 1.0_dp) <= eps_cf) exit
      end do
      value = h
   end function beta_continued_fraction

   subroutine invert_student_tail(target, lower_df, upper_df, df_value, status)
      real(dp), intent(in) :: target !! Target probability equal to `P(T_df <= -2)`.
      real(dp), intent(in) :: lower_df !! Lower degrees-of-freedom bracket, strictly above two.
      real(dp), intent(in) :: upper_df !! Upper degrees-of-freedom bracket greater than `lower_df`.
      real(dp), intent(out) :: df_value !! Bisection solution for the degrees of freedom.
      integer, intent(out) :: status !! Zero if the target is bracketed and the bisection converges.
      real(dp) :: f_lower, f_mid, f_upper, left, mid, right
      integer :: iter

      left = max(2.000001_dp, lower_df)
      right = max(left + 1.0e-6_dp, upper_df)
      f_lower = student_t_cdf(-2.0_dp, left) - target
      f_upper = student_t_cdf(-2.0_dp, right) - target
      if (f_lower * f_upper > 0.0_dp) then
         df_value = right
         status = ftsa_invalid_input
         return
      end if
      do iter = 1, 100
         mid = 0.5_dp * (left + right)
         f_mid = student_t_cdf(-2.0_dp, mid) - target
         if (abs(f_mid) < 1.0e-11_dp .or. right - left < 1.0e-8_dp) exit
         if (f_lower * f_mid <= 0.0_dp) then
            right = mid
            f_upper = f_mid
         else
            left = mid
            f_lower = f_mid
         end if
      end do
      df_value = 0.5_dp * (left + right)
      status = ftsa_ok
   end subroutine invert_student_tail

   subroutine select_var_aic(series, max_lag, selected_lag, coefficients, residuals, status)
      real(dp), intent(in) :: series(:, :) !! Time-by-variable transformed parameter series.
      integer, intent(in) :: max_lag !! Largest VAR lag considered, at least one.
      integer, intent(out) :: selected_lag !! Lag minimizing the Gaussian VAR AIC approximation.
      real(dp), allocatable, intent(out) :: coefficients(:, :) !! Coefficients of the selected constant-including VAR.
      real(dp), allocatable, intent(out) :: residuals(:, :) !! Residuals from the selected VAR fit.
      integer, intent(out) :: status !! Zero if at least one lag yields a nonsingular finite criterion.
      real(dp), allocatable :: candidate_coef(:, :), candidate_resid(:, :), covariance(:, :), eigenvalues(:), eigenvectors(:, :)
      real(dp) :: aic, best, logdet
      integer :: d, info, local_status, nobs, npar, p

      allocate(coefficients(0, 0), residuals(0, 0))
      selected_lag = 0
      d = size(series, 2)
      best = huge(1.0_dp)
      if (size(series, 1) < 4 .or. d < 1 .or. max_lag < 1) then
         status = ftsa_invalid_input
         return
      end if
      do p = 1, max_lag
         call fit_score_var(series, p, .true., candidate_coef, candidate_resid, local_status)
         if (local_status /= ftsa_ok .or. size(candidate_resid, 1) < 2) cycle
         allocate(covariance(d, d))
         covariance = matmul(transpose(candidate_resid), candidate_resid) / real(size(candidate_resid, 1), dp)
         call symmetric_eigen(0.5_dp * (covariance + transpose(covariance)), eigenvalues, eigenvectors, info)
         deallocate(covariance)
         if (info /= 0) cycle
         logdet = sum(log(max(eigenvalues, 1.0e-12_dp * max(1.0_dp, maxval(abs(eigenvalues))))))
         nobs = size(candidate_resid, 1)
         npar = d * (d * p + 1)
         aic = logdet + 2.0_dp * real(npar, dp) / real(nobs, dp)
         if (aic < best) then
            best = aic
            selected_lag = p
            if (allocated(coefficients)) deallocate(coefficients)
            if (allocated(residuals)) deallocate(residuals)
            allocate(coefficients(size(candidate_coef, 1), size(candidate_coef, 2)))
            allocate(residuals(size(candidate_resid, 1), size(candidate_resid, 2)))
            coefficients = candidate_coef
            residuals = candidate_resid
         end if
      end do
      if (selected_lag < 1) then
         status = ftsa_invalid_input
      else
         status = ftsa_ok
      end if
   end subroutine select_var_aic

   subroutine sort_simplex(simplex, values)
      real(dp), intent(inout) :: simplex(:, :) !! Nelder-Mead vertices, reordered from best to worst objective.
      real(dp), intent(inout) :: values(:) !! Objective values corresponding to `simplex` columns.
      real(dp) :: tmp_col(size(simplex, 1)), tmp_value
      integer :: i, j, best

      do i = 1, size(values) - 1
         best = i
         do j = i + 1, size(values)
            if (values(j) < values(best)) best = j
         end do
         if (best /= i) then
            tmp_value = values(i)
            values(i) = values(best)
            values(best) = tmp_value
            tmp_col = simplex(:, i)
            simplex(:, i) = simplex(:, best)
            simplex(:, best) = tmp_col
         end if
      end do
   end subroutine sort_simplex



   subroutine ftsa_coda_bayes_nw(densities, forecast, selected_k, cv, status, normalize)
      real(dp), intent(in) :: densities(:, :) !! Positive grid-by-time discretized densities/compositions used by the
      !! CoDa Nadaraya-Watson forecast.
      real(dp), allocatable, intent(out) :: forecast(:) !! One-step density/composition forecast on the original grid.
      integer, intent(out) :: selected_k !! Cross-validation-selected number of nearest functional neighbours.
      real(dp), intent(out) :: cv !! Mean squared leave-one-out functional-response cross-validation criterion.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid/nonpositive data or degenerate kernel weights.
      logical, intent(in), optional :: normalize !! Apply closure and centered-log-ratio normalization; defaults to true.
      real(dp), allocatable :: alpha(:), dat(:, :), dat_center(:, :), f(:, :), h(:, :), predicted(:)
      real(dp) :: constant
      integer :: g, i, n, t
      logical :: do_normalize

      allocate(forecast(0))
      selected_k = 0
      cv = huge(1.0_dp)
      g = size(densities, 1)
      n = size(densities, 2)
      if (g < 2 .or. n < 6 .or. any(densities <= 0.0_dp)) then
         status = ftsa_invalid_input
         return
      end if
      do_normalize = .true.
      if (present(normalize)) do_normalize = normalize
      constant = sum(densities(:, 1))
      allocate(dat(n, g), dat_center(n, g), alpha(g), f(n, g), h(n, g))
      dat = transpose(densities)
      dat_center = dat
      do t = 1, n
         dat_center(t, :) = dat_center(t, :) / sum(dat_center(t, :))
      end do
      do i = 1, g
         alpha(i) = exp(sum(log(dat_center(:, i))) / real(n, dp))
      end do
      do t = 1, n
         f(t, :) = dat(t, :) / alpha
         if (do_normalize) f(t, :) = f(t, :) / sum(f(t, :))
         if (do_normalize) then
            h(t, :) = log(f(t, :) / exp(sum(log(f(t, :))) / real(g, dp)))
         else
            h(t, :) = log(f(t, :))
         end if
      end do
      call functional_knn_gcv(h(1:n - 1, :), h(2:n, :), h(n, :), predicted, selected_k, cv, status)
      if (status /= ftsa_ok) return
      deallocate(forecast)
      allocate(forecast(g))
      forecast = exp(predicted - maxval(predicted))
      forecast = forecast / sum(forecast)
      forecast = forecast * alpha
      forecast = constant * forecast / sum(forecast)
      status = ftsa_ok
   end subroutine ftsa_coda_bayes_nw

   subroutine ftsa_coda_fpca(densities, varprop, forecast, ncomp, status, normalize, drift)
      real(dp), intent(in) :: densities(:, :) !! Positive grid-by-time discretized densities or compositions.
      real(dp), intent(in) :: varprop !! Cumulative singular-value variance fraction in `(0,1]` used for component selection.
      real(dp), allocatable, intent(out) :: forecast(:) !! One-step inverse-CLR density/composition forecast on the input grid.
      integer, intent(out) :: ncomp !! Number of retained log-ratio principal components.
      integer, intent(out) :: status !! Zero on success; nonzero for nonpositive inputs, invalid dimensions, or failed SVD.
      logical, intent(in), optional :: normalize !! Apply closure and centered-log-ratio normalization; defaults to true.
      logical, intent(in), optional :: drift !! Add the average score increment to the random-walk score forecast; default false.
      real(dp), allocatable :: alpha(:), basis(:, :), centered(:, :), dat(:, :), f(:, :), h(:, :), singular(:), u(:, :), vt(:, :)
      real(dp), allocatable :: scores(:, :), score_forecast(:), weights(:)
      real(dp) :: constant, cumulative, total
      integer :: g, i, info, n, r, t
      logical :: do_drift, do_normalize

      allocate(forecast(0))
      n = size(densities, 2)
      g = size(densities, 1)
      ncomp = 0
      if (g < 2 .or. n < 2 .or. varprop <= 0.0_dp .or. varprop > 1.0_dp .or. any(densities <= 0.0_dp)) then
         status = ftsa_invalid_input
         return
      end if
      do_normalize = .true.
      if (present(normalize)) do_normalize = normalize
      do_drift = .false.
      if (present(drift)) do_drift = drift
      constant = sum(densities(:, 1))

      allocate(centered(n, g), dat(n, g), alpha(g), f(n, g), h(n, g))
      dat = transpose(densities)
      centered = dat
      do t = 1, n
         centered(t, :) = centered(t, :) / sum(centered(t, :))
      end do
      do i = 1, g
         alpha(i) = exp(sum(log(centered(:, i))) / real(n, dp))
      end do
      do t = 1, n
         f(t, :) = dat(t, :) / alpha
         if (do_normalize) f(t, :) = f(t, :) / sum(f(t, :))
         if (do_normalize) then
            h(t, :) = log(f(t, :) / exp(sum(log(f(t, :))) / real(g, dp)))
         else
            h(t, :) = log(f(t, :))
         end if
      end do

      call thin_svd(h, u, singular, vt, info)
      if (info /= 0 .or. size(singular) < 1) then
         status = merge(info, ftsa_invalid_input, info /= 0)
         return
      end if
      total = sum(singular**2)
      if (total <= tiny(1.0_dp)) then
         status = ftsa_invalid_input
         return
      end if
      cumulative = 0.0_dp
      ncomp = size(singular)
      do r = 1, size(singular)
         cumulative = cumulative + singular(r)**2
         if (cumulative / total >= varprop) then
            ncomp = r
            exit
         end if
      end do
      allocate(basis(g, ncomp), scores(n, ncomp), score_forecast(ncomp))
      basis = transpose(vt(1:ncomp, :))
      scores = matmul(h, basis)
      do r = 1, ncomp
         score_forecast(r) = scores(n, r)
         if (do_drift .and. n > 1) score_forecast(r) = score_forecast(r) + &
            (scores(n, r) - scores(1, r)) / real(n - 1, dp)
      end do
      allocate(weights(g))
      weights = exp(matmul(basis, score_forecast) - maxval(matmul(basis, score_forecast)))
      weights = weights / sum(weights)
      deallocate(forecast)
      allocate(forecast(g))
      forecast = weights * alpha
      forecast = constant * forecast / sum(forecast)
      status = ftsa_ok
   end subroutine ftsa_coda_fpca

   subroutine ftsa_horta_ziegelmann_fpca(densities, du, lag_p, dimension, forecast, basis, scores, &
      eigenvalues, selected_var_lag, status, max_var_lag)
      real(dp), intent(in) :: densities(:, :) !! Grid-by-time density functions used by the dynamic FPCA estimator.
      real(dp), intent(in) :: du !! Positive integration-grid spacing for L2 inner products and density normalization.
      integer, intent(in) :: lag_p !! Positive lag count `p` entering the upstream K-star dynamic covariance construction.
      integer, intent(in) :: dimension !! Fixed number of dynamic components to retain; bootstrap selection is omitted here.
      real(dp), allocatable, intent(out) :: forecast(:) !! One-step nonnegative normalized density forecast.
      real(dp), allocatable, intent(out) :: basis(:, :) !! Grid-by-component dynamic eigenfunctions `psihat`.
      real(dp), allocatable, intent(out) :: scores(:, :) !! Component-by-time dynamic scores `etahat`.
      real(dp), allocatable, intent(out) :: eigenvalues(:) !! Retained K-star eigenvalues in descending order.
      integer, intent(out) :: selected_var_lag !! FPE-selected lag used for the score VAR forecast.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid dimensions or failed eigensystem/VAR fitting.
      integer, intent(in), optional :: max_var_lag !! Maximum score-VAR lag considered; default is three.
      real(dp), allocatable :: bhalf(:, :), binvhalf(:, :), core(:, :), cmat(:, :), eval_b(:), evec_b(:, :)
      real(dp), allocatable :: eval_c(:), evec_c(:, :), gamma(:, :), k0(:, :), ksum(:, :), root(:, :), ydev(:, :)
      real(dp), allocatable :: coefficients(:, :), residuals(:, :), score_forecasts(:, :), score_time(:, :)
      real(dp), allocatable :: work_forecast(:), ybar(:)
      real(dp) :: best_fpe, fpe, norm_value, scale, tol
      integer :: d, i, info, lag, local_status, n, nvmax, q, rank_b

      allocate(forecast(0), basis(0, 0), scores(0, 0), eigenvalues(0))
      selected_var_lag = 0
      q = size(densities, 1)
      n = size(densities, 2)
      if (q < 2 .or. n < 4 .or. du <= 0.0_dp .or. lag_p < 1 .or. lag_p >= n - 1 .or. &
         dimension < 1 .or. dimension > n - lag_p .or. any(densities < 0.0_dp)) then
         status = ftsa_invalid_input
         return
      end if
      nvmax = 3
      if (present(max_var_lag)) nvmax = max_var_lag
      if (nvmax < 0) then
         status = ftsa_invalid_input
         return
      end if

      allocate(ybar(q), ydev(q, n), core(n, n))
      ybar = sum(densities, dim=2) / real(n, dp)
      ydev = densities - spread(ybar, 2, n)
      core = du * matmul(transpose(ydev), ydev)
      d = n - lag_p
      allocate(k0(d, d), ksum(d, d))
      k0 = core(1:d, 1:d)
      ksum = 0.0_dp
      do lag = 1, lag_p
         ksum = ksum + core(lag + 1:lag + d, lag + 1:lag + d)
      end do
      call symmetric_eigen(k0, eval_b, evec_b, info, descending=.true.)
      if (info /= 0) then
         status = info
         return
      end if
      tol = epsilon(1.0_dp) * real(max(1, d), dp) * max(1.0_dp, maxval(abs(eval_b))) * 100.0_dp
      rank_b = count(eval_b > tol)
      if (rank_b < dimension) then
         status = ftsa_invalid_input
         return
      end if
      allocate(bhalf(d, d), binvhalf(d, d))
      bhalf = 0.0_dp
      binvhalf = 0.0_dp
      do i = 1, d
         if (eval_b(i) > tol) then
            root = spread(evec_b(:, i), 2, d) * spread(evec_b(:, i), 1, d)
            bhalf = bhalf + sqrt(eval_b(i)) * root
            binvhalf = binvhalf + root / sqrt(eval_b(i))
         end if
      end do
      allocate(cmat(d, d))
      cmat = matmul(bhalf, matmul(ksum, bhalf)) / real(d * d, dp)
      cmat = 0.5_dp * (cmat + transpose(cmat))
      call symmetric_eigen(cmat, eval_c, evec_c, info, descending=.true.)
      if (info /= 0) then
         status = info
         return
      end if
      if (count(eval_c > tol) < dimension) then
         status = ftsa_invalid_input
         return
      end if
      allocate(gamma(d, dimension))
      gamma = matmul(binvhalf, evec_c(:, 1:dimension))
      deallocate(basis, scores, eigenvalues)
      allocate(basis(q, dimension), scores(dimension, n), eigenvalues(dimension))
      basis = matmul(ydev(:, 1:d), gamma)
      do i = 1, dimension
         norm_value = sqrt(du * sum(basis(:, i)**2))
         if (norm_value <= tiny(1.0_dp)) then
            status = ftsa_invalid_input
            return
         end if
         basis(:, i) = basis(:, i) / norm_value
      end do
      scores = du * matmul(transpose(basis), ydev)
      eigenvalues = max(eval_c(1:dimension), 0.0_dp)

      allocate(score_time(n, dimension))
      score_time = transpose(scores)
      best_fpe = huge(1.0_dp)
      selected_var_lag = -1
      do lag = 0, min(nvmax, n - 2)
         call fit_score_var(score_time, lag, .true., coefficients, residuals, local_status)
         if (local_status /= ftsa_ok) cycle
         call fpe_trace(residuals, lag, fpe, local_status)
         if (local_status /= ftsa_ok) cycle
         if (fpe < best_fpe) then
            best_fpe = fpe
            selected_var_lag = lag
         end if
      end do
      if (selected_var_lag < 0) then
         status = ftsa_invalid_input
         return
      end if
      call fit_score_var(score_time, selected_var_lag, .true., coefficients, residuals, local_status)
      if (local_status /= ftsa_ok) then
         status = local_status
         return
      end if
      call forecast_score_var(score_time, coefficients, selected_var_lag, 1, .true., score_forecasts, local_status)
      if (local_status /= ftsa_ok) then
         status = local_status
         return
      end if
      allocate(work_forecast(q))
      work_forecast = ybar + matmul(basis, score_forecasts(1, :))
      work_forecast = max(work_forecast, 0.0_dp)
      scale = sum(work_forecast) * du
      if (scale <= tiny(1.0_dp)) then
         status = ftsa_invalid_input
         return
      end if
      work_forecast = work_forecast / scale
      deallocate(forecast)
      allocate(forecast(q))
      forecast = work_forecast
      status = ftsa_ok
   end subroutine ftsa_horta_ziegelmann_fpca

   subroutine ftsa_farforecast(y, horizon, dmax, pmax, forecasts, selected_p, selected_d, status, include_intercept)
      real(dp), intent(in) :: y(:, :) !! Regular-grid functional observations with grid points in rows and times in columns.
      integer, intent(in) :: horizon !! Number of future curves to forecast; must be positive.
      integer, intent(in) :: dmax !! Maximum functional principal-component dimension considered by FPE selection.
      integer, intent(in) :: pmax !! Maximum autoregressive lag considered by FPE selection; may be zero.
      real(dp), allocatable, intent(out) :: forecasts(:, :) !! Grid-by-horizon point forecasts from the selected score VAR.
      integer, intent(out) :: selected_p !! VAR lag minimizing the translated FPE criterion.
      integer, intent(out) :: selected_d !! Functional principal-component count minimizing the translated FPE criterion.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid dimensions or failed linear algebra.
      logical, intent(in), optional :: include_intercept !! Include a constant in score VAR fits; defaults to true, matching
      !! the upstream `var_type="const"` behavior.
      type(ftsm_fit) :: fit_all, fit_selected
      real(dp), allocatable :: coefficients(:, :), residuals(:, :), score_forecasts(:, :), scores(:, :)
      real(dp) :: best_value, explained, fpe, total_variance, value
      integer :: d, local_status, n, p
      logical :: with_intercept

      allocate(forecasts(0, 0))
      selected_p = 0
      selected_d = 0
      n = size(y, 2)
      if (size(y, 1) < 1 .or. n < 3 .or. horizon < 1 .or. dmax < 1 .or. pmax < 0) then
         status = ftsa_invalid_input
         return
      end if
      if (dmax > min(size(y, 1), n)) then
         status = ftsa_invalid_input
         return
      end if
      with_intercept = .true.
      if (present(include_intercept)) with_intercept = include_intercept

      call ftsa_ftsm(y, dmax, fit_all, local_status)
      if (local_status /= ftsa_ok) then
         status = local_status
         return
      end if
      allocate(scores(n, dmax))
      scores = transpose(fit_all%scores)
      total_variance = sum((y - spread(sum(y, dim=2) / real(n, dp), 2, n))**2) / real(n - 1, dp)
      best_value = huge(1.0_dp)
      do d = 1, dmax
         explained = sum(scores(:, 1:d)**2) / real(n - 1, dp)
         do p = 0, pmax
            call fit_score_var(scores(:, 1:d), p, with_intercept, coefficients, residuals, local_status)
            if (local_status /= ftsa_ok) cycle
            call fpe_trace(residuals, p, fpe, local_status)
            if (local_status /= ftsa_ok) cycle
            value = fpe + total_variance - explained
            if (value < best_value) then
               best_value = value
               selected_p = p
               selected_d = d
            end if
         end do
      end do
      if (selected_d < 1) then
         status = ftsa_invalid_input
         return
      end if

      call ftsa_ftsm(y, selected_d, fit_selected, local_status)
      if (local_status /= ftsa_ok) then
         status = local_status
         return
      end if
      deallocate(scores)
      allocate(scores(n, selected_d))
      scores = transpose(fit_selected%scores)
      call fit_score_var(scores, selected_p, with_intercept, coefficients, residuals, local_status)
      if (local_status /= ftsa_ok) then
         status = local_status
         return
      end if
      call forecast_score_var(scores, coefficients, selected_p, horizon, with_intercept, score_forecasts, local_status)
      if (local_status /= ftsa_ok) then
         status = local_status
         return
      end if
      deallocate(forecasts)
      allocate(forecasts(size(y, 1), horizon))
      forecasts = spread(fit_selected%mean_curve, 2, horizon) + matmul(fit_selected%basis, transpose(score_forecasts))
      status = ftsa_ok
   end subroutine ftsa_farforecast


   subroutine ftsa_lqdt_fpca(densities, grid, varprop, forecast, reconstruction, l2_diff, total_mass, &
      ncomp, status, lqd_size, drift)
      real(dp), intent(in) :: densities(:, :) !! Nonnegative grid-by-time density functions to transform and forecast.
      real(dp), intent(in) :: grid(:) !! Strictly increasing density support grid matching the first dimension of `densities`.
      real(dp), intent(in) :: varprop !! Target fraction in `(0,1]` for the upstream squared-variance component rule.
      real(dp), allocatable, intent(out) :: forecast(:) !! One-step density forecast interpolated to the input support grid.
      real(dp), allocatable, intent(out) :: reconstruction(:, :) !! Grid-by-time densities reconstructed through the LQD transform.
      real(dp), allocatable, intent(out) :: l2_diff(:) !! Trapezoidal L2 reconstruction error for each historical density.
      real(dp), allocatable, intent(out) :: total_mass(:) !! Reconstructed density mass on the original support grid.
      integer, intent(out) :: ncomp !! Number of retained LQD functional principal components.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid densities/grids or failed FTSM fitting.
      integer, intent(in), optional :: lqd_size !! Number of equally spaced `[0,1]` LQD grid points; default equals input grid size.
      logical, intent(in), optional :: drift !! Use random-walk-with-drift LQD score forecasts when true; default false.
      type(ftsm_fit) :: fit
      real(dp), allocatable :: c(:), dens(:, :), lqd(:, :), lqd_fore(:), one(:, :), t(:)
      real(dp) :: c_fore, cumulative, total, t0
      integer :: g, i, local_status, m, max_order, n, r
      logical :: use_drift

      allocate(forecast(0), reconstruction(0, 0), l2_diff(0), total_mass(0))
      ncomp = 0
      g = size(densities, 1)
      n = size(densities, 2)
      m = g
      if (present(lqd_size)) m = lqd_size
      use_drift = .false.
      if (present(drift)) use_drift = drift
      if (g < 3 .or. n < 3 .or. size(grid) /= g .or. m < 5 .or. varprop <= 0.0_dp .or. varprop > 1.0_dp) then
         status = ftsa_invalid_input
         return
      end if
      if (any(grid(2:g) <= grid(1:g - 1)) .or. any(densities < 0.0_dp)) then
         status = ftsa_invalid_input
         return
      end if
      allocate(t(m), dens(g, n), lqd(m, n), c(n))
      do i = 1, m
         t(i) = real(i - 1, dp) / real(m - 1, dp)
      end do
      t0 = grid(minloc(abs(grid), dim=1))
      do i = 1, n
         dens(:, i) = densities(:, i)
         call normalize_density(grid, dens(:, i), local_status)
         if (local_status /= ftsa_ok) then
            status = local_status
            return
         end if
         call density_to_lqd(dens(:, i), grid, t, t0, lqd(:, i), c(i), local_status)
         if (local_status /= ftsa_ok) then
            status = local_status
            return
         end if
      end do

      deallocate(reconstruction, l2_diff, total_mass)
      allocate(reconstruction(g, n), l2_diff(n), total_mass(n))
      do i = 1, n
         call lqd_to_density(lqd(:, i), t, t0, c(i), grid, reconstruction(:, i), local_status)
         if (local_status /= ftsa_ok) then
            status = local_status
            return
         end if
         total_mass(i) = trapz(grid, reconstruction(:, i))
         l2_diff(i) = sqrt(max(0.0_dp, trapz(grid, (dens(:, i) - reconstruction(:, i))**2)))
      end do

      max_order = min(10, min(m - 2, n))
      call ftsa_ftsm(lqd(2:m - 1, :), max_order, fit, local_status)
      if (local_status /= ftsa_ok) then
         status = local_status
         return
      end if
      total = sum(fit%varprop**2)
      if (total <= tiny(1.0_dp)) then
         status = ftsa_invalid_input
         return
      end if
      cumulative = 0.0_dp
      ncomp = max_order
      do r = 1, max_order
         cumulative = cumulative + fit%varprop(r)**2
         if (cumulative / total >= varprop) then
            ncomp = r
            exit
         end if
      end do
      call ftsa_ftsm(lqd(2:m - 1, :), ncomp, fit, local_status)
      if (local_status /= ftsa_ok) then
         status = local_status
         return
      end if
      call ftsa_forecast_ftsm(fit, 1, one, local_status, drift=use_drift)
      if (local_status /= ftsa_ok) then
         status = local_status
         return
      end if
      allocate(lqd_fore(m))
      lqd_fore(1) = 8.0_dp
      lqd_fore(m) = 8.0_dp
      lqd_fore(2:m - 1) = one(:, 1)
      call forecast_probability(c, c_fore, use_drift)
      deallocate(forecast)
      allocate(forecast(g))
      call lqd_to_density(lqd_fore, t, t0, c_fore, grid, forecast, local_status, boundary_cut=.true.)
      if (local_status /= ftsa_ok) then
         status = local_status
         return
      end if
      status = ftsa_ok
   end subroutine ftsa_lqdt_fpca

   subroutine ftsa_t_stationary(sample, statistic, p_value, ncomp, status, l_basis, j_terms, mc_rep, &
      cumulative_var, bandwidth, pivotal, kernel_kind)
      real(dp), intent(in) :: sample(:, :) !! Functional sample with grid points in rows and successive curves in columns.
      real(dp), intent(out) :: statistic !! Upstream stationarity statistic for the requested pivotal or nonpivotal form.
      real(dp), intent(out) :: p_value !! Monte Carlo upper-tail p-value based on the translated limiting distribution.
      integer, intent(out) :: ncomp !! Number of eigencomponents used by the statistic or Monte Carlo approximation.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid dimensions or degenerate covariance structure.
      integer, intent(in), optional :: l_basis !! Number of Fourier basis functions used in the covariance projection; default 49.
      integer, intent(in), optional :: j_terms !! Number of Karhunen-Loeve terms in each Monte Carlo draw; default 500.
      integer, intent(in), optional :: mc_rep !! Number of Monte Carlo replications; default 1000.
      real(dp), intent(in), optional :: cumulative_var !! Pivotal eigenvalue fraction used to select components; default 0.90.
      real(dp), intent(in), optional :: bandwidth !! Long-run covariance lag bandwidth; default is square root of curve count.
      logical, intent(in), optional :: pivotal !! Use the pivotal score statistic when true; default false.
      integer, intent(in), optional :: kernel_kind !! Lag kernel: 1 for upstream Ker1, 2 for upstream Ker2; default 2.
      real(dp), allocatable :: centered(:, :), covmat(:, :), dmat(:, :), eval_d(:), evec_d(:, :)
      real(dp), allocatable :: eval_z(:), evec_z(:, :), fourier(:, :), gamma(:, :), partial(:), sim_values(:)
      real(dp), allocatable :: score(:, :)
      real(dp) :: bw, cumulative, cv, pi, total_eval, total_score, z
      integer :: i, info, j, jj, k, l, lag, lb, mc, n, positive, t, tf
      logical :: use_pivotal

      statistic = 0.0_dp
      p_value = 0.0_dp
      ncomp = 0
      tf = size(sample, 1)
      n = size(sample, 2)
      if (tf < 2 .or. n < 3) then
         status = ftsa_invalid_input
         return
      end if
      lb = 49
      if (present(l_basis)) lb = l_basis
      jj = 500
      if (present(j_terms)) jj = j_terms
      mc = 1000
      if (present(mc_rep)) mc = mc_rep
      cv = 0.90_dp
      if (present(cumulative_var)) cv = cumulative_var
      bw = sqrt(real(n, dp))
      if (present(bandwidth)) bw = bandwidth
      use_pivotal = .false.
      if (present(pivotal)) use_pivotal = pivotal
      k = 2
      if (present(kernel_kind)) k = kernel_kind
      if (lb < 1 .or. jj < 1 .or. mc < 1 .or. cv <= 0.0_dp .or. cv > 1.0_dp .or. bw <= 0.0_dp) then
         status = ftsa_invalid_input
         return
      end if
      if (k /= 1 .and. k /= 2) then
         status = ftsa_invalid_input
         return
      end if

      allocate(centered(tf, n), covmat(tf, tf))
      centered = sample - spread(sum(sample, dim=2) / real(n, dp), 2, n)
      covmat = matmul(centered, transpose(centered)) / real(n, dp)
      do lag = 1, n - 1
         allocate(gamma(tf, tf))
         gamma = matmul(centered(:, lag + 1:n), transpose(centered(:, 1:n - lag))) / real(n, dp)
         covmat = covmat + lag_kernel(real(lag, dp) / bw, k) * (gamma + transpose(gamma))
         deallocate(gamma)
      end do

      allocate(fourier(tf, lb))
      call fourier_basis(tf, lb, fourier)
      allocate(dmat(lb, lb))
      dmat = matmul(transpose(fourier), matmul(covmat, fourier)) / real(tf * tf, dp)
      dmat = 0.5_dp * (dmat + transpose(dmat))
      call symmetric_eigen(dmat, eval_d, evec_d, info, descending=.true.)
      if (info /= 0) then
         status = info
         return
      end if
      call symmetric_eigen(covmat, eval_z, evec_z, info, descending=.true.)
      if (info /= 0) then
         status = info
         return
      end if
      eval_d = max(eval_d, 0.0_dp)
      positive = count(eval_d > sqrt(epsilon(1.0_dp)) * max(1.0_dp, maxval(eval_d)))
      if (positive < 1) then
         status = ftsa_invalid_input
         return
      end if
      pi = acos(-1.0_dp)

      if (use_pivotal) then
         total_eval = sum(eval_d)
         cumulative = 0.0_dp
         ncomp = positive
         do i = 1, positive
            cumulative = cumulative + eval_d(i)
            if (cumulative / total_eval >= cv) then
               ncomp = i
               exit
            end if
         end do
         allocate(score(ncomp, n))
         score = matmul(transpose(evec_z(:, 1:ncomp)), sample) / sqrt(real(tf, dp))
         statistic = 0.0_dp
         do i = 1, ncomp
            total_score = sum(score(i, :))
            do t = 1, n
               statistic = statistic + ((sum(score(i, 1:t)) - real(t, dp) / real(n, dp) * total_score) / &
                  sqrt(real(n, dp)))**2 / real(n, dp) / max(eval_d(i), tiny(1.0_dp))
            end do
         end do
         allocate(sim_values(mc))
         do i = 1, mc
            sim_values(i) = 0.0_dp
            do l = 1, ncomp
               do j = 1, jj
                  call standard_normal_local(z)
                  sim_values(i) = sim_values(i) + (z / (pi * real(j, dp)))**2
               end do
            end do
         end do
      else
         ncomp = min(positive, 15)
         allocate(partial(tf), sim_values(mc))
         statistic = 0.0_dp
         do t = 1, n
            partial = (sum(sample(:, 1:t), dim=2) - real(t, dp) / real(n, dp) * sum(sample, dim=2)) / sqrt(real(n, dp))
            statistic = statistic + sum(partial**2) / real(n * tf, dp)
         end do
         do i = 1, mc
            sim_values(i) = 0.0_dp
            do l = 1, ncomp
               do j = 1, jj
                  call standard_normal_local(z)
                  sim_values(i) = sim_values(i) + eval_d(l) * (z / (pi * real(j, dp)))**2
               end do
            end do
         end do
      end if
      p_value = real(count(sim_values >= statistic), dp) / real(mc, dp)
      status = ftsa_ok
   end subroutine ftsa_t_stationary



   subroutine density_to_lqd(density, support, lqd_support, t0, lqd, c_value, status)
      real(dp), intent(in) :: density(:) !! Normalized nonnegative density values on `support`.
      real(dp), intent(in) :: support(:) !! Strictly increasing density support.
      real(dp), intent(in) :: lqd_support(:) !! Target LQD probability grid, normally spanning zero to one.
      real(dp), intent(in) :: t0 !! Density-domain reference point whose CDF value is retained.
      real(dp), intent(out) :: lqd(size(lqd_support)) !! Interpolated log-quantile-density function.
      real(dp), intent(out) :: c_value !! CDF probability associated with the support point nearest `t0`.
      integer, intent(out) :: status !! Zero on success; nonzero when no positive density support remains.
      real(dp), allocatable :: cdf(:), d(:), x(:)
      integer :: first, i, last, nearest, n

      n = size(density)
      first = 0
      last = 0
      do i = 1, n
         if (density(i) > 0.0_dp) then
            if (first == 0) first = i
            last = i
         end if
      end do
      if (first == 0 .or. last - first < 1) then
         status = ftsa_invalid_input
         return
      end if
      allocate(d(last - first + 1), x(last - first + 1), cdf(last - first + 1))
      d = density(first:last)
      x = support(first:last)
      d = d / trapz(x, d)
      call cumulative_trapz(x, d, cdf)
      nearest = minloc(abs(x - t0), dim=1)
      c_value = cdf(nearest)
      call linear_interp(cdf, -log(max(d, tiny(1.0_dp))), lqd_support, lqd, .true.)
      status = ftsa_ok
   end subroutine density_to_lqd

   subroutine lqd_to_density(lqd, lqd_support, t0, c_value, output_grid, density, status, boundary_cut)
      real(dp), intent(in) :: lqd(:) !! Log-quantile-density values on `lqd_support`.
      real(dp), intent(in) :: lqd_support(:) !! Increasing probability support for the LQD function.
      real(dp), intent(in) :: t0 !! Density-domain point satisfying `F(t0)=c_value` approximately.
      real(dp), intent(in) :: c_value !! Reference CDF probability used to position the reconstructed support.
      real(dp), intent(in) :: output_grid(:) !! Density-domain grid on which the reconstructed density is returned.
      real(dp), intent(out) :: density(size(output_grid)) !! Reconstructed nonnegative density values.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid cuts or degenerate normalization.
      logical, intent(in), optional :: boundary_cut !! Drop leading/trailing LQD values above seven as in upstream forecasting.
      real(dp), allocatable :: dtemp(:), exp_lqd(:), ls(:), lt(:), raw(:), cumulative(:)
      real(dp) :: mass, offset
      integer :: cut_left, cut_right, i, i0, m, max_check, nuse
      logical :: do_cut

      m = size(lqd)
      if (m < 3 .or. size(lqd_support) /= m) then
         status = ftsa_invalid_input
         return
      end if
      do_cut = .false.
      if (present(boundary_cut)) do_cut = boundary_cut
      cut_left = 0
      cut_right = 0
      if (do_cut) then
         max_check = min(15, m / 2)
         do i = 1, max_check
            if (lqd(i) > 7.0_dp) cut_left = i
         end do
         do i = 1, max_check
            if (lqd(m - i + 1) > 7.0_dp) cut_right = i
         end do
      end if
      nuse = m - cut_left - cut_right
      if (nuse < 3) then
         status = ftsa_invalid_input
         return
      end if
      allocate(ls(nuse), lt(nuse), exp_lqd(nuse), cumulative(nuse), dtemp(nuse), raw(nuse))
      ls = lqd(cut_left + 1:m - cut_right)
      lt = lqd_support(cut_left + 1:m - cut_right)
      exp_lqd = exp(min(ls, log(huge(1.0_dp)) - 2.0_dp))
      call cumulative_trapz(lt, exp_lqd, cumulative)
      i0 = minloc(abs(lt - c_value), dim=1)
      offset = cumulative(i0)
      dtemp = t0 + cumulative - offset
      raw = exp(-ls)
      call linear_interp(dtemp, raw, output_grid, density, .false.)
      density = max(density, 0.0_dp)
      mass = trapz(output_grid, density)
      if (mass <= tiny(1.0_dp)) then
         status = ftsa_invalid_input
         return
      end if
      density = density / mass * (lt(nuse) - lt(1))
      status = ftsa_ok
   end subroutine lqd_to_density

   pure subroutine normalize_density(grid, density, status)
      real(dp), intent(in) :: grid(:) !! Strictly increasing density support grid.
      real(dp), intent(inout) :: density(:) !! Density vector normalized in place to unit trapezoidal mass.
      integer, intent(out) :: status !! Zero when the density has positive finite mass.
      real(dp) :: mass

      mass = trapz(grid, density)
      if (mass <= tiny(1.0_dp)) then
         status = ftsa_invalid_input
         return
      end if
      density = density / mass
      status = ftsa_ok
   end subroutine normalize_density

   pure real(dp) function trapz(x, y) result(value)
      real(dp), intent(in) :: x(:) !! Strictly increasing integration abscissae.
      real(dp), intent(in) :: y(:) !! Integrand values matching `x`.
      integer :: i

      value = 0.0_dp
      do i = 1, min(size(x), size(y)) - 1
         value = value + 0.5_dp * (x(i + 1) - x(i)) * (y(i + 1) + y(i))
      end do
   end function trapz

   pure subroutine cumulative_trapz(x, y, cumulative)
      real(dp), intent(in) :: x(:) !! Strictly increasing integration abscissae.
      real(dp), intent(in) :: y(:) !! Integrand values matching `x`.
      real(dp), intent(out) :: cumulative(size(x)) !! Cumulative trapezoidal integral beginning at zero.
      integer :: i

      cumulative(1) = 0.0_dp
      do i = 2, size(x)
         cumulative(i) = cumulative(i - 1) + 0.5_dp * (x(i) - x(i - 1)) * (y(i) + y(i - 1))
      end do
   end subroutine cumulative_trapz

   pure subroutine linear_interp(x, y, xout, yout, constant_extrapolation)
      real(dp), intent(in) :: x(:) !! Strictly increasing interpolation abscissae.
      real(dp), intent(in) :: y(:) !! Interpolation ordinates matching `x`.
      real(dp), intent(in) :: xout(:) !! Abscissae at which interpolated values are requested.
      real(dp), intent(out) :: yout(size(xout)) !! Piecewise-linear interpolated values.
      logical, intent(in) :: constant_extrapolation !! Use endpoint values outside the range when true, otherwise zero.
      real(dp) :: weight
      integer :: i, j, n

      n = size(x)
      do i = 1, size(xout)
         if (xout(i) <= x(1)) then
            if (constant_extrapolation) then
               yout(i) = y(1)
            else
               yout(i) = 0.0_dp
            end if
         else if (xout(i) >= x(n)) then
            if (constant_extrapolation) then
               yout(i) = y(n)
            else
               yout(i) = 0.0_dp
            end if
         else
            j = 1
            do while (j < n - 1 .and. x(j + 1) < xout(i))
               j = j + 1
            end do
            weight = (xout(i) - x(j)) / (x(j + 1) - x(j))
            yout(i) = (1.0_dp - weight) * y(j) + weight * y(j + 1)
         end if
      end do
   end subroutine linear_interp

   pure subroutine forecast_probability(probability, forecast, drift)
      real(dp), intent(in) :: probability(:) !! Historical probabilities in `[0,1]` transformed on the logit scale.
      real(dp), intent(out) :: forecast !! One-step probability forecast transformed back from the logit scale.
      logical, intent(in) :: drift !! Add the average historical logit increment when true.
      real(dp) :: eta, first_eta, last_eta
      integer :: n

      n = size(probability)
      last_eta = log(max(tiny(1.0_dp), min(1.0_dp - epsilon(1.0_dp), probability(n))) / &
         max(tiny(1.0_dp), 1.0_dp - min(1.0_dp - epsilon(1.0_dp), probability(n))))
      eta = last_eta
      if (drift .and. n > 1) then
         first_eta = log(max(tiny(1.0_dp), min(1.0_dp - epsilon(1.0_dp), probability(1))) / &
            max(tiny(1.0_dp), 1.0_dp - min(1.0_dp - epsilon(1.0_dp), probability(1))))
         eta = eta + (last_eta - first_eta) / real(n - 1, dp)
      end if
      forecast = 1.0_dp / (1.0_dp + exp(-eta))
   end subroutine forecast_probability

   subroutine functional_knn_gcv(curves, responses, pred, predicted, selected_k, cv, status)
      real(dp), intent(in) :: curves(:, :) !! Training predictor curves stored row by row.
      real(dp), intent(in) :: responses(:, :) !! Functional responses paired row by row with `curves`.
      real(dp), intent(in) :: pred(:) !! Predictor curve for which the functional response is forecast.
      real(dp), allocatable, intent(out) :: predicted(:) !! Kernel-weighted one-step functional response forecast.
      integer, intent(out) :: selected_k !! Nearest-neighbour count minimizing leave-one-out squared prediction error.
      real(dp), intent(out) :: cv !! Mean squared cross-validation criterion at `selected_k`.
      integer, intent(out) :: status !! Zero on success; nonzero for incompatible dimensions or zero kernel denominators.
      real(dp), allocatable :: dist(:), estimate(:), order_dist(:), weights(:)
      real(dp) :: bandwidth, criterion, denom, best
      integer, allocatable :: order(:)
      integer :: g, i, j, k, kmax, kmin, n

      allocate(predicted(0))
      selected_k = 0
      cv = huge(1.0_dp)
      n = size(curves, 1)
      g = size(curves, 2)
      if (n < 3 .or. g < 2 .or. size(responses, 1) /= n .or. size(responses, 2) /= g .or. size(pred) /= g) then
         status = ftsa_invalid_input
         return
      end if
      kmin = min(5, max(1, n / 3))
      kmax = max(kmin, min(n - 2, n / 2))
      if (kmax < 1) then
         status = ftsa_invalid_input
         return
      end if
      allocate(dist(n), order(n), order_dist(n), weights(n), estimate(g))
      best = huge(1.0_dp)
      do k = kmin, kmax
         criterion = 0.0_dp
         do i = 1, n
            do j = 1, n
               dist(j) = sqrt(sum((curves(j, :) - curves(i, :))**2) / real(max(1, g - 1), dp))
            end do
            call sort_with_index(dist, order_dist, order)
            bandwidth = 0.5_dp * (order_dist(k + 1) + order_dist(k + 2))
            if (bandwidth <= tiny(1.0_dp)) then
               estimate = responses(order(k + 1), :)
            else
               weights = 0.0_dp
               do j = 2, k + 1
                  weights(order(j)) = max(0.0_dp, 1.0_dp - (dist(order(j)) / bandwidth)**2)
               end do
               denom = sum(weights)
               if (denom <= tiny(1.0_dp)) then
                  estimate = responses(order(k + 1), :)
               else
                  estimate = matmul(weights, responses) / denom
               end if
            end if
            criterion = criterion + sum((estimate - responses(i, :))**2)
         end do
         if (criterion < best) then
            best = criterion
            selected_k = k
            cv = criterion / real(n * g, dp)
         end if
      end do
      if (selected_k < 1) then
         status = ftsa_invalid_input
         return
      end if

      do j = 1, n
         dist(j) = sqrt(sum((curves(j, :) - pred)**2) / real(max(1, g - 1), dp))
      end do
      call sort_with_index(dist, order_dist, order)
      bandwidth = 0.5_dp * (order_dist(selected_k) + order_dist(selected_k + 1))
      if (bandwidth <= tiny(1.0_dp)) then
         estimate = responses(order(1), :)
      else
         weights = 0.0_dp
         do j = 1, selected_k
            weights(order(j)) = max(0.0_dp, 1.0_dp - (dist(order(j)) / bandwidth)**2)
         end do
         denom = sum(weights)
         if (denom <= tiny(1.0_dp)) then
            estimate = responses(order(1), :)
         else
            estimate = matmul(weights, responses) / denom
         end if
      end if
      deallocate(predicted)
      allocate(predicted(g))
      predicted = estimate
      status = ftsa_ok
   end subroutine functional_knn_gcv

   pure subroutine sort_with_index(values, sorted, order)
      real(dp), intent(in) :: values(:) !! Values to sort in ascending order.
      real(dp), intent(out) :: sorted(size(values)) !! Ascending copy of `values`.
      integer, intent(out) :: order(size(values)) !! Original indices corresponding to `sorted`.
      real(dp) :: temp
      integer :: i, j, itemp

      sorted = values
      do i = 1, size(values)
         order(i) = i
      end do
      do i = 2, size(values)
         temp = sorted(i)
         itemp = order(i)
         j = i - 1
         do while (j >= 1)
            if (sorted(j) <= temp) exit
            sorted(j + 1) = sorted(j)
            order(j + 1) = order(j)
            j = j - 1
         end do
         sorted(j + 1) = temp
         order(j + 1) = itemp
      end do
   end subroutine sort_with_index

   subroutine fit_score_var(scores, lag_order, include_intercept, coefficients, residuals, status)
      real(dp), intent(in) :: scores(:, :) !! Time-by-component score matrix.
      integer, intent(in) :: lag_order !! VAR lag order; zero fits a constant or zero mean only.
      logical, intent(in) :: include_intercept !! Include a constant regression column when true.
      real(dp), allocatable, intent(out) :: coefficients(:, :) !! Regression coefficients, with intercept first when requested.
      real(dp), allocatable, intent(out) :: residuals(:, :) !! Effective-sample residual matrix.
      integer, intent(out) :: status !! Zero on success; nonzero when the regression is underdetermined or fails.
      real(dp), allocatable :: design(:, :), fitted(:, :), gram(:, :), response(:, :), rhs(:, :)
      real(dp) :: diagonal_scale, ridge
      integer :: col, d, info, lag, n, neff, npar

      n = size(scores, 1)
      d = size(scores, 2)
      neff = n - lag_order
      npar = d * lag_order
      if (include_intercept) npar = npar + 1
      allocate(coefficients(0, 0), residuals(0, 0))
      if (d < 1 .or. lag_order < 0 .or. neff < 2 .or. npar < 1 .or. neff < npar) then
         status = ftsa_invalid_input
         return
      end if
      allocate(design(neff, npar), response(neff, d), fitted(neff, d))
      col = 0
      if (include_intercept) then
         col = 1
         design(:, col) = 1.0_dp
      end if
      do lag = 1, lag_order
         design(:, col + 1:col + d) = scores(lag_order + 1 - lag:n - lag, :)
         col = col + d
      end do
      response = scores(lag_order + 1:n, :)
      if (lag_order == 0 .and. .not. include_intercept) then
         status = ftsa_invalid_input
         return
      end if
      allocate(gram(npar, npar), rhs(npar, d))
      gram = matmul(transpose(design), design)
      rhs = matmul(transpose(design), response)
      diagonal_scale = sum(abs([(gram(lag, lag), lag=1, npar)])) / real(npar, dp)
      ridge = sqrt(epsilon(1.0_dp)) * max(1.0_dp, diagonal_scale)
      do lag = 1, npar
         gram(lag, lag) = gram(lag, lag) + ridge
      end do
      deallocate(coefficients)
      allocate(coefficients(npar, d))
      call solve_spd(gram, rhs, coefficients, info)
      if (info /= 0) then
         status = info
         return
      end if
      fitted = matmul(design, coefficients)
      deallocate(residuals)
      allocate(residuals(neff, d))
      residuals = response - fitted
      status = ftsa_ok
   end subroutine fit_score_var

   subroutine forecast_score_var(scores, coefficients, lag_order, horizon, include_intercept, forecasts, status)
      real(dp), intent(in) :: scores(:, :) !! Historical time-by-component score matrix used for the fitted VAR.
      real(dp), intent(in) :: coefficients(:, :) !! Coefficients returned by `fit_score_var`.
      integer, intent(in) :: lag_order !! Fitted VAR lag order.
      integer, intent(in) :: horizon !! Number of recursive score forecasts required.
      logical, intent(in) :: include_intercept !! Whether the first coefficient row is an intercept.
      real(dp), allocatable, intent(out) :: forecasts(:, :) !! Horizon-by-component recursive score forecasts.
      integer, intent(out) :: status !! Zero on success; nonzero for incompatible dimensions.
      real(dp), allocatable :: history(:, :), next_history(:, :), row(:)
      integer :: col, d, h, lag, n, npar

      n = size(scores, 1)
      d = size(scores, 2)
      npar = d * lag_order
      if (include_intercept) npar = npar + 1
      allocate(forecasts(0, 0))
      if (horizon < 1 .or. d < 1 .or. n < max(1, lag_order) .or. size(coefficients, 1) /= npar .or. &
         size(coefficients, 2) /= d) then
         status = ftsa_invalid_input
         return
      end if
      deallocate(forecasts)
      allocate(history(n, d), forecasts(horizon, d), row(npar))
      history = scores
      do h = 1, horizon
         col = 0
         if (include_intercept) then
            col = 1
            row(col) = 1.0_dp
         end if
         do lag = 1, lag_order
            row(col + 1:col + d) = history(size(history, 1) + 1 - lag, :)
            col = col + d
         end do
         forecasts(h, :) = matmul(row, coefficients)
         allocate(next_history(size(history, 1) + 1, d))
         next_history(1:size(history, 1), :) = history
         next_history(size(history, 1) + 1, :) = forecasts(h, :)
         call move_alloc(next_history, history)
      end do
      status = ftsa_ok
   end subroutine forecast_score_var

   pure subroutine fpe_trace(residuals, lag_order, value, status)
      real(dp), intent(in) :: residuals(:, :) !! VAR or AR residual matrix with observations in rows.
      integer, intent(in) :: lag_order !! Autoregressive lag used to form the residuals.
      real(dp), intent(out) :: value !! Trace final-prediction-error criterion used by upstream `FPE.trace`.
      integer, intent(out) :: status !! Zero when the FPE denominator and residual dimensions are valid.
      real(dp) :: trace_cov
      real(dp) :: means(size(residuals, 2))
      integer :: d, n

      n = size(residuals, 1)
      d = size(residuals, 2)
      value = huge(1.0_dp)
      if (n < 2 .or. d < 1 .or. n <= lag_order * d) then
         status = ftsa_invalid_input
         return
      end if
      means = sum(residuals, dim=1) / real(n, dp)
      trace_cov = sum((residuals - spread(means, 1, n))**2) / real(n - 1, dp)
      value = real(lag_order * d + n, dp) / real(n - lag_order * d, dp) * trace_cov
      status = ftsa_ok
   end subroutine fpe_trace

   pure real(dp) function lag_kernel(x, kernel_kind) result(value)
      real(dp), intent(in) :: x !! Nonnegative lag divided by the bandwidth.
      integer, intent(in) :: kernel_kind !! One selects `min(1,max(1.1-|x|,0))`; two selects `min(1,max(2-2|x|,0))`.

      if (kernel_kind == 1) then
         value = min(1.0_dp, max(1.1_dp - abs(x), 0.0_dp))
      else
         value = min(1.0_dp, max(2.0_dp - 2.0_dp * abs(x), 0.0_dp))
      end if
   end function lag_kernel

   pure subroutine fourier_basis(n_grid, n_basis, basis)
      integer, intent(in) :: n_grid !! Number of equally spaced evaluation points on `(0,1]`.
      integer, intent(in) :: n_basis !! Number of orthonormal Fourier basis columns requested.
      real(dp), intent(out) :: basis(n_grid, n_basis) !! Analytic Fourier basis evaluated at `i/n_grid`.
      real(dp) :: pi, x
      integer :: i, j, harmonic

      pi = acos(-1.0_dp)
      basis(:, 1) = 1.0_dp
      do j = 2, n_basis
         harmonic = (j - 2) / 2 + 1
         do i = 1, n_grid
            x = real(i, dp) / real(n_grid, dp)
            if (mod(j, 2) == 0) then
               basis(i, j) = sqrt(2.0_dp) * sin(2.0_dp * pi * real(harmonic, dp) * x)
            else
               basis(i, j) = sqrt(2.0_dp) * cos(2.0_dp * pi * real(harmonic, dp) * x)
            end if
         end do
      end do
   end subroutine fourier_basis

   subroutine standard_normal_local(z)
      real(dp), intent(out) :: z !! One standard-normal random variate generated by Box-Muller transformation.
      real(dp) :: pi, u1, u2

      pi = acos(-1.0_dp)
      call random_number(u1)
      call random_number(u2)
      u1 = max(u1, tiny(1.0_dp))
      z = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * pi * u2)
   end subroutine standard_normal_local

end module ftsa_remaining
