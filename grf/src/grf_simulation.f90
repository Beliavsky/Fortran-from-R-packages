module grf_simulation
   use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
   use grf_kinds, only : dp
   use grf_rng, only : grf_rng_state
   implicit none
   private

   public :: generate_causal_data
   public :: generate_causal_survival_data

contains

   subroutine generate_causal_data(n, p, x, y, w, tau, m, propensity, info, seed, dgp, &
                                   sigma_m, sigma_tau, sigma_noise)
      integer, intent(in) :: n !! Number of simulated observations; must be positive.
      integer, intent(in) :: p !! Number of predictor columns; simple DGP requires at least three.
      real(dp), allocatable, intent(out) :: x(:,:) !! Simulated predictor matrix with n rows and p columns.
      real(dp), allocatable, intent(out) :: y(:) !! Simulated observed outcome m(X)+(W-e(X))*tau(X)+Gaussian noise.
      integer, allocatable, intent(out) :: w(:) !! Binary treatment assignment drawn from the propensity score.
      real(dp), allocatable, intent(out) :: tau(:) !! True conditional treatment effect for each simulated row.
      real(dp), allocatable, intent(out) :: m(:) !! Conditional outcome mean used by the selected benchmark DGP.
      real(dp), allocatable, intent(out) :: propensity(:) !! True treatment propensity e(X) for each simulated row.
      integer, intent(out) :: info !! Zero on success; negative for invalid dimensions or unsupported DGP selector.
      integer, intent(in), optional :: seed !! Deterministic pseudo-random seed; default 42.
      character(len=*), intent(in), optional :: dgp !! Benchmark selector; default simple. All twelve upstream causal designs are supported.
      real(dp), intent(in), optional :: sigma_m !! Target sample standard deviation of the baseline response; default one.
      real(dp), intent(in), optional :: sigma_tau !! Target sample standard deviation of nonconstant treatment effects; default 0.1.
      real(dp), intent(in), optional :: sigma_noise !! Conditional noise standard deviation; default one.
      type(grf_rng_state) :: rng
      character(len=:), allocatable :: name
      real(dp), allocatable :: beta(:)
      integer :: local_seed
      integer :: i
      integer :: j
      real(dp) :: baseline_scale
      real(dp) :: effect_scale
      real(dp) :: noise_scale
      real(dp) :: zeta1
      real(dp) :: zeta2
      real(dp) :: mu0
      real(dp) :: mu1
      real(dp) :: pi
      real(dp) :: innovation_scale

      info = 0
      name = 'simple'
      if (present(dgp)) name = trim(dgp)
      if (n < 1 .or. p < 1) then
         info = -1
         return
      end if
      select case (name)
      case ('simple', 'nw3')
         if (p < 3) then
            info = -2
            return
         end if
      case ('aw1', 'aw2', 'aw3', 'aw3reverse', 'ai1', 'kunzel')
         if (p < 2) then
            info = -2
            return
         end if
      case ('ai2')
         if (p < 6) then
            info = -2
            return
         end if
      case ('nw1', 'nw2', 'nw4')
         if (p < 5) then
            info = -2
            return
         end if
      case default
         info = -3
         return
      end select
      baseline_scale = 1.0_dp
      if (present(sigma_m)) baseline_scale = sigma_m
      effect_scale = 0.1_dp
      if (present(sigma_tau)) effect_scale = sigma_tau
      noise_scale = 1.0_dp
      if (present(sigma_noise)) noise_scale = sigma_noise
      if (baseline_scale < 0.0_dp .or. effect_scale < 0.0_dp .or. noise_scale < 0.0_dp) then
         info = -4
         return
      end if
      local_seed = 42
      if (present(seed)) local_seed = seed
      call rng%seed(local_seed)
      allocate(x(n,p), y(n), w(n), tau(n), m(n), propensity(n))
      allocate(beta(p))
      beta = 0.0_dp
      select case (name)
      case ('aw1', 'aw2', 'aw3', 'aw3reverse', 'nw1')
         do j = 1, p
            do i = 1, n
               x(i,j) = rng%uniform()
            end do
         end do
      case ('kunzel')
         innovation_scale = sqrt(0.75_dp)
         do i = 1, n
            x(i,1) = rng%normal()
            do j = 2, p
               x(i,j) = 0.5_dp * x(i,j - 1) + innovation_scale * rng%normal()
            end do
         end do
      case default
         do j = 1, p
            do i = 1, n
               x(i,j) = rng%normal()
            end do
         end do
      end select
      pi = acos(-1.0_dp)
      if (name == 'kunzel') then
         do j = 1, p
            beta(j) = 10.0_dp * rng%uniform() - 5.0_dp
         end do
      end if
      do i = 1, n
         select case (name)
         case ('simple')
            tau(i) = max(x(i,1), 0.0_dp)
            if (x(i,1) > 0.0_dp) then
               propensity(i) = 0.6_dp
            else
               propensity(i) = 0.4_dp
            end if
            m(i) = x(i,2) + min(x(i,3), 0.0_dp) + propensity(i) * tau(i)
         case ('aw1')
            tau(i) = 0.0_dp
            propensity(i) = 0.25_dp * (1.0_dp + beta_density_2_4(x(i,1)))
            m(i) = 2.0_dp * x(i,1) - 1.0_dp
         case ('aw2')
            zeta1 = 1.0_dp + logistic(20.0_dp * (x(i,1) - 1.0_dp / 3.0_dp))
            zeta2 = 1.0_dp + logistic(20.0_dp * (x(i,2) - 1.0_dp / 3.0_dp))
            tau(i) = zeta1 * zeta2
            propensity(i) = 0.5_dp
            m(i) = propensity(i) * tau(i)
         case ('aw3')
            zeta1 = 1.0_dp + logistic(20.0_dp * (x(i,1) - 1.0_dp / 3.0_dp))
            zeta2 = 1.0_dp + logistic(20.0_dp * (x(i,2) - 1.0_dp / 3.0_dp))
            tau(i) = zeta1 * zeta2
            propensity(i) = 0.25_dp * (1.0_dp + beta_density_2_4(x(i,1)))
            m(i) = 2.0_dp * x(i,1) - 1.0_dp + propensity(i) * tau(i)
         case ('aw3reverse')
            zeta1 = 1.0_dp + logistic(-20.0_dp * (x(i,1) - 1.0_dp / 3.0_dp))
            zeta2 = 1.0_dp + logistic(-20.0_dp * (x(i,2) - 1.0_dp / 3.0_dp))
            tau(i) = zeta1 * zeta2
            propensity(i) = 0.25_dp * (1.0_dp + beta_density_2_4(x(i,1)))
            m(i) = 2.0_dp * x(i,1) - 1.0_dp + propensity(i) * tau(i)
         case ('ai1')
            tau(i) = 0.25_dp * x(i,1)
            propensity(i) = 0.5_dp
            m(i) = 0.5_dp * x(i,1) + x(i,2) + propensity(i) * tau(i)
         case ('ai2')
            tau(i) = 0.5_dp * (max(x(i,1), 0.0_dp) + max(x(i,2), 0.0_dp))
            propensity(i) = 0.5_dp
            m(i) = 0.5_dp * x(i,1) + 0.5_dp * x(i,2) + sum(x(i,3:6)) + propensity(i) * tau(i)
         case ('kunzel')
            tau(i) = merge(8.0_dp, 0.0_dp, x(i,2) > 0.1_dp)
            propensity(i) = 0.01_dp
         case ('nw1')
            tau(i) = 0.5_dp * (x(i,1) + x(i,2))
            propensity(i) = max(0.1_dp, min(sin(pi * x(i,1) * x(i,2)), 0.9_dp))
            m(i) = sin(pi * x(i,1) * x(i,2)) + 2.0_dp * (x(i,3) - 0.5_dp)**2 + &
                   x(i,4) + 0.5_dp * x(i,5) + propensity(i) * tau(i)
         case ('nw2')
            tau(i) = x(i,1) + log(1.0_dp + exp(x(i,2)))
            propensity(i) = 0.5_dp
            m(i) = max(0.0_dp, x(i,1) + x(i,2), x(i,3)) + &
                   max(0.0_dp, x(i,4) + x(i,5)) + propensity(i) * tau(i)
         case ('nw3')
            tau(i) = 1.0_dp
            propensity(i) = 1.0_dp / (1.0_dp + exp(x(i,2) + x(i,3)))
            m(i) = 2.0_dp * log(1.0_dp + exp(x(i,1) + x(i,2) + x(i,3))) + propensity(i) * tau(i)
         case ('nw4')
            tau(i) = max(x(i,1) + x(i,2) + x(i,3), 0.0_dp) - max(x(i,4) + x(i,5), 0.0_dp)
            propensity(i) = 1.0_dp / (1.0_dp + exp(-x(i,1)) + exp(-x(i,2)))
            m(i) = 0.5_dp * (max(x(i,1) + x(i,2) + x(i,3), 0.0_dp) + &
                   max(x(i,4) + x(i,5), 0.0_dp)) + propensity(i) * tau(i)
         end select
         w(i) = rng%bernoulli(propensity(i))
         if (name == 'kunzel') then
            mu0 = dot_product(x(i,:), beta) + merge(5.0_dp, 0.0_dp, x(i,1) > 0.5_dp) + rng%normal()
            mu1 = mu0 + tau(i) + rng%normal()
            m(i) = real(w(i),dp) * mu1 + real(1 - w(i),dp) * mu0 - &
                   (real(w(i),dp) - propensity(i)) * tau(i)
         end if
      end do
      call scale_by_sample_sd(m, baseline_scale)
      call scale_by_sample_sd(tau, effect_scale)
      do i = 1, n
         y(i) = m(i) + (real(w(i),dp) - propensity(i)) * tau(i) + noise_scale * rng%normal()
      end do
   end subroutine generate_causal_data

   pure elemental real(dp) function beta_density_2_4(x) result(value)
      real(dp), intent(in) :: x !! Point at which the beta(2,4) density is evaluated.

      if (x < 0.0_dp .or. x > 1.0_dp) then
         value = 0.0_dp
      else
         value = 20.0_dp * x * (1.0_dp - x)**3
      end if
   end function beta_density_2_4

   pure elemental real(dp) function logistic(x) result(value)
      real(dp), intent(in) :: x !! Real-valued logistic argument.

      if (x >= 0.0_dp) then
         value = 1.0_dp / (1.0_dp + exp(-x))
      else
         value = exp(x) / (1.0_dp + exp(x))
      end if
   end function logistic

   pure subroutine scale_by_sample_sd(values, target_sd)
      real(dp), intent(inout) :: values(:) !! Values rescaled in place when their sample standard deviation is positive.
      real(dp), intent(in) :: target_sd !! Requested sample standard deviation after multiplicative rescaling.
      real(dp) :: average
      real(dp) :: standard_deviation

      if (size(values) < 2) return
      average = sum(values) / real(size(values), dp)
      standard_deviation = sqrt(sum((values - average)**2) / real(size(values) - 1, dp))
      if (standard_deviation > 0.0_dp) values = values * target_sd / standard_deviation
   end subroutine scale_by_sample_sd

   subroutine generate_causal_survival_data(n, p, x, y, w, event, cate, cate_probability, info, &
                                            seed, y_max, y0, dgp, rho, n_mc, x_input, cate_sign)
      integer, intent(in) :: n !! Number of simulated observations; must be positive.
      integer, intent(in) :: p !! Number of generated predictor columns; ignored when x_input is supplied.
      real(dp), allocatable, intent(out) :: x(:,:) !! Uniform-marginal predictor matrix used by the selected benchmark.
      real(dp), allocatable, intent(out) :: y(:) !! Observed minimum of failure and censoring times.
      integer, allocatable, intent(out) :: w(:) !! Randomized binary treatment with probability one half.
      integer, allocatable, intent(out) :: event(:) !! Failure indicator equal to one when failure precedes censoring.
      real(dp), allocatable, intent(out) :: cate(:) !! Analytic restricted-mean treatment effect at y_max for each row.
      real(dp), allocatable, intent(out) :: cate_probability(:) !! Analytic treatment effect on survival probability at y0.
      integer, intent(out) :: info !! Zero on success; negative values indicate invalid dimensions, controls, or selector.
      integer, intent(in), optional :: seed !! Deterministic pseudo-random seed; default 42.
      real(dp), intent(in), optional :: y_max !! Administrative maximum failure time; the default depends on dgp.
      real(dp), intent(in), optional :: y0 !! Survival-probability evaluation time; the default depends on dgp.
      character(len=*), intent(in), optional :: dgp !! Benchmark selector: simple1 or type1 through type5; default simple1.
      real(dp), intent(in), optional :: rho !! AR(1) latent-normal predictor correlation; default zero.
      integer, intent(in), optional :: n_mc !! Deterministic Monte Carlo draws for target effects; default 10000.
      real(dp), intent(in), optional :: x_input(:,:) !! Optional caller-supplied predictor matrix replacing generated covariates.
      real(dp), allocatable, intent(out), optional :: cate_sign(:) !! Sign of each target effect; type4 uses NaN below its active threshold.
      type(grf_rng_state) :: rng
      character(len=:), allocatable :: name
      real(dp), allocatable :: monte_carlo_draw(:)
      real(dp), allocatable :: sign_values(:)
      real(dp) :: max_time
      real(dp) :: probability_time
      real(dp) :: failure
      real(dp) :: censor
      real(dp) :: scale
      real(dp) :: correlation
      real(dp) :: innovation_scale
      real(dp) :: propensity
      real(dp) :: indicator
      real(dp) :: denominator
      real(dp) :: failure0
      real(dp) :: failure1
      real(dp) :: lambda0
      real(dp) :: lambda1
      real(dp) :: cate_sum
      real(dp) :: probability_sum
      integer :: local_seed
      integer :: monte_carlo_count
      integer :: n_obs
      integer :: p_obs
      integer :: i
      integer :: j
      integer :: k

      info = 0
      name = 'simple1'
      if (present(dgp)) name = trim(dgp)
      select case (name)
      case ('simple1')
      case ('type1', 'type2', 'type3', 'type4', 'type5')
      case default
         info = -3
         return
      end select
      if (present(x_input)) then
         n_obs = size(x_input,1)
         p_obs = size(x_input,2)
      else
         n_obs = n
         p_obs = p
      end if
      if (n_obs < 1 .or. p_obs < 1) then
         info = -1
         return
      end if
      if (name /= 'simple1' .and. p_obs < 5) then
         info = -2
         return
      end if
      max_time = 1.0_dp
      probability_time = 0.6_dp
      select case (name)
      case ('simple1')
      case ('type1')
         max_time = 1.5_dp
         probability_time = 0.8_dp
      case ('type2')
         max_time = 2.0_dp
         probability_time = 1.2_dp
      case ('type3')
         max_time = 15.0_dp
         probability_time = 10.0_dp
      case ('type4')
         max_time = 3.0_dp
         probability_time = 2.0_dp
      case ('type5')
         max_time = 2.0_dp
         probability_time = 0.17_dp
      end select
      if (present(y_max)) max_time = y_max
      if (present(y0)) probability_time = y0
      if (max_time <= 0.0_dp .or. probability_time < 0.0_dp) then
         info = -4
         return
      end if
      correlation = 0.0_dp
      if (present(rho)) correlation = rho
      if (abs(correlation) > 1.0_dp) then
         info = -5
         return
      end if
      monte_carlo_count = 10000
      if (present(n_mc)) monte_carlo_count = n_mc
      if (monte_carlo_count < 1) then
         info = -6
         return
      end if
      local_seed = 42
      if (present(seed)) local_seed = seed
      call rng%seed(local_seed)
      if (present(x_input)) then
         x = x_input
      else
         allocate(x(n_obs,p_obs))
         if (abs(correlation) <= epsilon(1.0_dp)) then
            do j = 1, p_obs
               do i = 1, n_obs
                  x(i,j) = rng%uniform()
               end do
            end do
         else
            innovation_scale = sqrt(max(0.0_dp, 1.0_dp - correlation**2))
            do i = 1, n_obs
               scale = rng%normal()
               x(i,1) = normal_cdf(scale)
               do j = 2, p_obs
                  scale = correlation * scale + innovation_scale * rng%normal()
                  x(i,j) = normal_cdf(scale)
               end do
            end do
         end if
      end if
      allocate(y(n_obs), w(n_obs), event(n_obs), cate(n_obs), cate_probability(n_obs), sign_values(n_obs))

      select case (name)
      case ('simple1')
         do i = 1, n_obs
            w(i) = rng%bernoulli(0.5_dp)
            scale = max(x(i,1), 1.0e-8_dp)
            failure = min(rng%exponential() * scale + real(w(i),dp), max_time)
            censor = 2.0_dp * rng%uniform()
            y(i) = min(failure, censor)
            event(i) = merge(1, 0, failure <= censor)
            cate(i) = restricted_mean_simple(scale, max_time, 1.0_dp) - &
                      restricted_mean_simple(scale, max_time, 0.0_dp)
            cate_probability(i) = exp(-max(probability_time - 1.0_dp, 0.0_dp) / scale) - &
                                  exp(-probability_time / scale)
            sign_values(i) = 1.0_dp
         end do
      case ('type1')
         do i = 1, n_obs
            indicator = merge(1.0_dp, 0.0_dp, x(i,1) < 0.5_dp)
            propensity = 0.25_dp * (1.0_dp + beta_density_2_4(x(i,1)))
            w(i) = rng%bernoulli(propensity)
            failure = exp(-1.85_dp - 0.8_dp * indicator + 0.7_dp * sqrt(x(i,2)) + &
                          0.2_dp * x(i,3) + (0.7_dp - 0.4_dp * indicator - &
                          0.4_dp * sqrt(x(i,2))) * real(w(i),dp) + rng%normal())
            failure = min(failure, max_time)
            denominator = exp(-1.75_dp - 0.5_dp * sqrt(x(i,2)) + 0.2_dp * x(i,3) + &
                              (1.15_dp + 0.5_dp * indicator - 0.3_dp * sqrt(x(i,2))) * real(w(i),dp))
            censor = sqrt(rng%exponential() / denominator)
            y(i) = min(failure, censor)
            event(i) = merge(1, 0, failure <= censor)
            sign_values(i) = signum(0.7_dp - 0.4_dp * indicator - 0.4_dp * sqrt(x(i,2)))
         end do
         allocate(monte_carlo_draw(monte_carlo_count))
         do k = 1, monte_carlo_count
            monte_carlo_draw(k) = rng%normal()
         end do
         do i = 1, n_obs
            indicator = merge(1.0_dp, 0.0_dp, x(i,1) < 0.5_dp)
            cate_sum = 0.0_dp
            probability_sum = 0.0_dp
            do k = 1, monte_carlo_count
               failure0 = exp(-1.85_dp - 0.8_dp * indicator + 0.7_dp * sqrt(x(i,2)) + &
                              0.2_dp * x(i,3) + monte_carlo_draw(k))
               failure1 = failure0 * exp(0.7_dp - 0.4_dp * indicator - 0.4_dp * sqrt(x(i,2)))
               cate_sum = cate_sum + min(failure1, max_time) - min(failure0, max_time)
               probability_sum = probability_sum + indicator_value(failure1 > probability_time) - &
                                 indicator_value(failure0 > probability_time)
            end do
            cate(i) = cate_sum / real(monte_carlo_count, dp)
            cate_probability(i) = probability_sum / real(monte_carlo_count, dp)
         end do
      case ('type2', 'type5')
         do i = 1, n_obs
            propensity = 0.25_dp * (1.0_dp + beta_density_2_4(x(i,1)))
            w(i) = rng%bernoulli(propensity)
            scale = merge(-0.4_dp, -0.5_dp, name == 'type5') + x(i,2)
            failure = (rng%exponential() / exp(x(i,1) + scale * real(w(i),dp)))**2
            failure = min(failure, max_time)
            if (name == 'type2') then
               censor = 3.0_dp * rng%uniform()
            else
               j = 1 + modulo(i - 1, p_obs)
               censor = exp(x(1,j) - x(i,3) * real(w(i),dp) + rng%normal())
            end if
            y(i) = min(failure, censor)
            event(i) = merge(1, 0, failure <= censor)
            sign_values(i) = -signum(scale)
         end do
         allocate(monte_carlo_draw(monte_carlo_count))
         do k = 1, monte_carlo_count
            monte_carlo_draw(k) = rng%exponential()
         end do
         do i = 1, n_obs
            scale = merge(-0.4_dp, -0.5_dp, name == 'type5') + x(i,2)
            cate_sum = 0.0_dp
            probability_sum = 0.0_dp
            do k = 1, monte_carlo_count
               failure0 = (monte_carlo_draw(k) / exp(x(i,1)))**2
               failure1 = (monte_carlo_draw(k) / exp(x(i,1) + scale))**2
               cate_sum = cate_sum + min(failure1, max_time) - min(failure0, max_time)
               probability_sum = probability_sum + indicator_value(failure1 > probability_time) - &
                                 indicator_value(failure0 > probability_time)
            end do
            cate(i) = cate_sum / real(monte_carlo_count, dp)
            cate_probability(i) = probability_sum / real(monte_carlo_count, dp)
         end do
      case ('type3', 'type4')
         do i = 1, n_obs
            if (name == 'type3') then
               propensity = 0.25_dp * (1.0_dp + beta_density_2_4(x(i,1)))
               lambda0 = x(i,2)**2 + x(i,3) + 6.0_dp
               lambda1 = lambda0 + 2.0_dp * (sqrt(x(i,1)) - 0.3_dp)
               denominator = 12.0_dp + log(1.0_dp + exp(x(i,3)))
               sign_values(i) = signum(sqrt(x(i,1)) - 0.3_dp)
            else
               propensity = logistic(x(i,1)) * logistic(x(i,2))
               lambda0 = x(i,2) + x(i,3)
               lambda1 = lambda0 + max(0.0_dp, x(i,1) - 0.3_dp)
               denominator = 1.0_dp + log(1.0_dp + exp(x(i,3)))
               sign_values(i) = signum(max(0.0_dp, x(i,1) - 0.3_dp))
               if (x(i,1) < 0.3_dp) sign_values(i) = ieee_value(0.0_dp, ieee_quiet_nan)
            end if
            w(i) = rng%bernoulli(propensity)
            scale = merge(lambda1, lambda0, w(i) == 1)
            failure = min(real(rng%poisson(scale),dp), max_time)
            censor = real(rng%poisson(denominator),dp)
            y(i) = min(failure, censor)
            event(i) = merge(1, 0, failure <= censor)
            cate_sum = 0.0_dp
            probability_sum = 0.0_dp
            do k = 1, monte_carlo_count
               failure0 = real(rng%poisson(lambda0),dp)
               failure1 = real(rng%poisson(lambda1),dp)
               cate_sum = cate_sum + min(failure1, max_time) - min(failure0, max_time)
               probability_sum = probability_sum + indicator_value(failure1 > probability_time) - &
                                 indicator_value(failure0 > probability_time)
            end do
            cate(i) = cate_sum / real(monte_carlo_count, dp)
            cate_probability(i) = probability_sum / real(monte_carlo_count, dp)
         end do
      end select
      if (present(cate_sign)) cate_sign = sign_values
   end subroutine generate_causal_survival_data

   pure elemental real(dp) function normal_cdf(x) result(value)
      real(dp), intent(in) :: x !! Standard-normal quantile.

      value = 0.5_dp * (1.0_dp + erf(x / sqrt(2.0_dp)))
   end function normal_cdf

   pure elemental real(dp) function indicator_value(condition) result(value)
      logical, intent(in) :: condition !! Condition converted to zero or one.

      value = merge(1.0_dp, 0.0_dp, condition)
   end function indicator_value

   pure elemental real(dp) function signum(x) result(value)
      real(dp), intent(in) :: x !! Value whose mathematical sign is returned.

      if (x > 0.0_dp) then
         value = 1.0_dp
      else if (x < 0.0_dp) then
         value = -1.0_dp
      else
         value = 0.0_dp
      end if
   end function signum

   pure elemental real(dp) function restricted_mean_simple(scale, horizon, treatment) result(value)
      real(dp), intent(in) :: scale !! Positive exponential scale multiplying the untreated failure time.
      real(dp), intent(in) :: horizon !! Administrative truncation horizon for restricted mean survival.
      real(dp), intent(in) :: treatment !! Deterministic treatment shift, zero or one in the simple1 benchmark.
      real(dp) :: remaining

      if (horizon <= treatment) then
         value = horizon
      else
         remaining = horizon - treatment
         value = treatment + scale * (1.0_dp - exp(-remaining / scale))
      end if
   end function restricted_mean_simple

end module grf_simulation
