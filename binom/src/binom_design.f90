! SPDX-License-Identifier: GPL-2.0-or-later
module binom_design
   use binom_intervals
   use binom_types, only: binom_ok, binom_domain, binom_convergence, binom_resource
   use binom_numerics
   use r_distributions, only: r_qnorm, r_pnorm, r_dnorm
   use r_transforms, only: r_log1p
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private
   integer, parameter, public :: binom_enumeration_limit = 10000
   public :: binom_coverage, binom_length, binom_power, cloglog_sample_size, binom_sim, binom_sim_counts
   public :: binom_design_result, binom_sim_result
   type :: binom_design_result
      real(dp) :: p_null = 0.5_dp
      real(dp) :: p_alt = 0.5_dp
      real(dp) :: delta = 0.0_dp
      real(dp) :: alpha = 0.05_dp
      real(dp) :: power = 0.8_dp
      real(dp) :: n = 0.0_dp
      real(dp) :: phi = 1.0_dp
      integer :: status = binom_ok
   end type binom_design_result
   type :: binom_sim_result
      type(binom_interval) :: capture_interval
      real(dp) :: width = 0.0_dp
      integer :: captured = 0
      integer :: replicates = 0
      integer :: status = binom_ok
   end type binom_sim_result
contains
   pure elemental real(dp) function binom_coverage(p, n, method, conf_level) result(value)
      real(dp), intent(in) :: p !! True binomial success probability in [0,1].
      integer, intent(in) :: n !! Positive trial count, no greater than binom_enumeration_limit.
      character(len=*), intent(in), optional :: method !! Interval method from binom_methods; default wilson.
      real(dp), intent(in), optional :: conf_level !! Confidence level strictly between zero and one; default 0.95.
      type(binom_interval) :: ci
      integer :: x
      real(dp) :: profile_step
      value = binom_nan()
      if (.not. ieee_is_finite(p) .or. p < 0.0_dp .or. p > 1.0_dp) return
      if (n < 1 .or. n > binom_enumeration_limit) return
      value = 0.0_dp
      profile_step = collection_profile_step(0, n, conf_level)
      do x = 0, n
         ci = collection_interval(x, n, method, conf_level, profile_step)
         if (ci%status /= binom_ok) then
            value = binom_nan()
            return
         end if
         if (p >= ci%lower .and. p <= ci%upper) value = value + binomial_mass(real(x, dp), real(n, dp), p)
      end do
   end function binom_coverage

   pure elemental real(dp) function binom_length(p, n, method, conf_level) result(value)
      real(dp), intent(in) :: p !! True binomial success probability in [0,1].
      integer, intent(in) :: n !! Positive trial count, no greater than binom_enumeration_limit.
      character(len=*), intent(in), optional :: method !! Interval method from binom_methods; default wilson.
      real(dp), intent(in), optional :: conf_level !! Confidence level strictly between zero and one; default 0.95.
      type(binom_interval) :: ci
      integer :: x
      real(dp) :: profile_step
      value = binom_nan()
      if (.not. ieee_is_finite(p) .or. p < 0.0_dp .or. p > 1.0_dp) return
      if (n < 1 .or. n > binom_enumeration_limit) return
      value = 0.0_dp
      profile_step = collection_profile_step(0, n, conf_level)
      do x = 0, n
         ci = collection_interval(x, n, method, conf_level, profile_step)
         if (ci%status /= binom_ok) then
            value = binom_nan()
            return
         end if
         value = value + (ci%upper - ci%lower) * binomial_mass(real(x, dp), real(n, dp), p)
      end do
   end function binom_length

   pure elemental real(dp) function binom_power(p_alt, n, p, alpha, phi, alternative, method) result(value)
      real(dp), intent(in) :: p_alt !! Alternative success probability, interior for transformed-normal methods.
      real(dp), intent(in), optional :: n !! Positive trial count; default 100.
      real(dp), intent(in), optional :: p !! Null probability, interior for transformed-normal methods; default 0.5.
      real(dp), intent(in), optional :: alpha !! Type-I error probability strictly between zero and one; default 0.05.
      real(dp), intent(in), optional :: phi !! Positive dispersion multiplier for normal approximations; default one.
      character(len=*), intent(in), optional :: alternative !! two.sided, greater or less; default two.sided.
      character(len=*), intent(in), optional :: method !! cloglog, logit, probit, asymp, lrt or exact; default cloglog.
      real(dp) :: trials, null_p, significance, dispersion, z, g0, g1, s0, s1, down, up, ci_alpha, log_scale
      character(len=16) :: which, selected
      type(binom_interval) :: ci
      trials = 100.0_dp
      if (present(n)) trials = n
      null_p = 0.5_dp
      if (present(p)) null_p = p
      significance = 0.05_dp
      if (present(alpha)) significance = alpha
      dispersion = 1.0_dp
      if (present(phi)) dispersion = phi
      which = 'two.sided'
      if (present(alternative)) which = alternative
      selected = 'cloglog'
      if (present(method)) selected = method
      value = binom_nan()
      if (.not. all(ieee_is_finite([p_alt, trials, null_p, significance, dispersion]))) return
      if (trials <= 0.0_dp .or. dispersion <= 0.0_dp .or. significance <= 0.0_dp .or. significance >= 1.0_dp) return
      if (p_alt < 0.0_dp .or. p_alt > 1.0_dp .or. null_p < 0.0_dp .or. null_p > 1.0_dp) return
      if (which /= 'two.sided' .and. which /= 'greater' .and. which /= 'less') return
      ci_alpha = significance
      if (which == 'two.sided') ci_alpha = 0.5_dp * significance
      z = r_qnorm(ci_alpha, lower_tail=.false.)
      if (selected == 'lrt' .or. selected == 'exact') then
         ci = binom_confint(trials * null_p, trials, selected, 1.0_dp - 2.0_dp * ci_alpha)
         if (ci%status /= binom_ok) return
         up = 0.0_dp
         down = 0.0_dp
         if (which /= 'less') up = beta_probability(ci%upper, trials * p_alt + 1.0_dp, &
            trials * (1.0_dp - p_alt), .true.)
         if (which /= 'greater') down = beta_probability(ci%lower, trials * p_alt, &
            trials * (1.0_dp - p_alt) + 1.0_dp)
         value = up + down
         return
      end if
      if (selected == 'asymp' .or. selected == 'asymptotic') then
         g0 = null_p
         g1 = p_alt
         s0 = sqrt(null_p * (1.0_dp - null_p))
         s1 = sqrt(p_alt * (1.0_dp - p_alt))
      else
         if (p_alt <= 0.0_dp .or. p_alt >= 1.0_dp .or. null_p <= 0.0_dp .or. null_p >= 1.0_dp) return
         call link_moments(null_p, selected, g0, s0)
         call link_moments(p_alt, selected, g1, s1)
      end if
      if (.not. all(ieee_is_finite([g0, g1, s0, s1]))) return
      if (p_alt == null_p) then
         value = significance
         return
      end if
      log_scale = 0.5_dp * (log(dispersion) - log(trials))
      if (max(s0, s1) > 0.0_dp) then
         if (log_scale + log(max(s0, s1)) >= log(huge(1.0_dp)) - log(max(1.0_dp, abs(z)))) return
      end if
      if (s0 > 0.0_dp) s0 = exp(log(s0) + log_scale)
      if (s1 > 0.0_dp) s1 = exp(log(s1) + log_scale)
      up = scaled_normal_probability(g1 - g0 - z * s0, s1)
      down = scaled_normal_probability(g0 - g1 - z * s0, s1)
      if (selected == 'cloglog') then
         value = up
         up = down
         down = value
      end if
      if (which == 'less') then
         value = down
      else if (which == 'greater') then
         value = up
      else
         value = up + down
      end if
   end function binom_power

   pure elemental real(dp) function scaled_normal_probability(numerator, scale) result(value)
      real(dp), intent(in) :: numerator !! Signed distance from a normal rejection threshold.
      real(dp), intent(in) :: scale !! Nonnegative alternative standard deviation; zero uses the degenerate limit.
      if (scale == 0.0_dp) then
         value = binom_nan()
         if (numerator > 0.0_dp) value = 1.0_dp
         if (numerator < 0.0_dp) value = 0.0_dp
      else
         value = r_pnorm(numerator / scale)
      end if
   end function scaled_normal_probability

   pure subroutine link_moments(p, method, gamma, sd)
      real(dp), intent(in) :: p !! Probability strictly between zero and one.
      character(len=*), intent(in) :: method !! Transformed-normal method identifier.
      real(dp), intent(out) :: gamma !! Link-transformed binomial probability.
      real(dp), intent(out) :: sd !! Standard deviation for a unit sample and unit dispersion.
      gamma = binom_nan()
      sd = binom_nan()
      select case (method)
      case ('cloglog')
         gamma = log(-log(p))
         sd = exp(0.5_dp * (r_log1p(-p) - log(p)) - log(abs(log(p))))
      case ('logit')
         gamma = logit(p)
         sd = exp(-0.5_dp * (log(p) + r_log1p(-p)))
      case ('probit')
         gamma = r_qnorm(p)
         sd = exp(0.5_dp * (log(p) + r_log1p(-p)) - r_dnorm(gamma, log_density=.true.))
      case ('asymp', 'asymptotic')
         gamma = p
         sd = sqrt(p * (1.0_dp - p))
      end select
   end subroutine link_moments

   pure elemental function cloglog_sample_size(p_alt, n, p, power, alpha, alternative, &
      exact_n, recompute_power, phi) result(design)
      real(dp), intent(in), optional :: p_alt !! Alternative probability; omit with n to solve for a detectable alternative.
      real(dp), intent(in), optional :: n !! Positive trial count; omit with p_alt supplied to solve for sample size.
      real(dp), intent(in), optional :: p !! Null probability strictly between zero and one; default 0.5.
      real(dp), intent(in), optional :: power !! Target power strictly between alpha and one; default 0.8.
      real(dp), intent(in), optional :: alpha !! Type-I error probability strictly between zero and one; default 0.05.
      character(len=*), intent(in), optional :: alternative !! two.sided (default), greater or less; two-sided inversion is upward.
      logical, intent(in), optional :: exact_n !! Return the fractional sample size, floored at one; default false rounds upward.
      logical, intent(in), optional :: recompute_power !! Recompute attained power after rounding sample size; default false.
      real(dp), intent(in), optional :: phi !! Positive dispersion multiplier; default one.
      type(binom_design_result) :: design
      real(dp) :: g0, g1, sd0, sd1, z, zp, low, high, midpoint, flow, fhigh, fmid, tail
      character(len=16) :: which
      logical :: fractional, recompute
      integer :: i
      design%status = binom_domain
      if (present(p)) design%p_null = p
      if (present(power)) design%power = power
      if (present(alpha)) design%alpha = alpha
      if (present(phi)) design%phi = phi
      which = 'two.sided'
      if (present(alternative)) which = alternative
      fractional = .false.
      if (present(exact_n)) fractional = exact_n
      recompute = .false.
      if (present(recompute_power)) recompute = recompute_power
      if (.not. all(ieee_is_finite([design%p_null, design%power, design%alpha, design%phi]))) return
      if (design%p_null <= 0.0_dp .or. design%p_null >= 1.0_dp .or. design%phi <= 0.0_dp) return
      if (design%alpha <= 0.0_dp .or. design%alpha >= 1.0_dp) return
      if (which /= 'two.sided' .and. which /= 'greater' .and. which /= 'less') return
      if (.not. present(p_alt) .and. .not. present(n)) return
      tail = design%alpha
      if (which == 'two.sided') tail = 0.5_dp * tail
      z = r_qnorm(tail, lower_tail=.false.)
      if (present(n)) then
         if (.not. ieee_is_finite(n) .or. n <= 0.0_dp) return
         design%n = n
      end if
      if (present(p_alt)) then
         if (.not. ieee_is_finite(p_alt) .or. p_alt <= 0.0_dp .or. p_alt >= 1.0_dp) return
         design%p_alt = p_alt
      end if
      if (present(n) .and. present(p_alt)) then
         design%power = binom_power(p_alt, n, design%p_null, design%alpha, design%phi, which, 'cloglog')
      else
         if (design%power <= design%alpha .or. design%power >= 1.0_dp) return
         if (.not. present(n)) then
            if (design%p_alt == design%p_null) return
            call link_moments(design%p_null, 'cloglog', g0, sd0)
            call link_moments(design%p_alt, 'cloglog', g1, sd1)
            zp = r_qnorm(design%power)
            design%n = design%phi * ((sd0 * z + sd1 * zp) / (g1 - g0)) ** 2
            if (.not. ieee_is_finite(design%n)) then
               design%status = binom_resource
               return
            end if
            if (fractional) then
               design%n = max(1.0_dp, design%n)
            else
               ! REAL ceiling avoids overflow in default INTEGER for very large designs.
               midpoint = aint(design%n)
               if (midpoint < design%n) midpoint = midpoint + 1.0_dp
               design%n = midpoint
            end if
            if (recompute .and. .not. fractional) design%power = &
               binom_power(design%p_alt, design%n, design%p_null, design%alpha, design%phi, which, 'cloglog')
         else
            low = design%p_null
            high = 1.0_dp - 1.0e-8_dp
            if (which == 'less') then
               low = 1.0e-8_dp
               high = design%p_null
            end if
            flow = binom_power(low, design%n, design%p_null, design%alpha, design%phi, which) - design%power
            fhigh = binom_power(high, design%n, design%p_null, design%alpha, design%phi, which) - design%power
            if (.not. all(ieee_is_finite([flow, fhigh]))) then
               design%status = binom_convergence
               return
            end if
            if ((flow > 0.0_dp .and. fhigh > 0.0_dp) .or. (flow < 0.0_dp .and. fhigh < 0.0_dp)) then
               design%status = binom_convergence
               return
            end if
            do i = 1, 100
               midpoint = low + 0.5_dp * (high - low)
               fmid = binom_power(midpoint, design%n, design%p_null, design%alpha, design%phi, which) - design%power
               if (.not. ieee_is_finite(fmid)) then
                  design%status = binom_convergence
                  return
               end if
               if ((fmid <= 0.0_dp .and. flow >= 0.0_dp) .or. (fmid >= 0.0_dp .and. flow <= 0.0_dp)) then
                  high = midpoint
               else
                  low = midpoint
                  flow = fmid
               end if
            end do
            design%p_alt = low + 0.5_dp * (high - low)
         end if
      end if
      design%delta = design%p_alt - design%p_null
      if (.not. all(ieee_is_finite([design%n, design%p_alt, design%power]))) then
         design%status = binom_convergence
         return
      end if
      design%status = binom_ok
   end function cloglog_sample_size

   pure function binom_sim_counts(counts, n, p, method, conf_level) result(simulation)
      integer, intent(in) :: counts(:) !! Caller-supplied binomial success counts for each Monte Carlo replicate, nonempty.
      integer, intent(in) :: n !! Positive trial count per replicate.
      real(dp), intent(in) :: p !! True success probability in [0,1] against which coverage is assessed.
      character(len=*), intent(in), optional :: method !! Confidence interval method; default wilson.
      real(dp), intent(in), optional :: conf_level !! Confidence level for both replicate and capture intervals; default 0.95.
      type(binom_sim_result) :: simulation
      type(binom_interval) :: ci
      integer :: i
      real(dp) :: profile_step
      simulation%status = binom_domain
      simulation%width = binom_nan()
      if (size(counts) == 0 .or. n < 1) return
      if (.not. ieee_is_finite(p) .or. p < 0.0_dp .or. p > 1.0_dp) return
      if (any(counts < 0) .or. any(counts > n)) return
      simulation%width = 0.0_dp
      simulation%replicates = size(counts)
      profile_step = collection_profile_step(counts(1), n, conf_level)
      do i = 1, size(counts)
         ci = collection_interval(counts(i), n, method, conf_level, profile_step)
         if (ci%status /= binom_ok) then
            simulation%status = ci%status
            return
         end if
         if (p >= ci%lower .and. p <= ci%upper) simulation%captured = simulation%captured + 1
         simulation%width = simulation%width + (ci%upper - ci%lower) / real(size(counts), dp)
      end do
      simulation%capture_interval = binom_confint(real(simulation%captured, dp), &
         real(size(counts), dp), method, conf_level)
      simulation%status = simulation%capture_interval%status
   end function binom_sim_counts

   pure function binom_sim(uniforms, n, p, method, conf_level) result(simulation)
      real(dp), intent(in) :: uniforms(:) !! Uniforms in [0,1), one per replicate; caller controls the RNG.
      integer, intent(in) :: n !! Trial count from one through binom_enumeration_limit.
      real(dp), intent(in) :: p !! True success probability in [0,1].
      character(len=*), intent(in), optional :: method !! Confidence interval method; default wilson.
      real(dp), intent(in), optional :: conf_level !! Confidence level for replicate and capture intervals; default 0.95.
      type(binom_sim_result) :: simulation
      integer, allocatable :: counts(:)
      integer :: i, low, high, middle
      real(dp) :: cumulative
      simulation%status = binom_domain
      simulation%width = binom_nan()
      if (size(uniforms) == 0 .or. n < 1) return
      if (n > binom_enumeration_limit) then
         simulation%status = binom_resource
         return
      end if
      if (.not. all(ieee_is_finite(uniforms)) .or. .not. ieee_is_finite(p)) return
      if (any(uniforms < 0.0_dp) .or. any(uniforms >= 1.0_dp) .or. p < 0.0_dp .or. p > 1.0_dp) return
      allocate(counts(size(uniforms)))
      do i = 1, size(uniforms)
         if (p == 0.0_dp) then
            counts(i) = 0
         else if (p == 1.0_dp) then
            counts(i) = n
         else
            low = 0
            high = n
            do while (low < high)
               middle = low + (high - low) / 2
               cumulative = beta_probability(p, real(middle + 1, dp), real(n - middle, dp), .true.)
               if (cumulative >= uniforms(i)) then
                  high = middle
               else
                  low = middle + 1
               end if
            end do
            counts(i) = low
         end if
      end do
      simulation = binom_sim_counts(counts, n, p, method, conf_level)
   end function binom_sim
   pure function collection_profile_step(first, n, conf_level) result(step)
      integer, intent(in) :: first !! First success count in a collection passed to R profile.
      integer, intent(in) :: n !! Common positive trial count for the collection.
      real(dp), intent(in), optional :: conf_level !! Confidence probability; default 0.95.
      real(dp) :: step, level, trials
      level = 0.95_dp
      if (present(conf_level)) level = conf_level
      step = binom_nan()
      if (.not. ieee_is_finite(level) .or. level <= 0.0_dp .or. level >= 1.0_dp) return
      trials = real(n, dp)
      if (first == 0 .or. first == n) trials = trials + 1.0_dp
      step = r_qnorm(0.5_dp * (1.0_dp - level) / trials, lower_tail=.false.) / 5.0_dp
   end function collection_profile_step

   pure function collection_interval(x, n, method, conf_level, profile_step) result(ci)
      integer, intent(in) :: x !! Current success count in a binomial collection.
      integer, intent(in) :: n !! Common positive trial count.
      character(len=*), intent(in), optional :: method !! Named confidence method; default wilson.
      real(dp), intent(in), optional :: conf_level !! Confidence probability; default 0.95.
      real(dp), intent(in) :: profile_step !! Step computed from the first observation, matching R lazy default evaluation.
      type(binom_interval) :: ci
      if (present(method)) then
         if (method == 'profile') then
            ci = binom_profile(real(x, dp), real(n, dp), conf_level, del=profile_step)
            return
         end if
      end if
      ci = binom_confint(real(x, dp), real(n, dp), method, conf_level)
   end function collection_interval
end module binom_design
