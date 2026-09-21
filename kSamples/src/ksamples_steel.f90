module ksamples_steel
   use, intrinsic :: ieee_arithmetic, only : ieee_negative_inf, ieee_positive_inf, ieee_value
   use r_kinds, only : dp
   use ksamples_types, only : steel_confint_result, steel_result
   use ksamples_utils, only : average_ranks, multinomial_count
   use ksamples_utils, only : initial_group_labels, next_group_labels, grouped_values
   use ksamples_utils, only : normal_cdf, normal_pdf, normal_quantile, shuffle_real
   use ksamples_utils, only : sort_real_in_place
   implicit none
   private

   public :: steel_test, steel_confint

contains

   subroutine steel_test(x, ns, alternative, method, nsim, keep_null, continuity_correction, result)
      real(dp), intent(in) :: x(:) !! Samples concatenated with the common control first, followed by treatment groups.
      integer, intent(in) :: ns(:) !! Positive sample sizes; ns(1) is the control size and sum(ns)=size(x).
      character(len=*), intent(in), optional :: alternative !! Direction: greater, less, or two-sided; default greater.
      character(len=*), intent(in), optional :: method !! Requested p-value method: asymptotic, simulated, or exact.
      integer, intent(in), optional :: nsim !! Number of simulations and exact-enumeration threshold; default 10000.
      logical, intent(in), optional :: keep_null !! Retain scalar Steel null statistics when true.
      logical, intent(in), optional :: continuity_correction !! Apply continuity corrections in normal p-values; default true.
      type(steel_result), intent(out) :: result !! Statistics, adjusted p-values, and optional null distribution.
      real(dp), allocatable :: ranks(:), mu(:), tau(:), sig(:), corr(:), null(:), work(:), grouped(:)
      integer, allocatable :: labels(:)
      real(dp) :: sig0, ncomb, observed
      integer :: alt_code, simulations, stored, i
      logical :: keep, correct, more
      character(len=10) :: chosen, direction
      if (size(ns) < 2 .or. sum(ns) /= size(x) .or. any(ns <= 0)) &
         error stop 'steel_test: invalid group sizes'
      direction = 'greater'
      if (present(alternative)) direction = adjustl(alternative)
      select case (trim(direction))
      case ('greater')
         alt_code = 1
      case ('less')
         alt_code = -1
      case ('two-sided', 'two.sided')
         alt_code = 0
         direction = 'two-sided'
      case default
         error stop 'steel_test: invalid alternative'
      end select
      result%alternative = direction
      call average_ranks(x, ranks)
      allocate(result%w(size(ns) - 1), result%w_standardized(size(ns) - 1))
      allocate(mu(size(ns) - 1), tau(size(ns) - 1), sig(size(ns) - 1), corr(size(ns) - 1))
      call steel_parameters(x, ranks, ns, mu, sig0, sig, tau, corr)
      call steel_w_values(ranks, ns, result%w)
      result%w_standardized = (result%w - mu)/tau
      observed = steel_scalar(result%w_standardized, alt_code)
      result%statistic = observed
      correct = .true.
      if (present(continuity_correction)) correct = continuity_correction
      call steel_normal_pvalues(mu, sig0, sig, tau, result%w, ns(2:), corr, &
         alt_code, correct, result%adjusted_asymptotic_p)
      result%asymptotic_p = minval(result%adjusted_asymptotic_p)
      chosen = 'asymptotic'
      if (present(method)) chosen = adjustl(method)
      simulations = 10000
      if (present(nsim)) simulations = max(1, nsim)
      keep = .false.
      if (present(keep_null)) keep = keep_null
      if (trim(chosen) == 'asymptotic') then
         result%method = 'asymptotic'
         return
      end if
      ncomb = multinomial_count(ns)
      if (trim(chosen) == 'exact' .and. ncomb <= real(simulations,dp) .and. ncomb <= real(huge(1),dp)) then
         allocate(null(nint(ncomb)))
         call initial_group_labels(ns, labels)
         stored = 0
         more = .true.
         do while (more)
            call grouped_values(ranks, labels, size(ns), grouped)
            stored = stored + 1
            null(stored) = steel_statistic_from_ranks(grouped, ns, mu, tau, alt_code)
            call next_group_labels(labels, more)
         end do
         result%method = 'exact'
      else
         allocate(null(simulations), work(size(ranks)))
         work = ranks
         do i = 1, simulations
            call shuffle_real(work)
            null(i) = steel_statistic_from_ranks(work, ns, mu, tau, alt_code)
         end do
         result%method = 'simulated'
      end if
      select case (alt_code)
      case (1)
         result%randomization_p = real(count(null >= observed),dp)/real(size(null),dp)
      case (-1)
         result%randomization_p = real(count(null <= observed),dp)/real(size(null),dp)
      case default
         result%randomization_p = real(count(null >= observed),dp)/real(size(null),dp)
      end select
      allocate(result%adjusted_randomization_p(size(ns) - 1))
      do i = 1, size(ns) - 1
         select case (alt_code)
         case (1)
            result%adjusted_randomization_p(i) = &
               real(count(null >= result%w_standardized(i)),dp)/real(size(null),dp)
         case (-1)
            result%adjusted_randomization_p(i) = &
               real(count(null <= result%w_standardized(i)),dp)/real(size(null),dp)
         case default
            result%adjusted_randomization_p(i) = &
               real(count(null >= abs(result%w_standardized(i))),dp)/real(size(null),dp)
         end select
      end do
      if (keep) then
         allocate(result%null_dist(size(null)))
         result%null_dist = null
      end if
   end subroutine steel_test

   pure subroutine steel_parameters(x, ranks, ns, mu, sig0, sig, tau, corr)
      real(dp), intent(in) :: x(:) !! Original pooled observations in group order.
      real(dp), intent(in) :: ranks(:) !! Average ranks of x in the same order.
      integer, intent(in) :: ns(:) !! Control and treatment sample sizes.
      real(dp), intent(out) :: mu(:) !! Null means of the treatment-control Mann-Whitney statistics.
      real(dp), intent(out) :: sig0 !! Common-control component standard deviation under the observed tie pattern.
      real(dp), intent(out) :: sig(:) !! Treatment-specific residual standard deviations under the observed tie pattern.
      real(dp), intent(out) :: tau(:) !! Standard deviations of each Mann-Whitney statistic under randomization.
      real(dp), intent(out) :: corr(:) !! Pairwise continuity corrections: 0.25 with ties and 0.5 without ties.
      real(dp), allocatable :: sorted(:)
      integer :: i, left, right, n, n0, offset
      real(dp) :: d2, d3, n2, n3, sig02
      n = size(x)
      n0 = ns(1)
      sorted = ranks
      call sort_real_in_place(sorted)
      d2 = 0.0_dp
      d3 = 0.0_dp
      left = 1
      do while (left <= n)
         right = left
         do while (right < n)
            if (sorted(right + 1) /= sorted(left)) exit
            right = right + 1
         end do
         d2 = d2 + real((right - left + 1)*(right - left),dp)
         d3 = d3 + real((right - left + 1)*(right - left)*(right - left - 1),dp)
         left = right + 1
      end do
      n2 = real(n*(n - 1),dp)
      if (n > 2) then
         n3 = real(n*(n - 1)*(n - 2),dp)
      else
         n3 = huge(1.0_dp)
      end if
      sig02 = real(n0,dp)/12.0_dp*(1.0_dp - d3/n3)
      sig0 = sqrt(max(0.0_dp, sig02))
      offset = n0
      do i = 1, size(ns) - 1
         mu(i) = 0.5_dp*real(n0*ns(i + 1),dp)
         sig(i) = sqrt(max(0.0_dp, real(n0*ns(i + 1),dp)/12.0_dp* &
            (real(n0 + 1,dp) - 3.0_dp*d2/n2 - real(n0 - 2,dp)*d3/n3)))
         tau(i) = sqrt(real(ns(i + 1)*ns(i + 1),dp)*sig02 + sig(i)*sig(i))
         if (pair_has_ties(x(1:n0), x(offset + 1:offset + ns(i + 1)))) then
            corr(i) = 0.25_dp
         else
            corr(i) = 0.5_dp
         end if
         offset = offset + ns(i + 1)
      end do
   end subroutine steel_parameters

   pure logical function pair_has_ties(control, treatment) result(has_tie)
      real(dp), intent(in) :: control(:) !! Common-control observations.
      real(dp), intent(in) :: treatment(:) !! One treatment sample.
      integer :: i, j
      has_tie = .false.
      do i = 1, size(control)
         do j = i + 1, size(control)
            if (control(i) == control(j)) then
               has_tie = .true.
               return
            end if
         end do
         do j = 1, size(treatment)
            if (control(i) == treatment(j)) then
               has_tie = .true.
               return
            end if
         end do
      end do
      do i = 1, size(treatment)
         do j = i + 1, size(treatment)
            if (treatment(i) == treatment(j)) then
               has_tie = .true.
               return
            end if
         end do
      end do
   end function pair_has_ties

   pure subroutine steel_w_values(ranks, ns, w)
      real(dp), intent(in) :: ranks(:) !! Pooled midranks arranged by control and treatment groups.
      integer, intent(in) :: ns(:) !! Control and treatment sample sizes.
      real(dp), intent(out) :: w(:) !! Mann-Whitney treatment-versus-control statistics.
      integer :: i, j, m, offset
      w = 0.0_dp
      offset = ns(1)
      do i = 1, size(ns) - 1
         do j = offset + 1, offset + ns(i + 1)
            do m = 1, ns(1)
               if (ranks(m) < ranks(j)) then
                  w(i) = w(i) + 1.0_dp
               else if (ranks(m) == ranks(j)) then
                  w(i) = w(i) + 0.5_dp
               end if
            end do
         end do
         offset = offset + ns(i + 1)
      end do
   end subroutine steel_w_values

   pure real(dp) function steel_statistic_from_ranks(ranks, ns, mu, tau, alt_code) result(statistic)
      real(dp), intent(in) :: ranks(:) !! Pooled midranks regrouped into control and treatment samples.
      integer, intent(in) :: ns(:) !! Control and treatment sample sizes.
      real(dp), intent(in) :: mu(:) !! Null means of pairwise Mann-Whitney statistics.
      real(dp), intent(in) :: tau(:) !! Pairwise randomization standard deviations.
      integer, intent(in) :: alt_code !! Direction code: 1 greater, -1 less, 0 two-sided.
      real(dp) :: w(size(ns) - 1), standardized(size(ns) - 1)
      call steel_w_values(ranks, ns, w)
      standardized = (w - mu)/tau
      statistic = steel_scalar(standardized, alt_code)
   end function steel_statistic_from_ranks

   pure real(dp) function steel_scalar(standardized, alt_code) result(statistic)
      real(dp), intent(in) :: standardized(:) !! Pairwise standardized Mann-Whitney statistics.
      integer, intent(in) :: alt_code !! Direction code: 1 greater, -1 less, 0 two-sided.
      select case (alt_code)
      case (1)
         statistic = maxval(standardized)
      case (-1)
         statistic = minval(standardized)
      case default
         statistic = maxval(abs(standardized))
      end select
   end function steel_scalar

   pure subroutine steel_normal_pvalues(mu, sig0, sig, tau, w, ni, corr, alt_code, correct, pvalues)
      real(dp), intent(in) :: mu(:) !! Null means of the pairwise Mann-Whitney statistics.
      real(dp), intent(in) :: sig0 !! Common-control normal component standard deviation.
      real(dp), intent(in) :: sig(:) !! Treatment-specific normal component standard deviations.
      real(dp), intent(in) :: tau(:) !! Marginal standard deviations of the pairwise statistics.
      real(dp), intent(in) :: w(:) !! Observed pairwise Mann-Whitney statistics.
      integer, intent(in) :: ni(:) !! Treatment sample sizes.
      real(dp), intent(in) :: corr(:) !! Pairwise continuity corrections.
      integer, intent(in) :: alt_code !! Direction code: 1 greater, -1 less, 0 two-sided.
      logical, intent(in) :: correct !! Whether continuity corrections are active.
      real(dp), allocatable, intent(out) :: pvalues(:) !! Steel-adjusted asymptotic p-values for each treatment comparison.
      integer, parameter :: n_integral = 1200
      integer :: j, i, iz, weight
      real(dp) :: cc, svalue, probability, h, z, factor, upper, lower, term
      allocate(pvalues(size(w)))
      h = 20.0_dp/real(n_integral, dp)
      do j = 1, size(w)
         cc = 0.0_dp
         if (correct) cc = corr(j)
         select case (alt_code)
         case (1)
            svalue = (w(j) - cc - mu(j))/tau(j)
         case (-1)
            svalue = (w(j) + cc - mu(j))/tau(j)
         case default
            svalue = abs(w(j) - cc - mu(j))/tau(j)
         end select
         probability = 0.0_dp
         do iz = 0, n_integral
            z = -10.0_dp + real(iz, dp)*h
            factor = 1.0_dp
            do i = 1, size(ni)
               upper = (svalue*tau(i) - real(ni(i), dp)*sig0*z)/sig(i)
               select case (alt_code)
               case (1)
                  factor = factor*normal_cdf(upper)
               case (-1)
                  factor = factor*(1.0_dp - normal_cdf(upper))
               case default
                  lower = (-svalue*tau(i) - real(ni(i), dp)*sig0*z)/sig(i)
                  factor = factor*(normal_cdf(upper) - normal_cdf(lower))
               end select
            end do
            term = normal_pdf(z)*factor
            if (iz == 0 .or. iz == n_integral) then
               weight = 1
            else if (mod(iz, 2) == 0) then
               weight = 2
            else
               weight = 4
            end if
            probability = probability + real(weight, dp)*term
         end do
         probability = probability*h/3.0_dp
         pvalues(j) = min(1.0_dp, max(0.0_dp, 1.0_dp - probability))
      end do
   end subroutine steel_normal_pvalues

   pure subroutine steel_confint(x, ns, conf_level, alternative, method, nsim, result)
      real(dp), intent(in) :: x(:) !! Samples concatenated with control first and treatments following in ns order.
      integer, intent(in) :: ns(:) !! Positive sample sizes; ns(1) is the common control.
      real(dp), intent(in), optional :: conf_level !! Simultaneous confidence level in (0,1); default 0.95.
      character(len=*), intent(in), optional :: alternative !! less, greater, or two-sided; default two-sided.
      character(len=*), intent(in), optional :: method !! Requested method; asymptotic is implemented for interval calibration.
      integer, intent(in), optional :: nsim !! Reserved simulation/exact calibration count; retained for API correspondence.
      type(steel_confint_result), intent(out) :: result !! Conservative and closest asymptotic simultaneous interval bounds.
      integer, allocatable :: ell_up_p(:), ell_up_c(:), ell_low_p(:), ell_low_c(:)
      integer :: i, offset, n0, nt, unused_nsim
      real(dp), allocatable :: mu(:), tau(:), ni_real(:), diffs(:)
      real(dp) :: level, gamma, cgamma
      real(dp) :: prob_up_p, prob_up_c, prob_low_p, prob_low_c
      real(dp) :: pinf, ninf
      character(len=10) :: direction, requested
      if (size(ns) < 2 .or. sum(ns) /= size(x) .or. any(ns <= 0)) &
         error stop 'steel_confint: invalid group sizes'
      level = 0.95_dp
      if (present(conf_level)) level = conf_level
      if (level <= 0.0_dp .or. level >= 1.0_dp) error stop 'steel_confint: conf_level must be in (0,1)'
      direction = 'two.sided'
      if (present(alternative)) direction = adjustl(alternative)
      if (trim(direction) == 'two-sided') direction = 'two.sided'
      if (trim(direction) /= 'less' .and. trim(direction) /= 'greater' .and. &
          trim(direction) /= 'two.sided') error stop 'steel_confint: invalid alternative'
      requested = 'asymptotic'
      if (present(method)) requested = adjustl(method)
      select case (trim(requested))
      case ('asymptotic', 'exact', 'simulated')
      case default
         error stop 'steel_confint: method must be asymptotic, exact, or simulated'
      end select
      unused_nsim = 10000
      if (present(nsim)) unused_nsim = nsim
      if (unused_nsim < 1) error stop 'steel_confint: nsim must be positive'
      result%method = 'asymptotic'
      result%alternative = direction
      n0 = ns(1)
      nt = size(ns) - 1
      allocate(mu(nt), tau(nt), ni_real(nt))
      ni_real = real(ns(2:),dp)
      mu = 0.5_dp*real(n0,dp)*ni_real
      tau = sqrt(real(n0,dp)*ni_real*(real(n0,dp) + ni_real + 1.0_dp)/12.0_dp)
      gamma = level
      if (trim(direction) == 'two.sided') gamma = 1.0_dp - 0.5_dp*(1.0_dp - level)
      cgamma = qmax_wilcox(gamma, ns)
      allocate(ell_up_p(nt), ell_up_c(nt), ell_low_p(nt), ell_low_c(nt))
      ell_up_p = huge(1)
      ell_up_c = huge(1)
      ell_low_p = -huge(1)
      ell_low_c = -huge(1)
      prob_up_p = 1.0_dp
      prob_up_c = 1.0_dp
      prob_low_p = 1.0_dp
      prob_low_c = 1.0_dp
      if (trim(direction) /= 'greater') &
         call choose_upper_indices(cgamma, gamma, n0, ns(2:), mu, tau, ell_up_p, ell_up_c, prob_up_p, prob_up_c)
      if (trim(direction) /= 'less') &
         call choose_lower_indices(cgamma, gamma, n0, ns(2:), mu, tau, ell_low_p, ell_low_c, prob_low_p, prob_low_c)
      allocate(result%lower_conservative(nt), result%upper_conservative(nt))
      allocate(result%lower_closest(nt), result%upper_closest(nt))
      pinf = ieee_value(0.0_dp, ieee_positive_inf)
      ninf = ieee_value(0.0_dp, ieee_negative_inf)
      result%lower_conservative = ninf
      result%lower_closest = ninf
      result%upper_conservative = pinf
      result%upper_closest = pinf
      offset = n0
      do i = 1, nt
         call pair_differences(x(1:n0), x(offset + 1:offset + ns(i + 1)), diffs)
         if (trim(direction) /= 'greater') then
            if (ell_up_p(i) <= size(diffs)) result%upper_conservative(i) = diffs(max(1,ell_up_p(i)))
            if (ell_up_c(i) <= size(diffs)) result%upper_closest(i) = diffs(max(1,ell_up_c(i)))
         end if
         if (trim(direction) /= 'less') then
            if (ell_low_p(i) >= 1) result%lower_conservative(i) = diffs(min(size(diffs),ell_low_p(i)))
            if (ell_low_c(i) >= 1) result%lower_closest(i) = diffs(min(size(diffs),ell_low_c(i)))
         end if
         offset = offset + ns(i + 1)
      end do
      select case (trim(direction))
      case ('less')
         result%achieved_conservative = prob_up_p
         result%achieved_closest = prob_up_c
      case ('greater')
         result%achieved_conservative = prob_low_p
         result%achieved_closest = prob_low_c
      case default
         result%achieved_conservative = prob_up_p + prob_low_p - 1.0_dp
         result%achieved_closest = prob_up_c + prob_low_c - 1.0_dp
      end select
   end subroutine steel_confint

   pure subroutine choose_upper_indices(cgamma, gamma, n0, ni, mu, tau, conservative, closest, pcons, pclose)
      real(dp), intent(in) :: cgamma !! Asymptotic quantile of the maximum standardized Wilcoxon statistic.
      real(dp), intent(in) :: gamma !! One-sided target joint probability.
      integer, intent(in) :: n0 !! Common-control sample size.
      integer, intent(in) :: ni(:) !! Treatment sample sizes.
      real(dp), intent(in) :: mu(:) !! Continuous-null Mann-Whitney means.
      real(dp), intent(in) :: tau(:) !! Continuous-null Mann-Whitney standard deviations.
      integer, intent(out) :: conservative(:) !! Upper order-statistic indices with joint probability at least gamma.
      integer, intent(out) :: closest(:) !! Upper indices with joint probability closest to gamma.
      real(dp), intent(out) :: pcons !! Achieved probability for conservative indices.
      real(dp), intent(out) :: pclose !! Achieved probability for closest indices.
      integer :: i, candidates(size(ni),3)
      real(dp) :: probs(3), u(size(ni))
      candidates(:,1) = ceiling(tau*cgamma + mu + 1.0_dp)
      candidates(:,2) = nint(tau*cgamma + mu + 1.0_dp)
      candidates(:,3) = floor(tau*cgamma + mu + 1.0_dp)
      do i = 1, 3
         call upper_probability_vector(candidates(:,i), n0, ni, mu, tau, u)
         probs(i) = prob_wilcox(u, n0, ni)
      end do
      call choose_candidate(candidates, probs, gamma, conservative, closest, pcons, pclose)
   end subroutine choose_upper_indices

   pure subroutine choose_lower_indices(cgamma, gamma, n0, ni, mu, tau, conservative, closest, pcons, pclose)
      real(dp), intent(in) :: cgamma !! Asymptotic quantile of the maximum standardized Wilcoxon statistic.
      real(dp), intent(in) :: gamma !! One-sided target joint probability.
      integer, intent(in) :: n0 !! Common-control sample size.
      integer, intent(in) :: ni(:) !! Treatment sample sizes.
      real(dp), intent(in) :: mu(:) !! Continuous-null Mann-Whitney means.
      real(dp), intent(in) :: tau(:) !! Continuous-null Mann-Whitney standard deviations.
      integer, intent(out) :: conservative(:) !! Lower order-statistic indices with joint probability at least gamma.
      integer, intent(out) :: closest(:) !! Lower indices with joint probability closest to gamma.
      real(dp), intent(out) :: pcons !! Achieved probability for conservative indices.
      real(dp), intent(out) :: pclose !! Achieved probability for closest indices.
      integer :: i, candidates(size(ni),3), ell(size(ni))
      real(dp) :: probs(3), u(size(ni))
      ell = ceiling(real(n0*ni,dp) - tau*cgamma - mu)
      candidates(:,1) = floor(real(n0*ni,dp) - tau*cgamma - mu)
      candidates(:,2) = nint(real(n0*ni,dp) - tau*cgamma - mu)
      candidates(:,3) = ell
      do i = 1, 3
         call lower_probability_vector(candidates(:,i), n0, ni, mu, tau, u)
         probs(i) = prob_wilcox(u, n0, ni)
      end do
      call choose_candidate(candidates, probs, gamma, conservative, closest, pcons, pclose)
   end subroutine choose_lower_indices

   pure subroutine choose_candidate(candidates, probabilities, target, conservative, closest, pcons, pclose)
      integer, intent(in) :: candidates(:, :) !! Three candidate order-index vectors in columns.
      real(dp), intent(in) :: probabilities(:) !! Joint probabilities corresponding to candidate columns.
      real(dp), intent(in) :: target !! Target joint probability.
      integer, intent(out) :: conservative(:) !! Candidate vector closest to target subject to probability >= target.
      integer, intent(out) :: closest(:) !! Candidate vector having probability closest to target without constraint.
      real(dp), intent(out) :: pcons !! Joint probability of conservative choice.
      real(dp), intent(out) :: pclose !! Joint probability of closest choice.
      integer :: i, icons, iclose
      icons = 1
      pcons = huge(1.0_dp)
      do i = 1, size(probabilities)
         if (probabilities(i) >= target .and. probabilities(i) < pcons) then
            pcons = probabilities(i)
            icons = i
         end if
      end do
      if (pcons == huge(1.0_dp)) then
         icons = maxloc(probabilities,dim=1)
         pcons = probabilities(icons)
      end if
      iclose = 1
      do i = 2, size(probabilities)
         if (abs(probabilities(i) - target) < abs(probabilities(iclose) - target)) iclose = i
      end do
      pclose = probabilities(iclose)
      conservative = candidates(:,icons)
      closest = candidates(:,iclose)
   end subroutine choose_candidate

   pure subroutine upper_probability_vector(indices, n0, ni, mu, tau, u)
      integer, intent(in) :: indices(:) !! Candidate upper order-statistic indices.
      integer, intent(in) :: n0 !! Common-control sample size.
      integer, intent(in) :: ni(:) !! Treatment sample sizes.
      real(dp), intent(in) :: mu(:) !! Mann-Whitney means.
      real(dp), intent(in) :: tau(:) !! Mann-Whitney standard deviations.
      real(dp), intent(out) :: u(:) !! Standardized thresholds for joint upper-bound probability.
      integer :: i
      do i = 1, size(ni)
         if (indices(i) > n0*ni(i)) then
            u(i) = huge(1.0_dp)
         else
            u(i) = (real(indices(i) - 1,dp) - mu(i))/tau(i)
         end if
      end do
   end subroutine upper_probability_vector

   pure subroutine lower_probability_vector(indices, n0, ni, mu, tau, u)
      integer, intent(in) :: indices(:) !! Candidate lower order-statistic indices.
      integer, intent(in) :: n0 !! Common-control sample size.
      integer, intent(in) :: ni(:) !! Treatment sample sizes.
      real(dp), intent(in) :: mu(:) !! Mann-Whitney means.
      real(dp), intent(in) :: tau(:) !! Mann-Whitney standard deviations.
      real(dp), intent(out) :: u(:) !! Reflected standardized thresholds used for joint lower-bound probability.
      integer :: i
      do i = 1, size(ni)
         if (indices(i) < 1) then
            u(i) = huge(1.0_dp)
         else
            u(i) = (real(n0*ni(i) - indices(i),dp) - mu(i))/tau(i)
         end if
      end do
   end subroutine lower_probability_vector

   pure real(dp) function qmax_wilcox(probability, ns) result(quantile)
      real(dp), intent(in) :: probability !! Target CDF probability in (0,1).
      integer, intent(in) :: ns(:) !! Common-control and treatment sample sizes.
      real(dp) :: lo, hi, mid
      integer :: iteration
      lo = normal_quantile(probability) - 2.0_dp
      hi = normal_quantile(probability) + 2.0_dp
      do while (pmax_wilcox(lo, ns) > probability)
         lo = lo - 1.0_dp
      end do
      do while (pmax_wilcox(hi, ns) < probability)
         hi = hi + 1.0_dp
      end do
      do iteration = 1, 80
         mid = 0.5_dp*(lo + hi)
         if (pmax_wilcox(mid, ns) < probability) then
            lo = mid
         else
            hi = mid
         end if
      end do
      quantile = 0.5_dp*(lo + hi)
   end function qmax_wilcox

   pure real(dp) function pmax_wilcox(uvalue, ns) result(probability)
      real(dp), intent(in) :: uvalue !! Common standardized upper threshold for all treatment comparisons.
      integer, intent(in) :: ns(:) !! Common-control and treatment sample sizes.
      integer, parameter :: n_integral = 1600
      real(dp) :: f1(size(ns) - 1), f2(size(ns) - 1), h, z, term, total
      integer :: n0, iz, weight
      n0 = ns(1)
      f1 = sqrt(1.0_dp + real(ns(2:), dp)/real(n0 + 1, dp))
      f2 = sqrt(real(ns(2:), dp)/real(n0 + 1, dp))
      h = 20.0_dp/real(n_integral, dp)
      total = 0.0_dp
      do iz = 0, n_integral
         z = -10.0_dp + real(iz, dp)*h
         term = normal_pdf(z)*product(normal_cdf(uvalue*f1 - z*f2))
         if (iz == 0 .or. iz == n_integral) then
            weight = 1
         else if (mod(iz, 2) == 0) then
            weight = 2
         else
            weight = 4
         end if
         total = total + real(weight, dp)*term
      end do
      probability = total*h/3.0_dp
   end function pmax_wilcox

   pure real(dp) function prob_wilcox(u, n0, ni) result(probability)
      real(dp), intent(in) :: u(:) !! Per-comparison standardized upper thresholds.
      integer, intent(in) :: n0 !! Common-control sample size.
      integer, intent(in) :: ni(:) !! Treatment sample sizes.
      integer, parameter :: n_integral = 1600
      real(dp) :: f1(size(ni)), f2(size(ni)), h, z, factor, term, total
      integer :: i, iz, weight
      f1 = sqrt(1.0_dp + real(ni, dp)/real(n0 + 1, dp))
      f2 = sqrt(real(ni, dp)/real(n0 + 1, dp))
      h = 20.0_dp/real(n_integral, dp)
      total = 0.0_dp
      do iz = 0, n_integral
         z = -10.0_dp + real(iz, dp)*h
         factor = 1.0_dp
         do i = 1, size(u)
            if (u(i) < 0.5_dp*huge(1.0_dp)) factor = factor*normal_cdf(u(i)*f1(i) - z*f2(i))
         end do
         term = normal_pdf(z)*factor
         if (iz == 0 .or. iz == n_integral) then
            weight = 1
         else if (mod(iz, 2) == 0) then
            weight = 2
         else
            weight = 4
         end if
         total = total + real(weight, dp)*term
      end do
      probability = total*h/3.0_dp
   end function prob_wilcox

   pure subroutine pair_differences(control, treatment, differences)
      real(dp), intent(in) :: control(:) !! Common-control sample.
      real(dp), intent(in) :: treatment(:) !! One treatment sample.
      real(dp), allocatable, intent(out) :: differences(:) !! Sorted treatment-minus-control pairwise differences.
      integer :: i, j, n
      allocate(differences(size(control)*size(treatment)))
      n = 0
      do i = 1, size(treatment)
         do j = 1, size(control)
            n = n + 1
            differences(n) = treatment(i) - control(j)
         end do
      end do
      call sort_real_in_place(differences)
   end subroutine pair_differences

end module ksamples_steel
