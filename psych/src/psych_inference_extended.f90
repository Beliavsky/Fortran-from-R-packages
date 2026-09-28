module psych_inference_extended
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use r_kinds, only : dp
   use psych_basic, only : fisher_z, fisher_r, d_to_t, t_to_d, winsor_variance
   use psych_correlations, only : pseudoinverse
   use psych_distributions, only : normal_cdf, normal_quantile, student_t_cdf, chi_square_cdf, f_cdf
   implicit none
   private
   public :: alpha_confidence_interval, effect_size_confidence_interval
   public :: robust_effect_size_two_groups, correlation_test_general
   public :: correlation_information, jennrich_correlation_test
   public :: correlation_matrix_test, correlation_normal_test

contains

   pure subroutine alpha_confidence_interval(alpha, n_obs, n_var, p_value, lower, upper, r_bar)
      real(dp), intent(in) :: alpha !! Observed coefficient alpha.
      integer, intent(in) :: n_obs !! Number of observations used to estimate alpha.
      integer, intent(in) :: n_var !! Number of items entering alpha.
      real(dp), intent(in) :: p_value !! Two-sided tail probability, usually 0.05.
      real(dp), intent(out) :: lower !! Lower confidence limit for alpha.
      real(dp), intent(out) :: upper !! Upper confidence limit for alpha.
      real(dp), intent(out) :: r_bar !! Average inter-item correlation implied by alpha.
      real(dp) :: q_low, q_high
      integer :: df1, df2

      df1 = max(1, n_obs - 1)
      df2 = max(1, (n_obs - 1) * max(1, n_var - 1))
      q_low = f_quantile(0.5_dp * p_value, real(df1, dp), real(df2, dp))
      q_high = f_quantile(1.0_dp - 0.5_dp * p_value, real(df1, dp), real(df2, dp))
      upper = 1.0_dp - (1.0_dp - alpha) * q_low
      lower = 1.0_dp - (1.0_dp - alpha) * q_high
      if (n_var > 0) then
         r_bar = alpha / (real(n_var, dp) - alpha * real(n_var - 1, dp))
      else
         r_bar = 0.0_dp
      end if
   end subroutine alpha_confidence_interval

   pure subroutine effect_size_confidence_interval(d, n1, n2, alpha, lower, upper)
      real(dp), intent(in) :: d !! Standardized mean difference.
      integer, intent(in) :: n1 !! First-group sample size.
      integer, intent(in) :: n2 !! Second-group sample size.
      real(dp), intent(in) :: alpha !! Two-sided tail probability, usually 0.05.
      real(dp), intent(out) :: lower !! Approximate lower confidence limit for d.
      real(dp), intent(out) :: upper !! Approximate upper confidence limit for d.
      real(dp) :: se, z
      integer :: ntotal

      ntotal = max(2, n1 + n2)
      if (n1 <= 0 .or. n2 <= 0) then
         lower = d
         upper = d
         return
      end if
      se = sqrt(real(ntotal, dp) / real(n1 * n2, dp) + d * d / (2.0_dp * real(max(1, ntotal - 2), dp)))
      z = normal_quantile(1.0_dp - 0.5_dp * alpha)
      lower = d - z * se
      upper = d + z * se
   end subroutine effect_size_confidence_interval

   pure subroutine robust_effect_size_two_groups(x, group, trim, means, variances, pooled_sd, robust_d, counts)
      real(dp), intent(in) :: x(:, :) !! Observation-by-variable matrix with IEEE NaNs allowed.
      integer, intent(in) :: group(:) !! Two-group labels; the two distinct integer values define the comparison.
      real(dp), intent(in) :: trim !! Fraction trimmed from each tail before computing group means.
      real(dp), allocatable, intent(out) :: means(:, :) !! Two-by-variable matrix of trimmed group means.
      real(dp), allocatable, intent(out) :: variances(:, :) !! Two-by-variable matrix of winsorized group variances.
      real(dp), allocatable, intent(out) :: pooled_sd(:) !! Per-variable pooled winsorized standard deviations.
      real(dp), allocatable, intent(out) :: robust_d(:) !! Algina-style robust standardized mean differences.
      integer, allocatable, intent(out) :: counts(:, :) !! Two-by-variable finite-observation counts.
      real(dp), allocatable :: values(:)
      integer :: labels(2), n_labels, i, j, g, n_use

      call first_two_labels(group, labels, n_labels)
      allocate(means(2, size(x, 2)), variances(2, size(x, 2)))
      allocate(pooled_sd(size(x, 2)), robust_d(size(x, 2)), counts(2, size(x, 2)))
      means = 0.0_dp
      variances = 0.0_dp
      pooled_sd = 0.0_dp
      robust_d = 0.0_dp
      counts = 0
      if (n_labels < 2) return

      do j = 1, size(x, 2)
         do g = 1, 2
            allocate(values(size(x, 1)))
            n_use = 0
            do i = 1, min(size(x, 1), size(group))
               if (group(i) /= labels(g)) cycle
               if (ieee_is_nan(x(i, j))) cycle
               n_use = n_use + 1
               values(n_use) = x(i, j)
            end do
            counts(g, j) = n_use
            if (n_use > 0) then
               call sort_real(values(1:n_use))
               means(g, j) = trimmed_mean_sorted(values(1:n_use), trim)
               variances(g, j) = winsor_variance(values(1:n_use), trim)
            end if
            deallocate(values)
         end do
         if (counts(1, j) + counts(2, j) > 2) then
            pooled_sd(j) = sqrt(max(0.0_dp, &
               (variances(1, j) * real(max(0, counts(1, j) - 1), dp) + &
                variances(2, j) * real(max(0, counts(2, j) - 1), dp)) / &
               real(max(1, counts(1, j) + counts(2, j) - 2), dp)))
         end if
         if (pooled_sd(j) > tiny(1.0_dp)) then
            robust_d(j) = 0.642_dp * (means(2, j) - means(1, j)) / pooled_sd(j)
         end if
      end do
   end subroutine robust_effect_size_two_groups

   pure subroutine correlation_test_general(n, r12, statistic, p_value, df, r34, r23, r13, r14, r24, n2, pooled, two_tailed)
      integer, intent(in) :: n !! Sample size for the first correlation or dependent-correlation set.
      real(dp), intent(in) :: r12 !! First correlation.
      real(dp), intent(out) :: statistic !! Student-t or normal-z statistic for the selected comparison.
      real(dp), intent(out) :: p_value !! Tail probability for the selected test.
      real(dp), intent(out) :: df !! Degrees of freedom for t cases; zero for z cases.
      real(dp), intent(in), optional :: r34 !! Second independent or nonoverlapping correlation.
      real(dp), intent(in), optional :: r23 !! Correlation linking variables 2 and 3 in dependent tests.
      real(dp), intent(in), optional :: r13 !! Alternative correlation sharing variable 1, or cross-correlation 1-3.
      real(dp), intent(in), optional :: r14 !! Cross-correlation 1-4 for nonoverlapping dependent tests.
      real(dp), intent(in), optional :: r24 !! Cross-correlation 2-4 for nonoverlapping dependent tests.
      integer, intent(in), optional :: n2 !! Sample size for an independent second correlation.
      logical, intent(in), optional :: pooled !! Use the pooled-correlation covariance approximation in four-variable tests.
      logical, intent(in), optional :: two_tailed !! Double the one-sided probability when true; default true.
      real(dp) :: av, cube, determin, diff, r1234, z12, z34, z1234, se
      real(dp) :: rr34, rr13, rr23, rr14, rr24, pooled_r
      integer :: second_n
      logical :: do_pooled, do_two

      do_two = .true.
      if (present(two_tailed)) do_two = two_tailed
      do_pooled = .true.
      if (present(pooled)) do_pooled = pooled
      df = 0.0_dp

      if (.not. present(r34) .and. .not. present(r13) .and. .not. present(r23)) then
         statistic = r12 * sqrt(real(max(1, n - 2), dp) / max(tiny(1.0_dp), 1.0_dp - r12 * r12))
         df = real(max(1, n - 2), dp)
         p_value = 1.0_dp - student_t_cdf(abs(statistic), df)
      else if (present(r34) .and. .not. present(r23)) then
         rr34 = r34
         second_n = n
         if (present(n2)) second_n = n2
         z12 = fisher_z(r12)
         z34 = fisher_z(rr34)
         se = sqrt(1.0_dp / real(max(1, n - 3), dp) + 1.0_dp / real(max(1, second_n - 3), dp))
         statistic = abs(z12 - z34) / max(tiny(1.0_dp), se)
         p_value = 1.0_dp - normal_cdf(abs(statistic))
      else if (present(r23) .and. present(r13) .and. .not. present(r14)) then
         rr23 = r23
         rr13 = r13
         diff = r12 - rr13
         determin = 1.0_dp - r12 * r12 - rr23 * rr23 - rr13 * rr13 + 2.0_dp * r12 * rr23 * rr13
         av = 0.5_dp * (r12 + rr13)
         cube = (1.0_dp - rr23) ** 3
         statistic = diff * sqrt(real(max(1, n - 1), dp) * (1.0_dp + rr23) / &
            max(tiny(1.0_dp), 2.0_dp * real(max(1, n - 1), dp) / real(max(1, n - 3), dp) * determin + av * av * cube))
         df = real(max(1, n - 3), dp)
         p_value = 1.0_dp - student_t_cdf(abs(statistic), df)
      else if (present(r34) .and. present(r23) .and. present(r13) .and. present(r14) .and. present(r24)) then
         rr34 = r34
         rr23 = r23
         rr13 = r13
         rr14 = r14
         rr24 = r24
         z12 = fisher_z(r12)
         z34 = fisher_z(rr34)
         if (do_pooled) then
            pooled_r = 0.5_dp * (r12 + rr34)
            r1234 = 0.5_dp * ((rr13 - pooled_r * rr23) * (rr24 - rr23 * pooled_r) + &
               (rr14 - rr13 * pooled_r) * (rr23 - pooled_r * rr13) + &
               (rr13 - rr14 * pooled_r) * (rr24 - pooled_r * rr14) + &
               (rr14 - pooled_r * rr24) * (rr23 - rr24 * pooled_r))
            z1234 = r1234 / max(tiny(1.0_dp), (1.0_dp - pooled_r * pooled_r) ** 2)
         else
            r1234 = 0.5_dp * ((rr13 - r12 * rr23) * (rr24 - rr23 * rr34) + &
               (rr14 - rr13 * rr34) * (rr23 - r12 * rr13) + &
               (rr13 - rr14 * rr34) * (rr24 - r12 * rr14) + &
               (rr14 - r12 * rr24) * (rr23 - rr24 * rr34))
            z1234 = r1234 / max(tiny(1.0_dp), (1.0_dp - r12 * r12) * (1.0_dp - rr34 * rr34))
         end if
         statistic = (z12 - z34) * sqrt(real(max(1, n - 3), dp) / max(tiny(1.0_dp), 2.0_dp * (1.0_dp - z1234)))
         p_value = 1.0_dp - normal_cdf(abs(statistic))
      else
         statistic = 0.0_dp
         p_value = 1.0_dp
      end if
      if (do_two) p_value = min(1.0_dp, 2.0_dp * p_value)
   end subroutine correlation_test_general

   pure subroutine correlation_information(r, n_pairs, fisher, information, chisq, df, fisher_information)
      real(dp), intent(in) :: r(:, :) !! Correlation matrix.
      real(dp), intent(in) :: n_pairs(:, :) !! Pairwise sample-size matrix with the same shape as `r`.
      logical, intent(in) :: fisher !! Apply Fisher-z transformation before accumulating information when true.
      real(dp), intent(out) :: information !! Sum of `(n-2)/(1-r^2)^2` over unique off-diagonal pairs.
      real(dp), intent(out) :: chisq !! Sum of `r^2*(n-3)` over unique off-diagonal pairs.
      real(dp), intent(out) :: df !! Number of unique off-diagonal correlations.
      real(dp), intent(out) :: fisher_information !! Sum of `n*(1+r^2)/(1-r^2)^2` over unique pairs.
      real(dp) :: value, denom, nij
      integer :: i, j, p

      p = min(size(r, 1), size(r, 2))
      information = 0.0_dp
      fisher_information = 0.0_dp
      chisq = 0.0_dp
      df = real(p * max(0, p - 1) / 2, dp)
      do j = 1, p
         do i = j + 1, p
            value = r(i, j)
            nij = n_pairs(i, j)
            if (ieee_is_nan(value) .or. ieee_is_nan(nij)) cycle
            if (fisher) value = fisher_z(max(-0.999999999_dp, min(0.999999999_dp, value)))
            denom = max(tiny(1.0_dp), (1.0_dp - value * value) ** 2)
            chisq = chisq + value * value * (nij - 3.0_dp)
            information = information + (nij - 2.0_dp) / denom
            fisher_information = fisher_information + nij * (1.0_dp + value * value) / denom
         end do
      end do
   end subroutine correlation_information

   subroutine jennrich_correlation_test(r1, r2, n1, n2, chisq, p_value, df, info)
      real(dp), intent(in) :: r1(:, :) !! First correlation matrix.
      real(dp), intent(in) :: r2(:, :) !! Second correlation matrix of matching size.
      integer, intent(in) :: n1 !! First sample size.
      integer, intent(in) :: n2 !! Second sample size.
      real(dp), intent(out) :: chisq !! Jennrich chi-square statistic.
      real(dp), intent(out) :: p_value !! Upper-tail chi-square probability.
      real(dp), intent(out) :: df !! Degrees of freedom `p*(p-1)/2`.
      integer, intent(out) :: info !! Zero on successful pseudoinverse operations.
      real(dp), allocatable :: r(:, :), s(:, :), rinv(:, :), sinv(:, :), z(:, :), d(:), work(:)
      real(dp) :: c
      integer :: i, p, info1, info2

      p = size(r1, 1)
      c = real(max(1, n1) * max(1, n2), dp) / real(max(1, n1 + n2), dp)
      allocate(r(p, p), s(p, p), z(p, p), d(p), work(p))
      r = (real(n1, dp) * r1 + real(n2, dp) * r2) / real(max(1, n1 + n2), dp)
      s = r * r
      call pseudoinverse(r, rinv, info=info1)
      call pseudoinverse(s, sinv, info=info2)
      info = max(info1, info2)
      if (info /= 0) then
         chisq = 0.0_dp
         p_value = 1.0_dp
         df = real(p * max(0, p - 1) / 2, dp)
         return
      end if
      z = sqrt(c) * matmul(rinv, r1 - r2)
      do i = 1, p
         d(i) = z(i, i)
      end do
      work = matmul(sinv, d)
      chisq = 0.5_dp * sum(z * z) - dot_product(d, work)
      df = real(p * max(0, p - 1) / 2, dp)
      p_value = 1.0_dp - chi_square_cdf(max(0.0_dp, chisq), df)
   end subroutine jennrich_correlation_test

   subroutine correlation_matrix_test(r1, n1, chisq, p_value, df, rmsea, srmr, r2, n2, info)
      real(dp), intent(in) :: r1(:, :) !! First correlation matrix.
      integer, intent(in) :: n1 !! First sample size.
      real(dp), intent(out) :: chisq !! Bartlett or two-matrix likelihood chi-square statistic.
      real(dp), intent(out) :: p_value !! Upper-tail chi-square probability.
      real(dp), intent(out) :: df !! Test degrees of freedom.
      real(dp), intent(out) :: rmsea !! Approximate RMSEA for the matrix comparison.
      real(dp), intent(out) :: srmr !! Root mean squared off-diagonal difference.
      real(dp), intent(in), optional :: r2(:, :) !! Optional second correlation matrix.
      integer, intent(in), optional :: n2 !! Second sample size; defaults to `n1`.
      integer, intent(out) :: info !! Zero when determinant and pseudoinverse computations succeed.
      real(dp), allocatable :: inv1(:, :), inv2(:, :), cross12(:, :), cross21(:, :), diff(:, :)
      real(dp) :: det1, det12, det21, e1, e2, z
      integer :: i, j, p, second_n, info1, info2

      p = size(r1, 1)
      df = real(p * max(0, p - 1) / 2, dp)
      rmsea = 0.0_dp
      srmr = 0.0_dp
      info = 0
      if (.not. present(r2)) then
         call determinant_local(r1, det1, info)
         if (info /= 0 .or. det1 <= 0.0_dp) then
            chisq = 0.0_dp
            p_value = 1.0_dp
            return
         end if
         chisq = -log(det1) * (real(max(1, n1 - 1), dp) - real(2 * p + 5, dp) / 6.0_dp)
         p_value = 1.0_dp - chi_square_cdf(max(0.0_dp, chisq), df)
         if (df > 0.0_dp .and. n1 > 1) rmsea = sqrt(max(chisq / (df * real(n1, dp)) - 1.0_dp / real(n1 - 1, dp), 0.0_dp))
         return
      end if

      second_n = n1
      if (present(n2)) second_n = n2
      call pseudoinverse(r1, inv1, info=info1)
      call pseudoinverse(r2, inv2, info=info2)
      info = max(info1, info2)
      if (info /= 0) then
         chisq = 0.0_dp
         p_value = 1.0_dp
         return
      end if
      cross12 = matmul(inv1, r2)
      cross21 = matmul(inv2, r1)
      call determinant_local(cross12, det12, info1)
      call determinant_local(cross21, det21, info2)
      info = max(info1, info2)
      if (info /= 0 .or. det12 <= 0.0_dp .or. det21 <= 0.0_dp) then
         chisq = 0.0_dp
         p_value = 1.0_dp
         return
      end if
      e1 = 0.5_dp * (trace_local(cross12) - log(det12) - real(p, dp))
      e2 = 0.5_dp * (trace_local(cross21) - log(det21) - real(p, dp))
      chisq = e1 * (real(max(1, n1 - 1), dp) - real(2 * p - 5, dp) / 6.0_dp) + &
         e2 * (real(max(1, second_n - 1), dp) - real(2 * p - 5, dp) / 6.0_dp)
      df = real(p * max(0, p - 1), dp)
      diff = r1 - r2
      z = 0.0_dp
      do j = 1, p
         do i = 1, p
            z = z + diff(i, j) * diff(i, j)
         end do
      end do
      if (df > 0.0_dp) srmr = sqrt(z / df)
      if (df > 0.0_dp .and. n1 > 1) rmsea = sqrt(max(chisq / (df * real(n1, dp)) - 1.0_dp / real(n1 - 1, dp), 0.0_dp))
      p_value = 1.0_dp - chi_square_cdf(max(0.0_dp, chisq), df)
   end subroutine correlation_matrix_test

   pure subroutine correlation_normal_test(r1, n1, fisher, chisq, p_value, df, rmsea, srmr, r2, n2)
      real(dp), intent(in) :: r1(:, :) !! First correlation matrix.
      integer, intent(in) :: n1 !! First sample size.
      logical, intent(in) :: fisher !! Compare Fisher-z transformed correlations when true.
      real(dp), intent(out) :: chisq !! Normal-theory chi-square statistic.
      real(dp), intent(out) :: p_value !! Upper-tail chi-square probability.
      real(dp), intent(out) :: df !! Number of unique correlations.
      real(dp), intent(out) :: rmsea !! Approximate RMSEA.
      real(dp), intent(out) :: srmr !! Root mean squared transformed-correlation magnitude or difference.
      real(dp), intent(in), optional :: r2(:, :) !! Optional second correlation matrix.
      integer, intent(in), optional :: n2 !! Second sample size, defaulting to `n1`.
      real(dp) :: a, b, e, n_eff, z
      integer :: i, j, p, second_n

      p = size(r1, 1)
      df = real(p * max(0, p - 1) / 2, dp)
      e = 0.0_dp
      second_n = n1
      if (present(n2)) second_n = n2
      do j = 1, p
         do i = j + 1, p
            a = r1(i, j)
            if (fisher) a = fisher_z(max(-0.999999999_dp, min(0.999999999_dp, a)))
            if (present(r2)) then
               b = r2(i, j)
               if (fisher) b = fisher_z(max(-0.999999999_dp, min(0.999999999_dp, b)))
               z = a - b
            else
               z = a
            end if
            e = e + z * z
         end do
      end do
      if (present(r2)) then
         n_eff = 2.0_dp * real(max(1, n1) * max(1, second_n), dp) / real(max(1, n1 + second_n), dp)
      else
         n_eff = real(max(1, n1), dp)
      end if
      chisq = e * max(0.0_dp, n_eff - 3.0_dp)
      if (df > 0.0_dp) then
         srmr = sqrt(e / df)
      else
         srmr = 0.0_dp
      end if
      if (df > 0.0_dp .and. n_eff > 1.0_dp) then
         rmsea = sqrt(max(chisq / (df * n_eff) - 1.0_dp / (n_eff - 1.0_dp), 0.0_dp))
      else
         rmsea = 0.0_dp
      end if
      p_value = 1.0_dp - chi_square_cdf(max(0.0_dp, chisq), df)
   end subroutine correlation_normal_test

   pure real(dp) function f_quantile(probability, df1, df2) result(value)
      real(dp), intent(in) :: probability !! Probability in the closed unit interval.
      real(dp), intent(in) :: df1 !! Positive numerator degrees of freedom.
      real(dp), intent(in) :: df2 !! Positive denominator degrees of freedom.
      real(dp) :: lo, hi, mid
      integer :: iter

      if (probability <= 0.0_dp) then
         value = 0.0_dp
         return
      end if
      if (probability >= 1.0_dp) then
         value = huge(1.0_dp)
         return
      end if
      lo = 0.0_dp
      hi = 1.0_dp
      do while (f_cdf(hi, df1, df2) < probability .and. hi < 1.0e12_dp)
         hi = 2.0_dp * hi
      end do
      do iter = 1, 120
         mid = 0.5_dp * (lo + hi)
         if (f_cdf(mid, df1, df2) < probability) then
            lo = mid
         else
            hi = mid
         end if
      end do
      value = 0.5_dp * (lo + hi)
   end function f_quantile

   pure subroutine first_two_labels(group, labels, n_labels)
      integer, intent(in) :: group(:) !! Integer group labels.
      integer, intent(out) :: labels(2) !! First two distinct labels encountered.
      integer, intent(out) :: n_labels !! Number of distinct labels found, capped at two.
      integer :: i

      labels = 0
      n_labels = 0
      do i = 1, size(group)
         if (n_labels == 0) then
            n_labels = 1
            labels(1) = group(i)
         else if (group(i) /= labels(1) .and. n_labels == 1) then
            n_labels = 2
            labels(2) = group(i)
            exit
         end if
      end do
   end subroutine first_two_labels

   pure subroutine sort_real(x)
      real(dp), intent(inout) :: x(:) !! Values sorted in ascending order in place.
      real(dp) :: value
      integer :: i, j

      do i = 2, size(x)
         value = x(i)
         j = i - 1
         do while (j >= 1)
            if (x(j) <= value) exit
            x(j + 1) = x(j)
            j = j - 1
         end do
         x(j + 1) = value
      end do
   end subroutine sort_real

   pure real(dp) function trimmed_mean_sorted(x, trim) result(value)
      real(dp), intent(in) :: x(:) !! Sorted finite values.
      real(dp), intent(in) :: trim !! Fraction removed from each tail.
      integer :: cut, first, last

      if (size(x) == 0) then
         value = 0.0_dp
         return
      end if
      cut = min(size(x) / 2, max(0, int(max(0.0_dp, min(0.499999_dp, trim)) * real(size(x), dp))))
      first = 1 + cut
      last = size(x) - cut
      if (last < first) then
         value = sum(x) / real(size(x), dp)
      else
         value = sum(x(first:last)) / real(last - first + 1, dp)
      end if
   end function trimmed_mean_sorted

   pure subroutine determinant_local(a, value, info)
      real(dp), intent(in) :: a(:, :) !! Square matrix whose determinant is required.
      real(dp), intent(out) :: value !! Matrix determinant.
      integer, intent(out) :: info !! Zero on success; nonzero for singular or nonsquare input.
      real(dp), allocatable :: work(:, :), row(:)
      real(dp) :: factor, pivot_value
      integer :: i, k, pivot, n

      n = size(a, 1)
      if (size(a, 2) /= n) then
         value = 0.0_dp
         info = 1
         return
      end if
      allocate(work(n, n), row(n))
      work = a
      value = 1.0_dp
      info = 0
      do k = 1, n
         pivot = k
         pivot_value = abs(work(k, k))
         do i = k + 1, n
            if (abs(work(i, k)) > pivot_value) then
               pivot = i
               pivot_value = abs(work(i, k))
            end if
         end do
         if (pivot_value <= epsilon(1.0_dp)) then
            value = 0.0_dp
            info = 1
            return
         end if
         if (pivot /= k) then
            row = work(k, :)
            work(k, :) = work(pivot, :)
            work(pivot, :) = row
            value = -value
         end if
         value = value * work(k, k)
         do i = k + 1, n
            factor = work(i, k) / work(k, k)
            work(i, k:n) = work(i, k:n) - factor * work(k, k:n)
         end do
      end do
   end subroutine determinant_local

   pure real(dp) function trace_local(a) result(value)
      real(dp), intent(in) :: a(:, :) !! Matrix whose diagonal is summed.
      integer :: i
      value = 0.0_dp
      do i = 1, min(size(a, 1), size(a, 2))
         value = value + a(i, i)
      end do
   end function trace_local

end module psych_inference_extended
