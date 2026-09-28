module psych_additional
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
   use r_kinds, only : dp
   use r_linalg, only : determinant, symmetric_eigen
   use psych_correlations, only : pseudoinverse
   use psych_basic, only : phi_coefficient, yule_q
   use psych_polychoric, only : tetrachoric_table
   use psych_descriptive, only : describe_vector, pairwise_count_matrix
   use psych_types, only : describe_result
   use psych_distributions, only : f_cdf, normal_cdf, normal_quantile, student_t_cdf
   use psych_stats, only : correlation_matrix, covariance_matrix, finite_correlation, finite_covariance, finite_quantile
   implicit none
   private
   public :: replication_probability, replication_probability_f, replication_probability_r
   public :: replication_probability_t, cohen_profile_similarity, factor_score_indeterminacy
   public :: cluster_fit_metrics, super_correlation_matrix, cor2_matrix
   public :: interbattery_loadings, block_random_design
   public :: scaling_fit, comorbidity_statistics
   public :: eigen_confidence_intervals, pairwise_sample_description

contains

   subroutine eigen_confidence_intervals(x, eigenvalues, normal_lower, normal_upper, bootstrap_lower, bootstrap_upper, &
      n_iter, alpha, info)
      real(dp), intent(in) :: x(:, :) !! Observation-by-variable numeric matrix used to estimate correlation eigenvalues.
      real(dp), allocatable, intent(out) :: eigenvalues(:) !! Eigenvalues of the observed correlation matrix, descending.
      real(dp), allocatable, intent(out) :: normal_lower(:) !! Normal-approximation lower confidence limits.
      real(dp), allocatable, intent(out) :: normal_upper(:) !! Normal-approximation upper confidence limits.
      real(dp), allocatable, intent(out) :: bootstrap_lower(:) !! Row-bootstrap lower quantiles for each eigenvalue.
      real(dp), allocatable, intent(out) :: bootstrap_upper(:) !! Row-bootstrap upper quantiles for each eigenvalue.
      integer, intent(in), optional :: n_iter !! Number of row-bootstrap replicates; defaults to 1000 and is at least one.
      real(dp), intent(in), optional :: alpha !! Two-sided interval error probability in `(0,1)`; defaults to `0.05`.
      integer, intent(out), optional :: info !! Zero on success or the first nonzero symmetric eigensolver status.
      real(dp), allocatable :: boot(:, :), cor(:, :), sample(:, :), values(:), vectors(:, :)
      real(dp) :: a, u, z
      integer :: b, i, j, local_info, n, nboot, p, row_index

      n = size(x, 1)
      p = size(x, 2)
      nboot = 1000
      if (present(n_iter)) nboot = max(1, n_iter)
      a = 0.05_dp
      if (present(alpha)) a = min(1.0_dp, max(0.0_dp, alpha))
      allocate(cor(p, p))
      call correlation_matrix(x, cor)
      call symmetric_eigen(cor, eigenvalues, vectors, local_info, descending=.true.)
      if (present(info)) info = local_info
      if (local_info /= 0) then
         allocate(normal_lower(0), normal_upper(0), bootstrap_lower(0), bootstrap_upper(0))
         return
      end if

      z = normal_quantile(0.5_dp * a)
      allocate(normal_lower(p), normal_upper(p), bootstrap_lower(p), bootstrap_upper(p), boot(nboot, p), sample(n, p))
      normal_lower = eigenvalues + z * abs(eigenvalues) * sqrt(2.0_dp / real(max(1, n), dp))
      normal_upper = eigenvalues - z * abs(eigenvalues) * sqrt(2.0_dp / real(max(1, n), dp))

      do b = 1, nboot
         do i = 1, n
            call random_number(u)
            row_index = min(n, int(u * real(n, dp)) + 1)
            sample(i, :) = x(row_index, :)
         end do
         call correlation_matrix(sample, cor)
         call symmetric_eigen(cor, values, vectors, local_info, descending=.true.)
         if (local_info /= 0) then
            if (present(info)) info = local_info
            return
         end if
         boot(b, :) = values
      end do
      do j = 1, p
         bootstrap_lower(j) = finite_quantile(boot(:, j), 0.5_dp * a)
         bootstrap_upper(j) = finite_quantile(boot(:, j), 1.0_dp - 0.5_dp * a)
      end do
   end subroutine eigen_confidence_intervals

   subroutine pairwise_sample_description(x, result, quantiles, sample_size, diagonal)
      real(dp), intent(in) :: x(:, :) !! Observation-by-item matrix whose NaNs mark missing responses.
      type(describe_result), intent(out) :: result !! Descriptive summary of sampled pairwise observation counts.
      real(dp), allocatable, intent(out) :: quantiles(:) !! Count quantiles at probabilities 0.05, 0.25, 0.50, 0.75, 0.95.
      integer, intent(in), optional :: sample_size !! Number of columns sampled without replacement; default 100 or all items.
      logical, intent(in), optional :: diagonal !! Include diagonal counts when true; default false.
      integer, allocatable :: counts(:, :), order(:)
      real(dp), allocatable :: sampled(:, :), values(:)
      real(dp) :: u
      integer :: i, j, k, n_items, n_selected, n_values, swap_index, temp
      logical :: include_diagonal

      n_items = size(x, 2)
      n_selected = min(100, n_items)
      if (present(sample_size)) n_selected = min(max(0, sample_size), n_items)
      include_diagonal = .false.
      if (present(diagonal)) include_diagonal = diagonal
      allocate(order(n_items))
      order = [(i, i = 1, n_items)]
      do i = n_items, 2, -1
         call random_number(u)
         swap_index = min(i, int(u * real(i, dp)) + 1)
         temp = order(i)
         order(i) = order(swap_index)
         order(swap_index) = temp
      end do

      allocate(sampled(size(x, 1), n_selected), counts(n_selected, n_selected))
      if (n_selected > 0) sampled = x(:, order(1:n_selected))
      call pairwise_count_matrix(sampled, counts)
      n_values = n_selected * (n_selected - 1) / 2
      if (include_diagonal) n_values = n_values + n_selected
      allocate(values(n_values))
      k = 0
      do j = 1, n_selected
         do i = j + merge(0, 1, include_diagonal), n_selected
            k = k + 1
            values(k) = real(counts(i, j), dp)
         end do
      end do
      call describe_vector(values, result)
      allocate(quantiles(5))
      quantiles(1) = finite_quantile(values, 0.05_dp)
      quantiles(2) = finite_quantile(values, 0.25_dp)
      quantiles(3) = finite_quantile(values, 0.50_dp)
      quantiles(4) = finite_quantile(values, 0.75_dp)
      quantiles(5) = finite_quantile(values, 0.95_dp)
   end subroutine pairwise_sample_description


   pure subroutine scaling_fit(model, data, goodness_of_fit, original_ss, residual_ss, residual, test, digits, rowwise, &
      probability)
      real(dp), intent(in) :: model(:) !! Scale locations used to form all pairwise choice probabilities.
      real(dp), intent(in) :: data(:, :) !! Observed square pairwise-choice matrix corresponding to `model`.
      real(dp), intent(out) :: goodness_of_fit !! Upstream fit coefficient `1 - residual_ss / original_ss`.
      real(dp), intent(out) :: original_ss !! Sum of squared finite observed values.
      real(dp), intent(out) :: residual_ss !! Sum of squared finite residuals.
      real(dp), allocatable, intent(out) :: residual(:, :) !! Rounded observed-minus-predicted residual matrix.
      character(len=*), intent(in), optional :: test !! Probability link: `logit`, `normal`, or `choice`; default `logit`.
      integer, intent(in), optional :: digits !! Decimal places retained in the returned residual matrix; default two.
      logical, intent(in), optional :: rowwise !! Reverse pairwise probabilities as in upstream rowwise mode; default true.
      real(dp), allocatable, intent(out), optional :: probability(:, :) !! Optional unrounded pairwise probability matrix.
      real(dp), allocatable :: adjusted(:), error(:, :), prob(:, :)
      real(dp) :: denom, scale
      integer :: i, j, n, ndigits
      logical :: reverse_rows
      character(len=16) :: link

      n = size(model)
      allocate(adjusted(n), error(n, n), prob(n, n), residual(n, n))
      adjusted = model
      link = 'logit'
      if (present(test)) link = trim(test)
      reverse_rows = .true.
      if (present(rowwise)) reverse_rows = rowwise
      ndigits = 2
      if (present(digits)) ndigits = digits
      if (trim(link) == 'choice' .and. minval(adjusted) <= 0.0_dp) adjusted = adjusted - minval(adjusted)

      do j = 1, n
         do i = 1, n
            select case (trim(link))
            case ('choice')
               denom = adjusted(i) + adjusted(j)
               if (abs(denom) <= tiny(1.0_dp)) then
                  prob(i, j) = ieee_value(0.0_dp, ieee_quiet_nan)
               else
                  prob(i, j) = adjusted(i) / denom
               end if
            case ('normal')
               prob(i, j) = normal_cdf(model(i) - model(j))
            case default
               prob(i, j) = 1.0_dp / (1.0_dp + exp(-(model(i) - model(j))))
            end select
            if (reverse_rows .and. .not. ieee_is_nan(prob(i, j))) prob(i, j) = 1.0_dp - prob(i, j)
            error(i, j) = data(i, j) - prob(i, j)
         end do
      end do

      original_ss = 0.0_dp
      residual_ss = 0.0_dp
      do j = 1, min(size(data, 2), n)
         do i = 1, min(size(data, 1), n)
            if (.not. ieee_is_nan(data(i, j))) original_ss = original_ss + data(i, j) * data(i, j)
            if (.not. ieee_is_nan(error(i, j))) residual_ss = residual_ss + error(i, j) * error(i, j)
         end do
      end do
      if (original_ss > 0.0_dp) then
         goodness_of_fit = 1.0_dp - residual_ss / original_ss
      else
         goodness_of_fit = ieee_value(0.0_dp, ieee_quiet_nan)
      end if

      scale = 10.0_dp ** real(ndigits, dp)
      residual = error
      do j = 1, n
         do i = 1, n
            if (.not. ieee_is_nan(residual(i, j))) residual(i, j) = anint(residual(i, j) * scale) / scale
         end do
      end do
      if (present(probability)) probability = prob
   end subroutine scaling_fit

   pure subroutine comorbidity_statistics(d1, d2, comorbidity, table2x2, phi, yule, tetrachoric)
      real(dp), intent(in) :: d1 !! Marginal prevalence or probability of the first dichotomous condition.
      real(dp), intent(in) :: d2 !! Marginal prevalence or probability of the second dichotomous condition.
      real(dp), intent(in) :: comorbidity !! Joint prevalence or probability that both conditions occur.
      real(dp), intent(out) :: table2x2(2, 2) !! Joint probability table in the upstream `comorbidity` arrangement.
      real(dp), intent(out) :: phi !! Pearson phi coefficient of the resulting 2-by-2 table.
      real(dp), intent(out) :: yule !! Yule Q association coefficient of the resulting table.
      real(dp), intent(out) :: tetrachoric !! Tetrachoric correlation estimated from the resulting table.

      table2x2 = reshape([comorbidity, d1 - comorbidity, d2 - comorbidity, &
         1.0_dp - d1 - d2 + comorbidity], [2, 2])
      phi = phi_coefficient(table2x2)
      yule = yule_q(table2x2)
      call tetrachoric_table(table2x2, tetrachoric)
   end subroutine comorbidity_statistics


   pure real(dp) function student_t_quantile_local(probability, df) result(value)
      real(dp), intent(in) :: probability !! Lower-tail probability in `[0,1]`.
      real(dp), intent(in) :: df !! Positive Student-t degrees of freedom.
      real(dp) :: lo, hi, mid
      integer :: iter

      if (probability <= 0.0_dp) then
         value = -huge(1.0_dp)
         return
      end if
      if (probability >= 1.0_dp) then
         value = huge(1.0_dp)
         return
      end if
      if (df <= 0.0_dp) then
         value = 0.0_dp
         return
      end if
      if (abs(probability - 0.5_dp) <= epsilon(1.0_dp)) then
         value = 0.0_dp
         return
      end if

      lo = -1.0_dp
      hi = 1.0_dp
      do while (student_t_cdf(lo, df) > probability .and. abs(lo) < 1.0e12_dp)
         lo = 2.0_dp * lo
      end do
      do while (student_t_cdf(hi, df) < probability .and. hi < 1.0e12_dp)
         hi = 2.0_dp * hi
      end do
      do iter = 1, 100
         mid = 0.5_dp * (lo + hi)
         if (student_t_cdf(mid, df) < probability) then
            lo = mid
         else
            hi = mid
         end if
      end do
      value = 0.5_dp * (lo + hi)
   end function student_t_quantile_local

   pure subroutine replication_probability(p, p_rep, n, twotailed, d_prime, r_equiv)
      real(dp), intent(in) :: p !! Tail probability transformed to a replication probability.
      real(dp), intent(out) :: p_rep !! Killeen-style replication probability corresponding to `p`.
      integer, intent(in), optional :: n !! Optional sample size used to derive equivalent `d` and `r`; must exceed two.
      logical, intent(in), optional :: twotailed !! When true, doubles `p` before the upstream transformation; default false.
      real(dp), intent(out), optional :: d_prime !! Equivalent standardized mean-difference measure when `n` is supplied.
      real(dp), intent(out), optional :: r_equiv !! Equivalent correlation magnitude when `n` is supplied.
      real(dp) :: probability, df, t_value
      logical :: two

      two = .false.
      if (present(twotailed)) two = twotailed
      probability = p
      if (two) probability = 2.0_dp * probability
      probability = min(1.0_dp, max(0.0_dp, probability))
      p_rep = normal_cdf(normal_quantile(1.0_dp - probability) / sqrt(2.0_dp))

      if (present(d_prime)) d_prime = 0.0_dp
      if (present(r_equiv)) r_equiv = 0.0_dp
      if (.not. present(n)) return
      if (n <= 2) return

      df = real(n - 2, dp)
      t_value = -student_t_quantile_local(0.5_dp * probability, df)
      if (present(d_prime)) d_prime = 2.0_dp * t_value / sqrt(df)
      if (present(r_equiv)) r_equiv = sqrt(t_value * t_value / (t_value * t_value + df))
   end subroutine replication_probability

   pure subroutine replication_probability_f(f_statistic, df2, p_rep, d_prime, probability, r_equiv, twotailed)
      real(dp), intent(in) :: f_statistic !! Observed F statistic with one numerator degree of freedom.
      real(dp), intent(in) :: df2 !! Positive denominator degrees of freedom.
      real(dp), intent(out) :: p_rep !! Replication probability from the upstream `p.rep.f` transformation.
      real(dp), intent(out) :: d_prime !! Standardized effect-size approximation `sqrt(4F/df2)`.
      real(dp), intent(out) :: probability !! Tail probability returned by upstream `p.rep.f`.
      real(dp), intent(out) :: r_equiv !! Equivalent correlation magnitude.
      logical, intent(in), optional :: twotailed !! Apply the upstream two-tailed probability adjustment when true; default false.
      real(dp) :: work_p
      logical :: two

      two = .false.
      if (present(twotailed)) two = twotailed
      if (df2 <= 0.0_dp .or. f_statistic < 0.0_dp) then
         p_rep = 0.0_dp
         d_prime = 0.0_dp
         probability = 1.0_dp
         r_equiv = 0.0_dp
         return
      end if

      work_p = 1.0_dp - f_cdf(f_statistic, 1.0_dp, df2)
      d_prime = sqrt(4.0_dp * f_statistic / df2)
      if (two) work_p = 2.0_dp * work_p
      work_p = min(1.0_dp, max(0.0_dp, work_p))
      p_rep = normal_cdf(normal_quantile(1.0_dp - work_p) / sqrt(2.0_dp))
      probability = work_p
      if (two) probability = 2.0_dp * probability
      r_equiv = sqrt(f_statistic / (f_statistic + df2))
   end subroutine replication_probability_f

   pure subroutine replication_probability_r(r, n, p_rep, d_prime, probability, twotailed)
      real(dp), intent(in) :: r !! Observed correlation, expected strictly between `-1` and `1`.
      integer, intent(in) :: n !! Sample size; upstream approximation requires `n > 4`.
      real(dp), intent(out) :: p_rep !! Replication probability from the upstream `p.rep.r` transformation.
      real(dp), intent(out) :: d_prime !! Correlation converted to the upstream `d-prime` scale.
      real(dp), intent(out) :: probability !! One- or two-tailed normal-approximation probability.
      logical, intent(in), optional :: twotailed !! Double the returned probability when true; default true.
      real(dp) :: sigma_d, z_value, one_tail
      logical :: two

      two = .true.
      if (present(twotailed)) two = twotailed
      if (n <= 4 .or. abs(r) >= 1.0_dp) then
         p_rep = 0.0_dp
         d_prime = 0.0_dp
         probability = 1.0_dp
         return
      end if

      d_prime = 2.0_dp * r / sqrt(1.0_dp - r * r)
      sigma_d = sqrt(4.0_dp / real(n - 4, dp))
      z_value = d_prime / sigma_d
      one_tail = 1.0_dp - normal_cdf(z_value)
      p_rep = normal_cdf(normal_quantile(1.0_dp - one_tail) / sqrt(2.0_dp))
      probability = one_tail
      if (two) probability = 2.0_dp * probability
   end subroutine replication_probability_r

   pure subroutine replication_probability_t(t_statistic, df, p_rep, d_prime, probability, r_equiv, df2, twotailed)
      real(dp), intent(in) :: t_statistic !! Observed Student-t statistic.
      real(dp), intent(in) :: df !! Degrees of freedom for a one-sample/paired test, or first-group `n-1` when `df2` is supplied.
      real(dp), intent(out) :: p_rep !! Replication probability from the upstream `p.rep.t` transformation.
      real(dp), intent(out) :: d_prime !! Standardized effect-size approximation used by upstream `p.rep.t`.
      real(dp), intent(out) :: probability !! One- or two-tailed Student-t probability.
      real(dp), intent(out) :: r_equiv !! Equivalent correlation magnitude.
      real(dp), intent(in), optional :: df2 !! Optional second-group `n-1`; when present, the two-group upstream formula is used.
      logical, intent(in), optional :: twotailed !! Double the tail probability when true; default true.
      real(dp) :: effective_df, n1, n2, nc
      logical :: two

      two = .true.
      if (present(twotailed)) two = twotailed
      effective_df = df
      if (present(df2)) then
         n1 = df + 1.0_dp
         n2 = df2 + 1.0_dp
         effective_df = df + df2
         if (n1 <= 0.0_dp .or. n2 <= 0.0_dp .or. effective_df <= 0.0_dp) then
            p_rep = 0.0_dp
            d_prime = 0.0_dp
            probability = 1.0_dp
            r_equiv = 0.0_dp
            return
         end if
         nc = ((n1 + n2) / 2.0_dp) / ((2.0_dp * n1 * n2) / (n1 + n2))
         d_prime = (2.0_dp * t_statistic / sqrt(effective_df)) * sqrt(nc)
      else
         if (effective_df <= 0.0_dp) then
            p_rep = 0.0_dp
            d_prime = 0.0_dp
            probability = 1.0_dp
            r_equiv = 0.0_dp
            return
         end if
         d_prime = 2.0_dp * t_statistic / sqrt(effective_df)
      end if

      probability = 1.0_dp - student_t_cdf(t_statistic, effective_df)
      r_equiv = sqrt(t_statistic * t_statistic / (t_statistic * t_statistic + effective_df))
      if (two) probability = 2.0_dp * probability
      probability = min(1.0_dp, max(0.0_dp, probability))
      p_rep = normal_cdf(normal_quantile(1.0_dp - probability) / sqrt(2.0_dp))
   end subroutine replication_probability_t

   pure subroutine cohen_profile_similarity(x, similarity, y, midpoint)
      real(dp), intent(in) :: x(:, :) !! First matrix, with profiles in columns and common scale values in rows.
      real(dp), allocatable, intent(out) :: similarity(:, :) !! Profile congruence; rows index `y` and columns index `x`.
      real(dp), intent(in), optional :: y(:, :) !! Optional second profile matrix; defaults to `x`.
      real(dp), intent(in), optional :: midpoint !! Optional scale midpoint; default is half the finite input range.
      real(dp), allocatable :: y_work(:, :)
      real(dp) :: center, min_value, max_value, denom, sx, sy
      integer :: i, j, k
      logical :: found

      if (present(y)) then
         allocate(y_work(size(y, 1), size(y, 2)))
         y_work = y
      else
         allocate(y_work(size(x, 1), size(x, 2)))
         y_work = x
      end if
      allocate(similarity(size(y_work, 2), size(x, 2)))
      similarity = 0.0_dp
      if (size(y_work, 1) /= size(x, 1)) return

      if (present(midpoint)) then
         center = midpoint
      else
         min_value = huge(1.0_dp)
         max_value = -huge(1.0_dp)
         found = .false.
         do j = 1, size(x, 2)
            do i = 1, size(x, 1)
               if (.not. ieee_is_nan(x(i, j))) then
                  min_value = min(min_value, x(i, j))
                  max_value = max(max_value, x(i, j))
                  found = .true.
               end if
            end do
         end do
         do j = 1, size(y_work, 2)
            do i = 1, size(y_work, 1)
               if (.not. ieee_is_nan(y_work(i, j))) then
                  min_value = min(min_value, y_work(i, j))
                  max_value = max(max_value, y_work(i, j))
                  found = .true.
               end if
            end do
         end do
         if (found) then
            center = 0.5_dp * (min_value + max_value)
         else
            center = 0.0_dp
         end if
      end if

      do i = 1, size(x, 2)
         do j = 1, size(y_work, 2)
            sx = 0.0_dp
            sy = 0.0_dp
            denom = 0.0_dp
            do k = 1, size(x, 1)
               if (.not. ieee_is_nan(x(k, i)) .and. .not. ieee_is_nan(y_work(k, j))) then
                  sx = sx + (x(k, i) - center) ** 2
                  sy = sy + (y_work(k, j) - center) ** 2
                  denom = denom + (x(k, i) - center) * (y_work(k, j) - center)
               end if
            end do
            if (sx > 0.0_dp .and. sy > 0.0_dp) then
               similarity(j, i) = denom / sqrt(sx * sy)
            end if
         end do
      end do
   end subroutine cohen_profile_similarity

   subroutine factor_score_indeterminacy(f, r2, info, phi, r, grice, h)
      real(dp), intent(in) :: f(:, :) !! Item-by-factor loading matrix.
      real(dp), allocatable, intent(out) :: r2(:) !! Factor-score determinacy-squared values.
      integer, intent(out) :: info !! Zero on success; negative for invalid dimensions or missing Grice input.
      real(dp), intent(in), optional :: phi(:, :) !! Optional factor-correlation matrix; identity is used when omitted.
      real(dp), intent(in), optional :: r(:, :) !! Optional observed correlation matrix required by the Grice branch.
      logical, intent(in), optional :: grice !! Use the Grice regression-score formula when true; default false.
      real(dp), allocatable, intent(out), optional :: h(:, :) !! Optional leading factor block of the non-Grice projection.
      real(dp), allocatable :: factor_phi(:, :), s(:, :), weights(:, :), pinv(:, :)
      real(dp), allocatable :: communalities(:), lmat(:, :), projection(:, :)
      integer :: i, n_items, n_factors, pinv_info
      logical :: use_grice

      n_items = size(f, 1)
      n_factors = size(f, 2)
      info = 0
      use_grice = .false.
      if (present(grice)) use_grice = grice
      allocate(factor_phi(n_factors, n_factors))
      factor_phi = 0.0_dp
      do i = 1, n_factors
         factor_phi(i, i) = 1.0_dp
      end do
      if (present(phi)) then
         if (size(phi, 1) /= n_factors .or. size(phi, 2) /= n_factors) then
            info = -1
            allocate(r2(0))
            if (present(h)) allocate(h(0, 0))
            return
         end if
         factor_phi = phi
      end if

      if (use_grice) then
         if (.not. present(r)) then
            info = -2
            allocate(r2(0))
            if (present(h)) allocate(h(0, 0))
            return
         end if
         if (size(r, 1) /= n_items .or. size(r, 2) /= n_items) then
            info = -3
            allocate(r2(0))
            if (present(h)) allocate(h(0, 0))
            return
         end if
         allocate(s(n_items, n_factors))
         s = matmul(f, factor_phi)
         call pseudoinverse(r, pinv, info=pinv_info)
         if (pinv_info /= 0) then
            info = pinv_info
            allocate(r2(n_factors))
            r2 = 0.0_dp
            if (present(h)) allocate(h(0, 0))
            return
         end if
         allocate(weights(n_items, n_factors), r2(n_factors))
         weights = matmul(pinv, s)
         do i = 1, n_factors
            r2(i) = min(1.0_dp, dot_product(weights(:, i), s(:, i)))
         end do
         if (present(h)) allocate(h(0, 0))
         return
      end if

      allocate(communalities(n_items))
      communalities = sum(f * f, dim=2)
      communalities = min(1.0_dp, max(0.0_dp, communalities))
      allocate(lmat(n_items, n_factors + n_items))
      lmat = 0.0_dp
      if (present(phi)) then
         lmat(:, 1:n_factors) = matmul(f, factor_phi)
      else
         lmat(:, 1:n_factors) = f
      end if
      do i = 1, n_items
         lmat(i, n_factors + i) = sqrt(max(0.0_dp, 1.0_dp - communalities(i)))
      end do
      call pseudoinverse(lmat, pinv, info=pinv_info)
      if (pinv_info /= 0) then
         info = pinv_info
         allocate(r2(n_factors))
         r2 = 0.0_dp
         if (present(h)) allocate(h(0, 0))
         return
      end if
      projection = matmul(pinv, lmat)
      allocate(r2(n_factors))
      do i = 1, n_factors
         r2(i) = projection(i, i)
      end do
      if (present(h)) then
         allocate(h(n_factors, n_factors))
         h = projection(1:n_factors, 1:n_factors)
      end if
   end subroutine factor_score_indeterminacy

   subroutine cluster_fit_metrics(original, loadings, clusters, cluster_fit, structure_fit, pattern_fit, &
      cluster_rmse, pattern_rmse, pattern_available, info, diagonal)
      real(dp), intent(in) :: original(:, :) !! Observed square correlation or covariance matrix.
      real(dp), intent(in) :: loadings(:, :) !! Structure/loadings matrix with one row per observed variable.
      real(dp), intent(in) :: clusters(:, :) !! Item-by-cluster key matrix conformable with `loadings`.
      real(dp), intent(out) :: cluster_fit !! Fit of the cluster-key constrained model.
      real(dp), intent(out) :: structure_fit !! Fit of the raw loading structure model.
      real(dp), intent(out) :: pattern_fit !! Fit of the pattern model, or zero when the cluster correlation matrix is singular.
      real(dp), intent(out) :: cluster_rmse !! RMSE of the cluster-key constrained residual matrix.
      real(dp), intent(out) :: pattern_rmse !! RMSE of the pattern residual matrix, or zero when unavailable.
      logical, intent(out) :: pattern_available !! True when the cluster correlation matrix is nonsingular enough to invert.
      integer, intent(out) :: info !! Zero on success; negative for invalid dimensions, positive for linear-algebra failures.
      logical, intent(in), optional :: diagonal !! Match upstream diagonal accounting when true; default false.
      real(dp), allocatable :: model(:, :), residual(:, :), covar(:, :), phi(:, :), phi_inv(:, :), pattern(:, :), keyed(:, :)
      real(dp) :: total_original, total_residual, df, det_phi
      integer :: i, j, p, k, det_info, inv_info
      logical :: include_diagonal_flag

      p = size(original, 1)
      k = size(loadings, 2)
      cluster_fit = 0.0_dp
      structure_fit = 0.0_dp
      pattern_fit = 0.0_dp
      cluster_rmse = 0.0_dp
      pattern_rmse = 0.0_dp
      pattern_available = .false.
      info = 0
      if (size(original, 2) /= p .or. size(loadings, 1) /= p .or. size(clusters, 1) /= p .or. size(clusters, 2) /= k) then
         info = -1
         return
      end if
      if (p < 2) then
         info = -2
         return
      end if
      include_diagonal_flag = .false.
      if (present(diagonal)) include_diagonal_flag = diagonal
      df = real(p * (p - 1), dp) / 2.0_dp

      total_original = sum(original * original)
      if (include_diagonal_flag) then
         do i = 1, p
            total_original = total_original - original(i, i) * original(i, i)
         end do
      end if
      if (total_original <= tiny(1.0_dp)) then
         info = -3
         return
      end if

      allocate(model(p, p), residual(p, p))
      model = matmul(loadings, transpose(loadings))
      residual = original - model
      total_residual = sum(residual * residual)
      if (include_diagonal_flag) then
         do i = 1, p
            total_residual = total_residual - residual(i, i) * residual(i, i)
         end do
      end if
      structure_fit = 1.0_dp - total_residual / total_original

      allocate(covar(k, k), phi(k, k))
      covar = matmul(transpose(clusters), matmul(original, clusters))
      phi = 0.0_dp
      do j = 1, k
         do i = 1, k
            if (covar(i, i) > 0.0_dp .and. covar(j, j) > 0.0_dp) then
               phi(i, j) = covar(i, j) / sqrt(covar(i, i) * covar(j, j))
            end if
         end do
      end do
      call determinant(phi, det_phi, det_info)
      if (det_info == 0 .and. abs(det_phi) > 100.0_dp * epsilon(1.0_dp)) then
         call pseudoinverse(phi, phi_inv, info=inv_info)
         if (inv_info == 0) then
            pattern_available = .true.
            allocate(pattern(p, k))
            pattern = matmul(loadings, phi_inv)
            model = matmul(pattern, transpose(loadings))
            residual = original - model
            total_residual = sum(residual * residual)
            if (.not. include_diagonal_flag) then
               do i = 1, p
                  total_residual = total_residual - residual(i, i) * residual(i, i)
               end do
            end if
            pattern_rmse = sqrt(max(0.0_dp, total_residual / (2.0_dp * df)))
            pattern_fit = 1.0_dp - total_residual / total_original
         else
            info = inv_info
         end if
      end if

      allocate(keyed(p, k))
      keyed = loadings * abs(clusters)
      model = matmul(matmul(keyed, phi), transpose(keyed))
      residual = original - model
      total_residual = sum(residual * residual)
      if (include_diagonal_flag) then
         do i = 1, p
            total_residual = total_residual - residual(i, i) * residual(i, i)
         end do
      end if
      cluster_fit = 1.0_dp - total_residual / total_original
      cluster_rmse = sqrt(max(0.0_dp, total_residual / (2.0_dp * df)))

   end subroutine cluster_fit_metrics

   pure subroutine super_correlation_matrix(x, xy, result, y, info)
      real(dp), intent(in) :: x(:, :) !! Score-by-score correlation block for the first score set.
      real(dp), intent(in) :: xy(:, :) !! Item/second-score by first-score correlation block.
      real(dp), allocatable, intent(out) :: result(:, :) !! Combined block matrix matching upstream `superCor` numeric assembly.
      real(dp), intent(in), optional :: y(:, :) !! Optional second score-by-score correlation block.
      integer, intent(out) :: info !! Zero on conformable inputs; negative on dimension mismatch.
      integer :: nx, ny

      nx = size(x, 2)
      ny = size(xy, 1)
      info = 0
      if (size(x, 1) /= nx .or. size(xy, 2) /= nx) then
         info = -1
         allocate(result(0, 0))
         return
      end if
      if (present(y)) then
         if (size(y, 1) /= ny .or. size(y, 2) /= ny) then
            info = -2
            allocate(result(0, 0))
            return
         end if
         allocate(result(nx + ny, nx + ny))
         result = 0.0_dp
         result(1:nx, 1:nx) = x
         result(nx + 1:nx + ny, 1:nx) = xy
         result(1:nx, nx + 1:nx + ny) = transpose(xy)
         result(nx + 1:nx + ny, nx + 1:nx + ny) = y
      else
         allocate(result(nx + ny, nx))
         result(1:nx, 1:nx) = x
         result(nx + 1:nx + ny, 1:nx) = xy
      end if
   end subroutine super_correlation_matrix

   pure subroutine cor2_matrix(x, result, y, covariance_mode)
      real(dp), intent(in) :: x(:, :) !! First observation-by-variable numeric matrix.
      real(dp), allocatable, intent(out) :: result(:, :) !! Pearson correlation or covariance matrix/block.
      real(dp), intent(in), optional :: y(:, :) !! Optional second observation-by-variable matrix for cross-correlations.
      logical, intent(in), optional :: covariance_mode !! If true, return covariance rather than Pearson correlation.
      integer :: i, j
      logical :: use_covariance

      use_covariance = .false.
      if (present(covariance_mode)) use_covariance = covariance_mode
      if (.not. present(y)) then
         allocate(result(size(x, 2), size(x, 2)))
         if (use_covariance) then
            call covariance_matrix(x, result)
         else
            call correlation_matrix(x, result)
         end if
         return
      end if

      allocate(result(size(x, 2), size(y, 2)))
      if (size(x, 1) /= size(y, 1)) then
         result = 0.0_dp
         return
      end if
      do j = 1, size(y, 2)
         do i = 1, size(x, 2)
            if (use_covariance) then
               result(i, j) = finite_covariance(x(:, i), y(:, j))
            else
               result(i, j) = finite_correlation(x(:, i), y(:, j))
            end if
         end do
      end do
   end subroutine cor2_matrix

   subroutine interbattery_loadings(r12, n_factors, a1, a2, info)
      real(dp), intent(in) :: r12(:, :) !! Cross-correlation matrix between the X and Y variable sets.
      integer, intent(in) :: n_factors !! Number of leading interbattery dimensions to retain.
      real(dp), allocatable, intent(out) :: a1(:, :) !! Interbattery loading matrix for the X-variable set.
      real(dp), allocatable, intent(out) :: a2(:, :) !! Interbattery loading matrix for the Y-variable set.
      integer, intent(out) :: info !! Zero on successful eigendecomposition; negative for invalid dimensions.
      real(dp), allocatable :: h1(:, :), values(:), vectors(:, :), w1(:, :), w2(:, :), gamma(:)
      real(dp) :: scale_factor, sign_sum
      integer :: i, j, k

      k = min(max(1, n_factors), min(size(r12, 1), size(r12, 2)))
      if (size(r12, 1) < 1 .or. size(r12, 2) < 1) then
         info = -1
         allocate(a1(0, 0), a2(0, 0))
         return
      end if
      allocate(h1(size(r12, 1), size(r12, 1)))
      h1 = matmul(r12, transpose(r12))
      call symmetric_eigen(h1, values, vectors, info, descending=.true.)
      if (info /= 0) then
         allocate(a1(0, 0), a2(0, 0))
         return
      end if
      allocate(w1(size(r12, 1), k), w2(size(r12, 2), k), gamma(k))
      allocate(a1(size(r12, 1), k), a2(size(r12, 2), k))
      w1 = vectors(:, 1:k)
      gamma = sqrt(max(0.0_dp, values(1:k)))
      w2 = 0.0_dp
      a1 = 0.0_dp
      a2 = 0.0_dp
      do j = 1, k
         if (gamma(j) > sqrt(epsilon(1.0_dp))) then
            scale_factor = sqrt(gamma(j))
            a1(:, j) = w1(:, j) * scale_factor
            w2(:, j) = matmul(transpose(r12), w1(:, j)) / gamma(j)
            a2(:, j) = w2(:, j) * scale_factor
         end if
         sign_sum = 0.0_dp
         do i = 1, size(a1, 1)
            if (a1(i, j) > 0.0_dp) sign_sum = sign_sum + 1.0_dp
            if (a1(i, j) < 0.0_dp) sign_sum = sign_sum - 1.0_dp
         end do
         if (sign_sum < 0.0_dp) a1(:, j) = -a1(:, j)
         sign_sum = 0.0_dp
         do i = 1, size(a2, 1)
            if (a2(i, j) > 0.0_dp) sign_sum = sign_sum + 1.0_dp
            if (a2(i, j) < 0.0_dp) sign_sum = sign_sum - 1.0_dp
         end do
         if (sign_sum < 0.0_dp) a2(:, j) = -a2(:, j)
      end do
   end subroutine interbattery_loadings

   subroutine block_random_design(n, blocks, info, conditions)
      integer, intent(in) :: n !! Number of subjects; must be a multiple of the product of condition counts.
      integer, allocatable, intent(out) :: blocks(:, :) !! Assignments: block id followed by one column per independent variable.
      integer, intent(out) :: info !! Zero on success; negative for invalid condition counts or incompatible subject count.
      integer, intent(in), optional :: conditions(:) !! Levels per factor; default is one two-level factor.
      integer, allocatable :: cond(:), permutation(:)
      real(dp) :: u
      integer :: block, i, j, total_conditions, n_blocks, tmp, previous_product, cumulative_product, code

      if (present(conditions)) then
         allocate(cond(size(conditions)))
         cond = conditions
      else
         allocate(cond(1))
         cond = [2]
      end if
      info = 0
      if (size(cond) == 0 .or. any(cond <= 0)) then
         info = -1
         allocate(blocks(0, 0))
         return
      end if
      total_conditions = product(cond)
      if (n <= 0 .or. mod(n, total_conditions) /= 0) then
         info = -2
         allocate(blocks(0, 0))
         return
      end if
      n_blocks = n / total_conditions
      allocate(blocks(n, size(cond) + 1), permutation(total_conditions))
      do block = 1, n_blocks
         do i = 1, total_conditions
            permutation(i) = i
         end do
         do i = total_conditions, 2, -1
            call random_number(u)
            j = 1 + int(u * real(i, dp))
            j = min(i, max(1, j))
            tmp = permutation(i)
            permutation(i) = permutation(j)
            permutation(j) = tmp
         end do
         do i = 1, total_conditions
            blocks((block - 1) * total_conditions + i, 1) = block
            code = permutation(i)
            previous_product = 1
            do j = 1, size(cond)
               cumulative_product = previous_product * cond(j)
               blocks((block - 1) * total_conditions + i, j + 1) = &
                  (mod(code, cumulative_product) + previous_product) / previous_product
               previous_product = cumulative_product
            end do
         end do
      end do
   end subroutine block_random_design

end module psych_additional
