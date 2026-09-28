module psych_scoring_extended
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
   use r_kinds, only : dp
   use psych_types, only : describe_result
   use psych_stats, only : covariance_matrix, correlation_matrix, finite_mean
   use psych_stats, only : finite_sd, finite_median
   use psych_correlations, only : squared_multiple_correlations, pseudoinverse, correct_correlation
   use psych_descriptive, only : describe_vector
   use psych_scoring, only : score_key_matrix
   implicit none
   private
   public :: cluster_assignments_to_keys, cluster_correlation_summary, cluster_loadings_summary
   public :: score_items_summary, score_overlap_summary, grouped_describe_matrix
   public :: pairwise_describe_counts, pairwise_impute_by_keys, pairwise_count_blocked
   public :: item_sort_indices, key_sort_columns

contains

   pure subroutine cluster_assignments_to_keys(clusters, keys)
      integer, intent(in) :: clusters(:) !! Positive one-based cluster assignment for each item.
      real(dp), allocatable, intent(out) :: keys(:, :) !! Item-by-cluster zero/one key matrix.
      integer :: i, n_clusters

      if (size(clusters) == 0) then
         allocate(keys(0, 0))
         return
      end if
      n_clusters = max(0, maxval(clusters))
      allocate(keys(size(clusters), n_clusters))
      keys = 0.0_dp
      do i = 1, size(clusters)
         if (clusters(i) >= 1 .and. clusters(i) <= n_clusters) keys(i, clusters(i)) = 1.0_dp
      end do
   end subroutine cluster_assignments_to_keys

   subroutine cluster_correlation_summary(keys, r, cluster_cor, alpha, average_r, g6, scale_sd, corrected, info)
      real(dp), intent(in) :: keys(:, :) !! Item-by-scale signed key matrix.
      real(dp), intent(in) :: r(:, :) !! Item correlation matrix.
      real(dp), allocatable, intent(out) :: cluster_cor(:, :) !! Raw correlations among keyed scale sums.
      real(dp), allocatable, intent(out) :: alpha(:) !! Standardized alpha estimate for each keyed scale.
      real(dp), allocatable, intent(out) :: average_r(:) !! Average inter-item correlation implied by each alpha.
      real(dp), allocatable, intent(out) :: g6(:) !! Guttman lambda-6 reliability estimate for each scale.
      real(dp), allocatable, intent(out) :: scale_sd(:) !! Standard deviation of each keyed scale sum.
      real(dp), allocatable, intent(out) :: corrected(:, :) !! Scale correlations corrected for attenuation.
      integer, intent(out) :: info !! Zero when SMC calculation succeeds; nonzero otherwise.
      real(dp), allocatable :: covar(:, :), smc(:), sizes(:), smc_sum(:)
      real(dp) :: denom
      integer :: i, j, k

      k = size(keys, 2)
      covar = matmul(transpose(keys), matmul(r, keys))
      allocate(cluster_cor(k, k), alpha(k), average_r(k), g6(k), scale_sd(k), corrected(k, k))
      allocate(sizes(k), smc(k), smc_sum(k))
      do i = 1, k
         scale_sd(i) = sqrt(max(0.0_dp, covar(i, i)))
         sizes(i) = sum(abs(keys(:, i)) * abs(keys(:, i)))
      end do
      cluster_cor = 0.0_dp
      do j = 1, k
         do i = 1, k
            denom = scale_sd(i) * scale_sd(j)
            if (denom > tiny(1.0_dp)) cluster_cor(i, j) = covar(i, j) / denom
         end do
      end do
      call squared_multiple_correlations(r, smc, info)
      if (info /= 0) smc = 1.0_dp
      smc_sum = matmul(transpose(keys), smc)
      do i = 1, k
         if (sizes(i) > 1.0_dp .and. covar(i, i) > tiny(1.0_dp)) then
            alpha(i) = ((covar(i, i) - sizes(i)) / covar(i, i)) * sizes(i) / (sizes(i) - 1.0_dp)
         else
            alpha(i) = 1.0_dp
         end if
         average_r(i) = alpha(i) / max(tiny(1.0_dp), sizes(i) - alpha(i) * (sizes(i) - 1.0_dp))
         if (covar(i, i) > tiny(1.0_dp)) then
            g6(i) = (covar(i, i) - sizes(i) + smc_sum(i)) / covar(i, i)
         else
            g6(i) = 0.0_dp
         end if
      end do
      call correct_correlation(cluster_cor, alpha, corrected)
   end subroutine cluster_correlation_summary

   subroutine cluster_loadings_summary(keys, r, loadings, pattern, cluster_cor, corrected, alpha, average_r, g6, info)
      real(dp), intent(in) :: keys(:, :) !! Item-by-cluster signed key matrix.
      real(dp), intent(in) :: r(:, :) !! Item correlation matrix.
      real(dp), allocatable, intent(out) :: loadings(:, :) !! Item correlations with keyed scale composites.
      real(dp), allocatable, intent(out) :: pattern(:, :) !! Item pattern coefficients after scale-correlation adjustment.
      real(dp), allocatable, intent(out) :: cluster_cor(:, :) !! Raw correlations among keyed scales.
      real(dp), allocatable, intent(out) :: corrected(:, :) !! Attenuation-corrected scale correlation matrix.
      real(dp), allocatable, intent(out) :: alpha(:) !! Alpha reliability for each scale.
      real(dp), allocatable, intent(out) :: average_r(:) !! Average inter-item correlation for each scale.
      real(dp), allocatable, intent(out) :: g6(:) !! Lambda-6 reliability for each scale.
      integer, intent(out) :: info !! Zero when all inversions and SMC computations succeed.
      real(dp), allocatable :: covar(:, :), item_covar(:, :), sd(:), inv_cluster(:, :), smc(:), sizes(:), corrected_r(:, :)
      real(dp) :: denom
      integer :: i, j, k, smc_info

      call cluster_correlation_summary(keys, r, cluster_cor, alpha, average_r, g6, sd, corrected, info)
      k = size(keys, 2)
      covar = matmul(transpose(keys), matmul(r, keys))
      item_covar = matmul(r, keys)
      allocate(loadings(size(r, 1), k), pattern(size(r, 1), k), smc(size(r, 1)), sizes(k), corrected_r(size(r, 1), size(r, 2)))
      corrected_r = r
      call squared_multiple_correlations(r, smc, smc_info)
      if (smc_info /= 0) smc = 1.0_dp
      do i = 1, size(r, 1)
         corrected_r(i, i) = max(tiny(1.0_dp), smc(i))
      end do
      item_covar = matmul(corrected_r, keys)
      do j = 1, k
         sizes(j) = sum(abs(keys(:, j)) * abs(keys(:, j)))
         do i = 1, size(r, 1)
            denom = sqrt(max(tiny(1.0_dp), corrected_r(i, i) * covar(j, j)))
            loadings(i, j) = item_covar(i, j) / denom
         end do
      end do
      call pseudoinverse(cluster_cor, inv_cluster, info=info)
      if (info == 0) then
         pattern = matmul(loadings, inv_cluster)
      else
         pattern = loadings
      end if
   end subroutine cluster_loadings_summary

   subroutine score_items_summary(keys, items, scores, alpha, average_r, g6, scale_cor, item_cor, info, totals, impute_mode)
      real(dp), intent(in) :: keys(:, :) !! Item-by-scale signed key matrix.
      real(dp), intent(in) :: items(:, :) !! Subject-by-item response matrix with IEEE NaNs allowed.
      real(dp), allocatable, intent(out) :: scores(:, :) !! Subject-by-scale mean or total scores.
      real(dp), allocatable, intent(out) :: alpha(:) !! Reliability alpha for each scale.
      real(dp), allocatable, intent(out) :: average_r(:) !! Average inter-item correlation for each scale.
      real(dp), allocatable, intent(out) :: g6(:) !! Lambda-6 reliability for each scale.
      real(dp), allocatable, intent(out) :: scale_cor(:, :) !! Correlations among computed scale scores.
      real(dp), allocatable, intent(out) :: item_cor(:, :) !! Item-by-scale correlations.
      integer, intent(out) :: info !! Zero when the reliability calculations succeed.
      logical, intent(in), optional :: totals !! Return total rather than mean scale scores when true.
      integer, intent(in), optional :: impute_mode !! Missing handling passed to `score_key_matrix`: 0 none, 1 mean, 2 median.
      real(dp), allocatable :: item_cov(:, :), score_cov(:, :), score_sd(:), item_sd(:), smc(:), sizes(:)
      real(dp), allocatable :: minima(:), maxima(:), smc_sum(:)
      real(dp) :: denom
      integer :: i, j, k, mode, smc_info
      logical :: use_totals

      k = size(keys, 2)
      mode = 2
      if (present(impute_mode)) mode = impute_mode
      use_totals = .false.
      if (present(totals)) use_totals = totals
      allocate(minima(size(items, 2)), maxima(size(items, 2)))
      do j = 1, size(items, 2)
         minima(j) = finite_min(items(:, j))
         maxima(j) = finite_max(items(:, j))
      end do
      call score_key_matrix(keys, items, scores, totals=use_totals, impute_mode=mode, minima=minima, maxima=maxima)
      allocate(item_cov(size(items, 2), size(items, 2)), score_cov(k, k))
      call covariance_matrix(items, item_cov)
      call covariance_matrix(scores, score_cov)
      allocate(scale_cor(k, k), item_cor(size(items, 2), k), alpha(k), average_r(k), g6(k))
      allocate(score_sd(k), item_sd(size(items, 2)), smc(size(items, 2)), sizes(k), smc_sum(k))
      do j = 1, k
         score_sd(j) = sqrt(max(0.0_dp, score_cov(j, j)))
      end do
      do i = 1, size(items, 2)
         item_sd(i) = sqrt(max(0.0_dp, item_cov(i, i)))
      end do
      do j = 1, k
         do i = 1, k
            denom = score_sd(i) * score_sd(j)
            if (denom > tiny(1.0_dp)) then
               scale_cor(i, j) = score_cov(i, j) / denom
            else
               scale_cor(i, j) = 0.0_dp
            end if
         end do
      end do
      do j = 1, k
         do i = 1, size(items, 2)
            denom = item_sd(i) * score_sd(j)
            if (denom > tiny(1.0_dp)) then
               item_cor(i, j) = finite_covariance_local(items(:, i), scores(:, j)) / denom
            else
               item_cor(i, j) = 0.0_dp
            end if
         end do
      end do
      call squared_multiple_correlations(covariance_to_correlation_local(item_cov), smc, smc_info)
      info = smc_info
      if (smc_info /= 0) smc = 1.0_dp
      smc_sum = matmul(transpose(abs(keys)), smc)
      do j = 1, k
         sizes(j) = sum(abs(keys(:, j)) * abs(keys(:, j)))
         if (sizes(j) > 1.0_dp .and. score_cov(j, j) > tiny(1.0_dp)) then
            alpha(j) = (score_cov(j, j) - dot_product(diagonal_local(item_cov), abs(keys(:, j)))) * &
               sizes(j) / ((sizes(j) - 1.0_dp) * score_cov(j, j))
         else
            alpha(j) = 1.0_dp
         end if
         average_r(j) = alpha(j) / max(tiny(1.0_dp), sizes(j) - alpha(j) * (sizes(j) - 1.0_dp))
         if (score_cov(j, j) > tiny(1.0_dp)) then
            g6(j) = (score_cov(j, j) - dot_product(diagonal_local(item_cov), abs(keys(:, j))) + smc_sum(j)) / score_cov(j, j)
         else
            g6(j) = 0.0_dp
         end if
      end do
   end subroutine score_items_summary

   subroutine score_overlap_summary(keys, r, adjusted_cor, alpha, average_r, g6, item_cor, info)
      real(dp), intent(in) :: keys(:, :) !! Item-by-scale signed key matrix with possible overlap.
      real(dp), intent(in) :: r(:, :) !! Item correlation matrix.
      real(dp), allocatable, intent(out) :: adjusted_cor(:, :) !! Scale correlations after correcting overlapping items.
      real(dp), allocatable, intent(out) :: alpha(:) !! Scale alpha reliability estimates.
      real(dp), allocatable, intent(out) :: average_r(:) !! Scale average inter-item correlations.
      real(dp), allocatable, intent(out) :: g6(:) !! Lambda-6 reliability estimates.
      real(dp), allocatable, intent(out) :: item_cor(:, :) !! Corrected item-by-scale correlations.
      integer, intent(out) :: info !! Zero when SMC and matrix operations succeed.
      real(dp), allocatable :: raw_cov(:, :), adj_cov(:, :), sizes(:), sd(:), smc(:), smc_sum(:), rwork(:, :)
      real(dp) :: overlap, replacement, denom
      integer :: i, j, k, smc_info

      k = size(keys, 2)
      raw_cov = matmul(transpose(keys), matmul(r, keys))
      adj_cov = raw_cov
      allocate(adjusted_cor(k, k), alpha(k), average_r(k), g6(k), sizes(k), sd(k))
      allocate(smc(size(r, 1)), smc_sum(k), rwork(size(r, 1), size(r, 2)), item_cor(size(r, 1), k))
      call squared_multiple_correlations(r, smc, smc_info)
      info = smc_info
      if (smc_info /= 0) smc = 1.0_dp
      do i = 1, k
         sizes(i) = sum(abs(keys(:, i)) * abs(keys(:, i)))
         sd(i) = sqrt(max(0.0_dp, raw_cov(i, i)))
         if (sizes(i) > 1.0_dp .and. raw_cov(i, i) > tiny(1.0_dp)) then
            alpha(i) = ((raw_cov(i, i) - sizes(i)) / raw_cov(i, i)) * sizes(i) / (sizes(i) - 1.0_dp)
         else
            alpha(i) = 1.0_dp
         end if
         average_r(i) = alpha(i) / max(tiny(1.0_dp), sizes(i) - alpha(i) * (sizes(i) - 1.0_dp))
      end do
      do j = 1, k
         do i = 1, j
            overlap = sum(keys(:, i) * keys(:, j))
            replacement = sum(keys(:, i) * keys(:, j)) * sqrt(max(0.0_dp, average_r(i) * average_r(j)))
            adj_cov(i, j) = raw_cov(i, j) - overlap + replacement
            adj_cov(j, i) = adj_cov(i, j)
         end do
      end do
      do i = 1, k
         adj_cov(i, i) = raw_cov(i, i)
      end do
      do j = 1, k
         do i = 1, k
            denom = sd(i) * sd(j)
            if (denom > tiny(1.0_dp)) adjusted_cor(i, j) = adj_cov(i, j) / denom
         end do
      end do
      smc_sum = matmul(transpose(abs(keys)), smc)
      do i = 1, k
         if (raw_cov(i, i) > tiny(1.0_dp)) then
            g6(i) = (raw_cov(i, i) - sizes(i) + smc_sum(i)) / raw_cov(i, i)
         else
            g6(i) = 0.0_dp
         end if
      end do
      rwork = r
      do i = 1, size(r, 1)
         rwork(i, i) = smc(i)
      end do
      item_cor = matmul(rwork, keys)
      do j = 1, k
         denom = sqrt(max(tiny(1.0_dp), g6(j) * raw_cov(j, j)))
         item_cor(:, j) = item_cor(:, j) / denom
      end do
   end subroutine score_overlap_summary

   pure subroutine grouped_describe_matrix(x, groups, group_values, means, sds, counts)
      real(dp), intent(in) :: x(:, :) !! Observation-by-variable data matrix with IEEE NaNs allowed.
      integer, intent(in) :: groups(:) !! Group label for each observation.
      integer, allocatable, intent(out) :: group_values(:) !! Distinct group labels in first-observed order.
      real(dp), allocatable, intent(out) :: means(:, :) !! Group-by-variable finite means.
      real(dp), allocatable, intent(out) :: sds(:, :) !! Group-by-variable sample standard deviations.
      integer, allocatable, intent(out) :: counts(:, :) !! Group-by-variable finite counts.
      real(dp), allocatable :: values(:)
      integer, allocatable :: unique_work(:)
      integer :: g, i, j, n_use, n_group
      logical :: seen

      allocate(unique_work(size(groups)))
      n_group = 0
      do i = 1, size(groups)
         seen = .false.
         do g = 1, n_group
            if (unique_work(g) == groups(i)) then
               seen = .true.
               exit
            end if
         end do
         if (.not. seen) then
            n_group = n_group + 1
            unique_work(n_group) = groups(i)
         end if
      end do
      allocate(group_values(n_group))
      if (n_group > 0) group_values = unique_work(1:n_group)
      allocate(means(n_group, size(x, 2)), sds(n_group, size(x, 2)))
      allocate(counts(n_group, size(x, 2)))
      means = 0.0_dp
      sds = 0.0_dp
      counts = 0
      do g = 1, size(group_values)
         do j = 1, size(x, 2)
            allocate(values(size(x, 1)))
            n_use = 0
            do i = 1, min(size(x, 1), size(groups))
               if (groups(i) /= group_values(g)) cycle
               if (ieee_is_nan(x(i, j))) cycle
               n_use = n_use + 1
               values(n_use) = x(i, j)
            end do
            counts(g, j) = n_use
            if (n_use > 0) means(g, j) = finite_mean(values(1:n_use))
            if (n_use > 1) sds(g, j) = finite_sd(values(1:n_use))
            deallocate(values)
         end do
      end do
   end subroutine grouped_describe_matrix

   pure subroutine pairwise_describe_counts(x, diagonal, summary)
      real(dp), intent(in) :: x(:, :) !! Observation-by-variable data matrix with IEEE NaNs indicating missingness.
      logical, intent(in) :: diagonal !! Include diagonal pair counts in the descriptive summary when true.
      type(describe_result), intent(out) :: summary !! Descriptive statistics of the selected pairwise counts.
      real(dp), allocatable :: counts(:, :), values(:)
      integer :: i, j, k, p

      p = size(x, 2)
      allocate(counts(p, p))
      call pairwise_count_local(x, counts)
      allocate(values(p * p))
      k = 0
      do j = 1, p
         do i = 1, p
            if (.not. diagonal .and. i == j) cycle
            k = k + 1
            values(k) = counts(i, j)
         end do
      end do
      call describe_vector(values(1:k), summary)
   end subroutine pairwise_describe_counts

   pure subroutine pairwise_impute_by_keys(keys, r, average_r, finite_counts, finite_fraction, imputed, fix)
      real(dp), intent(in) :: keys(:, :) !! Item-by-scale key matrix; positive entries select an item for a scale.
      real(dp), intent(in) :: r(:, :) !! Correlation matrix with IEEE NaNs allowed off diagonal.
      real(dp), allocatable, intent(out) :: average_r(:, :) !! Mean available item correlation for each pair of scales.
      integer, allocatable, intent(out) :: finite_counts(:, :) !! Number of finite selected correlations per scale pair.
      real(dp), allocatable, intent(out) :: finite_fraction(:, :) !! Fraction of selected correlations that are finite.
      real(dp), allocatable, intent(out) :: imputed(:, :) !! Correlation matrix with selected missing cells optionally imputed.
      logical, intent(in), optional :: fix !! Replace selected missing correlations by scale-pair means when true.
      real(dp) :: total, count_total
      integer :: i, j, a, b, k, n_selected
      logical :: do_fix

      k = size(keys, 2)
      allocate(average_r(k, k), finite_counts(k, k), finite_fraction(k, k), imputed(size(r, 1), size(r, 2)))
      average_r = 0.0_dp
      finite_counts = 0
      finite_fraction = 0.0_dp
      imputed = r
      do j = 1, k
         do i = 1, k
            total = 0.0_dp
            n_selected = 0
            do b = 1, size(r, 2)
               if (keys(b, j) <= 0.0_dp) cycle
               do a = 1, size(r, 1)
                  if (keys(a, i) <= 0.0_dp) cycle
                  n_selected = n_selected + 1
                  if (.not. ieee_is_nan(r(a, b))) then
                     finite_counts(i, j) = finite_counts(i, j) + 1
                     total = total + r(a, b)
                  end if
               end do
            end do
            if (finite_counts(i, j) > 0) average_r(i, j) = total / real(finite_counts(i, j), dp)
            count_total = real(max(1, n_selected), dp)
            finite_fraction(i, j) = real(finite_counts(i, j), dp) / count_total
         end do
      end do
      do_fix = .false.
      if (present(fix)) do_fix = fix
      if (do_fix) then
         do j = 1, k
            do i = 1, k
               do b = 1, size(r, 2)
                  if (keys(b, j) <= 0.0_dp) cycle
                  do a = 1, size(r, 1)
                     if (keys(a, i) <= 0.0_dp) cycle
                     if (ieee_is_nan(imputed(a, b))) imputed(a, b) = average_r(i, j)
                  end do
               end do
            end do
         end do
         do i = 1, min(size(imputed, 1), size(imputed, 2))
            imputed(i, i) = 1.0_dp
         end do
      end if
   end subroutine pairwise_impute_by_keys

   pure subroutine pairwise_count_blocked(x, block_size, counts)
      real(dp), intent(in) :: x(:, :) !! Observation-by-variable matrix with IEEE NaNs indicating missing values.
      integer, intent(in) :: block_size !! Number of variables processed per output block; values below one use all variables.
      integer, allocatable, intent(out) :: counts(:, :) !! Pairwise finite-observation counts.
      integer :: a, b, i, j, n, p, block

      n = size(x, 1)
      p = size(x, 2)
      block = max(1, block_size)
      allocate(counts(p, p))
      counts = 0
      do b = 1, p, block
         do a = 1, p, block
            do j = b, min(p, b + block - 1)
               do i = a, min(p, a + block - 1)
                  counts(i, j) = 0
                  do n = 1, size(x, 1)
                     if (.not. ieee_is_nan(x(n, i)) .and. .not. ieee_is_nan(x(n, j))) counts(i, j) = counts(i, j) + 1
                  end do
               end do
            end do
         end do
      end do
   end subroutine pairwise_count_blocked

   pure subroutine item_sort_indices(values, ascending, order)
      real(dp), intent(in) :: values(:) !! Item summary values used as the sorting key.
      logical, intent(in) :: ascending !! Sort from smallest to largest when true.
      integer, allocatable, intent(out) :: order(:) !! One-based permutation of items.
      integer :: i, j, idx, tmp

      allocate(order(size(values)))
      do i = 1, size(values)
         order(i) = i
      end do
      do i = 1, size(values) - 1
         idx = i
         do j = i + 1, size(values)
            if (ascending) then
               if (values(order(j)) < values(order(idx))) idx = j
            else
               if (values(order(j)) > values(order(idx))) idx = j
            end if
         end do
         if (idx /= i) then
            tmp = order(i)
            order(i) = order(idx)
            order(idx) = tmp
         end if
      end do
   end subroutine item_sort_indices

   pure subroutine key_sort_columns(keys, sorted, order)
      real(dp), intent(in) :: keys(:, :) !! Item-by-scale key matrix.
      real(dp), allocatable, intent(out) :: sorted(:, :) !! Key matrix with columns ordered by weighted item positions.
      integer, allocatable, intent(out) :: order(:) !! Original one-based column indices in sorted order.
      real(dp), allocatable :: weights(:)
      integer :: i, j, idx, tmp

      allocate(weights(size(keys, 2)), order(size(keys, 2)), sorted(size(keys, 1), size(keys, 2)))
      do j = 1, size(keys, 2)
         weights(j) = 0.0_dp
         do i = 1, size(keys, 1)
            weights(j) = weights(j) + real(i, dp) * abs(keys(i, j))
         end do
         order(j) = j
      end do
      do i = 1, size(order) - 1
         idx = i
         do j = i + 1, size(order)
            if (weights(order(j)) < weights(order(idx))) idx = j
         end do
         if (idx /= i) then
            tmp = order(i)
            order(i) = order(idx)
            order(idx) = tmp
         end if
      end do
      do j = 1, size(order)
         sorted(:, j) = keys(:, order(j))
      end do
   end subroutine key_sort_columns

   pure real(dp) function finite_min(x) result(value)
      real(dp), intent(in) :: x(:) !! Values whose finite minimum is required.
      integer :: i
      value = huge(1.0_dp)
      do i = 1, size(x)
         if (.not. ieee_is_nan(x(i))) value = min(value, x(i))
      end do
      if (value > 0.5_dp * huge(1.0_dp)) value = 0.0_dp
   end function finite_min

   pure real(dp) function finite_max(x) result(value)
      real(dp), intent(in) :: x(:) !! Values whose finite maximum is required.
      integer :: i
      value = -huge(1.0_dp)
      do i = 1, size(x)
         if (.not. ieee_is_nan(x(i))) value = max(value, x(i))
      end do
      if (value < -0.5_dp * huge(1.0_dp)) value = 0.0_dp
   end function finite_max

   pure real(dp) function finite_covariance_local(x, y) result(value)
      real(dp), intent(in) :: x(:) !! First variable with IEEE NaNs allowed.
      real(dp), intent(in) :: y(:) !! Second variable with IEEE NaNs allowed.
      real(dp) :: mx, my
      integer :: i, n

      mx = 0.0_dp
      my = 0.0_dp
      n = 0
      do i = 1, min(size(x), size(y))
         if (ieee_is_nan(x(i)) .or. ieee_is_nan(y(i))) cycle
         n = n + 1
         mx = mx + x(i)
         my = my + y(i)
      end do
      if (n <= 1) then
         value = 0.0_dp
         return
      end if
      mx = mx / real(n, dp)
      my = my / real(n, dp)
      value = 0.0_dp
      do i = 1, min(size(x), size(y))
         if (ieee_is_nan(x(i)) .or. ieee_is_nan(y(i))) cycle
         value = value + (x(i) - mx) * (y(i) - my)
      end do
      value = value / real(n - 1, dp)
   end function finite_covariance_local

   pure function covariance_to_correlation_local(cov) result(cor)
      real(dp), intent(in) :: cov(:, :) !! Covariance matrix.
      real(dp), allocatable :: cor(:, :) !! Correlation matrix standardized from `cov`.
      real(dp), allocatable :: sd(:)
      integer :: i, j

      allocate(cor(size(cov, 1), size(cov, 2)), sd(size(cov, 1)))
      do i = 1, size(cov, 1)
         sd(i) = sqrt(max(0.0_dp, cov(i, i)))
      end do
      cor = 0.0_dp
      do j = 1, size(cov, 2)
         do i = 1, size(cov, 1)
            if (sd(i) * sd(j) > tiny(1.0_dp)) cor(i, j) = cov(i, j) / (sd(i) * sd(j))
         end do
      end do
   end function covariance_to_correlation_local

   pure function diagonal_local(a) result(d)
      real(dp), intent(in) :: a(:, :) !! Matrix whose diagonal is returned.
      real(dp), allocatable :: d(:) !! Diagonal vector.
      integer :: i
      allocate(d(min(size(a, 1), size(a, 2))))
      do i = 1, size(d)
         d(i) = a(i, i)
      end do
   end function diagonal_local

   pure subroutine pairwise_count_local(x, counts)
      real(dp), intent(in) :: x(:, :) !! Observation-by-variable matrix with IEEE NaNs indicating missingness.
      real(dp), intent(out) :: counts(:, :) !! Pairwise finite-observation counts as reals.
      integer :: i, j, row

      counts = 0.0_dp
      do j = 1, size(x, 2)
         do i = 1, size(x, 2)
            do row = 1, size(x, 1)
               if (.not. ieee_is_nan(x(row, i)) .and. .not. ieee_is_nan(x(row, j))) counts(i, j) = counts(i, j) + 1.0_dp
            end do
         end do
      end do
   end subroutine pairwise_count_local

end module psych_scoring_extended
