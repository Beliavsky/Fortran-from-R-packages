module psych_iclust_extended
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use r_kinds, only : dp
   implicit none
   private
   public :: iclust_hierarchical, iclust_sort_loadings

contains

   pure subroutine iclust_hierarchical(r, n_clusters, clusters, merge_left, merge_right, merge_sign, &
      merge_similarity, merge_alpha, merge_beta, merge_size, n_merges, info, reverse, beta_min, &
      alpha_rule, beta_rule, alpha_size, beta_size, correct_similarity)
      real(dp), intent(in) :: r(:, :) !! Item correlation matrix used to form homogeneous hierarchical composites.
      integer, intent(in) :: n_clusters !! Requested final cluster count; zero applies reliability stopping criteria.
      real(dp), allocatable, intent(out) :: clusters(:, :) !! Final signed item-by-cluster key matrix.
      integer, allocatable, intent(out) :: merge_left(:) !! Left active-cluster slot for each accepted merge.
      integer, allocatable, intent(out) :: merge_right(:) !! Right active-cluster slot for each accepted merge.
      integer, allocatable, intent(out) :: merge_sign(:) !! Sign applied to the right cluster at each merge, plus or minus one.
      real(dp), allocatable, intent(out) :: merge_similarity(:) !! Reliability-corrected similarity selected at each merge.
      real(dp), allocatable, intent(out) :: merge_alpha(:) !! Coefficient alpha of each newly formed cluster.
      real(dp), allocatable, intent(out) :: merge_beta(:) !! Weighted beta of each newly formed cluster.
      integer, allocatable, intent(out) :: merge_size(:) !! Number of items in each newly formed cluster.
      integer, intent(out) :: n_merges !! Number of accepted hierarchical merges.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid dimensions or missing correlations.
      logical, intent(in), optional :: reverse !! Permit sign reversal when the largest absolute similarity is negative;
                                             !! default true.
      real(dp), intent(in), optional :: beta_min !! Minimum accepted beta when no fixed cluster count is requested; default 0.5.
      integer, intent(in), optional :: alpha_rule !! Alpha rule: 0 none, 1 minimum, 2 mean, 3 maximum; default 3.
      integer, intent(in), optional :: beta_rule !! Beta rule: 0 none, 1 minimum, 2 mean, 3 maximum; default 1.
      integer, intent(in), optional :: alpha_size !! Apply the alpha rule only when both child clusters exceed this size; default 3.
      integer, intent(in), optional :: beta_size !! Apply the beta rule only when both child clusters exceed this size; default 4.
      logical, intent(in), optional :: correct_similarity !! Correct composite correlations for reliability; default true.
      real(dp), allocatable :: keys(:, :), variance(:), alpha(:), beta(:), reliability(:), similarity(:, :)
      real(dp), allocatable :: covariance(:, :)
      integer, allocatable :: sizes(:)
      logical, allocatable :: active(:), rejected(:, :)
      real(dp) :: selected_similarity, raw_correlation, cov12, r12, new_variance, new_alpha, new_beta
      real(dp) :: alpha_reference, beta_reference, threshold_beta, denom
      integer :: p, i, j, left, right, sign_right, active_count, target, arule, brule, asize, bsize
      integer :: local_merges, out_col
      logical :: allow_reverse, use_correction, acceptable, found_pair

      p = size(r, 1)
      info = 0
      n_merges = 0
      if (p < 3 .or. size(r, 2) /= p) info = 1
      if (n_clusters < 0 .or. n_clusters > p) info = 2
      do j = 1, size(r, 2)
         do i = 1, size(r, 1)
            if (ieee_is_nan(r(i, j))) info = 3
         end do
      end do
      allocate(merge_left(max(0, p - 1)), merge_right(max(0, p - 1)), merge_sign(max(0, p - 1)))
      allocate(merge_similarity(max(0, p - 1)), merge_alpha(max(0, p - 1)), merge_beta(max(0, p - 1)))
      allocate(merge_size(max(0, p - 1)))
      merge_left = 0
      merge_right = 0
      merge_sign = 0
      merge_similarity = 0.0_dp
      merge_alpha = 0.0_dp
      merge_beta = 0.0_dp
      merge_size = 0
      if (info /= 0) then
         allocate(clusters(p, 0))
         return
      end if

      allow_reverse = .true.
      if (present(reverse)) allow_reverse = reverse
      threshold_beta = 0.5_dp
      if (present(beta_min)) threshold_beta = beta_min
      arule = 3
      if (present(alpha_rule)) arule = max(0, min(3, alpha_rule))
      brule = 1
      if (present(beta_rule)) brule = max(0, min(3, beta_rule))
      asize = 3
      if (present(alpha_size)) asize = max(0, alpha_size)
      bsize = 4
      if (present(beta_size)) bsize = max(0, beta_size)
      use_correction = .true.
      if (present(correct_similarity)) use_correction = correct_similarity
      target = max(1, n_clusters)

      allocate(keys(p, p), variance(p), alpha(p), beta(p), reliability(p), sizes(p), active(p))
      allocate(similarity(p, p), covariance(p, p), rejected(p, p))
      keys = 0.0_dp
      do i = 1, p
         keys(i, i) = 1.0_dp
      end do
      active = .true.
      sizes = 1
      variance = 1.0_dp
      alpha = 0.0_dp
      beta = 0.0_dp
      do i = 1, p
         reliability(i) = max_abs_off_diagonal(r, i)
         alpha(i) = reliability(i)
         beta(i) = reliability(i)
      end do
      active_count = p
      local_merges = 0

      do while (active_count > 1)
         if (n_clusters > 0 .and. active_count <= target) exit
         covariance = 0.0_dp
         similarity = -huge(1.0_dp)
         rejected = .false.
         do j = 1, p
            if (.not. active(j)) cycle
            variance(j) = dot_product(keys(:, j), matmul(r, keys(:, j)))
         end do
         do j = 2, p
            if (.not. active(j)) cycle
            do i = 1, j - 1
               if (.not. active(i)) cycle
               covariance(i, j) = dot_product(keys(:, i), matmul(r, keys(:, j)))
               covariance(j, i) = covariance(i, j)
               denom = sqrt(max(tiny(1.0_dp), variance(i) * variance(j)))
               raw_correlation = covariance(i, j) / denom
               if (use_correction) then
                  denom = sqrt(max(tiny(1.0_dp), reliability(i) * reliability(j)))
                  similarity(i, j) = raw_correlation / denom
               else
                  similarity(i, j) = raw_correlation
               end if
               similarity(j, i) = similarity(i, j)
            end do
         end do

         found_pair = .false.
         do
            call select_similarity_pair(similarity, active, rejected, allow_reverse, left, right, sign_right, &
               selected_similarity, found_pair)
            if (.not. found_pair) exit
            cov12 = real(sign_right, dp) * covariance(left, right)
            r12 = cov12 / real(max(1, sizes(left) * sizes(right)), dp)
            new_variance = variance(left) + variance(right) + 2.0_dp * cov12
            if (new_variance <= tiny(1.0_dp)) then
               rejected(left, right) = .true.
               rejected(right, left) = .true.
               cycle
            end if
            new_alpha = alpha_from_variance(new_variance, sizes(left) + sizes(right))
            new_beta = real((sizes(left) + sizes(right)) ** 2, dp) * r12 / new_variance
            acceptable = .true.
            if (n_clusters == 0) then
               if (arule > 0 .and. min(sizes(left), sizes(right)) > asize) then
                  alpha_reference = rule_reference(alpha(left), alpha(right), arule)
                  if (new_alpha < alpha_reference) acceptable = .false.
               end if
               if (brule > 0 .and. min(sizes(left), sizes(right)) > bsize) then
                  beta_reference = rule_reference(beta(left), beta(right), brule)
                  if (new_beta < beta_reference) acceptable = .false.
               end if
               if (new_beta < threshold_beta) acceptable = .false.
            end if
            if (acceptable) exit
            rejected(left, right) = .true.
            rejected(right, left) = .true.
         end do
         if (.not. found_pair) exit

         local_merges = local_merges + 1
         merge_left(local_merges) = left
         merge_right(local_merges) = right
         merge_sign(local_merges) = sign_right
         merge_similarity(local_merges) = selected_similarity
         merge_alpha(local_merges) = new_alpha
         merge_beta(local_merges) = new_beta
         merge_size(local_merges) = sizes(left) + sizes(right)
         keys(:, left) = keys(:, left) + real(sign_right, dp) * keys(:, right)
         keys(:, right) = 0.0_dp
         sizes(left) = sizes(left) + sizes(right)
         sizes(right) = 0
         variance(left) = new_variance
         alpha(left) = new_alpha
         beta(left) = new_beta
         reliability(left) = max(tiny(1.0_dp), new_alpha)
         active(right) = .false.
         active_count = active_count - 1
      end do

      allocate(clusters(p, active_count))
      out_col = 0
      do i = 1, p
         if (.not. active(i)) cycle
         out_col = out_col + 1
         clusters(:, out_col) = keys(:, i)
      end do
      n_merges = local_merges
   end subroutine iclust_hierarchical

   pure subroutine iclust_sort_loadings(loadings, cut, cluster_order, assignments, item_order, sorted_loadings, &
      pattern, cluster_sort, info)
      real(dp), intent(in) :: loadings(:, :) !! Item-by-cluster loading matrix to classify and sort.
      real(dp), intent(in) :: cut !! Minimum absolute primary loading required for cluster assignment.
      integer, allocatable, intent(out) :: cluster_order(:) !! Cluster-column order, descending by eigenvalue-like strength.
      integer, allocatable, intent(out) :: assignments(:) !! Assigned cluster for each original item; zero means below `cut`.
      integer, allocatable, intent(out) :: item_order(:) !! Item permutation grouped by cluster and descending primary loading.
      real(dp), allocatable, intent(out) :: sorted_loadings(:, :) !! Loading matrix reordered by `item_order` and `cluster_order`.
      real(dp), intent(in), optional :: pattern(:, :) !! Optional pattern matrix used with loadings for cluster strength.
      logical, intent(in), optional :: cluster_sort !! Reorder cluster columns by strength when true; default true.
      integer, intent(out), optional :: info !! Zero on success; nonzero for incompatible optional pattern dimensions.
      real(dp), allocatable :: strength(:), ordered(:, :)
      integer, allocatable :: work_items(:)
      real(dp) :: best, key_strength
      integer :: p, k, i, j, c, key_cluster, key_item, pos, start_pos, end_pos
      logical :: do_cluster_sort

      p = size(loadings, 1)
      k = size(loadings, 2)
      allocate(cluster_order(k), assignments(p), item_order(p), sorted_loadings(p, k), strength(k), ordered(p, k))
      if (present(pattern)) then
         if (size(pattern, 1) /= p .or. size(pattern, 2) /= k) then
            cluster_order = 0
            assignments = 0
            item_order = 0
            sorted_loadings = 0.0_dp
            if (present(info)) info = 1
            return
         end if
         do c = 1, k
            strength(c) = dot_product(pattern(:, c), loadings(:, c))
         end do
      else
         do c = 1, k
            strength(c) = sum(loadings(:, c) * loadings(:, c))
         end do
      end if
      do c = 1, k
         cluster_order(c) = c
      end do
      do_cluster_sort = .true.
      if (present(cluster_sort)) do_cluster_sort = cluster_sort
      if (do_cluster_sort) then
         do i = 2, k
            key_cluster = cluster_order(i)
            key_strength = strength(key_cluster)
            j = i - 1
            do while (j >= 1)
               if (strength(cluster_order(j)) >= key_strength) exit
               cluster_order(j + 1) = cluster_order(j)
               j = j - 1
            end do
            cluster_order(j + 1) = key_cluster
         end do
      end if
      ordered = loadings(:, cluster_order)

      assignments = 0
      do i = 1, p
         best = 0.0_dp
         do c = 1, k
            if (abs(ordered(i, c)) > best) then
               best = abs(ordered(i, c))
               assignments(i) = c
            end if
         end do
         if (best < cut) assignments(i) = 0
      end do
      allocate(work_items(p))
      pos = 0
      do c = 1, k
         start_pos = pos + 1
         do i = 1, p
            if (assignments(i) == c) then
               pos = pos + 1
               work_items(pos) = i
            end if
         end do
         end_pos = pos
         do i = start_pos + 1, end_pos
            key_item = work_items(i)
            j = i - 1
            do while (j >= start_pos)
               if (abs(ordered(work_items(j), c)) >= abs(ordered(key_item, c))) exit
               work_items(j + 1) = work_items(j)
               j = j - 1
            end do
            work_items(j + 1) = key_item
         end do
      end do
      do i = 1, p
         if (assignments(i) == 0) then
            pos = pos + 1
            work_items(pos) = i
         end if
      end do
      item_order = work_items
      sorted_loadings = ordered(item_order, :)
      if (present(info)) info = 0
   end subroutine iclust_sort_loadings

   pure real(dp) function max_abs_off_diagonal(r, row) result(value)
      real(dp), intent(in) :: r(:, :) !! Correlation matrix searched for an item's largest off-diagonal relationship.
      integer, intent(in) :: row !! One-based row/item index.
      integer :: j

      value = tiny(1.0_dp)
      do j = 1, size(r, 2)
         if (j == row) cycle
         value = max(value, abs(r(row, j)))
      end do
   end function max_abs_off_diagonal

   pure real(dp) function alpha_from_variance(variance, size_cluster) result(value)
      real(dp), intent(in) :: variance !! Variance of the signed unit-weight cluster composite.
      integer, intent(in) :: size_cluster !! Number of items in the cluster.

      if (size_cluster > 1 .and. variance > tiny(1.0_dp)) then
         value = ((variance - real(size_cluster, dp)) / variance) * real(size_cluster, dp) / &
            real(size_cluster - 1, dp)
      else
         value = 0.0_dp
      end if
   end function alpha_from_variance

   pure real(dp) function rule_reference(first, second, rule) result(value)
      real(dp), intent(in) :: first !! Reliability of the first child cluster.
      real(dp), intent(in) :: second !! Reliability of the second child cluster.
      integer, intent(in) :: rule !! Rule code: 1 minimum, 2 mean, 3 maximum.

      select case (rule)
      case (1)
         value = min(first, second)
      case (2)
         value = 0.5_dp * (first + second)
      case default
         value = max(first, second)
      end select
   end function rule_reference

   pure subroutine select_similarity_pair(similarity, active, rejected, allow_reverse, left, right, sign_right, &
      selected_similarity, found)
      real(dp), intent(in) :: similarity(:, :) !! Symmetric active-cluster similarity matrix.
      logical, intent(in) :: active(:) !! Active-cluster flags.
      logical, intent(in) :: rejected(:, :) !! Pair flags excluded by reliability criteria in the current iteration.
      logical, intent(in) :: allow_reverse !! Select by absolute similarity and reverse negative pairs when true.
      integer, intent(out) :: left !! Left selected cluster slot.
      integer, intent(out) :: right !! Right selected cluster slot.
      integer, intent(out) :: sign_right !! Sign applied to the right cluster, plus or minus one.
      real(dp), intent(out) :: selected_similarity !! Similarity value used to select the pair.
      logical, intent(out) :: found !! True when at least one candidate pair remains.
      real(dp) :: criterion, best_criterion
      integer :: i, j

      found = .false.
      left = 0
      right = 0
      sign_right = 1
      selected_similarity = 0.0_dp
      best_criterion = -huge(1.0_dp)
      do j = 2, size(active)
         if (.not. active(j)) cycle
         do i = 1, j - 1
            if (.not. active(i) .or. rejected(i, j)) cycle
            if (allow_reverse) then
               criterion = abs(similarity(i, j))
            else
               criterion = similarity(i, j)
            end if
            if (criterion > best_criterion) then
               best_criterion = criterion
               left = i
               right = j
               selected_similarity = similarity(i, j)
               sign_right = 1
               if (allow_reverse .and. similarity(i, j) < 0.0_dp) sign_right = -1
               found = .true.
            end if
         end do
      end do
   end subroutine select_similarity_pair

end module psych_iclust_extended
