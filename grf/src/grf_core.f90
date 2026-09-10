module grf_core
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
   use grf_kinds, only : dp
   use grf_rng, only : grf_rng_state
   use grf_types, only : grf_options, grf_tree, grf_forest
   use grf_types, only : grf_regression, grf_multi_regression, grf_probability, grf_quantile
   use grf_types, only : grf_causal, grf_instrumental, grf_survival, grf_causal_survival, grf_lm
   use grf_utils, only : invert_matrix, sort_indices_by_values, weighted_quantile, solve_linear_system
   implicit none
   private

   public :: train_forest_core
   public :: find_leaf_node
   public :: compute_forest_weights

contains

   subroutine train_forest_core(forest, options, info)
      type(grf_forest), intent(inout) :: forest !! Forest object containing training data, optional cluster metadata, and family information; trees are replaced.
      type(grf_options), intent(in) :: options !! Tree count, cluster/row sampling, honesty, split-balance, grouped-tree, and random-seed settings.
      integer, intent(out) :: info !! Zero on success; negative values indicate invalid training dimensions or options.
      type(grf_rng_state) :: rng
      integer, allocatable :: all_index(:)
      integer, allocatable :: half_sample(:)
      integer :: group_first
      integer :: num_trees
      integer :: group_last
      integer :: half_size
      integer :: population_size
      integer :: i
      integer :: t

      info = validate_training_data(forest, options)
      if (info /= 0) return
      forest%options = options
      num_trees = options%num_trees
      if (options%ci_group_size > 1) then
         if (mod(num_trees, options%ci_group_size) /= 0) then
            num_trees = num_trees + options%ci_group_size - mod(num_trees, options%ci_group_size)
         end if
      end if
      forest%options%num_trees = num_trees
      if (allocated(forest%trees)) deallocate(forest%trees)
      allocate(forest%trees(num_trees))
      call rng%seed(options%seed)
      if (options%ci_group_size == 1) then
         do t = 1, num_trees
            call train_one_tree(forest, options, rng, forest%trees(t), info)
            if (info /= 0) return
         end do
         return
      end if

      if (allocated(forest%clusters)) then
         population_size = forest%n_clusters
      else
         population_size = size(forest%x, 1)
      end if
      half_size = max(1, int(0.5_dp * real(population_size, dp)))
      allocate(all_index(population_size), half_sample(half_size))
      all_index = [(i, i = 1, population_size)]
      group_first = 1
      do while (group_first <= num_trees)
         call rng%shuffle(all_index)
         half_sample = all_index(1:half_size)
         group_last = min(num_trees, group_first + options%ci_group_size - 1)
         do t = group_first, group_last
            call train_one_tree(forest, options, rng, forest%trees(t), info, half_sample)
            if (info /= 0) return
         end do
         group_first = group_last + 1
      end do
   end subroutine train_forest_core

   pure integer function validate_training_data(forest, options) result(info)
      type(grf_forest), intent(in) :: forest !! Candidate forest object whose training arrays are checked for consistency.
      type(grf_options), intent(in) :: options !! Candidate forest options checked for valid ranges.

      info = 0
      if (.not. allocated(forest%x)) then
         info = -1
         return
      end if
      if (size(forest%x,1) < 2 .or. size(forest%x,2) < 1) then
         info = -2
         return
      end if
      if (options%num_trees < 1 .or. options%min_node_size < 1 .or. options%ci_group_size < 1) then
         info = -3
         return
      end if
      if (options%sample_fraction <= 0.0_dp .or. options%sample_fraction > 1.0_dp) then
         info = -4
         return
      end if
      if (options%ci_group_size > 1 .and. options%sample_fraction > 0.5_dp) then
         info = -8
         return
      end if
      if (options%honesty) then
         if (options%honesty_fraction <= 0.0_dp .or. options%honesty_fraction >= 1.0_dp) then
            info = -5
            return
         end if
      end if
      if (.not. allocated(forest%sample_weights)) then
         info = -6
         return
      end if
      if (size(forest%sample_weights) /= size(forest%x,1)) then
         info = -7
         return
      end if
   end function validate_training_data

   subroutine train_one_tree(forest, options, rng, tree, info, candidate_pool)
      type(grf_forest), intent(in) :: forest !! Training data, optional cluster assignments, and forest-family specification used by this tree.
      type(grf_options), intent(in) :: options !! Sampling and splitting options applied to this tree.
      type(grf_rng_state), intent(inout) :: rng !! Mutable deterministic generator used for subsampling and variable selection.
      type(grf_tree), intent(out) :: tree !! Trained tree with split arrays, honest leaf membership, and cluster-aware in-bag mask.
      integer, intent(out) :: info !! Zero on success; negative values indicate an impossible row or cluster subsample.
      integer, intent(in), optional :: candidate_pool(:) !! Optional grouped half-sample pool: row IDs without clustering or cluster IDs with clustering.
      integer, allocatable :: sampled(:)
      integer, allocatable :: grow_samples(:)
      integer, allocatable :: estimate_samples(:)
      integer, allocatable :: inbag_samples(:)
      integer, allocatable :: leaf_counts(:)
      integer :: n
      integer :: next_node
      integer :: max_nodes
      integer :: pruned_root

      n = size(forest%x, 1)
      if (allocated(forest%clusters)) then
         call prepare_cluster_tree_samples(forest, options, rng, sampled, grow_samples, estimate_samples, &
                                           inbag_samples, info, candidate_pool)
      else
         call prepare_row_tree_samples(n, options, rng, sampled, grow_samples, estimate_samples, info, candidate_pool)
         if (info == 0) inbag_samples = sampled
      end if
      if (info /= 0) return
      if (size(grow_samples) < 1 .or. size(estimate_samples) < 1) then
         info = -20
         return
      end if

      max_nodes = max(3, 2 * size(grow_samples) + 1)
      call allocate_tree(tree, max_nodes, n)
      tree%inbag = .false.
      tree%inbag(inbag_samples) = .true.
      tree%n_nodes = 1
      tree%root = 1
      tree%depth(1) = 0
      next_node = 1
      call grow_node(forest, options, rng, tree, grow_samples, 1, next_node)
      tree%n_nodes = next_node

      allocate(leaf_counts(tree%n_nodes))
      call count_leaf_members(tree, forest%x, estimate_samples, leaf_counts)
      if (options%honesty .and. options%honesty_prune_leaves) then
         pruned_root = prune_empty_subtrees(tree, tree%root, leaf_counts)
         tree%root = pruned_root
      end if
      call populate_leaf_members(tree, forest%x, estimate_samples)
      info = 0
   end subroutine train_one_tree

   subroutine prepare_row_tree_samples(n, options, rng, sampled, grow_samples, estimate_samples, info, candidate_pool)
      integer, intent(in) :: n !! Total number of training rows available when cluster sampling is disabled.
      type(grf_options), intent(in) :: options !! Row-sampling fraction and honesty controls for the current tree.
      type(grf_rng_state), intent(inout) :: rng !! Mutable deterministic generator used to shuffle row identifiers.
      integer, allocatable, intent(out) :: sampled(:) !! Rows selected for the tree before the honesty split.
      integer, allocatable, intent(out) :: grow_samples(:) !! Rows used to choose tree splits.
      integer, allocatable, intent(out) :: estimate_samples(:) !! Rows used to populate honest terminal leaves.
      integer, intent(out) :: info !! Zero on success or a negative code when too few rows are available.
      integer, intent(in), optional :: candidate_pool(:) !! Optional row IDs in the group's half-sample.
      integer, allocatable :: pool(:)
      integer :: n_sample
      integer :: n_grow
      real(dp) :: fraction
      integer :: i

      info = 0
      if (present(candidate_pool)) then
         pool = candidate_pool
         fraction = min(1.0_dp, 2.0_dp * options%sample_fraction)
         n_sample = ceiling(fraction * real(size(pool), dp))
      else
         allocate(pool(n))
         pool = [(i, i = 1, n)]
         n_sample = ceiling(options%sample_fraction * real(n, dp))
      end if
      n_sample = max(1, min(size(pool), n_sample))
      if (options%honesty .and. n_sample < 2) then
         info = -21
         return
      end if
      call rng%shuffle(pool)
      sampled = pool(1:n_sample)
      if (options%honesty) then
         n_grow = max(1, min(n_sample - 1, ceiling(options%honesty_fraction * real(n_sample, dp))))
         grow_samples = sampled(1:n_grow)
         estimate_samples = sampled(n_grow+1:n_sample)
      else
         grow_samples = sampled
         estimate_samples = sampled
      end if
   end subroutine prepare_row_tree_samples

   subroutine prepare_cluster_tree_samples(forest, options, rng, sampled_clusters, grow_samples, estimate_samples, &
                                           inbag_samples, info, candidate_pool)
      type(grf_forest), intent(in) :: forest !! Training forest with one-based canonical cluster IDs and samples-per-cluster limit.
      type(grf_options), intent(in) :: options !! Cluster-sampling fraction and honesty controls for the current tree.
      type(grf_rng_state), intent(inout) :: rng !! Mutable deterministic generator used for cluster and within-cluster sampling.
      integer, allocatable, intent(out) :: sampled_clusters(:) !! Cluster IDs selected for the tree before the honesty split.
      integer, allocatable, intent(out) :: grow_samples(:) !! Within-cluster sampled rows used to choose tree splits.
      integer, allocatable, intent(out) :: estimate_samples(:) !! Independently sampled rows from honesty estimate clusters.
      integer, allocatable, intent(out) :: inbag_samples(:) !! All rows belonging to selected clusters, matching upstream cluster-aware OOB status.
      integer, intent(out) :: info !! Zero on success or negative when too few clusters are available.
      integer, intent(in), optional :: candidate_pool(:) !! Optional cluster IDs from grouped half-sample.
      integer, allocatable :: pool(:)
      integer, allocatable :: grow_clusters(:)
      integer, allocatable :: estimate_clusters(:)
      integer :: n_sample
      integer :: n_grow
      real(dp) :: fraction
      integer :: i

      info = 0
      if (present(candidate_pool)) then
         pool = candidate_pool
         fraction = min(1.0_dp, 2.0_dp * options%sample_fraction)
         n_sample = ceiling(fraction * real(size(pool), dp))
      else
         allocate(pool(forest%n_clusters))
         pool = [(i, i = 1, forest%n_clusters)]
         n_sample = ceiling(options%sample_fraction * real(forest%n_clusters, dp))
      end if
      n_sample = max(1, min(size(pool), n_sample))
      if (options%honesty .and. n_sample < 2) then
         info = -22
         return
      end if
      call rng%shuffle(pool)
      sampled_clusters = pool(1:n_sample)
      call rows_in_clusters(forest%clusters, sampled_clusters, inbag_samples)
      if (options%honesty) then
         n_grow = max(1, min(n_sample - 1, ceiling(options%honesty_fraction * real(n_sample, dp))))
         grow_clusters = sampled_clusters(1:n_grow)
         estimate_clusters = sampled_clusters(n_grow+1:n_sample)
         call sample_rows_from_clusters(forest, grow_clusters, rng, grow_samples)
         call sample_rows_from_clusters(forest, estimate_clusters, rng, estimate_samples)
      else
         call sample_rows_from_clusters(forest, sampled_clusters, rng, grow_samples)
         estimate_samples = grow_samples
      end if
   end subroutine prepare_cluster_tree_samples

   subroutine sample_rows_from_clusters(forest, selected_clusters, rng, samples)
      type(grf_forest), intent(in) :: forest !! Forest containing canonical cluster membership and per-cluster sampling cap.
      integer, intent(in) :: selected_clusters(:) !! Cluster IDs to sample from independently.
      type(grf_rng_state), intent(inout) :: rng !! Mutable generator used to randomize rows within each selected cluster.
      integer, allocatable, intent(out) :: samples(:) !! Concatenated selected row IDs from all requested clusters.
      integer, allocatable :: work(:)
      integer, allocatable :: members(:)
      integer :: max_rows
      integer :: n_samples
      integer :: n_members
      integer :: take
      integer :: c
      integer :: i

      max_rows = max(1, forest%samples_per_cluster) * size(selected_clusters)
      allocate(work(max_rows))
      n_samples = 0
      do c = 1, size(selected_clusters)
         n_members = count(forest%clusters == selected_clusters(c))
         allocate(members(n_members))
         n_members = 0
         do i = 1, forest%n_train
            if (forest%clusters(i) == selected_clusters(c)) then
               n_members = n_members + 1
               members(n_members) = i
            end if
         end do
         call rng%shuffle(members)
         take = min(size(members), forest%samples_per_cluster)
         if (take > 0) then
            work(n_samples+1:n_samples+take) = members(1:take)
            n_samples = n_samples + take
         end if
         deallocate(members)
      end do
      allocate(samples(n_samples))
      if (n_samples > 0) samples = work(1:n_samples)
   end subroutine sample_rows_from_clusters

   pure subroutine rows_in_clusters(clusters, selected_clusters, rows)
      integer, intent(in) :: clusters(:) !! One-based cluster identifier for every training observation.
      integer, intent(in) :: selected_clusters(:) !! Cluster IDs whose complete observation sets are requested.
      integer, allocatable, intent(out) :: rows(:) !! Row IDs belonging to any selected cluster.
      integer :: n_rows
      integer :: i
      integer :: j
      logical :: selected

      n_rows = 0
      do i = 1, size(clusters)
         selected = .false.
         do j = 1, size(selected_clusters)
            if (clusters(i) == selected_clusters(j)) then
               selected = .true.
               exit
            end if
         end do
         if (selected) n_rows = n_rows + 1
      end do
      allocate(rows(n_rows))
      n_rows = 0
      do i = 1, size(clusters)
         selected = .false.
         do j = 1, size(selected_clusters)
            if (clusters(i) == selected_clusters(j)) then
               selected = .true.
               exit
            end if
         end do
         if (selected) then
            n_rows = n_rows + 1
            rows(n_rows) = i
         end if
      end do
   end subroutine rows_in_clusters

   pure subroutine allocate_tree(tree, max_nodes, n_train)
      type(grf_tree), intent(out) :: tree !! Tree object whose allocatable storage is initialized to empty-node defaults.
      integer, intent(in) :: max_nodes !! Maximum number of nodes permitted for the tree.
      integer, intent(in) :: n_train !! Number of training observations represented by the in-bag mask.

      allocate(tree%left(max_nodes), tree%right(max_nodes), tree%split_var(max_nodes))
      allocate(tree%split_value(max_nodes), tree%missing_left(max_nodes), tree%depth(max_nodes))
      allocate(tree%leaf_start(max_nodes), tree%leaf_count(max_nodes), tree%inbag(n_train))
      tree%left = 0
      tree%right = 0
      tree%split_var = 0
      tree%split_value = 0.0_dp
      tree%missing_left = .true.
      tree%depth = 0
      tree%leaf_start = 0
      tree%leaf_count = 0
   end subroutine allocate_tree

   recursive subroutine grow_node(forest, options, rng, tree, samples, node, next_node)
      type(grf_forest), intent(in) :: forest !! Training data and family metadata used to form node pseudo-responses.
      type(grf_options), intent(in) :: options !! Split-balance, minimum-node, and mtry settings for recursion.
      type(grf_rng_state), intent(inout) :: rng !! Mutable generator used to choose the candidate split variables.
      type(grf_tree), intent(inout) :: tree !! Tree arrays updated with the selected split and descendants.
      integer, intent(in) :: samples(:) !! One-based training-row indices currently assigned to the node.
      integer, intent(in) :: node !! One-based node index being split or finalized as a leaf.
      integer, intent(inout) :: next_node !! Largest node index allocated so far; increased when children are created.
      real(dp), allocatable :: responses(:,:)
      integer, allocatable :: left_samples(:)
      integer, allocatable :: right_samples(:)
      integer :: best_var
      real(dp) :: best_value
      logical :: best_missing_left
      logical :: found
      logical :: valid_response
      integer :: left_node
      integer :: right_node

      if (size(samples) <= options%min_node_size) return
      call build_node_responses(forest, samples, responses, valid_response)
      if (.not. valid_response) return
      if (maxval(maxval(responses, dim=1) - minval(responses, dim=1)) <= &
          100.0_dp * epsilon(1.0_dp)) return
      call find_best_split(forest, options, rng, samples, responses, best_var, best_value, &
                           best_missing_left, found)
      if (.not. found) return
      call partition_samples(forest%x(:,best_var), samples, best_value, best_missing_left, left_samples, right_samples)
      if (size(left_samples) == 0 .or. size(right_samples) == 0) return
      if (next_node + 2 > size(tree%left)) return

      left_node = next_node + 1
      right_node = next_node + 2
      next_node = right_node
      tree%left(node) = left_node
      tree%right(node) = right_node
      tree%split_var(node) = best_var
      tree%split_value(node) = best_value
      tree%missing_left(node) = best_missing_left
      tree%depth(left_node) = tree%depth(node) + 1
      tree%depth(right_node) = tree%depth(node) + 1
      call grow_node(forest, options, rng, tree, left_samples, left_node, next_node)
      call grow_node(forest, options, rng, tree, right_samples, right_node, next_node)
   end subroutine grow_node

   pure subroutine build_node_responses(forest, samples, responses, valid)
      type(grf_forest), intent(in) :: forest !! Forest family and training outcomes/treatments used for relabeling.
      integer, intent(in) :: samples(:) !! One-based row indices in the current tree-growing node.
      real(dp), allocatable, intent(out) :: responses(:,:) !! Node-specific GRF pseudo-responses, one row per sample.
      logical, intent(out) :: valid !! False when a causal local moment matrix is singular or the node has no survival signal.
      integer :: i
      integer :: j
      integer :: k
      integer :: n
      integer :: nt
      integer :: no
      real(dp), allocatable :: yy(:,:)
      real(dp), allocatable :: ww(:,:)
      real(dp), allocatable :: zz(:)
      real(dp), allocatable :: sw(:)
      real(dp), allocatable :: wc(:,:)
      real(dp), allocatable :: yc(:,:)
      real(dp), allocatable :: gram(:,:)
      real(dp), allocatable :: gram_inv(:,:)
      real(dp), allocatable :: beta(:,:)
      real(dp), allocatable :: residual(:,:)
      real(dp), allocatable :: rho(:,:)
      real(dp) :: ymean
      real(dp) :: wmean
      real(dp) :: zmean
      real(dp) :: zreg_mean
      real(dp) :: numerator
      real(dp) :: denominator
      real(dp) :: tau
      real(dp) :: zreg
      real(dp) :: total_weight
      real(dp) :: qvalue
      integer :: info
      integer :: cls

      n = size(samples)
      valid = .true.
      allocate(sw(n))
      sw = forest%sample_weights(samples)
      total_weight = sum(sw)
      if (total_weight <= 100.0_dp * epsilon(1.0_dp)) then
         allocate(responses(n,1))
         responses = 0.0_dp
         valid = .false.
         return
      end if

      select case (forest%family)
      case (grf_regression, grf_multi_regression)
         if (forest%family == grf_regression .and. forest%ll_split_enabled) then
            call build_ll_split_responses(forest, samples, responses, valid)
            return
         end if
         allocate(responses(n, forest%n_outputs))
         responses = forest%y(samples,1:forest%n_outputs)

      case (grf_probability)
         allocate(responses(n, forest%n_classes))
         responses = 0.0_dp
         do i = 1, n
            cls = forest%classes(samples(i))
            if (cls >= 1 .and. cls <= forest%n_classes) responses(i,cls) = 1.0_dp
         end do

      case (grf_quantile)
         allocate(responses(n, max(2, size(forest%quantiles) + 1)))
         responses = 0.0_dp
         do i = 1, n
            cls = 1
            do j = 1, size(forest%quantiles)
               call weighted_quantile(forest%y(samples,1), sw, forest%quantiles(j), qvalue)
               if (forest%y(samples(i),1) > qvalue) cls = cls + 1
            end do
            responses(i,min(cls,size(responses,2))) = 1.0_dp
         end do

      case (grf_causal, grf_instrumental)
         allocate(yy(n,1), ww(n,1), zz(n), responses(n,1))
         yy(:,1) = forest%y(samples,1) - forest%y_hat(samples,1)
         ww(:,1) = forest%w(samples,1) - forest%w_hat(samples,1)
         if (forest%family == grf_instrumental) then
            zz = forest%z(samples) - forest%z_hat(samples)
         else
            zz = ww(:,1)
         end if
         ymean = sum(sw * yy(:,1)) / total_weight
         wmean = sum(sw * ww(:,1)) / total_weight
         zmean = sum(sw * zz) / total_weight
         zreg_mean = (1.0_dp - forest%reduced_form_weight) * zmean + forest%reduced_form_weight * wmean
         numerator = 0.0_dp
         denominator = 0.0_dp
         do i = 1, n
            zreg = (1.0_dp - forest%reduced_form_weight) * zz(i) + forest%reduced_form_weight * ww(i,1)
            numerator = numerator + sw(i) * (zreg - zreg_mean) * (yy(i,1) - ymean)
            denominator = denominator + sw(i) * (zreg - zreg_mean) * (ww(i,1) - wmean)
         end do
         if (abs(denominator) <= 1.0e-10_dp * max(1.0_dp, abs(numerator))) then
            responses = 0.0_dp
            valid = .false.
            return
         end if
         tau = numerator / denominator
         do i = 1, n
            zreg = (1.0_dp - forest%reduced_form_weight) * zz(i) + forest%reduced_form_weight * ww(i,1)
            responses(i,1) = (zreg - zreg_mean) * ((yy(i,1) - ymean) - tau * (ww(i,1) - wmean))
         end do

      case (grf_causal_survival)
         allocate(responses(n,1))
         numerator = sum(sw * forest%causal_survival_numerator(samples))
         denominator = sum(sw * forest%causal_survival_denominator(samples))
         if (abs(denominator) <= 1.0e-10_dp .or. total_weight <= 1.0e-16_dp) then
            responses = 0.0_dp
            valid = .false.
            return
         end if
         tau = numerator / denominator
         do i = 1, n
            responses(i,1) = (forest%causal_survival_numerator(samples(i)) - &
               forest%causal_survival_denominator(samples(i)) * tau) / denominator
         end do

      case (grf_lm)
         nt = forest%n_treatments
         no = forest%n_outputs
         if (n <= nt) then
            allocate(responses(n,max(1,nt*no)))
            responses = 0.0_dp
            valid = .false.
            return
         end if
         allocate(yy(n,no), ww(n,nt), wc(n,nt), yc(n,no))
         yy = forest%y(samples,1:no) - forest%y_hat(samples,1:no)
         ww = forest%w(samples,1:nt) - forest%w_hat(samples,1:nt)
         do j = 1, no
            yc(:,j) = yy(:,j) - sum(sw * yy(:,j)) / total_weight
         end do
         do j = 1, nt
            wc(:,j) = ww(:,j) - sum(sw * ww(:,j)) / total_weight
         end do
         allocate(gram(nt,nt), gram_inv(nt,nt), beta(nt,no), residual(n,no), rho(n,nt))
         gram = 0.0_dp
         do i = 1, n
            do j = 1, nt
               do k = 1, nt
                  gram(j,k) = gram(j,k) + sw(i) * wc(i,j) * wc(i,k)
               end do
            end do
         end do
         call invert_matrix(gram, gram_inv, info)
         if (info /= 0) then
            allocate(responses(n,max(1,nt*no)))
            responses = 0.0_dp
            valid = .false.
            return
         end if
         beta = 0.0_dp
         do j = 1, nt
            do k = 1, no
               beta(j,k) = sum(sw * wc(:,j) * yc(:,k))
            end do
         end do
         beta = matmul(gram_inv, beta)
         residual = yc - matmul(wc, beta)
         rho = matmul(wc, transpose(gram_inv))
         allocate(responses(n,nt*no))
         do i = 1, n
            cls = 0
            do k = 1, no
               do j = 1, nt
                  cls = cls + 1
                  responses(i,cls) = rho(i,j) * residual(i,k)
               end do
            end do
         end do

      case (grf_survival)
         allocate(responses(n,1))
         responses(:,1) = real(forest%survival_time_index(samples), dp)
         valid = survival_node_has_split_signal(forest, samples)

      case default
         allocate(responses(n,1))
         responses = 0.0_dp
         valid = .false.
      end select
   end subroutine build_node_responses

   pure subroutine build_ll_split_responses(forest, samples, responses, valid)
      type(grf_forest), intent(in) :: forest !! Local-linear forest carrying split ridge controls and optional overall coefficients.
      integer, intent(in) :: samples(:) !! One-based training rows in the current tree-growing node.
      real(dp), allocatable, intent(out) :: responses(:,:) !! Ridge residual relabeling values used by the ordinary regression split rule.
      logical, intent(out) :: valid !! True when local or fallback ridge coefficients produce finite residual responses.
      real(dp), allocatable :: design(:,:)
      real(dp), allocatable :: gram(:,:)
      real(dp), allocatable :: rhs(:)
      real(dp), allocatable :: beta(:)
      real(dp) :: normalization
      integer :: i
      integer :: j
      integer :: solve_info

      allocate(design(size(samples),size(forest%ll_split_variables)+1))
      design(:,1) = 1.0_dp
      do j = 1, size(forest%ll_split_variables)
         do i = 1, size(samples)
            if (ieee_is_nan(forest%x(samples(i),forest%ll_split_variables(j)))) then
               design(i,j+1) = 0.0_dp
            else
               design(i,j+1) = forest%x(samples(i),forest%ll_split_variables(j))
            end if
         end do
      end do
      allocate(beta(size(design,2)))
      if (size(samples) < forest%ll_split_cutoff .and. size(forest%ll_split_overall_beta) == size(beta)) then
         beta = forest%ll_split_overall_beta
      else
         allocate(gram(size(design,2),size(design,2)), rhs(size(design,2)))
         gram = matmul(transpose(design), design)
         rhs = matmul(transpose(design), forest%y(samples,1))
         if (forest%ll_split_weight_penalty) then
            do j = 2, size(gram,1)
               gram(j,j) = gram(j,j) + forest%ll_split_lambda * gram(j,j)
            end do
         else
            normalization = 0.0_dp
            do j = 1, size(gram,1)
               normalization = normalization + gram(j,j)
            end do
            normalization = normalization / real(size(gram,1), dp)
            do j = 2, size(gram,1)
               gram(j,j) = gram(j,j) + forest%ll_split_lambda * normalization
            end do
         end if
         call solve_linear_system(gram, rhs, beta, solve_info)
         if (solve_info /= 0) then
            allocate(responses(size(samples),1))
            responses = 0.0_dp
            valid = .false.
            return
         end if
      end if
      allocate(responses(size(samples),1))
      responses(:,1) = matmul(design, beta) - forest%y(samples,1)
      valid = all(.not. ieee_is_nan(responses(:,1)))
   end subroutine build_ll_split_responses

   subroutine find_best_split(forest, options, rng, samples, responses, best_var, best_value, best_missing_left, found)
      type(grf_forest), intent(in) :: forest !! Training predictors, treatment/event metadata, and sample weights used to score splits.
      type(grf_options), intent(in) :: options !! Candidate-variable count, balance fraction, and imbalance penalty.
      type(grf_rng_state), intent(inout) :: rng !! Mutable generator used for Poisson mtry and candidate-variable shuffling.
      integer, intent(in) :: samples(:) !! One-based training rows in the current node.
      real(dp), intent(in) :: responses(:,:) !! Node pseudo-response matrix produced by the family-specific relabeling step.
      integer, intent(out) :: best_var !! One-based predictor index of the selected split, or zero when no split is found.
      real(dp), intent(out) :: best_value !! Numeric split threshold; may be NaN for a missing-versus-observed split.
      logical, intent(out) :: best_missing_left !! Whether missing predictor values are routed to the selected left child.
      logical, intent(out) :: found !! True when a valid positive-score split satisfying balance constraints is found.
      integer, allocatable :: candidates(:)
      integer, allocatable :: order(:)
      real(dp), allocatable :: node_values(:)
      integer :: p
      integer :: mtry_base
      integer :: mtry_draw
      integer :: c
      integer :: pos
      integer :: send_case
      integer :: var
      integer :: min_child
      real(dp) :: v1
      real(dp) :: v2
      real(dp) :: threshold
      real(dp) :: score
      real(dp) :: best_score
      logical :: send_missing_left
      logical :: valid
      logical :: has_missing
      logical :: has_finite

      p = size(forest%x,2)
      mtry_base = options%mtry
      if (mtry_base <= 0) mtry_base = min(p, ceiling(sqrt(real(p,dp)) + 20.0_dp))
      mtry_draw = max(1, min(p, rng%poisson(real(mtry_base,dp))))
      allocate(candidates(p))
      candidates = [(c, c = 1, p)]
      call rng%shuffle(candidates)
      min_child = max(1, ceiling(options%alpha * real(size(samples),dp)))
      best_var = 0
      best_value = 0.0_dp
      best_missing_left = .true.
      best_score = -huge(1.0_dp)
      found = .false.

      do c = 1, mtry_draw
         var = candidates(c)
         allocate(node_values(size(samples)), order(size(samples)))
         node_values = forest%x(samples,var)
         call sort_indices_by_values(node_values, order)
         has_missing = any([(ieee_is_nan(node_values(pos)), pos=1,size(node_values))])
         has_finite = any([(.not. ieee_is_nan(node_values(pos)), pos=1,size(node_values))])
         if (.not. has_finite) then
            deallocate(node_values, order)
            cycle
         end if

         if (has_missing .and. count(.not. ieee_is_nan(node_values)) > 0) then
            threshold = ieee_value(0.0_dp, ieee_quiet_nan)
            call evaluate_split(forest, options, samples, responses, var, threshold, .true., min_child, score, valid)
            if (valid .and. score > best_score) then
               best_score = score
               best_var = var
               best_value = threshold
               best_missing_left = .true.
               found = .true.
            end if
         end if

         do pos = 1, size(order) - 1
            v1 = node_values(order(pos))
            v2 = node_values(order(pos+1))
            if (ieee_is_nan(v1) .or. ieee_is_nan(v2)) cycle
            if (abs(v2 - v1) <= 10.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(v1), abs(v2))) cycle
            threshold = v1 + 0.5_dp * (v2 - v1)
            do send_case = 1, merge(2,1,has_missing)
               send_missing_left = send_case == 1
               call evaluate_split(forest, options, samples, responses, var, threshold, send_missing_left, &
                                   min_child, score, valid)
               if (valid .and. score > best_score) then
                  best_score = score
                  best_var = var
                  best_value = threshold
                  best_missing_left = send_missing_left
                  found = .true.
               end if
            end do
         end do
         deallocate(node_values, order)
      end do
      if (best_score <= 0.0_dp) found = .false.
   end subroutine find_best_split

   pure subroutine evaluate_split(forest, options, samples, responses, var, threshold, missing_left, min_child, score, valid)
      type(grf_forest), intent(in) :: forest !! Training data supplying predictor values, weights, treatment, and event indicators.
      type(grf_options), intent(in) :: options !! Imbalance penalty and minimum node size used in split validation.
      integer, intent(in) :: samples(:) !! One-based training rows in the current node.
      real(dp), intent(in) :: responses(:,:) !! Pseudo-responses whose child sums define the generalized split score.
      integer, intent(in) :: var !! One-based predictor column evaluated at this candidate split.
      real(dp), intent(in) :: threshold !! Numeric threshold, or NaN to split missing values from all finite values.
      logical, intent(in) :: missing_left !! Whether missing predictor values join the left child for finite thresholds.
      integer, intent(in) :: min_child !! Minimum raw observation count required in each child.
      real(dp), intent(out) :: score !! Generalized child-sum score minus the requested family-specific imbalance penalty.
      logical, intent(out) :: valid !! True when both children satisfy count and active family-specific stabilization constraints.
      real(dp), allocatable :: sum_left(:)
      real(dp), allocatable :: sum_right(:)
      real(dp) :: weight_left
      real(dp) :: weight_right
      real(dp) :: value
      real(dp) :: penalty
      real(dp) :: stabilize_size_left
      real(dp) :: stabilize_size_right
      integer :: i
      integer :: j
      integer :: n_left
      integer :: n_right
      integer :: failures_left
      integer :: failures_right
      integer :: min_failures
      logical :: goes_left
      logical :: scalar_stabilized
      logical :: multi_stabilized

      if (forest%family == grf_survival) then
         call evaluate_survival_split(forest, options, samples, var, threshold, missing_left, score, valid)
         return
      end if

      allocate(sum_left(size(responses,2)), sum_right(size(responses,2)))
      sum_left = 0.0_dp
      sum_right = 0.0_dp
      weight_left = 0.0_dp
      weight_right = 0.0_dp
      n_left = 0
      n_right = 0
      failures_left = 0
      failures_right = 0
      do i = 1, size(samples)
         value = forest%x(samples(i),var)
         if (ieee_is_nan(threshold)) then
            goes_left = ieee_is_nan(value)
         else if (ieee_is_nan(value)) then
            goes_left = missing_left
         else
            goes_left = value <= threshold
         end if
         if (goes_left) then
            n_left = n_left + 1
            weight_left = weight_left + forest%sample_weights(samples(i))
            do j = 1, size(responses,2)
               sum_left(j) = sum_left(j) + forest%sample_weights(samples(i)) * responses(i,j)
            end do
            if (allocated(forest%event)) failures_left = failures_left + forest%event(samples(i))
         else
            n_right = n_right + 1
            weight_right = weight_right + forest%sample_weights(samples(i))
            do j = 1, size(responses,2)
               sum_right(j) = sum_right(j) + forest%sample_weights(samples(i)) * responses(i,j)
            end do
            if (allocated(forest%event)) failures_right = failures_right + forest%event(samples(i))
         end if
      end do
      valid = n_left >= min_child .and. n_right >= min_child .and. weight_left > 0.0_dp .and. weight_right > 0.0_dp
      if (.not. valid) then
         score = -huge(1.0_dp)
         return
      end if

      scalar_stabilized = options%stabilize_splits .and. &
         (forest%family == grf_causal .or. forest%family == grf_instrumental .or. &
          forest%family == grf_causal_survival)
      multi_stabilized = options%stabilize_splits .and. forest%family == grf_lm

      if (scalar_stabilized) then
         if (forest%family == grf_causal_survival) then
            min_failures = max(1, ceiling(options%alpha * real(size(samples), dp)))
            if (failures_left < min_failures .or. failures_right < min_failures) then
               valid = .false.
               score = -huge(1.0_dp)
               return
            end if
         end if
         call scalar_stabilized_split_ok(forest, options, samples, var, threshold, missing_left, &
                                         stabilize_size_left, stabilize_size_right, valid)
         if (.not. valid) then
            score = -huge(1.0_dp)
            return
         end if
      else if (multi_stabilized) then
         call multi_stabilized_split_ok(forest, options, samples, var, threshold, missing_left, valid)
         if (.not. valid) then
            score = -huge(1.0_dp)
            return
         end if
      end if

      score = sum(sum_left * sum_left) / weight_left + sum(sum_right * sum_right) / weight_right
      if (scalar_stabilized) then
         penalty = options%imbalance_penalty * (1.0_dp / stabilize_size_left + 1.0_dp / stabilize_size_right)
      else
         penalty = options%imbalance_penalty * (1.0_dp / real(n_left,dp) + 1.0_dp / real(n_right,dp))
      end if
      score = score - penalty
   end subroutine evaluate_split

   pure subroutine scalar_stabilized_split_ok(forest, options, samples, var, threshold, missing_left, &
                                               size_left, size_right, valid)
      type(grf_forest), intent(in) :: forest !! Scalar-treatment causal/IV forest supplying residualized treatment or instrument values.
      type(grf_options), intent(in) :: options !! Stabilized-split alpha, minimum-node-size, and imbalance settings.
      integer, intent(in) :: samples(:) !! One-based rows in the parent node used to define parent treatment moments.
      integer, intent(in) :: var !! Predictor column defining the candidate split.
      real(dp), intent(in) :: threshold !! Candidate finite or NaN split threshold.
      logical, intent(in) :: missing_left !! Missing-value routing rule for finite thresholds.
      real(dp), intent(out) :: size_left !! Weighted child treatment/instrument centered sum of squares used by the GRF penalty.
      real(dp), intent(out) :: size_right !! Weighted right-child treatment/instrument centered sum of squares used by the GRF penalty.
      logical, intent(out) :: valid !! True when both children meet parent-mean balance and parent-variance alpha constraints.
      real(dp) :: weight_node
      real(dp) :: sum_node_z
      real(dp) :: sum_node_z_squared
      real(dp) :: mean_node_z
      real(dp) :: size_node
      real(dp) :: min_child_size
      real(dp) :: weight_left
      real(dp) :: weight_right
      real(dp) :: sum_left_z
      real(dp) :: sum_left_z_squared
      real(dp) :: sum_right_z
      real(dp) :: sum_right_z_squared
      real(dp) :: z_value
      real(dp) :: x_value
      integer :: i
      integer :: row
      integer :: n_left
      integer :: n_right
      integer :: num_node_small_z
      integer :: num_left_small_z
      integer :: num_right_small_z
      integer :: num_left_large_z
      integer :: num_right_large_z
      logical :: goes_left

      weight_node = 0.0_dp
      sum_node_z = 0.0_dp
      sum_node_z_squared = 0.0_dp
      num_node_small_z = 0
      do i = 1, size(samples)
         row = samples(i)
         z_value = split_stabilizer_value(forest, row, 1)
         weight_node = weight_node + forest%sample_weights(row)
         sum_node_z = sum_node_z + forest%sample_weights(row) * z_value
         sum_node_z_squared = sum_node_z_squared + forest%sample_weights(row) * z_value * z_value
      end do
      if (weight_node <= 0.0_dp) then
         valid = .false.
         size_left = 0.0_dp
         size_right = 0.0_dp
         return
      end if
      mean_node_z = sum_node_z / weight_node
      do i = 1, size(samples)
         if (split_stabilizer_value(forest, samples(i), 1) < mean_node_z) num_node_small_z = num_node_small_z + 1
      end do
      size_node = max(0.0_dp, sum_node_z_squared - sum_node_z * sum_node_z / weight_node)
      min_child_size = options%alpha * size_node

      weight_left = 0.0_dp
      sum_left_z = 0.0_dp
      sum_left_z_squared = 0.0_dp
      n_left = 0
      num_left_small_z = 0
      do i = 1, size(samples)
         row = samples(i)
         x_value = forest%x(row,var)
         if (ieee_is_nan(threshold)) then
            goes_left = ieee_is_nan(x_value)
         else if (ieee_is_nan(x_value)) then
            goes_left = missing_left
         else
            goes_left = x_value <= threshold
         end if
         if (.not. goes_left) cycle
         z_value = split_stabilizer_value(forest, row, 1)
         n_left = n_left + 1
         weight_left = weight_left + forest%sample_weights(row)
         sum_left_z = sum_left_z + forest%sample_weights(row) * z_value
         sum_left_z_squared = sum_left_z_squared + forest%sample_weights(row) * z_value * z_value
         if (z_value < mean_node_z) num_left_small_z = num_left_small_z + 1
      end do
      n_right = size(samples) - n_left
      num_right_small_z = num_node_small_z - num_left_small_z
      num_left_large_z = n_left - num_left_small_z
      num_right_large_z = n_right - num_right_small_z
      weight_right = weight_node - weight_left
      sum_right_z = sum_node_z - sum_left_z
      sum_right_z_squared = sum_node_z_squared - sum_left_z_squared
      if (weight_left <= 0.0_dp .or. weight_right <= 0.0_dp) then
         valid = .false.
         size_left = 0.0_dp
         size_right = 0.0_dp
         return
      end if
      size_left = max(0.0_dp, sum_left_z_squared - sum_left_z * sum_left_z / weight_left)
      size_right = max(0.0_dp, sum_right_z_squared - sum_right_z * sum_right_z / weight_right)
      valid = num_left_small_z >= options%min_node_size .and. num_left_large_z >= options%min_node_size .and. &
              num_right_small_z >= options%min_node_size .and. num_right_large_z >= options%min_node_size .and. &
              size_left >= min_child_size .and. size_right >= min_child_size
      if (options%imbalance_penalty > 0.0_dp) then
         valid = valid .and. size_left > 0.0_dp .and. size_right > 0.0_dp
      end if
   end subroutine scalar_stabilized_split_ok

   pure subroutine multi_stabilized_split_ok(forest, options, samples, var, threshold, missing_left, valid)
      type(grf_forest), intent(in) :: forest !! Local-coefficient/multi-causal forest supplying centered treatment columns.
      type(grf_options), intent(in) :: options !! Stabilized-split alpha and minimum-node-size controls.
      integer, intent(in) :: samples(:) !! One-based rows in the parent node used to define treatment-wise parent moments.
      integer, intent(in) :: var !! Predictor column defining the candidate split.
      real(dp), intent(in) :: threshold !! Candidate finite or NaN split threshold.
      logical, intent(in) :: missing_left !! Missing-value routing rule for finite thresholds.
      logical, intent(out) :: valid !! True when every treatment dimension meets parent-mean balance and variance constraints in both children.
      real(dp), allocatable :: sum_node_w(:)
      real(dp), allocatable :: sum_node_w_squared(:)
      real(dp), allocatable :: mean_node_w(:)
      real(dp), allocatable :: min_child_size(:)
      real(dp), allocatable :: sum_left_w(:)
      real(dp), allocatable :: sum_left_w_squared(:)
      real(dp), allocatable :: size_left(:)
      real(dp), allocatable :: size_right(:)
      integer, allocatable :: num_node_small_w(:)
      integer, allocatable :: num_left_small_w(:)
      integer :: i
      integer :: k
      integer :: row
      integer :: n_left
      integer :: n_right
      real(dp) :: weight_node
      real(dp) :: weight_left
      real(dp) :: weight_right
      real(dp) :: w_value
      real(dp) :: x_value
      logical :: goes_left

      allocate(sum_node_w(forest%n_treatments), sum_node_w_squared(forest%n_treatments))
      allocate(mean_node_w(forest%n_treatments), min_child_size(forest%n_treatments))
      allocate(sum_left_w(forest%n_treatments), sum_left_w_squared(forest%n_treatments))
      allocate(size_left(forest%n_treatments), size_right(forest%n_treatments))
      allocate(num_node_small_w(forest%n_treatments), num_left_small_w(forest%n_treatments))
      sum_node_w = 0.0_dp
      sum_node_w_squared = 0.0_dp
      num_node_small_w = 0
      weight_node = 0.0_dp
      do i = 1, size(samples)
         row = samples(i)
         weight_node = weight_node + forest%sample_weights(row)
         do k = 1, forest%n_treatments
            w_value = split_stabilizer_value(forest, row, k)
            sum_node_w(k) = sum_node_w(k) + forest%sample_weights(row) * w_value
            sum_node_w_squared(k) = sum_node_w_squared(k) + forest%sample_weights(row) * w_value * w_value
         end do
      end do
      if (weight_node <= 0.0_dp) then
         valid = .false.
         return
      end if
      mean_node_w = sum_node_w / weight_node
      min_child_size = options%alpha * max(0.0_dp, sum_node_w_squared - sum_node_w * sum_node_w / weight_node)
      do i = 1, size(samples)
         row = samples(i)
         do k = 1, forest%n_treatments
            if (split_stabilizer_value(forest, row, k) < mean_node_w(k)) then
               num_node_small_w(k) = num_node_small_w(k) + 1
            end if
         end do
      end do

      sum_left_w = 0.0_dp
      sum_left_w_squared = 0.0_dp
      num_left_small_w = 0
      weight_left = 0.0_dp
      n_left = 0
      do i = 1, size(samples)
         row = samples(i)
         x_value = forest%x(row,var)
         if (ieee_is_nan(threshold)) then
            goes_left = ieee_is_nan(x_value)
         else if (ieee_is_nan(x_value)) then
            goes_left = missing_left
         else
            goes_left = x_value <= threshold
         end if
         if (.not. goes_left) cycle
         n_left = n_left + 1
         weight_left = weight_left + forest%sample_weights(row)
         do k = 1, forest%n_treatments
            w_value = split_stabilizer_value(forest, row, k)
            sum_left_w(k) = sum_left_w(k) + forest%sample_weights(row) * w_value
            sum_left_w_squared(k) = sum_left_w_squared(k) + forest%sample_weights(row) * w_value * w_value
            if (w_value < mean_node_w(k)) num_left_small_w(k) = num_left_small_w(k) + 1
         end do
      end do
      n_right = size(samples) - n_left
      weight_right = weight_node - weight_left
      if (weight_left <= 0.0_dp .or. weight_right <= 0.0_dp) then
         valid = .false.
         return
      end if
      size_left = max(0.0_dp, sum_left_w_squared - sum_left_w * sum_left_w / weight_left)
      size_right = max(0.0_dp, (sum_node_w_squared - sum_left_w_squared) - &
                       (sum_node_w - sum_left_w) * (sum_node_w - sum_left_w) / weight_right)
      valid = .true.
      do k = 1, forest%n_treatments
         if (num_left_small_w(k) < options%min_node_size .or. &
             n_left - num_left_small_w(k) < options%min_node_size .or. &
             num_node_small_w(k) - num_left_small_w(k) < options%min_node_size .or. &
             n_right - num_node_small_w(k) + num_left_small_w(k) < options%min_node_size .or. &
             size_left(k) < min_child_size(k) .or. size_right(k) < min_child_size(k)) then
            valid = .false.
            return
         end if
      end do
      if (options%imbalance_penalty > 0.0_dp) then
         if (all(size_left <= 100.0_dp * epsilon(1.0_dp)) .or. &
             all(size_right <= 100.0_dp * epsilon(1.0_dp))) valid = .false.
      end if
   end subroutine multi_stabilized_split_ok

   pure real(dp) function split_stabilizer_value(forest, row, treatment_index) result(value)
      type(grf_forest), intent(in) :: forest !! Forest supplying centered treatment or instrument quantities for stabilized splitting.
      integer, intent(in) :: row !! One-based training-row index whose stabilization quantity is requested.
      integer, intent(in) :: treatment_index !! One-based treatment dimension; scalar causal/IV families use one.

      if (forest%family == grf_instrumental) then
         value = forest%z(row) - forest%z_hat(row)
      else
         value = forest%w(row,treatment_index) - forest%w_hat(row,treatment_index)
      end if
   end function split_stabilizer_value

   pure logical function survival_node_has_split_signal(forest, samples) result(valid)
      type(grf_forest), intent(in) :: forest !! Survival forest containing event indicators and relabeled failure-time indices.
      integer, intent(in) :: samples(:) !! One-based rows in the candidate parent node.
      integer :: first_time
      integer :: i

      valid = .false.
      first_time = -huge(1)
      do i = 1, size(samples)
         if (forest%event(samples(i)) /= 1) cycle
         if (first_time == -huge(1)) then
            first_time = forest%survival_time_index(samples(i))
         else if (forest%survival_time_index(samples(i)) /= first_time) then
            valid = .true.
            return
         end if
      end do
   end function survival_node_has_split_signal

   pure subroutine evaluate_survival_split(forest, options, samples, var, threshold, missing_left, score, valid)
      type(grf_forest), intent(in) :: forest !! Survival forest supplying times, events, weights, and relabeled failure-time indices.
      type(grf_options), intent(in) :: options !! Survival split controls; fast_logrank selects the accelerated approximation.
      integer, intent(in) :: samples(:) !! One-based rows in the parent node.
      integer, intent(in) :: var !! Predictor column defining the candidate split.
      real(dp), intent(in) :: threshold !! Candidate split threshold, or NaN for missing-versus-finite splitting.
      logical, intent(in) :: missing_left !! Missing-value routing rule for finite thresholds.
      real(dp), intent(out) :: score !! Nonnegative log-rank split statistic minus the imbalance penalty.
      logical, intent(out) :: valid !! True when both children contain enough failures and the statistic is defined.
      logical, allocatable :: left(:)
      real(dp) :: numerator
      real(dp) :: denominator
      real(dp) :: risk_total
      real(dp) :: risk_left
      real(dp) :: fail_total
      real(dp) :: fail_left
      real(dp) :: event_time
      real(dp) :: value
      real(dp) :: expected_left
      real(dp) :: fraction_left
      real(dp) :: penalty
      integer :: failures_left
      integer :: failures_right
      integer :: min_failures
      integer :: n_left
      integer :: n_right
      integer :: i
      integer :: j
      logical :: goes_left

      allocate(left(size(samples)))
      n_left = 0
      n_right = 0
      failures_left = 0
      failures_right = 0
      do i = 1, size(samples)
         value = forest%x(samples(i),var)
         if (ieee_is_nan(threshold)) then
            goes_left = ieee_is_nan(value)
         else if (ieee_is_nan(value)) then
            goes_left = missing_left
         else
            goes_left = value <= threshold
         end if
         left(i) = goes_left
         if (goes_left) then
            n_left = n_left + 1
            failures_left = failures_left + forest%event(samples(i))
         else
            n_right = n_right + 1
            failures_right = failures_right + forest%event(samples(i))
         end if
      end do
      min_failures = max(1, ceiling(options%alpha * real(size(samples),dp)))
      valid = n_left >= 1 .and. n_right >= 1 .and. failures_left >= min_failures .and. failures_right >= min_failures
      if (.not. valid) then
         score = -huge(1.0_dp)
         return
      end if

      numerator = 0.0_dp
      denominator = 0.0_dp
      if (options%fast_logrank) then
         do j = 1, forest%n_times
            event_time = forest%failure_times(j)
            risk_total = 0.0_dp
            fail_total = 0.0_dp
            risk_left = 0.0_dp
            do i = 1, size(samples)
               if (forest%y(samples(i),1) >= event_time) then
                  risk_total = risk_total + forest%sample_weights(samples(i))
                  if (left(i)) risk_left = risk_left + forest%sample_weights(samples(i))
               end if
               if (forest%event(samples(i)) == 1 .and. forest%survival_time_index(samples(i)) == j) then
                  fail_total = fail_total + forest%sample_weights(samples(i))
               end if
            end do
            if (risk_total > 100.0_dp * epsilon(1.0_dp)) then
               numerator = numerator + fail_total * risk_left / risk_total
            end if
         end do
         fail_left = 0.0_dp
         fail_total = 0.0_dp
         do i = 1, size(samples)
            if (forest%event(samples(i)) == 1) then
               fail_total = fail_total + forest%sample_weights(samples(i))
               if (left(i)) fail_left = fail_left + forest%sample_weights(samples(i))
            end if
         end do
         numerator = fail_left - numerator
         fraction_left = real(n_left,dp) / real(size(samples),dp)
         denominator = max(epsilon(1.0_dp), fail_total * fraction_left * (1.0_dp - fraction_left))
      else
         do j = 1, forest%n_times
            event_time = forest%failure_times(j)
            risk_total = 0.0_dp
            risk_left = 0.0_dp
            fail_total = 0.0_dp
            fail_left = 0.0_dp
            do i = 1, size(samples)
               if (forest%y(samples(i),1) >= event_time) then
                  risk_total = risk_total + forest%sample_weights(samples(i))
                  if (left(i)) risk_left = risk_left + forest%sample_weights(samples(i))
               end if
               if (forest%event(samples(i)) == 1 .and. forest%survival_time_index(samples(i)) == j) then
                  fail_total = fail_total + forest%sample_weights(samples(i))
                  if (left(i)) fail_left = fail_left + forest%sample_weights(samples(i))
               end if
            end do
            if (risk_total <= 100.0_dp * epsilon(1.0_dp)) cycle
            expected_left = fail_total * risk_left / risk_total
            numerator = numerator + fail_left - expected_left
            if (risk_total > 1.0_dp) then
               denominator = denominator + fail_total * (risk_left / risk_total) * (1.0_dp - risk_left / risk_total) * &
                  max(0.0_dp, (risk_total - fail_total) / max(risk_total - 1.0_dp, epsilon(1.0_dp)))
            end if
         end do
      end if
      if (denominator <= 100.0_dp * epsilon(1.0_dp)) then
         valid = .false.
         score = -huge(1.0_dp)
         return
      end if
      score = numerator * numerator / denominator
      penalty = options%imbalance_penalty * (1.0_dp / real(n_left,dp) + 1.0_dp / real(n_right,dp))
      score = score - penalty
   end subroutine evaluate_survival_split

   pure subroutine partition_samples(values, samples, threshold, missing_left, left_samples, right_samples)
      real(dp), intent(in) :: values(:) !! Full training predictor column indexed by the sample identifiers.
      integer, intent(in) :: samples(:) !! One-based row identifiers to partition.
      real(dp), intent(in) :: threshold !! Numeric threshold, or NaN for missing-versus-finite partitioning.
      logical, intent(in) :: missing_left !! Missing-value routing rule used with a finite threshold.
      integer, allocatable, intent(out) :: left_samples(:) !! Sample identifiers routed to the left child.
      integer, allocatable, intent(out) :: right_samples(:) !! Sample identifiers routed to the right child.
      integer, allocatable :: left_work(:)
      integer, allocatable :: right_work(:)
      integer :: nl
      integer :: nr
      integer :: i
      real(dp) :: value
      logical :: goes_left

      allocate(left_work(size(samples)), right_work(size(samples)))
      nl = 0
      nr = 0
      do i = 1, size(samples)
         value = values(samples(i))
         if (ieee_is_nan(threshold)) then
            goes_left = ieee_is_nan(value)
         else if (ieee_is_nan(value)) then
            goes_left = missing_left
         else
            goes_left = value <= threshold
         end if
         if (goes_left) then
            nl = nl + 1
            left_work(nl) = samples(i)
         else
            nr = nr + 1
            right_work(nr) = samples(i)
         end if
      end do
      allocate(left_samples(nl), right_samples(nr))
      if (nl > 0) left_samples = left_work(1:nl)
      if (nr > 0) right_samples = right_work(1:nr)
   end subroutine partition_samples

   pure integer function find_leaf_node(tree, row) result(node)
      type(grf_tree), intent(in) :: tree !! Trained tree whose split arrays define the traversal path.
      real(dp), intent(in) :: row(:) !! Predictor row ordered as in the tree's training matrix.
      integer :: var
      real(dp) :: value
      real(dp) :: threshold

      node = tree%root
      do while (tree%left(node) /= 0 .or. tree%right(node) /= 0)
         var = tree%split_var(node)
         if (var < 1 .or. var > size(row)) return
         value = row(var)
         threshold = tree%split_value(node)
         if (ieee_is_nan(threshold)) then
            if (ieee_is_nan(value)) then
               node = tree%left(node)
            else
               node = tree%right(node)
            end if
         else if (ieee_is_nan(value)) then
            if (tree%missing_left(node)) then
               node = tree%left(node)
            else
               node = tree%right(node)
            end if
         else if (value <= threshold) then
            node = tree%left(node)
         else
            node = tree%right(node)
         end if
         if (node == 0) exit
      end do
   end function find_leaf_node

   pure subroutine count_leaf_members(tree, x, samples, counts)
      type(grf_tree), intent(in) :: tree !! Tree structure used to route the supplied training rows.
      real(dp), intent(in) :: x(:,:) !! Full training predictor matrix.
      integer, intent(in) :: samples(:) !! One-based rows whose terminal-node membership is counted.
      integer, intent(out) :: counts(:) !! Per-node membership counts; only terminal nodes receive nonzero counts.
      integer :: i
      integer :: node

      counts = 0
      do i = 1, size(samples)
         node = find_leaf_node(tree, x(samples(i),:))
         if (node >= 1 .and. node <= size(counts)) counts(node) = counts(node) + 1
      end do
   end subroutine count_leaf_members

   recursive integer function prune_empty_subtrees(tree, node, counts) result(new_node)
      type(grf_tree), intent(inout) :: tree !! Tree whose empty honest leaves are pruned by child promotion.
      integer, intent(in) :: node !! Current node index in the recursive pruning walk.
      integer, intent(in) :: counts(:) !! Honest estimation-sample counts for terminal nodes before pruning.
      integer :: left_node
      integer :: right_node

      if (tree%left(node) == 0 .and. tree%right(node) == 0) then
         new_node = node
         return
      end if
      left_node = prune_empty_subtrees(tree, tree%left(node), counts)
      right_node = prune_empty_subtrees(tree, tree%right(node), counts)
      tree%left(node) = left_node
      tree%right(node) = right_node
      if (is_empty_leaf(tree, left_node, counts)) then
         new_node = right_node
      else if (is_empty_leaf(tree, right_node, counts)) then
         new_node = left_node
      else
         new_node = node
      end if
   end function prune_empty_subtrees

   pure logical function is_empty_leaf(tree, node, counts) result(empty)
      type(grf_tree), intent(in) :: tree !! Tree structure used to determine whether node is terminal.
      integer, intent(in) :: node !! Node index to test for terminal emptiness.
      integer, intent(in) :: counts(:) !! Honest estimation-sample counts corresponding to tree nodes.

      empty = tree%left(node) == 0 .and. tree%right(node) == 0 .and. counts(node) == 0
   end function is_empty_leaf

   pure subroutine populate_leaf_members(tree, x, samples)
      type(grf_tree), intent(inout) :: tree !! Tree receiving flattened honest estimation-sample memberships.
      real(dp), intent(in) :: x(:,:) !! Full training predictor matrix used to route estimation rows.
      integer, intent(in) :: samples(:) !! One-based honest estimation rows assigned to terminal nodes.
      integer, allocatable :: counts(:)
      integer, allocatable :: cursor(:)
      integer :: node
      integer :: i
      integer :: total

      allocate(counts(tree%n_nodes), cursor(tree%n_nodes))
      call count_leaf_members(tree, x, samples, counts)
      tree%leaf_start = 0
      tree%leaf_count = 0
      total = 0
      do node = 1, tree%n_nodes
         if (tree%left(node) == 0 .and. tree%right(node) == 0 .and. counts(node) > 0) then
            tree%leaf_start(node) = total + 1
            tree%leaf_count(node) = counts(node)
            total = total + counts(node)
         end if
      end do
      if (allocated(tree%leaf_samples)) deallocate(tree%leaf_samples)
      allocate(tree%leaf_samples(total))
      cursor = tree%leaf_start
      do i = 1, size(samples)
         node = find_leaf_node(tree, x(samples(i),:))
         if (node < 1 .or. node > tree%n_nodes) cycle
         if (tree%leaf_count(node) <= 0) cycle
         tree%leaf_samples(cursor(node)) = samples(i)
         cursor(node) = cursor(node) + 1
      end do
   end subroutine populate_leaf_members

   pure subroutine compute_forest_weights(forest, x_new, kernel_weights, oob)
      type(grf_forest), intent(in) :: forest !! Trained forest whose honest terminal memberships define GRF kernel weights.
      real(dp), intent(in) :: x_new(:,:) !! Query predictor matrix with columns ordered as in forest%x.
      real(dp), allocatable, intent(out) :: kernel_weights(:,:) !! Dense query-by-training matrix of normalized forest kernel weights.
      logical, intent(in), optional :: oob !! If true, query row i skips trees in which training observation i was in-bag.
      logical :: use_oob
      integer :: i
      integer :: t
      integer :: node
      integer :: start
      integer :: count_leaf
      integer :: j
      integer :: sample
      integer :: trees_used
      real(dp) :: leaf_weight_sum

      allocate(kernel_weights(size(x_new,1), forest%n_train))
      kernel_weights = 0.0_dp
      use_oob = .false.
      if (present(oob)) use_oob = oob
      do i = 1, size(x_new,1)
         trees_used = 0
         do t = 1, size(forest%trees)
            if (use_oob .and. i <= forest%n_train) then
               if (forest%trees(t)%inbag(i)) cycle
            end if
            node = find_leaf_node(forest%trees(t), x_new(i,:))
            if (node < 1 .or. node > forest%trees(t)%n_nodes) cycle
            count_leaf = forest%trees(t)%leaf_count(node)
            if (count_leaf <= 0) cycle
            start = forest%trees(t)%leaf_start(node)
            leaf_weight_sum = 0.0_dp
            do j = start, start + count_leaf - 1
               sample = forest%trees(t)%leaf_samples(j)
               leaf_weight_sum = leaf_weight_sum + forest%sample_weights(sample)
            end do
            if (leaf_weight_sum <= 0.0_dp) cycle
            trees_used = trees_used + 1
            do j = start, start + count_leaf - 1
               sample = forest%trees(t)%leaf_samples(j)
               kernel_weights(i,sample) = kernel_weights(i,sample) + forest%sample_weights(sample) / leaf_weight_sum
            end do
         end do
         if (trees_used > 0) kernel_weights(i,:) = kernel_weights(i,:) / real(trees_used,dp)
      end do
   end subroutine compute_forest_weights

end module grf_core
