! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
!
! Computational translations from ape R/gammaStat.R, R/yule.R,
! R/coalescent.intervals.R, R/ltt.plot.R, R/MoranI.R, R/mst.R, and src/delta_plot.c.
module ape_statistics
   use ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
   use r_kinds, only : dp
   use ape_types, only : phylo_tree, child_counts
   use ape_tree_algorithms, only : branching_times, node_depth_edgelength
   use ape_topology, only : is_binary_tree, is_ultrametric_tree
   use ape_tree_edit, only : has_singles, collapse_singles, multi2di
   implicit none
   private

   type, public :: moran_result
      real(dp) :: observed = 0.0_dp
      real(dp) :: expected = 0.0_dp
      real(dp) :: standard_deviation = 0.0_dp
      real(dp) :: p_value = 0.0_dp
   end type moran_result

   type, public :: yule_result
      real(dp) :: lambda = 0.0_dp
      real(dp) :: standard_error = 0.0_dp
      real(dp) :: log_likelihood = 0.0_dp
   end type yule_result

   public :: gamma_stat
   public :: yule_fit
   public :: coalescent_intervals
   public :: ltt_coordinates
   public :: moran_i
   public :: minimum_spanning_tree
   public :: delta_plot_statistics

contains

   pure subroutine gamma_stat(tree, gamma, info)
      !! Computes the Pybus-Harvey gamma statistic from ape-style branching times.
      type(phylo_tree), intent(in) :: tree !! Rooted tree with branch lengths; a fully bifurcating ultrametric tree is expected.
      real(dp), intent(out) :: gamma !! Pybus-Harvey gamma statistic.
      integer, intent(out) :: info !! Status code: zero on success, nonzero when the statistic is undefined.
      real(dp), allocatable :: bt(:)
      real(dp), allocatable :: g(:)
      real(dp) :: cumulative
      real(dp) :: mean_stat
      real(dp) :: scale
      real(dp) :: stat
      real(dp) :: st
      integer :: i
      integer :: n

      call branching_times(tree, bt, info)
      gamma = 0.0_dp
      if (info /= 0) return
      n = tree%n_tip
      if (n < 3 .or. size(bt) /= n - 1) then
         info = 1
         return
      end if
      call sort_ascending(bt)
      allocate(g(n - 1))
      g(n - 1) = bt(1)
      do i = 2, n - 1
         g(n - i) = bt(i) - bt(i - 1)
      end do
      st = 0.0_dp
      do i = 1, n - 1
         st = st + real(i + 1, dp) * g(i)
      end do
      if (n == 3) then
         stat = 2.0_dp * g(1)
      else
         cumulative = 0.0_dp
         stat = 0.0_dp
         do i = 1, n - 2
            cumulative = cumulative + real(i + 1, dp) * g(i)
            stat = stat + cumulative
         end do
         stat = stat / real(n - 2, dp)
      end if
      mean_stat = 0.5_dp * st
      if (n <= 2 .or. abs(st) <= tiny(1.0_dp)) then
         info = 2
         return
      end if
      scale = st * sqrt(1.0_dp / (12.0_dp * real(n - 2, dp)))
      if (abs(scale) <= tiny(1.0_dp)) then
         info = 3
         return
      end if
      gamma = (stat - mean_stat) / scale
   end subroutine gamma_stat

   pure subroutine yule_fit(tree, result, info, root_edge, use_root_edge)
      !! Fits ape's constant-rate Yule model by its closed-form maximum likelihood estimate.
      type(phylo_tree), intent(in) :: tree !! Rooted dichotomous tree with branch lengths.
      type(yule_result), intent(out) :: result !! Estimated rate, standard error, and log-likelihood.
      integer, intent(out) :: info !! Status code: zero on success, nonzero for invalid topology or branch lengths.
      real(dp), intent(in), optional :: root_edge !! Optional root-edge length to include when `use_root_edge` is true.
      logical, intent(in), optional :: use_root_edge !! Whether the supplied root edge contributes to exposure and event count.
      integer, allocatable :: nchild(:)
      integer :: events
      integer :: node
      real(dp) :: exposure
      logical :: include_root

      result = yule_result()
      info = 0
      if (.not. tree%valid() .or. .not. tree%has_lengths()) then
         info = 1
         return
      end if
      nchild = child_counts(tree)
      do node = tree%n_tip + 1, tree%total_nodes()
         if (nchild(node) /= 2) then
            info = 2
            return
         end if
      end do
      exposure = sum(tree%edge_length)
      events = tree%n_node
      include_root = .false.
      if (present(use_root_edge)) include_root = use_root_edge
      if (include_root .and. present(root_edge)) then
         if (root_edge < 0.0_dp) then
            info = 3
            return
         end if
         exposure = exposure + root_edge
      else
         events = events - 1
      end if
      if (events <= 0 .or. exposure <= 0.0_dp) then
         info = 4
         return
      end if
      result%lambda = real(events, dp) / exposure
      result%standard_error = result%lambda / sqrt(real(events, dp))
      result%log_likelihood = -result%lambda * exposure + log_gamma(real(tree%n_node + 1, dp)) &
         + real(events, dp) * log(result%lambda)
   end subroutine yule_fit

   pure subroutine coalescent_intervals(tree, lineages, widths, total_depth, info)
      !! Constructs ape coalescent interval widths and lineage counts from a binary tree.
      type(phylo_tree), intent(in) :: tree !! Rooted binary tree with branch lengths.
      integer, allocatable, intent(out) :: lineages(:) !! Number of lineages at the beginning of each interval.
      real(dp), allocatable, intent(out) :: widths(:) !! Successive ordered branching-time interval widths.
      real(dp), intent(out) :: total_depth !! Sum of all returned interval widths.
      integer, intent(out) :: info !! Status code: zero on success, nonzero if branching times are unavailable.
      real(dp), allocatable :: times(:)
      integer :: i
      integer :: n

      call branching_times(tree, times, info)
      n = size(times)
      allocate(lineages(n), widths(n))
      widths = 0.0_dp
      lineages = 0
      total_depth = 0.0_dp
      if (info /= 0) return
      call sort_ascending(times)
      if (n == 0) return
      widths(1) = times(1)
      do i = 2, n
         widths(i) = times(i) - times(i - 1)
      end do
      do i = 1, n
         lineages(i) = n + 2 - i
      end do
      total_depth = sum(widths)
   end subroutine coalescent_intervals

   pure subroutine ltt_coordinates(tree, time, lineages, info, backward, tolerance, step_type)
      !! Computes the numerical coordinates used by ape `ltt.plot.coords` without plotting.
      type(phylo_tree), intent(in) :: tree !! Rooted tree with lengths; singleton nodes are collapsed first.
      real(dp), allocatable, intent(out) :: time(:) !! Event times backward from latest tip, or shifted to start at zero.
      integer, allocatable, intent(out) :: lineages(:) !! Number of lineages associated with each returned time coordinate.
      integer, intent(out) :: info !! Status code: zero on success; nonzero for invalid trees, missing lengths, or invalid options.
      logical, intent(in), optional :: backward !! If true (default), latest-tip time is zero and older events are negative.
      real(dp), intent(in), optional :: tolerance !! Nonnegative ultrametric/present-event tolerance; default `1e-6`.
      character(len=*), intent(in), optional :: step_type !! Step convention `S` (default) or `s`, matching ape coordinate shifting.
      type(phylo_tree) :: work
      type(phylo_tree) :: resolved
      real(dp), allocatable :: branching(:)
      real(dp), allocatable :: depth(:)
      real(dp), allocatable :: event_time(:)
      real(dp), allocatable :: kept_time(:)
      integer, allocatable :: counts(:)
      integer, allocatable :: event(:)
      integer, allocatable :: kept_event(:)
      integer, allocatable :: order(:)
      character(len=1) :: convention
      integer :: i
      integer :: j
      integer :: k
      integer :: node
      integer :: root_node
      integer :: status
      integer :: tmp_event
      integer :: tmp_order
      real(dp) :: present_time
      real(dp) :: tmp_time
      real(dp) :: tol
      logical :: use_backward

      info = 0
      allocate(time(0), lineages(0))
      if (.not. tree%valid() .or. .not. tree%has_lengths() .or. tree%root() == 0) then
         info = 1
         return
      end if
      tol = 1.0e-6_dp
      if (present(tolerance)) tol = tolerance
      if (tol < 0.0_dp) then
         info = 2
         return
      end if
      use_backward = .true.
      if (present(backward)) use_backward = backward
      convention = 'S'
      if (present(step_type)) then
         if (len_trim(step_type) /= 1) then
            info = 2
            return
         end if
         convention = step_type(1:1)
      end if
      if (convention /= 'S' .and. convention /= 's') then
         info = 2
         return
      end if

      work = tree
      if (has_singles(work)) then
         call collapse_singles(work, resolved, status)
         if (status /= 0) then
            info = 3
            return
         end if
         work = resolved
      end if

      if (is_ultrametric_tree(work, tol)) then
         call branching_times(work, branching, status)
         if (status /= 0) then
            info = 3
            return
         end if
         allocate(order(work%n_node), counts(work%total_nodes()))
         counts = child_counts(work)
         do i = 1, work%n_node
            order(i) = i
         end do
         do i = 2, work%n_node
            tmp_time = branching(i)
            tmp_order = order(i)
            j = i - 1
            do while (j >= 1)
               if (branching(j) >= tmp_time) exit
               branching(j + 1) = branching(j)
               order(j + 1) = order(j)
               j = j - 1
            end do
            branching(j + 1) = tmp_time
            order(j + 1) = tmp_order
         end do
         deallocate(time, lineages)
         allocate(time(work%n_node + 1), lineages(work%n_node + 1))
         do i = 1, work%n_node
            time(i) = -branching(i)
         end do
         time(work%n_node + 1) = 0.0_dp
         lineages(1) = 1
         do i = 1, work%n_node
            node = work%n_tip + order(i)
            lineages(i + 1) = lineages(i) + counts(node) - 1
         end do
      else
         if (.not. is_binary_tree(work)) then
            call multi2di(work, resolved, status)
            if (status /= 0) then
               info = 3
               return
            end if
            work = resolved
         end if
         call node_depth_edgelength(work, depth, status)
         if (status /= 0) then
            info = 3
            return
         end if
         present_time = maxval(depth(1:work%n_tip))
         root_node = work%root()
         allocate(event_time(work%total_nodes()), event(work%total_nodes()))
         event_time = depth
         event = 1
         event(1:work%n_tip) = -1
         event(root_node) = 1
         k = count(present_time - event_time > tol)
         allocate(kept_time(k), kept_event(k))
         j = 0
         do i = 1, work%total_nodes()
            if (present_time - event_time(i) <= tol) cycle
            j = j + 1
            kept_time(j) = event_time(i)
            kept_event(j) = event(i)
         end do
         do i = 2, k
            tmp_time = kept_time(i)
            tmp_event = kept_event(i)
            j = i - 1
            do while (j >= 1)
               if (kept_time(j) <= tmp_time) exit
               kept_time(j + 1) = kept_time(j)
               kept_event(j + 1) = kept_event(j)
               j = j - 1
            end do
            kept_time(j + 1) = tmp_time
            kept_event(j + 1) = tmp_event
         end do
         deallocate(time, lineages)
         allocate(time(k + 1), lineages(k + 1))
         if (k > 0) time(1:k) = kept_time - present_time
         time(k + 1) = 0.0_dp
         lineages(1) = 1
         do i = 1, k
            lineages(i + 1) = lineages(i) + kept_event(i)
         end do
      end if

      if (.not. use_backward) time = time - time(1)
      if (convention == 's' .and. size(lineages) > 1) then
         lineages(1:size(lineages) - 1) = lineages(2:size(lineages))
      end if
   end subroutine ltt_coordinates

   pure subroutine moran_i(x, weight, result, info, scaled, remove_nan, alternative)
      !! Computes Moran's I, its randomization standard deviation, and a normal-approximation p-value.
      real(dp), intent(in) :: x(:) !! Observation vector; NaNs may be removed when `remove_nan` is true.
      real(dp), intent(in) :: weight(:, :) !! Square nonnegative spatial or taxonomic weight matrix.
      type(moran_result), intent(out) :: result !! Observed/expected Moran's I, standard deviation, and p-value.
      integer, intent(out) :: info !! Status code: zero on success, nonzero for invalid dimensions or degenerate data.
      logical, intent(in), optional :: scaled !! If true, apply ape's optional upper-bound scaling to the observed index.
      logical, intent(in), optional :: remove_nan !! If true, remove observations with NaN values and matching weight rows/columns.
      character(len=*), intent(in), optional :: alternative !! Tail choice: `two.sided`, `less`, or `greater`.
      real(dp), allocatable :: xx(:)
      real(dp), allocatable :: w(:, :)
      real(dp), allocatable :: row_sum(:)
      real(dp), allocatable :: col_sum(:)
      real(dp), allocatable :: y(:)
      integer, allocatable :: keep(:)
      integer :: i
      integer :: j
      integer :: kkeep
      integer :: n
      real(dp) :: cv
      real(dp) :: kurtosis_factor
      real(dp) :: mean_x
      real(dp) :: numerator
      real(dp) :: s
      real(dp) :: s1
      real(dp) :: s2
      real(dp) :: s_sq
      real(dp) :: variance_sum
      real(dp) :: z_cdf
      real(dp) :: i_max
      logical :: do_remove
      logical :: do_scale
      character(len=16) :: tail

      result = moran_result()
      info = 0
      n = size(x)
      if (size(weight, 1) /= n .or. size(weight, 2) /= n .or. n < 4) then
         info = 1
         return
      end if
      do_remove = .false.
      if (present(remove_nan)) do_remove = remove_nan
      kkeep = 0
      allocate(keep(n))
      do i = 1, n
         if (ieee_is_nan(x(i))) then
            if (.not. do_remove) then
               result%expected = -1.0_dp / real(n - 1, dp)
               result%observed = ieee_value(0.0_dp, ieee_quiet_nan)
               result%standard_deviation = ieee_value(0.0_dp, ieee_quiet_nan)
               result%p_value = ieee_value(0.0_dp, ieee_quiet_nan)
               info = 2
               return
            end if
         else
            kkeep = kkeep + 1
            keep(kkeep) = i
         end if
      end do
      if (kkeep < 4) then
         info = 3
         return
      end if
      allocate(xx(kkeep), w(kkeep, kkeep), row_sum(kkeep), col_sum(kkeep), y(kkeep))
      do i = 1, kkeep
         xx(i) = x(keep(i))
         do j = 1, kkeep
            w(i, j) = weight(keep(i), keep(j))
         end do
      end do
      n = kkeep
      result%expected = -1.0_dp / real(n - 1, dp)
      do i = 1, n
         row_sum(i) = sum(w(i, :))
         if (abs(row_sum(i)) <= tiny(1.0_dp)) row_sum(i) = 1.0_dp
         w(i, :) = w(i, :) / row_sum(i)
      end do
      row_sum = sum(w, dim=2)
      col_sum = sum(w, dim=1)
      s = sum(w)
      mean_x = sum(xx) / real(n, dp)
      y = xx - mean_x
      variance_sum = sum(y * y)
      if (abs(s) <= tiny(1.0_dp) .or. variance_sum <= tiny(1.0_dp)) then
         info = 4
         return
      end if
      cv = 0.0_dp
      do i = 1, n
         do j = 1, n
            cv = cv + w(i, j) * y(i) * y(j)
         end do
      end do
      result%observed = real(n, dp) * cv / (s * variance_sum)
      do_scale = .false.
      if (present(scaled)) do_scale = scaled
      if (do_scale) then
         i_max = real(n, dp) / s * sample_sd(row_sum * y) / sqrt(variance_sum / real(n - 1, dp))
         if (abs(i_max) <= tiny(1.0_dp)) then
            info = 5
            return
         end if
         result%observed = result%observed / i_max
      end if
      s1 = 0.0_dp
      do i = 1, n
         do j = 1, n
            s1 = s1 + (w(i, j) + w(j, i))**2
         end do
      end do
      s1 = 0.5_dp * s1
      s2 = sum((row_sum + col_sum)**2)
      s_sq = s * s
      kurtosis_factor = (sum(y**4) / real(n, dp)) / (variance_sum / real(n, dp))**2
      numerator = real(n, dp) * (real(n * n - 3 * n + 3, dp) * s1 - real(n, dp) * s2 + 3.0_dp * s_sq) &
         - kurtosis_factor * (real(n * (n - 1), dp) * s1 - 2.0_dp * real(n, dp) * s2 + 6.0_dp * s_sq)
      numerator = numerator / (real((n - 1) * (n - 2) * (n - 3), dp) * s_sq) &
         - 1.0_dp / real((n - 1)**2, dp)
      if (numerator <= 0.0_dp) then
         info = 6
         return
      end if
      result%standard_deviation = sqrt(numerator)
      z_cdf = normal_cdf((result%observed - result%expected) / result%standard_deviation)
      tail = 'two.sided'
      if (present(alternative)) tail = trim(adjustl(alternative))
      select case (trim(tail))
      case ('two.sided')
         if (result%observed <= result%expected) then
            result%p_value = 2.0_dp * z_cdf
         else
            result%p_value = 2.0_dp * (1.0_dp - z_cdf)
         end if
         result%p_value = min(1.0_dp, result%p_value)
      case ('less')
         result%p_value = z_cdf
      case ('greater')
         result%p_value = 1.0_dp - z_cdf
      case default
         info = 7
      end select
   end subroutine moran_i

   pure subroutine minimum_spanning_tree(distance, adjacency, total_weight, info)
      !! Builds a minimum spanning tree by Prim's algorithm.
      real(dp), intent(in) :: distance(:, :) !! Symmetric finite distance matrix; diagonal entries are ignored.
      integer, allocatable, intent(out) :: adjacency(:, :) !! Symmetric zero-one adjacency matrix of the spanning tree.
      real(dp), intent(out) :: total_weight !! Sum of selected edge weights.
      integer, intent(out) :: info !! Status code: zero on success, nonzero for invalid or disconnected input.
      logical, allocatable :: in_tree(:)
      integer :: best_i
      integer :: best_j
      integer :: i
      integer :: j
      integer :: n
      integer :: step
      real(dp) :: best

      n = size(distance, 1)
      allocate(adjacency(n, size(distance, 2)))
      adjacency = 0
      total_weight = 0.0_dp
      info = 0
      if (n < 1 .or. size(distance, 2) /= n) then
         info = 1
         return
      end if
      allocate(in_tree(n))
      in_tree = .false.
      in_tree(1) = .true.
      do step = 1, n - 1
         best = huge(1.0_dp)
         best_i = 0
         best_j = 0
         do i = 1, n
            if (.not. in_tree(i)) cycle
            do j = 1, n
               if (in_tree(j) .or. i == j) cycle
               if (ieee_is_nan(distance(i, j))) cycle
               if (distance(i, j) < best) then
                  best = distance(i, j)
                  best_i = i
                  best_j = j
               end if
            end do
         end do
         if (best_i == 0) then
            info = 2
            return
         end if
         adjacency(best_i, best_j) = 1
         adjacency(best_j, best_i) = 1
         total_weight = total_weight + best
         in_tree(best_j) = .true.
      end do
   end subroutine minimum_spanning_tree

   pure subroutine delta_plot_statistics(distance, nbins, counts, delta_bar, info)
      !! Computes Holland et al. quartet delta histogram counts and per-taxon mean delta values.
      real(dp), intent(in) :: distance(:, :) !! Symmetric taxon-distance matrix.
      integer, intent(in) :: nbins !! Number of equal-width histogram bins over delta in `[0,1]`; must be positive.
      integer, allocatable, intent(out) :: counts(:) !! Histogram counts; delta exactly one is included in the final bin.
      real(dp), allocatable, intent(out) :: delta_bar(:) !! Mean quartet delta involving each taxon.
      integer, intent(out) :: info !! Status code: zero on success, nonzero for invalid dimensions or too few taxa.
      integer, allocatable :: involvements(:)
      integer :: bin
      integer :: n
      integer :: u
      integer :: v
      integer :: x
      integer :: y
      real(dp) :: a
      real(dp) :: b
      real(dp) :: c
      real(dp) :: delta

      n = size(distance, 1)
      allocate(counts(max(0, nbins)), delta_bar(n), involvements(n))
      counts = 0
      delta_bar = 0.0_dp
      involvements = 0
      info = 0
      if (size(distance, 2) /= n .or. n < 4 .or. nbins <= 0) then
         info = 1
         return
      end if
      do x = 1, n - 3
         do y = x + 1, n - 2
            do u = y + 1, n - 1
               do v = u + 1, n
                  a = distance(x, v) + distance(y, u)
                  b = distance(x, u) + distance(y, v)
                  c = distance(x, y) + distance(u, v)
                  delta = quartet_delta(a, b, c)
                  bin = min(nbins, int(delta * real(nbins, dp)) + 1)
                  counts(bin) = counts(bin) + 1
                  delta_bar(x) = delta_bar(x) + delta
                  delta_bar(y) = delta_bar(y) + delta
                  delta_bar(u) = delta_bar(u) + delta
                  delta_bar(v) = delta_bar(v) + delta
                  involvements(x) = involvements(x) + 1
                  involvements(y) = involvements(y) + 1
                  involvements(u) = involvements(u) + 1
                  involvements(v) = involvements(v) + 1
               end do
            end do
         end do
      end do
      do x = 1, n
         if (involvements(x) > 0) delta_bar(x) = delta_bar(x) / real(involvements(x), dp)
      end do
   end subroutine delta_plot_statistics

   pure elemental real(dp) function quartet_delta(a, b, c) result(delta)
      real(dp), intent(in) :: a !! First four-point sum for a quartet.
      real(dp), intent(in) :: b !! Second four-point sum for a quartet.
      real(dp), intent(in) :: c !! Third four-point sum for a quartet.
      real(dp) :: values(3)
      real(dp) :: denominator

      values = [a, b, c]
      call sort_three(values)
      denominator = values(3) - values(1)
      if (abs(denominator) <= epsilon(1.0_dp) * max(1.0_dp, maxval(abs(values)))) then
         delta = 0.0_dp
      else
         delta = (values(3) - values(2)) / denominator
      end if
   end function quartet_delta

   pure subroutine sort_three(values)
      real(dp), intent(inout) :: values(3) !! Three values sorted in ascending order in place.
      real(dp) :: tmp

      if (values(1) > values(2)) then
         tmp = values(1)
         values(1) = values(2)
         values(2) = tmp
      end if
      if (values(2) > values(3)) then
         tmp = values(2)
         values(2) = values(3)
         values(3) = tmp
      end if
      if (values(1) > values(2)) then
         tmp = values(1)
         values(1) = values(2)
         values(2) = tmp
      end if
   end subroutine sort_three

   pure subroutine sort_ascending(values)
      real(dp), intent(inout) :: values(:) !! Real vector sorted in ascending order in place.
      integer :: i
      integer :: j
      real(dp) :: key

      do i = 2, size(values)
         key = values(i)
         j = i - 1
         do while (j >= 1)
            if (values(j) <= key) exit
            values(j + 1) = values(j)
            j = j - 1
         end do
         values(j + 1) = key
      end do
   end subroutine sort_ascending

   pure real(dp) function sample_sd(values) result(sd)
      real(dp), intent(in) :: values(:) !! Sample whose standard deviation uses denominator `n-1`.
      real(dp) :: mean_value

      if (size(values) <= 1) then
         sd = 0.0_dp
         return
      end if
      mean_value = sum(values) / real(size(values), dp)
      sd = sqrt(sum((values - mean_value)**2) / real(size(values) - 1, dp))
   end function sample_sd

   pure elemental real(dp) function normal_cdf(z) result(value)
      real(dp), intent(in) :: z !! Standard-normal variate at which the cumulative distribution is evaluated.

      value = 0.5_dp * erfc(-z / sqrt(2.0_dp))
   end function normal_cdf

end module ape_statistics
