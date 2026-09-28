module pdfcluster_graph
   use pdfcluster_kinds, only : dp
   use geometry, only : delaunayn_points, delaunay_result
   use pdfcluster_types, only : kepdf_result, density_profile_result
   use pdfcluster_stats, only : sample_sd_columns, quantile_type7
   use pdfcluster_kde, only : kepdf_estimate
   implicit none
   private

   public :: graph_unidimensional, graph_delaunay, graph_pairs
   public :: valley_measure, connected_components, density_profile, core_labels_from_profile

contains

   pure subroutine graph_unidimensional(x, adjacency)
      real(dp), intent(in) :: x(:) !! Scalar observations whose sorted nearest neighbours define the graph.
      logical, allocatable, intent(out) :: adjacency(:, :) !! Symmetric n-by-n adjacency matrix including self edges.
      integer, allocatable :: order(:)
      integer :: i
      integer :: j
      integer :: key
      integer :: n

      n = size(x)
      allocate(adjacency(n, n), order(n))
      adjacency = .false.
      order = [(i, i=1, n)]
      do i = 2, n
         key = order(i)
         j = i - 1
         do while (j >= 1)
            if (x(order(j)) <= x(key)) exit
            order(j + 1) = order(j)
            j = j - 1
         end do
         order(j + 1) = key
      end do
      do i = 1, n
         adjacency(i, i) = .true.
      end do
      do i = 1, n - 1
         adjacency(order(i), order(i + 1)) = .true.
         adjacency(order(i + 1), order(i)) = .true.
      end do
   end subroutine graph_unidimensional

   pure subroutine graph_delaunay(x, adjacency, success)
      real(dp), intent(in) :: x(:, :) !! Multivariate observations used to build the standardized Delaunay graph.
      logical, allocatable, intent(out) :: adjacency(:, :) !! Symmetric neighbour graph induced by Delaunay simplices.
      logical, intent(out) :: success !! True when a nonempty Delaunay triangulation is available.
      type(delaunay_result) :: tri
      real(dp), allocatable :: xs(:, :)
      real(dp), allocatable :: sd(:)
      real(dp) :: mu
      integer :: a
      integer :: b
      integer :: d
      integer :: i
      integer :: j
      integer :: n
      integer :: t

      n = size(x, 1)
      d = size(x, 2)
      allocate(adjacency(n, n))
      adjacency = .false.
      success = .false.
      if (n <= 0 .or. d <= 0 .or. n < d + 1) return
      allocate(xs(n, d), sd(d))
      call sample_sd_columns(x, sd)
      if (any(sd <= 0.0_dp)) return
      do j = 1, d
         mu = sum(x(:, j))/real(n, dp)
         xs(:, j) = (x(:, j) - mu)/sd(j)
      end do
      call delaunayn_points(xs, tri)
      if (.not. tri%success) return
      do i = 1, n
         adjacency(i, i) = .true.
      end do
      do t = 1, size(tri%tri, 1)
         do a = 1, size(tri%tri, 2)
            i = tri%tri(t, a)
            do b = 1, size(tri%tri, 2)
               j = tri%tri(t, b)
               adjacency(i, j) = .true.
            end do
         end do
      end do
      success = .true.
   end subroutine graph_delaunay

   pure subroutine valley_measure(f, area, filled)
      real(dp), intent(in) :: f(:) !! Density profile sampled at equally spaced points along a segment.
      real(dp), intent(out) :: area !! Upstream valley index: largest filled valley mass divided by final filled mass.
      real(dp), allocatable, intent(out), optional :: filled(:) !! Final profile after iterative valley filling.
      real(dp), allocatable :: f0(:)
      real(dp), allocatable :: ftemp(:)
      real(dp) :: best
      real(dp) :: inc
      real(dp) :: level
      integer, allocatable :: antimodes(:)
      integer, allocatable :: modes(:)
      integer :: am
      integer :: i
      integer :: left
      integer :: na
      integer :: nm
      integer :: right
      integer :: u
      logical :: changed

      area = 0.0_dp
      if (size(f) == 0) then
         if (present(filled)) allocate(filled(0))
         return
      end if
      allocate(f0(size(f)), ftemp(size(f)), antimodes(size(f)), modes(size(f)))
      f0 = f
      best = 0.0_dp
      do
         ftemp = f0
         na = 0
         nm = 1
         modes(1) = 1
         do i = 1, size(f) - 2
            if ((ftemp(i + 1) - ftemp(i))*(ftemp(i + 2) - ftemp(i + 1)) <= 0.0_dp) then
               if (ftemp(i + 1) - ftemp(i) <= 0.0_dp) then
                  na = na + 1
                  antimodes(na) = i + 1
               else
                  nm = nm + 1
                  modes(nm) = i + 1
               end if
            end if
         end do
         if (size(f) > 1) then
            nm = nm + 1
            modes(nm) = size(f)
         end if
         do am = 1, na
            left = modes(1)
            do i = 1, nm
               if (modes(i) < antimodes(am)) left = max(left, modes(i))
            end do
            right = modes(nm)
            do i = nm, 1, -1
               if (modes(i) > antimodes(am)) right = min(right, modes(i))
            end do
            level = min(ftemp(left), ftemp(right))
            inc = 0.0_dp
            do u = left, right - 1
               if (ftemp(u) <= level) then
                  f0(u) = level
                  inc = inc + f0(u) - f(u)
               end if
            end do
            best = max(best, inc)
         end do
         changed = any(abs(f0 - ftemp) > 0.0_dp)
         if (.not. changed) exit
      end do
      if (sum(f0) > 0.0_dp) area = best/sum(f0)
      if (present(filled)) then
         allocate(filled(size(f0)))
         filled = f0
      end if
   end subroutine valley_measure

   pure subroutine graph_pairs(x, pdf, lambda, grid_pairs, adjacency, pair_areas, pair_indices, success)
      real(dp), intent(in) :: x(:, :) !! Observations whose pairwise density valleys determine graph edges.
      type(kepdf_result), intent(in) :: pdf !! Density model fitted to all observations.
      real(dp), intent(in) :: lambda !! Maximum valley index for retaining an edge between a pair of observations.
      integer, intent(in) :: grid_pairs !! Number of equally spaced profile points sampled along each pairwise segment.
      logical, allocatable, intent(out) :: adjacency(:, :) !! Symmetric pair graph including self edges.
      real(dp), allocatable, intent(out) :: pair_areas(:) !! Valley index for every unordered observation pair.
      integer, allocatable, intent(out) :: pair_indices(:, :) !! Two-by-Npair matrix of one-based observation indices.
      logical, intent(out) :: success !! True when every segment density profile is evaluated successfully.
      real(dp), allocatable :: eval_points(:, :)
      real(dp), allocatable :: f(:)
      real(dp) :: alpha
      integer :: d
      integer :: g
      integer :: i
      integer :: j
      integer :: k
      integer :: n
      integer :: np
      logical :: ok

      success = .false.
      n = size(x, 1)
      d = size(x, 2)
      np = n*(n - 1)/2
      allocate(adjacency(n, n), pair_areas(np), pair_indices(2, np))
      adjacency = .false.
      pair_areas = 0.0_dp
      do i = 1, n
         adjacency(i, i) = .true.
      end do
      if (grid_pairs < 2 .or. d <= 0) return
      allocate(eval_points(grid_pairs, d))
      k = 0
      do i = 1, n - 1
         do j = i + 1, n
            k = k + 1
            pair_indices(:, k) = [i, j]
            do g = 1, grid_pairs
               alpha = real(g - 1, dp)/real(grid_pairs - 1, dp)
               eval_points(g, :) = alpha*x(i, :) + (1.0_dp - alpha)*x(j, :)
            end do
            if (trim(pdf%bwtype) == "adaptive") then
               call kepdf_estimate(x, eval_points, hx=pdf%hx, kernel=pdf%kernel, estimate=f, success=ok)
            else
               call kepdf_estimate(x, eval_points, h=pdf%h, kernel=pdf%kernel, estimate=f, success=ok)
            end if
            if (.not. ok) return
            call valley_measure(f, pair_areas(k))
            if (pair_areas(k) <= lambda) then
               adjacency(i, j) = .true.
               adjacency(j, i) = .true.
            end if
         end do
      end do
      success = .true.
   end subroutine graph_pairs

   pure subroutine connected_components(adjacency, active, comp, ncomp)
      logical, intent(in) :: adjacency(:, :) !! Symmetric graph adjacency matrix.
      logical, intent(in) :: active(:) !! Mask of vertices included in the induced graph.
      integer, allocatable, intent(out) :: comp(:) !! Component label for each vertex; zero denotes an inactive vertex.
      integer, intent(out) :: ncomp !! Number of connected components among active vertices.
      integer, allocatable :: stack(:)
      integer :: i
      integer :: j
      integer :: sp
      integer :: v

      allocate(comp(size(active)), stack(size(active)))
      comp = 0
      ncomp = 0
      if (size(adjacency, 1) /= size(active) .or. size(adjacency, 2) /= size(active)) return
      do i = 1, size(active)
         if (.not. active(i) .or. comp(i) /= 0) cycle
         ncomp = ncomp + 1
         sp = 1
         stack(1) = i
         comp(i) = ncomp
         do while (sp > 0)
            v = stack(sp)
            sp = sp - 1
            do j = 1, size(active)
               if (.not. active(j) .or. comp(j) /= 0) cycle
               if (adjacency(v, j)) then
                  comp(j) = ncomp
                  sp = sp + 1
                  stack(sp) = j
               end if
            end do
         end do
      end do
   end subroutine connected_components

   pure subroutine density_profile(adjacency, density, profile_grid, correct, result)
      logical, intent(in) :: adjacency(:, :) !! Graph used to define connected density level sets.
      real(dp), intent(in) :: density(:) !! Density estimate at each graph vertex.
      integer, intent(in) :: profile_grid !! Number of interior probability levels between zero and one.
      logical, intent(in) :: correct !! Whether singleton connected components are removed as in upstream pdfCluster.
      type(density_profile_result), intent(out) :: result !! Component-count and membership profile over density levels.
      logical, allocatable :: active(:)
      integer, allocatable :: comp(:)
      integer, allocatable :: count(:)
      integer :: c
      integer :: i
      integer :: k
      integer :: n
      integer :: ncomp
      integer :: nlev
      integer :: newlabel

      n = size(density)
      nlev = max(0, profile_grid) + 2
      allocate(result%p(nlev), result%q(nlev), result%nc(nlev), result%ids(n, nlev))
      allocate(active(n), count(max(1, n)))
      result%ids = -1
      result%nc = 0
      do k = 1, nlev
         result%p(k) = real(k - 1, dp)/real(max(1, nlev - 1), dp)
      end do
      if (n == 0) then
         result%q = 0.0_dp
         return
      end if
      result%q(1) = maxval(density)
      result%q(nlev) = minval(density)
      result%ids(:, 1) = -1
      result%ids(:, nlev) = 1
      result%nc(1) = 0
      result%nc(nlev) = 1
      do k = 2, nlev - 1
         result%q(k) = quantile_type7(density, 1.0_dp - result%p(k))
         active = density >= result%q(k)
         call connected_components(adjacency, active, comp, ncomp)
         count = 0
         do i = 1, n
            if (comp(i) > 0) count(comp(i)) = count(comp(i)) + 1
         end do
         result%ids(:, k) = -1
         newlabel = 0
         do c = 1, ncomp
            if (correct .and. count(c) <= 1) cycle
            newlabel = newlabel + 1
            do i = 1, n
               if (comp(i) == c) result%ids(i, k) = newlabel
            end do
         end do
         result%nc(k) = newlabel
      end do
   end subroutine density_profile

   pure subroutine core_labels_from_profile(profile, density, cores, noc)
      type(density_profile_result), intent(in) :: profile !! Connected-component profile over increasing density mass.
      real(dp), intent(in) :: density(:) !! Density values used to order the resulting cluster-core labels.
      integer, allocatable, intent(out) :: cores(:) !! Core membership labels; zero denotes not yet allocated.
      integer, intent(out) :: noc !! Number of distinct cluster cores.
      integer, allocatable :: jump_index(:)
      integer, allocatable :: labels(:)
      integer, allocatable :: remap(:)
      real(dp), allocatable :: peak(:)
      integer :: c
      integer :: i
      integer :: j
      integer :: k
      integer :: lab
      integer :: n
      integer :: njump
      integer :: nseen
      integer :: tmpi
      real(dp) :: tmpr
      logical :: has_new

      n = size(density)
      allocate(cores(n), labels(max(1, n)), remap(max(1, n)), peak(max(1, n)))
      cores = 0
      noc = 0
      if (n == 0 .or. size(profile%nc) < 2) return

      njump = count(profile%nc(2:) /= profile%nc(:size(profile%nc) - 1))
      if (njump <= 1) then
         noc = 1
         cores = 1
         return
      end if
      allocate(jump_index(njump))
      j = 0
      do i = 1, size(profile%nc) - 1
         if (profile%nc(i + 1) == profile%nc(i)) cycle
         j = j + 1
         jump_index(j) = i
      end do

      ! con2tree() allocates cores only at density-profile jump levels, and skips
      ! the first jump.  New observations are attached to an existing core only
      ! while their connected component contains at most one previously created core.
      do j = 2, njump
         k = jump_index(j)
         do c = 1, profile%nc(k)
            labels = 0
            nseen = 0
            has_new = .false.
            do i = 1, n
               if (profile%ids(i, k) /= c) cycle
               if (cores(i) == 0) has_new = .true.
               if (cores(i) <= 0) cycle
               if (nseen == 0 .or. .not. any(labels(1:nseen) == cores(i))) then
                  nseen = nseen + 1
                  labels(nseen) = cores(i)
               end if
            end do
            if (.not. has_new) cycle
            if (nseen == 0) then
               noc = noc + 1
               do i = 1, n
                  if (profile%ids(i, k) == c) cores(i) = noc
               end do
            else if (nseen == 1) then
               lab = labels(1)
               do i = 1, n
                  if (profile%ids(i, k) == c) cores(i) = lab
               end do
            end if
         end do
      end do

      if (noc == 0) then
         noc = 1
         cores = 1
         return
      end if

      ! Match con2tree()'s final relabeling: cores are ordered by the highest
      ! observed density attained by each core, from highest to lowest.
      peak(1:noc) = -huge(1.0_dp)
      remap(1:noc) = [(i, i=1, noc)]
      do lab = 1, noc
         do i = 1, n
            if (cores(i) == lab) peak(lab) = max(peak(lab), density(i))
         end do
      end do
      do i = 2, noc
         tmpi = remap(i)
         tmpr = peak(tmpi)
         k = i - 1
         do while (k >= 1)
            if (peak(remap(k)) >= tmpr) exit
            remap(k + 1) = remap(k)
            k = k - 1
         end do
         remap(k + 1) = tmpi
      end do
      labels = 0
      do i = 1, noc
         labels(remap(i)) = i
      end do
      do i = 1, n
         if (cores(i) > 0) cores(i) = labels(cores(i))
      end do
   end subroutine core_labels_from_profile

end module pdfcluster_graph
