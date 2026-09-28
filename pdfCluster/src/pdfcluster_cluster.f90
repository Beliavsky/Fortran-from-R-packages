module pdfcluster_cluster
   use pdfcluster_kinds, only : dp
   use pdfcluster_types, only : pdfcluster_result, density_profile_result
   use pdfcluster_kde, only : h_norm, hprop2f, kepdf_fit, kepdf_estimate
   use pdfcluster_graph, only : graph_unidimensional, graph_delaunay, graph_pairs
   use pdfcluster_graph, only : density_profile, core_labels_from_profile
   use pdfcluster_stats, only : quantile_type7, top_two_indices
   implicit none
   private

   public :: pdfcluster_fit, pdfcluster_regraph_pairs, pdf_classification, groups_at_stage

contains

   pure subroutine pdfcluster_fit(x, result, success, graphtype, hmult, kernel, bwtype, h, hx, alpha, lambda, &
                                  grid_pairs, n_grid, n_stage, se, hcores)
      real(dp), intent(in) :: x(:, :) !! Continuous numeric sample to cluster, with observations in rows.
      type(pdfcluster_result), intent(out) :: result !! Density-clustering result including graph, cores and stages.
      logical, intent(out) :: success !! True when density estimation, graph construction and core extraction succeed.
      character(len=*), intent(in), optional :: graphtype !! Graph mode; default depends on the data dimension.
      real(dp), intent(in), optional :: hmult !! Multiplicative factor applied to the fixed/pilot normal-reference bandwidth.
      character(len=*), intent(in), optional :: kernel !! Density kernel: gaussian or t7.
      character(len=*), intent(in), optional :: bwtype !! Bandwidth mode: fixed or adaptive.
      real(dp), intent(in), optional :: h(:) !! User-supplied fixed or pilot bandwidth before hmult scaling.
      real(dp), intent(in), optional :: hx(:, :) !! User-supplied adaptive row-specific bandwidth matrix.
      real(dp), intent(in), optional :: alpha !! Adaptive-bandwidth exponent; default is 0.5.
      real(dp), intent(in), optional :: lambda !! Pair-graph valley threshold; default is 0.10.
      integer, intent(in), optional :: grid_pairs !! Segment grid size for pair-graph valley calculations; default is 10.
      integer, intent(in), optional :: n_grid !! Total density-profile grid size including endpoints; default follows upstream rule.
      integer, intent(in), optional :: n_stage !! Number of staged classification passes; default is 5, zero skips classification.
      logical, intent(in), optional :: se !! Whether density log-ratios are standardized by their approximate standard error.
      logical, intent(in), optional :: hcores !! Whether classification bandwidths expand directly from cluster-core bandwidths.
      type(density_profile_result) :: profile
      real(dp), allocatable :: h0(:)
      real(dp) :: hm
      real(dp) :: lam
      character(len=16) :: graph_name
      character(len=12) :: kern
      character(len=12) :: bw
      integer :: d
      integer :: gp
      integer :: ng
      integer :: ns
      logical :: ok
      logical :: use_se
      logical :: use_hcores

      success = .false.
      if (size(x, 1) <= 0 .or. size(x, 2) <= 0) return
      d = size(x, 2)
      hm = merge(1.0_dp, 0.75_dp, d > 6)
      if (present(hmult)) hm = hmult
      kern = "gaussian"
      bw = "fixed"
      if (present(kernel)) kern = adjustl(kernel)
      if (present(bwtype)) bw = adjustl(bwtype)
      lam = 0.10_dp
      if (present(lambda)) lam = lambda
      gp = 10
      if (present(grid_pairs)) gp = grid_pairs
      ng = min(nint((5.0_dp + sqrt(real(size(x, 1), dp)))*4.0_dp), size(x, 1))
      if (present(n_grid)) ng = min(max(2, n_grid), size(x, 1))
      ns = 5
      if (present(n_stage)) ns = max(0, n_stage)
      use_se = .true.
      use_hcores = .false.
      if (present(se)) use_se = se
      if (present(hcores)) use_hcores = hcores
      if (present(graphtype)) then
         graph_name = adjustl(graphtype)
      else if (d == 1) then
         graph_name = "unidimensional"
      else if (d <= 6) then
         graph_name = "delaunay"
      else
         graph_name = "pairs"
      end if

      if (present(h)) then
         allocate(h0(size(h)))
         h0 = hm*h
      else
         call h_norm(x, h0)
         h0 = hm*h0
      end if
      if (any(h0 <= 0.0_dp)) return
      if (present(hx)) then
         call kepdf_fit(x, x, result%pdf, h=h0, hx=hx, kernel=kern, bwtype=bw, alpha=alpha, success=ok)
      else
         call kepdf_fit(x, x, result%pdf, h=h0, kernel=kern, bwtype=bw, alpha=alpha, success=ok)
      end if
      if (.not. ok) return
      result%x = x
      result%graphtype = graph_name
      result%lambda = lam

      select case (trim(graph_name))
      case ("unidimensional")
         if (d /= 1) return
         call graph_unidimensional(x(:, 1), result%adjacency)
         allocate(result%pair_valley_areas(0), result%pair_indices(2, 0))
      case ("delaunay")
         call graph_delaunay(x, result%adjacency, ok)
         if (.not. ok) return
         allocate(result%pair_valley_areas(0), result%pair_indices(2, 0))
      case ("pairs")
         call graph_pairs(x, result%pdf, lam, gp, result%adjacency, result%pair_valley_areas, result%pair_indices, ok)
         if (.not. ok) return
      case default
         return
      end select

      call density_profile(result%adjacency, result%pdf%estimate, max(0, ng - 2), .true., profile)
      call core_labels_from_profile(profile, result%pdf%estimate, result%cluster_cores, result%noc)
      allocate(result%stages(size(x, 1), 0), result%clusters(size(x, 1)))
      result%clusters = result%cluster_cores
      if (ns > 0) call pdf_classification(result, ns, use_se, use_hcores, ok)
      if (ns > 0 .and. .not. ok) return
      success = .true.
   end subroutine pdfcluster_fit

   pure subroutine pdfcluster_regraph_pairs(obj, result, success, lambda, n_grid, n_stage, se, hcores)
      type(pdfcluster_result), intent(in) :: obj !! Existing pair-graph clustering result whose valley threshold is changed.
      type(pdfcluster_result), intent(out) :: result !! Recomputed pair-graph clustering result using the stored valley measures.
      logical, intent(out) :: success !! True when stored pair information is valid and core/classification recomputation succeeds.
      real(dp), intent(in), optional :: lambda !! Replacement pair-graph valley threshold; default is 0.10.
      integer, intent(in), optional :: n_grid !! Total density-profile grid size including endpoints; default follows upstream rule.
      integer, intent(in), optional :: n_stage !! Number of classification stages; default is 5 and zero skips classification.
      logical, intent(in), optional :: se !! Standardize density log-ratios by their approximate standard error when true.
      logical, intent(in), optional :: hcores !! Expand bandwidths from core bandwidths directly when true.
      type(density_profile_result) :: profile
      real(dp) :: lam
      integer :: i
      integer :: k
      integer :: n
      integer :: ng
      integer :: ns
      logical :: ok
      logical :: use_hcores
      logical :: use_se

      success = .false.
      if (trim(obj%graphtype) /= "pairs") return
      if (.not. allocated(obj%x) .or. .not. allocated(obj%pdf%estimate)) return
      if (.not. allocated(obj%pair_valley_areas) .or. .not. allocated(obj%pair_indices)) return
      if (size(obj%pair_indices, 1) /= 2) return
      if (size(obj%pair_indices, 2) /= size(obj%pair_valley_areas)) return
      n = size(obj%x, 1)
      if (size(obj%pdf%estimate) /= n) return

      result = obj
      lam = 0.10_dp
      if (present(lambda)) lam = lambda
      result%lambda = lam
      if (allocated(result%adjacency)) deallocate(result%adjacency)
      allocate(result%adjacency(n, n))
      result%adjacency = .false.
      do i = 1, n
         result%adjacency(i, i) = .true.
      end do
      do k = 1, size(result%pair_valley_areas)
         if (result%pair_valley_areas(k) > lam) cycle
         i = result%pair_indices(1, k)
         if (i < 1 .or. i > n) return
         if (result%pair_indices(2, k) < 1 .or. result%pair_indices(2, k) > n) return
         result%adjacency(i, result%pair_indices(2, k)) = .true.
         result%adjacency(result%pair_indices(2, k), i) = .true.
      end do

      ng = min(nint((5.0_dp + sqrt(real(n, dp)))*4.0_dp), n)
      if (present(n_grid)) ng = min(max(2, n_grid), n)
      call density_profile(result%adjacency, result%pdf%estimate, max(0, ng - 2), .true., profile)
      call core_labels_from_profile(profile, result%pdf%estimate, result%cluster_cores, result%noc)

      ns = 5
      if (present(n_stage)) ns = max(0, n_stage)
      use_se = .true.
      use_hcores = .false.
      if (present(se)) use_se = se
      if (present(hcores)) use_hcores = hcores
      if (allocated(result%stages)) deallocate(result%stages)
      if (allocated(result%clusters)) deallocate(result%clusters)
      allocate(result%stages(n, 0), result%clusters(n))
      result%clusters = result%cluster_cores
      if (ns > 0) then
         call pdf_classification(result, ns, use_se, use_hcores, ok)
         if (.not. ok) return
      end if
      success = .true.
   end subroutine pdfcluster_regraph_pairs

   pure subroutine pdf_classification(obj, n_stage, se, hcores, success)
      type(pdfcluster_result), intent(inout) :: obj !! Existing clustering result whose unallocated points are classified.
      integer, intent(in) :: n_stage !! Positive number of staged classification passes.
      logical, intent(in) :: se !! Standardize log-density ratios using the upstream approximate variance formula when true.
      logical, intent(in) :: hcores !! Expand core bandwidths directly with stage number when true.
      logical, intent(out) :: success !! True when all stage calculations complete successfully.
      real(dp), allocatable :: f(:, :)
      real(dp), allocatable :: h(:)
      real(dp), allocatable :: h1(:)
      real(dp), allocatable :: h1x(:, :)
      real(dp), allocatable :: hxmed(:, :)
      real(dp), allocatable :: score(:)
      real(dp), allocatable :: valid_score(:)
      real(dp), allocatable :: se2(:, :)
      real(dp), allocatable :: values(:)
      real(dp) :: alp
      real(dp) :: h1scalar
      real(dp) :: threshold
      integer, allocatable :: g(:)
      integer, allocatable :: unallocated(:)
      logical, allocatable :: valid(:)
      integer :: count_un
      integer :: first
      integer :: i
      integer :: idx
      integer :: m
      integer :: nm
      integer :: nvalid
      integer :: second
      integer :: stage
      logical :: ok

      success = .false.
      if (n_stage <= 0 .or. obj%noc <= 0) return
      allocate(g(size(obj%x, 1)))
      g = obj%cluster_cores
      if (allocated(obj%stages)) deallocate(obj%stages)
      if (obj%noc == 1) then
         g = 1
         allocate(obj%stages(size(g), n_stage))
         do stage = 1, n_stage
            obj%stages(:, stage) = g
         end do
         obj%clusters = g
         obj%cluster_cores = merge(obj%cluster_cores, 1, obj%cluster_cores > 0)
         success = .true.
         return
      end if
      allocate(obj%stages(size(g), n_stage))
      obj%stages = 0
      do stage = 1, n_stage
         count_un = count(g == 0)
         if (count_un == 0) then
            obj%stages(:, stage:n_stage) = spread(g, dim=2, ncopies=n_stage - stage + 1)
            exit
         end if
         allocate(unallocated(count_un), f(count_un, obj%noc), se2(count_un, obj%noc))
         idx = 0
         do i = 1, size(g)
            if (g(i) == 0) then
               idx = idx + 1
               unallocated(idx) = i
            end if
         end do
         do m = 1, obj%noc
            nm = count(g == m)
            if (nm <= 0) return
            alp = real(nm, dp)/real(count(g /= 0), dp)
            if (trim(obj%pdf%bwtype) == "fixed") then
               allocate(h(size(obj%x, 2)))
               if (hcores) then
                  h = obj%pdf%h*(1.0_dp + real(stage, dp)/real(n_stage, dp))
               else
                  call h_norm(pack_rows(obj%x, g == m), h1)
                  h1scalar = max(1.0e-4_dp, maxval(h1))
                  h = exp((1.0_dp - alp)*log(obj%pdf%h) + alp*log(h1scalar))
               end if
               call kepdf_estimate(pack_rows(obj%x, g == m), obj%x(unallocated, :), h=h, kernel=obj%pdf%kernel, &
                                   estimate=values, success=ok)
               if (.not. ok) return
               f(:, m) = values
               se2(:, m) = 1.0_dp/max(values*real(nm, dp)*product(h), tiny(1.0_dp))
               deallocate(h)
               if (allocated(h1)) deallocate(h1)
            else
               if (hcores) then
                  hxmed = pack_rows(obj%pdf%hx, g == m)*(1.0_dp + real(stage, dp)/real(n_stage, dp))
               else
                  call h_norm(obj%x, h)
                  call hprop2f(pack_rows(obj%x, g == m), h, h1x, alpha=obj%pdf%alpha, kernel=obj%pdf%kernel, success=ok)
                  if (.not. ok) return
                  hxmed = exp((1.0_dp - alp)*log(pack_rows(obj%pdf%hx, g == m)) + alp*log(h1x))
                  deallocate(h)
               end if
               call kepdf_estimate(pack_rows(obj%x, g == m), obj%x(unallocated, :), hx=hxmed, kernel=obj%pdf%kernel, &
                                   estimate=values, success=ok)
               if (.not. ok) return
               f(:, m) = values
               se2(:, m) = 1.0_dp/max(values*real(nm, dp)*product(obj%pdf%h), tiny(1.0_dp))
               if (allocated(hxmed)) deallocate(hxmed)
               if (allocated(h1x)) deallocate(h1x)
            end if
         end do
         allocate(score(count_un), valid(count_un))
         valid = sum(f, dim=2) > 0.0_dp
         score = -huge(1.0_dp)
         do i = 1, count_un
            if (.not. valid(i)) cycle
            call top_two_indices(f(i, :), first, second)
            if (f(i, first) <= 0.0_dp) then
               valid(i) = .false.
            else if (f(i, second) <= 0.0_dp) then
               if (se) then
                  valid(i) = .false.
               else
                  score(i) = huge(1.0_dp)
               end if
            else if (se) then
               score(i) = (log(f(i, first)) - log(f(i, second)))/sqrt(se2(i, first) + se2(i, second))
            else
               score(i) = log(f(i, first)) - log(f(i, second))
            end if
         end do
         nvalid = count(valid)
         if (nvalid > 0) then
            allocate(valid_score(nvalid))
            valid_score = pack(score, valid)
            if (stage < n_stage) then
               threshold = quantile_type7(valid_score, real(n_stage - stage, dp)/real(n_stage, dp))
            else
               threshold = -huge(1.0_dp)
            end if
            do i = 1, count_un
               if (.not. valid(i) .or. score(i) < threshold) cycle
               call top_two_indices(f(i, :), first, second)
               g(unallocated(i)) = first
            end do
            deallocate(valid_score)
         end if
         obj%stages(:, stage) = g
         deallocate(unallocated, f, se2, score, valid)
         if (allocated(values)) deallocate(values)
      end do
      do stage = 2, n_stage
         if (all(obj%stages(:, stage) == 0)) obj%stages(:, stage) = obj%stages(:, stage - 1)
      end do
      obj%clusters = g
      success = .true.
   end subroutine pdf_classification

   pure subroutine groups_at_stage(obj, stage, groups, success)
      type(pdfcluster_result), intent(in) :: obj !! Clustering result from which core or staged group labels are requested.
      integer, intent(in), optional :: stage !! Stage number; zero requests cluster cores and omission requests the last stage.
      integer, allocatable, intent(out) :: groups(:) !! Group labels with zero denoting unallocated observations.
      logical, intent(out) :: success !! True when the requested stage exists.
      integer :: s

      success = .false.
      if (present(stage)) then
         s = stage
      else if (allocated(obj%stages)) then
         s = size(obj%stages, 2)
      else
         s = 0
      end if
      if (s == 0) then
         if (.not. allocated(obj%cluster_cores)) return
         groups = obj%cluster_cores
         success = .true.
      else if (allocated(obj%stages)) then
         if (s < 1 .or. s > size(obj%stages, 2)) return
         groups = obj%stages(:, s)
         success = .true.
      else if (.not. present(stage) .and. allocated(obj%clusters)) then
         groups = obj%clusters
         success = .true.
      end if
   end subroutine groups_at_stage

   pure function pack_rows(x, mask) result(y)
      real(dp), intent(in) :: x(:, :) !! Matrix from which selected rows are copied.
      logical, intent(in) :: mask(:) !! Row-selection mask with length equal to size(x,1).
      real(dp), allocatable :: y(:, :)
      integer :: i
      integer :: k

      allocate(y(count(mask), size(x, 2)))
      k = 0
      do i = 1, size(x, 1)
         if (.not. mask(i)) cycle
         k = k + 1
         y(k, :) = x(i, :)
      end do
   end function pack_rows

end module pdfcluster_cluster
