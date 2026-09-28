module pdfcluster_dbs
   use pdfcluster_kinds, only : dp
   use pdfcluster_types, only : pdfcluster_result, dbs_result
   use pdfcluster_kde, only : h_norm, kepdf_estimate
   use pdfcluster_cluster, only : groups_at_stage
   use pdfcluster_stats, only : top_two_indices
   implicit none
   private

   public :: dbs_from_labels, dbs_from_pdfcluster

contains

   pure subroutine dbs_from_labels(x, clusters, result, success, prior, hmult, kernel)
      real(dp), intent(in) :: x(:, :) !! Sample used to estimate cluster-specific densities at every observation.
      integer, intent(in) :: clusters(:) !! Positive cluster label for every row of x.
      type(dbs_result), intent(out) :: result !! Normalized density-based silhouette values and associated metadata.
      logical, intent(out) :: success !! True when at least two valid groups are fitted and all densities are finite.
      real(dp), intent(in), optional :: prior(:) !! Positive prior probability for each group; equal priors are the default.
      real(dp), intent(in), optional :: hmult !! Multiplier applied to each cluster's normal-reference bandwidth.
      character(len=*), intent(in), optional :: kernel !! Kernel name, gaussian or t7; Gaussian is the default.
      real(dp), allocatable :: f(:)
      real(dp), allocatable :: h(:)
      real(dp), allocatable :: joint(:, :)
      real(dp), allocatable :: tau(:, :)
      real(dp) :: hm
      real(dp) :: maxfinite
      real(dp) :: denom
      character(len=12) :: kern
      integer :: first
      integer :: i
      integer :: m
      integer :: noc
      integer :: second
      logical :: ok

      success = .false.
      if (size(clusters) /= size(x, 1) .or. any(clusters <= 0)) return
      noc = maxval(clusters)
      if (noc < 2) return
      result%noc = noc
      result%clusters = clusters
      allocate(result%prior(noc), result%values(size(clusters)), joint(size(x, 1), noc), tau(size(x, 1), noc))
      if (present(prior)) then
         if (size(prior) /= noc .or. any(prior < 0.0_dp) .or. sum(prior) <= 0.0_dp) return
         result%prior = prior/sum(prior)
      else
         result%prior = 1.0_dp/real(noc, dp)
      end if
      hm = 1.0_dp
      if (present(hmult)) hm = hmult
      kern = "gaussian"
      if (present(kernel)) kern = adjustl(kernel)
      do m = 1, noc
         call h_norm(pack_rows(x, clusters == m), h)
         h = max(h*hm, 1.0e-10_dp)
         call kepdf_estimate(pack_rows(x, clusters == m), x, h=h, kernel=kern, estimate=f, success=ok)
         if (.not. ok) return
         joint(:, m) = f*result%prior(m)
         deallocate(h, f)
      end do
      do i = 1, size(x, 1)
         denom = sum(joint(i, :))
         if (denom <= 0.0_dp) return
         tau(i, :) = joint(i, :)/denom
      end do
      do i = 1, size(x, 1)
         call top_two_indices(tau(i, :), first, second)
         if (clusters(i) == first) then
            result%values(i) = log(max(tau(i, first), tiny(1.0_dp))/max(joint(i, second), tiny(1.0_dp)))
         else
            result%values(i) = log(max(tau(i, clusters(i)), tiny(1.0_dp))/max(tau(i, first), tiny(1.0_dp)))
         end if
      end do
      maxfinite = maxval(result%values)
      if (maxfinite > 0.0_dp) result%values = result%values/maxfinite
      success = .true.
   end subroutine dbs_from_labels

   pure subroutine dbs_from_pdfcluster(obj, result, success, stage, prior, hmult, kernel)
      type(pdfcluster_result), intent(in) :: obj !! Existing pdfCluster result whose allocated observations are scored.
      type(dbs_result), intent(out) :: result !! Density-based silhouette values aligned to all original observations.
      logical, intent(out) :: success !! True when the requested stage has at least two groups and scoring succeeds.
      integer, intent(in), optional :: stage !! Stage to score; zero means core labels and omission means the final stage.
      real(dp), intent(in), optional :: prior(:) !! Optional group prior probabilities for the selected observations.
      real(dp), intent(in), optional :: hmult !! Multiplier applied to cluster-specific normal-reference bandwidths.
      character(len=*), intent(in), optional :: kernel !! Kernel name, Gaussian by default.
      type(dbs_result) :: tmp
      integer, allocatable :: g(:)
      logical, allocatable :: keep(:)
      logical :: ok

      success = .false.
      call groups_at_stage(obj, stage, g, ok)
      if (.not. ok) return
      allocate(keep(size(g)))
      keep = g > 0
      if (count(keep) == 0) return
      call dbs_from_labels(pack_rows(obj%x, keep), pack_int(g, keep), tmp, ok, prior=prior, hmult=hmult, kernel=kernel)
      if (.not. ok) return
      result%noc = tmp%noc
      result%prior = tmp%prior
      allocate(result%values(size(g)), result%clusters(size(g)))
      result%values = 0.0_dp
      result%clusters = 0
      result%values = unpack(tmp%values, keep, 0.0_dp)
      result%clusters = unpack(tmp%clusters, keep, 0)
      success = .true.
   end subroutine dbs_from_pdfcluster

   pure function pack_rows(x, mask) result(y)
      real(dp), intent(in) :: x(:, :) !! Matrix from which selected rows are copied.
      logical, intent(in) :: mask(:) !! Row-selection mask with one entry for each row of x.
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

   pure function pack_int(x, mask) result(y)
      integer, intent(in) :: x(:) !! Integer vector from which selected entries are copied.
      logical, intent(in) :: mask(:) !! Selection mask with the same length as x.
      integer, allocatable :: y(:)
      integer :: i
      integer :: k

      allocate(y(count(mask)))
      k = 0
      do i = 1, size(x)
         if (.not. mask(i)) cycle
         k = k + 1
         y(k) = x(i)
      end do
   end function pack_int

end module pdfcluster_dbs
