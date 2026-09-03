! Translation of R/genZcor.R and R/fixed2Zcor.R from geepack 1.3-13.
! Upstream license: GPL (>= 3). See LICENSE, NOTICE.md, and PROVENANCE.md.
! Exact upstream copyright-holder attribution is preserved in NOTICE.md and upstream/DESCRIPTION.
module geepack_design
   use r_kinds, only : dp
   use geepack_correlations, only : COR_INDEPENDENCE, COR_EXCHANGEABLE, COR_AR1, &
      COR_UNSTRUCTURED, COR_USERDEFINED
   use geepack_status, only : GEE_OK, GEE_ERR_ARGUMENT
   implicit none
   private

   public :: gen_zcor, gen_zodds, fixed_to_zcor

contains

   subroutine gen_zcor(cluster_sizes, waves, corstr, zcor, status)
      integer, intent(in) :: cluster_sizes(:) !! Number of observations in each contiguous cluster.
      integer, intent(in) :: waves(:) !! One-based wave indices for all observations in cluster order.
      integer, intent(in) :: corstr !! Correlation structure identifier COR_*.
      real(dp), allocatable, intent(out) :: zcor(:, :) !! Correlation design matrix.
      integer, intent(out) :: status !! GEE_OK on success or an error code.
      integer :: maxwave
      integer :: npair_total
      integer :: ncol
      integer :: row
      integer :: pos
      integer :: g
      integer :: i
      integer :: j
      integer :: idx

      status = GEE_OK
      if (sum(cluster_sizes) /= size(waves) .or. any(cluster_sizes < 1)) then
         allocate(zcor(0, 0))
         status = GEE_ERR_ARGUMENT
         return
      end if
      select case (corstr)
      case (COR_INDEPENDENCE)
         allocate(zcor(0, 0))
      case (COR_EXCHANGEABLE, COR_AR1)
         allocate(zcor(size(cluster_sizes), 1))
         zcor = 1.0_dp
      case (COR_UNSTRUCTURED, COR_USERDEFINED)
         maxwave = maxval(cluster_sizes)
         if (any(waves > maxwave)) then
            allocate(zcor(0, 0))
            status = GEE_ERR_ARGUMENT
            return
         end if
         ncol = maxwave * (maxwave - 1) / 2
         npair_total = sum(cluster_sizes * (cluster_sizes - 1) / 2)
         allocate(zcor(npair_total, ncol))
         zcor = 0.0_dp
         row = 0
         pos = 0
         do g = 1, size(cluster_sizes)
            do i = 1, cluster_sizes(g) - 1
               do j = i + 1, cluster_sizes(g)
                  row = row + 1
                  idx = pair_index(min(waves(pos + i), waves(pos + j)), &
                     max(waves(pos + i), waves(pos + j)), maxwave)
                  if (idx >= 1 .and. idx <= ncol) zcor(row, idx) = 1.0_dp
               end do
            end do
            pos = pos + cluster_sizes(g)
         end do
      case default
         allocate(zcor(0, 0))
         status = GEE_ERR_ARGUMENT
      end select
   end subroutine gen_zcor

   subroutine gen_zodds(cluster_sizes, waves, corstr, ncat, z, status)
      integer, intent(in) :: cluster_sizes(:) !! Number of ordinal observations in each cluster.
      integer, intent(in) :: waves(:) !! One-based wave indices for ordinal observations.
      integer, intent(in) :: corstr !! Ordinal association structure identifier.
      integer, intent(in) :: ncat !! Number of cumulative indicators, response levels minus one.
      real(dp), allocatable, intent(out) :: z(:, :) !! Odds-ratio design matrix.
      integer, intent(out) :: status !! GEE_OK on success or an error code.
      real(dp), allocatable :: base(:, :)
      integer :: c2
      integer :: row
      integer :: pair
      integer :: k

      c2 = ncat * ncat
      if (ncat < 1) then
         allocate(z(0, 0))
         status = GEE_ERR_ARGUMENT
         return
      end if
      if (corstr == COR_INDEPENDENCE) then
         allocate(z(0, 0))
         status = GEE_OK
         return
      end if
      if (corstr == COR_EXCHANGEABLE .or. corstr == COR_AR1) then
         allocate(z(sum(cluster_sizes * (cluster_sizes - 1) / 2) * c2, 1))
         z = 1.0_dp
         status = GEE_OK
         return
      end if
      call gen_zcor(cluster_sizes, waves, COR_UNSTRUCTURED, base, status)
      if (status /= GEE_OK) then
         allocate(z(0, 0))
         return
      end if
      allocate(z(size(base, 1) * c2, size(base, 2)))
      row = 0
      do pair = 1, size(base, 1)
         do k = 1, c2
            row = row + 1
            z(row, :) = base(pair, :)
         end do
      end do
   end subroutine gen_zodds

   subroutine fixed_to_zcor(correlation, cluster_sizes, waves, zcor, status)
      real(dp), intent(in) :: correlation(:, :) !! Fixed square correlation matrix indexed by wave.
      integer, intent(in) :: cluster_sizes(:) !! Number of observations in each cluster.
      integer, intent(in) :: waves(:) !! One-based wave indices for all observations.
      real(dp), allocatable, intent(out) :: zcor(:) !! Pair correlations in cluster upper-triangle order.
      integer, intent(out) :: status !! GEE_OK on success or an error code.
      integer :: g
      integer :: i
      integer :: j
      integer :: pos
      integer :: row

      status = GEE_OK
      if (size(correlation, 1) /= size(correlation, 2) .or. sum(cluster_sizes) /= size(waves)) then
         allocate(zcor(0))
         status = GEE_ERR_ARGUMENT
         return
      end if
      allocate(zcor(sum(cluster_sizes * (cluster_sizes - 1) / 2)))
      pos = 0
      row = 0
      do g = 1, size(cluster_sizes)
         do i = 1, cluster_sizes(g) - 1
            do j = i + 1, cluster_sizes(g)
               row = row + 1
               if (waves(pos + i) < 1 .or. waves(pos + j) < 1 .or. &
                   waves(pos + i) > size(correlation, 1) .or. &
                   waves(pos + j) > size(correlation, 1)) then
                  status = GEE_ERR_ARGUMENT
                  zcor = 0.0_dp
                  return
               end if
               zcor(row) = correlation(waves(pos + j), waves(pos + i))
            end do
         end do
         pos = pos + cluster_sizes(g)
      end do
   end subroutine fixed_to_zcor

   pure integer function pair_index(i, j, n) result(idx)
      integer, intent(in) :: i !! Smaller one-based wave index.
      integer, intent(in) :: j !! Larger one-based wave index.
      integer, intent(in) :: n !! Maximum wave index.
      integer :: a

      idx = 0
      do a = 1, i - 1
         idx = idx + n - a
      end do
      idx = idx + j - i
   end function pair_index

end module geepack_design
