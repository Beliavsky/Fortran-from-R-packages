! SPDX-License-Identifier: GPL-2.0-or-later
! Computational translation derived from R package spdep 1.4-2.
! Translated/modified 2026-08-30 for modern free-form Fortran.
module spdep_additional
   use spdep_kinds, only : dp
   use spdep_types, only : int_vector, real_vector, neighbor_list, spatial_weights
   use spdep_math, only : variance_dp, safe_nan
   use spdep_graph, only : dnearneigh, nbdists
   use spdep_weights, only : nb2listw, lag_listw
   implicit none
   private

   public :: rotation
   public :: complement_nb
   public :: nblag_cumul
   public :: nb2blocknb
   public :: nb2listwdist
   public :: autocov_dist
   public :: local_geary

contains

   pure function rotation(xy, angle) result(rotated)
      real(dp), intent(in) :: xy(:, :) !! Planar coordinate matrix; the first two columns are rotated about the origin.
      real(dp), intent(in) :: angle !! Counter-clockwise rotation angle in radians.
      real(dp), allocatable :: rotated(:, :)
      real(dp) :: cos_angle
      real(dp) :: sin_angle
      integer :: n

      n = size(xy, 1)
      allocate(rotated(n, 2))
      if (size(xy, 2) < 2) then
         rotated = safe_nan()
         return
      end if
      cos_angle = cos(angle)
      sin_angle = sin(angle)
      rotated(:, 1) = xy(:, 1) * cos_angle - xy(:, 2) * sin_angle
      rotated(:, 2) = xy(:, 1) * sin_angle + xy(:, 2) * cos_angle
   end function rotation

   pure function complement_nb(nb) result(out)
      type(neighbor_list), intent(in) :: nb !! Graph whose per-region set complement is returned; matches upstream complement.nb.
      type(neighbor_list) :: out
      logical, allocatable :: keep(:)
      integer :: n
      integer :: i
      integer :: j
      integer :: k
      integer :: m

      n = nb%size()
      allocate(out%neighbors(n), keep(n))
      do i = 1, n
         keep = .true.
         do k = 1, size(nb%neighbors(i)%values)
            j = nb%neighbors(i)%values(k)
            if (j >= 1 .and. j <= n) keep(j) = .false.
         end do
         m = count(keep)
         allocate(out%neighbors(i)%values(m))
         m = 0
         do j = 1, n
            if (.not. keep(j)) cycle
            m = m + 1
            out%neighbors(i)%values(m) = j
         end do
      end do
      out%self_included = .not. nb%self_included
   end function complement_nb

   pure function nblag_cumul(lags) result(out)
      type(neighbor_list), intent(in) :: lags(:) !! Exact graph-lag lists to union; all must describe the same regions.
      type(neighbor_list) :: out
      logical, allocatable :: adj(:, :)
      integer :: n
      integer :: ell
      integer :: i
      integer :: j
      integer :: k
      integer :: m

      if (size(lags) == 0) then
         allocate(out%neighbors(0))
         return
      end if
      n = lags(1)%size()
      if (any([(lags(ell)%size() /= n, ell = 1, size(lags))])) then
         allocate(out%neighbors(0))
         return
      end if
      allocate(adj(n, n))
      adj = .false.
      do ell = 1, size(lags)
         do i = 1, n
            do k = 1, size(lags(ell)%neighbors(i)%values)
               j = lags(ell)%neighbors(i)%values(k)
               if (j >= 1 .and. j <= n) adj(i, j) = .true.
            end do
         end do
      end do
      allocate(out%neighbors(n))
      do i = 1, n
         m = count(adj(i, :))
         allocate(out%neighbors(i)%values(m))
         m = 0
         do j = 1, n
            if (.not. adj(i, j)) cycle
            m = m + 1
            out%neighbors(i)%values(m) = j
         end do
      end do
      out%self_included = any([(adj(i, i), i = 1, n)])
   end function nblag_cumul

   pure function nb2blocknb(nb, id) result(out)
      type(neighbor_list), intent(in) :: nb !! Neighbor graph among blocks numbered from 1 through nb%size().
      integer, intent(in) :: id(:) !! Block number for each expanded entity; every value must be a valid block index.
      type(neighbor_list) :: out
      logical, allocatable :: allowed(:)
      logical, allocatable :: adj(:, :)
      integer :: n
      integer :: nblock
      integer :: i
      integer :: j
      integer :: k
      integer :: block
      integer :: other_block
      integer :: m

      n = size(id)
      nblock = nb%size()
      if (n == 0 .or. nblock == 0 .or. any(id < 1) .or. any(id > nblock)) then
         allocate(out%neighbors(0))
         return
      end if
      allocate(adj(n, n), allowed(nblock))
      adj = .false.
      do i = 1, n
         block = id(i)
         allowed = .false.
         allowed(block) = .true.
         do k = 1, size(nb%neighbors(block)%values)
            other_block = nb%neighbors(block)%values(k)
            if (other_block >= 1 .and. other_block <= nblock) allowed(other_block) = .true.
         end do
         do j = 1, n
            if (j /= i .and. allowed(id(j))) adj(i, j) = .true.
         end do
      end do
      allocate(out%neighbors(n))
      do i = 1, n
         m = count(adj(i, :))
         allocate(out%neighbors(i)%values(m))
         m = 0
         do j = 1, n
            if (.not. adj(i, j)) cycle
            m = m + 1
            out%neighbors(i)%values(m) = j
         end do
      end do
      out%self_included = .false.
   end function nb2blocknb

   pure function nb2listwdist(nb, coords, weight_type, style, alpha, dmax, longlat) result(listw)
      type(neighbor_list), intent(in) :: nb !! Neighbor topology whose links receive distance-decay weights.
      real(dp), intent(in) :: coords(:, :) !! Coordinate matrix with one row per region and at least two columns.
      character(len=*), intent(in), optional :: weight_type !! Distance weighting: idw, exp, or dpd; default is idw.
      character(len=*), intent(in), optional :: style !! Weight coding style raw, B, W, C, U, S, or minmax; default is raw.
      real(dp), intent(in), optional :: alpha !! Positive distance-decay exponent or rate; default is 1.
      real(dp), intent(in), optional :: dmax !! Positive truncation distance; required by dpd and optional for idw/exp.
      logical, intent(in), optional :: longlat !! If true, coordinates are longitude/latitude and dmax uses metres.
      type(spatial_weights) :: listw
      type(real_vector), allocatable :: distances(:)
      type(real_vector), allocatable :: glist(:)
      character(len=16) :: use_type
      character(len=8) :: use_style
      logical :: use_longlat
      real(dp) :: use_alpha
      real(dp) :: use_dmax
      real(dp) :: distance_value
      real(dp) :: weight_value
      real(dp) :: max_finite
      integer :: n
      integer :: i
      integer :: k

      n = nb%size()
      use_type = "idw"
      if (present(weight_type)) use_type = lower_ascii(trim(adjustl(weight_type)))
      use_style = "raw"
      if (present(style)) use_style = lower_ascii(trim(adjustl(style)))
      use_alpha = 1.0_dp
      if (present(alpha)) use_alpha = alpha
      use_longlat = .false.
      if (present(longlat)) use_longlat = longlat
      use_dmax = -1.0_dp
      if (present(dmax)) use_dmax = dmax

      distances = nbdists(nb, coords, use_longlat)
      allocate(glist(n))
      max_finite = 0.0_dp
      do i = 1, n
         allocate(glist(i)%values(size(distances(i)%values)))
         do k = 1, size(distances(i)%values)
            distance_value = distances(i)%values(k)
            if (use_longlat) distance_value = 1000.0_dp * distance_value
            select case (trim(use_type))
            case ("idw")
               if (distance_value > 0.0_dp) then
                  weight_value = distance_value ** (-use_alpha)
                  if (use_dmax > 0.0_dp .and. distance_value > use_dmax) weight_value = 0.0_dp
                  max_finite = max(max_finite, weight_value)
               else
                  weight_value = -1.0_dp
               end if
            case ("exp")
               weight_value = exp(-use_alpha * distance_value)
               if (use_dmax > 0.0_dp .and. distance_value > use_dmax) weight_value = 0.0_dp
            case ("dpd")
               if (use_dmax <= 0.0_dp .or. distance_value >= use_dmax) then
                  weight_value = 0.0_dp
               else
                  weight_value = (1.0_dp - (distance_value / use_dmax) ** use_alpha) ** use_alpha
               end if
            case default
               weight_value = safe_nan()
            end select
            glist(i)%values(k) = weight_value
         end do
      end do

      if (trim(use_type) == "idw") then
         if (max_finite <= 0.0_dp) max_finite = 1.0_dp
         do i = 1, n
            where (glist(i)%values < 0.0_dp)
               glist(i)%values = max_finite
            end where
         end do
      end if

      if (trim(use_style) == "raw") then
         listw%nb = nb
         allocate(listw%weights(n))
         do i = 1, n
            listw%weights(i)%values = glist(i)%values
         end do
         listw%style = "raw"
         listw%zero_policy = .true.
      else
         listw = nb2listw(nb, trim(use_style), glist)
      end if
   end function nb2listwdist

   pure function autocov_dist(z, coords, nbs, weight_type, style, longlat) result(value)
      real(dp), intent(in) :: z(:) !! Response values whose distance-weighted spatial autocovariate is required.
      real(dp), intent(in) :: coords(:, :) !! Coordinate matrix with one row per response and at least two columns.
      real(dp), intent(in), optional :: nbs !! Neighborhood radius in coordinate units or km for longlat; default is 1.
      character(len=*), intent(in), optional :: weight_type !! Weighting one, inverse, or inverse.squared; default is inverse.
      character(len=*), intent(in), optional :: style !! Weight coding style passed to nb2listw; default is B.
      logical, intent(in), optional :: longlat !! If true, use upstream-compatible WGS84 distances on degree coordinates.
      real(dp), allocatable :: value(:)
      type(neighbor_list) :: nb
      type(real_vector), allocatable :: distances(:)
      type(real_vector), allocatable :: glist(:)
      type(spatial_weights) :: listw
      character(len=16) :: use_type
      character(len=8) :: use_style
      logical :: use_longlat
      real(dp) :: radius
      integer :: exponent
      integer :: n
      integer :: i
      integer :: k

      n = size(z)
      allocate(value(n))
      if (size(coords, 1) /= n .or. size(coords, 2) < 2) then
         value = safe_nan()
         return
      end if
      radius = 1.0_dp
      if (present(nbs)) radius = nbs
      use_type = "inverse"
      if (present(weight_type)) use_type = lower_ascii(trim(adjustl(weight_type)))
      use_style = "B"
      if (present(style)) use_style = trim(adjustl(style))
      use_longlat = .false.
      if (present(longlat)) use_longlat = longlat
      select case (trim(use_type))
      case ("one")
         exponent = 0
      case ("inverse")
         exponent = 1
      case ("inverse.squared")
         exponent = 2
      case default
         value = safe_nan()
         return
      end select
      nb = dnearneigh(coords, 0.0_dp, radius, use_longlat)
      if (exponent == 0) then
         listw = nb2listw(nb, trim(use_style))
      else
         distances = nbdists(nb, coords, use_longlat)
         allocate(glist(n))
         do i = 1, n
            allocate(glist(i)%values(size(distances(i)%values)))
            do k = 1, size(distances(i)%values)
               if (distances(i)%values(k) > 0.0_dp) then
                  glist(i)%values(k) = 1.0_dp / distances(i)%values(k) ** exponent
               else
                  glist(i)%values(k) = 0.0_dp
               end if
            end do
         end do
         listw = nb2listw(nb, trim(use_style), glist)
      end if
      value = lag_listw(listw, z)
   end function autocov_dist

   pure function local_geary(x, listw, standardize) result(value)
      real(dp), intent(in) :: x(:) !! Numeric observations with one value per region.
      type(spatial_weights), intent(in) :: listw !! Spatial weights defining each local squared-difference contribution.
      logical, intent(in), optional :: standardize !! If true, sample-standardize x as upstream localC; default is true.
      real(dp), allocatable :: value(:)
      real(dp), allocatable :: z(:)
      logical :: use_standardize
      real(dp) :: mu
      real(dp) :: sd
      integer :: n
      integer :: i
      integer :: j
      integer :: k

      n = listw%size()
      allocate(value(n), z(n))
      if (size(x) /= n) then
         value = safe_nan()
         return
      end if
      use_standardize = .true.
      if (present(standardize)) use_standardize = standardize
      if (use_standardize) then
         if (n < 2) then
            value = safe_nan()
            return
         end if
         mu = sum(x) / real(n, dp)
         sd = sqrt(variance_dp(x, .true.))
         if (sd <= 0.0_dp) then
            value = safe_nan()
            return
         end if
         z = (x - mu) / sd
      else
         z = x
      end if
      value = 0.0_dp
      do i = 1, n
         do k = 1, size(listw%nb%neighbors(i)%values)
            j = listw%nb%neighbors(i)%values(k)
            if (k <= size(listw%weights(i)%values) .and. j >= 1 .and. j <= n) then
               value(i) = value(i) + listw%weights(i)%values(k) * (z(j) - z(i)) ** 2
            end if
         end do
      end do
   end function local_geary

   pure function lower_ascii(text) result(lower)
      character(len=*), intent(in) :: text !! ASCII text converted to lowercase for case-insensitive option matching.
      character(len=len(text)) :: lower
      integer :: i
      integer :: code

      lower = text
      do i = 1, len(text)
         code = iachar(lower(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code + 32)
      end do
   end function lower_ascii

end module spdep_additional
