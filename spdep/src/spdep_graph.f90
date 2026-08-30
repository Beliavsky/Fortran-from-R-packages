! SPDX-License-Identifier: GPL-2.0-or-later
! Computational translation derived from R package spdep 1.4-2.
! Translated/modified 2026-08-30 for modern free-form Fortran.
module spdep_graph
   use spdep_kinds, only : dp
   use spdep_types, only : int_vector, real_vector, neighbor_list, knn_result
   use spdep_math, only : distance_metric, safe_nan
   implicit none
   private

   real(dp), parameter :: earth_a_km = 6378.137_dp
   real(dp), parameter :: earth_f = 1.0_dp / 298.257223563_dp
   real(dp), parameter :: deg_to_rad = acos(-1.0_dp) / 180.0_dp

   public :: cell2nb
   public :: dnearneigh
   public :: knearneigh
   public :: knn2nb
   public :: nbdists
   public :: gabrielneigh
   public :: relative_neighborhood
   public :: tri2nb
   public :: connected_components
   public :: is_symmetric_nb
   public :: make_symmetric_nb
   public :: include_self
   public :: remove_self
   public :: nb_union
   public :: nb_intersection
   public :: nb_difference
   public :: droplinks
   public :: addlinks
   public :: nblag
   public :: nb_adjacency_matrix
   public :: graph_distance_matrix
   public :: euclidean_distance
   public :: great_circle_distance

contains

   pure elemental real(dp) function euclidean_distance(x1, y1, x2, y2) result(d)
      real(dp), intent(in) :: x1 !! First point x coordinate in the coordinate system of the input data.
      real(dp), intent(in) :: y1 !! First point y coordinate in the coordinate system of the input data.
      real(dp), intent(in) :: x2 !! Second point x coordinate in the coordinate system of the input data.
      real(dp), intent(in) :: y2 !! Second point y coordinate in the coordinate system of the input data.

      d = hypot(x1 - x2, y1 - y2)
   end function euclidean_distance

   pure elemental real(dp) function great_circle_distance(lon1, lat1, lon2, lat2) result(d)
      real(dp), intent(in) :: lon1 !! Longitude of the first point in degrees east.
      real(dp), intent(in) :: lat1 !! Latitude of the first point in degrees north.
      real(dp), intent(in) :: lon2 !! Longitude of the second point in degrees east.
      real(dp), intent(in) :: lat2 !! Latitude of the second point in degrees north.
      real(dp) :: fmid
      real(dp) :: ghalf
      real(dp) :: lhalf
      real(dp) :: s
      real(dp) :: c
      real(dp) :: w
      real(dp) :: r
      real(dp) :: h1
      real(dp) :: h2
      real(dp) :: base

      if (lon1 == lon2 .and. lat1 == lat2) then
         d = 0.0_dp
         return
      end if
      fmid = 0.5_dp * (lat1 + lat2) * deg_to_rad
      ghalf = 0.5_dp * (lat1 - lat2) * deg_to_rad
      lhalf = 0.5_dp * (lon1 - lon2) * deg_to_rad
      s = sin(ghalf) ** 2 * cos(lhalf) ** 2 + cos(fmid) ** 2 * sin(lhalf) ** 2
      c = cos(ghalf) ** 2 * cos(lhalf) ** 2 + sin(fmid) ** 2 * sin(lhalf) ** 2
      if (s <= tiny(1.0_dp)) then
         d = 0.0_dp
         return
      end if
      if (c <= tiny(1.0_dp)) then
         d = acos(-1.0_dp) * earth_a_km
         return
      end if
      w = atan(sqrt(s / c))
      if (abs(w) <= tiny(1.0_dp)) then
         d = 0.0_dp
         return
      end if
      r = sqrt(s * c) / w
      base = 2.0_dp * w * earth_a_km
      h1 = (3.0_dp * r - 1.0_dp) / (2.0_dp * c)
      h2 = (3.0_dp * r + 1.0_dp) / (2.0_dp * s)
      d = base * (1.0_dp + earth_f * h1 * sin(fmid) ** 2 * cos(ghalf) ** 2 &
         - earth_f * h2 * cos(fmid) ** 2 * sin(ghalf) ** 2)
   end function great_circle_distance

   pure function cell2nb(nrow, ncol, queen, torus) result(nb)
      integer, intent(in) :: nrow !! Number of grid rows; must be positive.
      integer, intent(in) :: ncol !! Number of grid columns; must be positive.
      logical, intent(in), optional :: queen !! If true, include diagonal neighbors; default is rook adjacency only.
      logical, intent(in), optional :: torus !! If true, wrap both grid dimensions periodically; default is false.
      type(neighbor_list) :: nb
      logical :: use_queen
      logical :: use_torus
      logical, allocatable :: adj(:, :)
      integer :: r
      integer :: c
      integer :: dr
      integer :: dc
      integer :: rr
      integer :: cc
      integer :: i
      integer :: j

      use_queen = .false.
      if (present(queen)) use_queen = queen
      use_torus = .false.
      if (present(torus)) use_torus = torus
      if (nrow <= 0 .or. ncol <= 0) then
         allocate(nb%neighbors(0))
         return
      end if
      allocate(adj(nrow * ncol, nrow * ncol))
      adj = .false.
      do c = 1, ncol
         do r = 1, nrow
            i = r + (c - 1) * nrow
            do dc = -1, 1
               do dr = -1, 1
                  if (dr == 0 .and. dc == 0) cycle
                  if (.not. use_queen .and. abs(dr) + abs(dc) /= 1) cycle
                  rr = r + dr
                  cc = c + dc
                  if (use_torus) then
                     rr = 1 + modulo(rr - 1, nrow)
                     cc = 1 + modulo(cc - 1, ncol)
                  else
                     if (rr < 1 .or. rr > nrow .or. cc < 1 .or. cc > ncol) cycle
                  end if
                  j = rr + (cc - 1) * nrow
                  if (j /= i) adj(i, j) = .true.
               end do
            end do
         end do
      end do
      nb = matrix_to_nb(adj)
   end function cell2nb

   pure function dnearneigh(coords, d1, d2, longlat) result(nb)
      real(dp), intent(in) :: coords(:, :) !! Point coordinates with shape (n,2); columns are x/y or longitude/latitude.
      real(dp), intent(in) :: d1 !! Inclusive lower distance bound in coordinate units or kilometres for longlat data.
      real(dp), intent(in) :: d2 !! Inclusive upper distance bound in coordinate units or kilometres for longlat data.
      logical, intent(in), optional :: longlat !! If true, use the upstream WGS84 distance approximation on degree coordinates.
      type(neighbor_list) :: nb
      logical :: use_longlat
      logical, allocatable :: adj(:, :)
      integer :: i
      integer :: j
      integer :: n
      real(dp) :: d

      n = size(coords, 1)
      use_longlat = .false.
      if (present(longlat)) use_longlat = longlat
      if (size(coords, 2) < 2 .or. d1 < 0.0_dp .or. d2 < d1) then
         allocate(nb%neighbors(0))
         return
      end if
      allocate(adj(n, n))
      adj = .false.
      do i = 1, n - 1
         do j = i + 1, n
            if (use_longlat) then
               d = great_circle_distance(coords(i, 1), coords(i, 2), coords(j, 1), coords(j, 2))
            else
               d = euclidean_distance(coords(i, 1), coords(i, 2), coords(j, 1), coords(j, 2))
            end if
            if (d >= d1 .and. d <= d2) then
               adj(i, j) = .true.
               adj(j, i) = .true.
            end if
         end do
      end do
      nb = matrix_to_nb(adj)
   end function dnearneigh

   pure function knearneigh(coords, k, longlat) result(knn)
      real(dp), intent(in) :: coords(:, :) !! Point coordinates with shape (n,2); columns are x/y or longitude/latitude.
      integer, intent(in) :: k !! Number of nearest neighbors requested for each point; must lie between 1 and n-1.
      logical, intent(in), optional :: longlat !! If true, use the upstream WGS84 distance approximation on degree coordinates.
      type(knn_result) :: knn
      logical :: use_longlat
      real(dp), allocatable :: d(:)
      integer, allocatable :: idx(:)
      integer :: i
      integer :: j
      integer :: m
      integer :: pos
      integer :: n
      integer :: tmpi
      real(dp) :: tmpd

      n = size(coords, 1)
      use_longlat = .false.
      if (present(longlat)) use_longlat = longlat
      if (size(coords, 2) < 2 .or. k < 1 .or. k >= n) then
         allocate(knn%index(0, 0), knn%distance(0, 0))
         return
      end if
      allocate(knn%index(n, k), knn%distance(n, k), d(n - 1), idx(n - 1))
      do i = 1, n
         m = 0
         do j = 1, n
            if (j == i) cycle
            m = m + 1
            idx(m) = j
            if (use_longlat) then
               d(m) = great_circle_distance(coords(i, 1), coords(i, 2), coords(j, 1), coords(j, 2))
            else
               d(m) = euclidean_distance(coords(i, 1), coords(i, 2), coords(j, 1), coords(j, 2))
            end if
         end do
         do j = 2, m
            tmpd = d(j)
            tmpi = idx(j)
            pos = j - 1
            do while (pos >= 1)
               if (d(pos) < tmpd) exit
               if (d(pos) == tmpd .and. idx(pos) < tmpi) exit
               d(pos + 1) = d(pos)
               idx(pos + 1) = idx(pos)
               pos = pos - 1
            end do
            d(pos + 1) = tmpd
            idx(pos + 1) = tmpi
         end do
         knn%index(i, :) = idx(1:k)
         knn%distance(i, :) = d(1:k)
      end do
   end function knearneigh

   pure function knn2nb(knn, sym) result(nb)
      type(knn_result), intent(in) :: knn !! K-nearest-neighbor indices with one row per observation.
      logical, intent(in), optional :: sym !! If true, return the union of directed KNN links and their reverses.
      type(neighbor_list) :: nb
      logical :: use_sym
      logical, allocatable :: adj(:, :)
      integer :: i
      integer :: j
      integer :: k
      integer :: n

      n = size(knn%index, 1)
      use_sym = .false.
      if (present(sym)) use_sym = sym
      allocate(adj(n, n))
      adj = .false.
      do i = 1, n
         do k = 1, size(knn%index, 2)
            j = knn%index(i, k)
            if (j < 1 .or. j > n .or. j == i) cycle
            adj(i, j) = .true.
            if (use_sym) adj(j, i) = .true.
         end do
      end do
      nb = matrix_to_nb(adj)
   end function knn2nb

   pure function nbdists(nb, coords, longlat) result(distances)
      type(neighbor_list), intent(in) :: nb !! Neighbor list whose link distances are required.
      real(dp), intent(in) :: coords(:, :) !! Point coordinates with one row per region and at least two columns.
      logical, intent(in), optional :: longlat !! If true, return WGS84 distances in kilometres rather than planar distances.
      type(real_vector), allocatable :: distances(:)
      logical :: use_longlat
      integer :: i
      integer :: j
      integer :: k
      integer :: n

      n = nb%size()
      use_longlat = .false.
      if (present(longlat)) use_longlat = longlat
      allocate(distances(n))
      if (size(coords, 1) /= n .or. size(coords, 2) < 2) then
         do i = 1, n
            allocate(distances(i)%values(0))
         end do
         return
      end if
      do i = 1, n
         allocate(distances(i)%values(size(nb%neighbors(i)%values)))
         do k = 1, size(nb%neighbors(i)%values)
            j = nb%neighbors(i)%values(k)
            if (use_longlat) then
               distances(i)%values(k) = great_circle_distance(&
                  coords(i, 1), coords(i, 2), coords(j, 1), coords(j, 2))
            else
               distances(i)%values(k) = euclidean_distance(&
                  coords(i, 1), coords(i, 2), coords(j, 1), coords(j, 2))
            end if
         end do
      end do
   end function nbdists

   pure function gabrielneigh(coords) result(nb)
      real(dp), intent(in) :: coords(:, :) !! Planar point coordinates with shape (n,2).
      type(neighbor_list) :: nb
      logical, allocatable :: adj(:, :)
      integer :: i
      integer :: j
      integer :: l
      integer :: n
      real(dp) :: mx
      real(dp) :: my
      real(dp) :: radius
      real(dp) :: dl
      logical :: blocked

      n = size(coords, 1)
      if (size(coords, 2) < 2) then
         allocate(nb%neighbors(0))
         return
      end if
      allocate(adj(n, n))
      adj = .false.
      do i = 1, n - 1
         do j = i + 1, n
            mx = 0.5_dp * (coords(i, 1) + coords(j, 1))
            my = 0.5_dp * (coords(i, 2) + coords(j, 2))
            radius = euclidean_distance(coords(i, 1), coords(i, 2), mx, my)
            blocked = .false.
            do l = 1, n
               if (l == i .or. l == j) cycle
               dl = euclidean_distance(coords(l, 1), coords(l, 2), mx, my)
               if (dl < radius) then
                  blocked = .true.
                  exit
               end if
            end do
            if (.not. blocked) then
               adj(i, j) = .true.
               adj(j, i) = .true.
            end if
         end do
      end do
      nb = matrix_to_nb(adj)
   end function gabrielneigh

   pure function relative_neighborhood(coords) result(nb)
      real(dp), intent(in) :: coords(:, :) !! Planar point coordinates with shape (n,2).
      type(neighbor_list) :: nb
      logical, allocatable :: adj(:, :)
      integer :: i
      integer :: j
      integer :: l
      integer :: n
      real(dp) :: dij
      real(dp) :: dil
      real(dp) :: djl
      logical :: blocked

      n = size(coords, 1)
      if (size(coords, 2) < 2) then
         allocate(nb%neighbors(0))
         return
      end if
      allocate(adj(n, n))
      adj = .false.
      do i = 1, n - 1
         do j = i + 1, n
            dij = euclidean_distance(coords(i, 1), coords(i, 2), coords(j, 1), coords(j, 2))
            blocked = .false.
            do l = 1, n
               if (l == i .or. l == j) cycle
               dil = euclidean_distance(coords(i, 1), coords(i, 2), coords(l, 1), coords(l, 2))
               djl = euclidean_distance(coords(j, 1), coords(j, 2), coords(l, 1), coords(l, 2))
               if (dil < dij .and. djl < dij) then
                  blocked = .true.
                  exit
               end if
            end do
            if (.not. blocked) then
               adj(i, j) = .true.
               adj(j, i) = .true.
            end if
         end do
      end do
      nb = matrix_to_nb(adj)
   end function relative_neighborhood

   pure function tri2nb(coords, tolerance) result(nb)
      real(dp), intent(in) :: coords(:, :) !! Planar point coordinates with shape (n,2) used for Delaunay-neighbor construction.
      real(dp), intent(in), optional :: tolerance !! Relative tolerance for collinearity and empty-circumcircle comparisons.
      type(neighbor_list) :: nb
      logical, allocatable :: adj(:, :)
      integer :: i
      integer :: j
      integer :: k
      integer :: l
      integer :: n
      real(dp) :: tol
      real(dp) :: ax
      real(dp) :: ay
      real(dp) :: bx
      real(dp) :: by
      real(dp) :: cx
      real(dp) :: cy
      real(dp) :: det
      real(dp) :: ux
      real(dp) :: uy
      real(dp) :: r2
      real(dp) :: dl2
      logical :: empty_circle

      n = size(coords, 1)
      tol = 100.0_dp * epsilon(1.0_dp)
      if (present(tolerance)) tol = max(0.0_dp, tolerance)
      if (size(coords, 2) < 2) then
         allocate(nb%neighbors(0))
         return
      end if
      allocate(adj(n, n))
      adj = .false.
      do i = 1, n - 2
         ax = coords(i, 1)
         ay = coords(i, 2)
         do j = i + 1, n - 1
            bx = coords(j, 1)
            by = coords(j, 2)
            do k = j + 1, n
               cx = coords(k, 1)
               cy = coords(k, 2)
               det = 2.0_dp * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
               if (abs(det) <= tol * max(1.0_dp, abs(ax) + abs(ay) + abs(bx) + abs(by) + abs(cx) + abs(cy))) cycle
               ux = ((ax * ax + ay * ay) * (by - cy) + &
                  (bx * bx + by * by) * (cy - ay) + &
                  (cx * cx + cy * cy) * (ay - by)) / det
               uy = ((ax * ax + ay * ay) * (cx - bx) + &
                  (bx * bx + by * by) * (ax - cx) + &
                  (cx * cx + cy * cy) * (bx - ax)) / det
               r2 = (ax - ux) ** 2 + (ay - uy) ** 2
               empty_circle = .true.
               do l = 1, n
                  if (l == i .or. l == j .or. l == k) cycle
                  dl2 = (coords(l, 1) - ux) ** 2 + (coords(l, 2) - uy) ** 2
                  if (dl2 < r2 - tol * max(1.0_dp, r2)) then
                     empty_circle = .false.
                     exit
                  end if
               end do
               if (empty_circle) then
                  adj(i, j) = .true.
                  adj(j, i) = .true.
                  adj(i, k) = .true.
                  adj(k, i) = .true.
                  adj(j, k) = .true.
                  adj(k, j) = .true.
               end if
            end do
         end do
      end do
      nb = matrix_to_nb(adj)
   end function tri2nb

   pure function connected_components(nb) result(component)
      type(neighbor_list), intent(in) :: nb !! Neighbor graph whose weakly connected component labels are required.
      integer, allocatable :: component(:)
      logical, allocatable :: adj(:, :)
      integer, allocatable :: stack(:)
      logical, allocatable :: seen(:)
      integer :: n
      integer :: i
      integer :: j
      integer :: top
      integer :: v
      integer :: label

      n = nb%size()
      allocate(component(n), adj(n, n), stack(max(1, n)), seen(n))
      component = 0
      seen = .false.
      adj = nb_adjacency_matrix(nb)
      adj = adj .or. transpose(adj)
      label = 0
      do i = 1, n
         if (seen(i)) cycle
         label = label + 1
         top = 1
         stack(1) = i
         seen(i) = .true.
         do while (top > 0)
            v = stack(top)
            top = top - 1
            component(v) = label
            do j = 1, n
               if (adj(v, j) .and. .not. seen(j)) then
                  top = top + 1
                  stack(top) = j
                  seen(j) = .true.
               end if
            end do
         end do
      end do
   end function connected_components

   pure logical function is_symmetric_nb(nb) result(ok)
      type(neighbor_list), intent(in) :: nb !! Neighbor graph tested for reciprocal links.
      logical, allocatable :: adj(:, :)

      adj = nb_adjacency_matrix(nb)
      ok = all(adj .eqv. transpose(adj))
   end function is_symmetric_nb

   pure function make_symmetric_nb(nb) result(out)
      type(neighbor_list), intent(in) :: nb !! Neighbor graph whose directed links are converted to an undirected union.
      type(neighbor_list) :: out
      logical, allocatable :: adj(:, :)

      adj = nb_adjacency_matrix(nb)
      out = matrix_to_nb(adj .or. transpose(adj))
      out%self_included = nb%self_included
   end function make_symmetric_nb

   pure function include_self(nb) result(out)
      type(neighbor_list), intent(in) :: nb !! Neighbor graph to augment with one self-link per region.
      type(neighbor_list) :: out
      logical, allocatable :: adj(:, :)
      integer :: i

      adj = nb_adjacency_matrix(nb)
      do i = 1, size(adj, 1)
         adj(i, i) = .true.
      end do
      out = matrix_to_nb(adj)
      out%self_included = .true.
   end function include_self

   pure function remove_self(nb) result(out)
      type(neighbor_list), intent(in) :: nb !! Neighbor graph from which all self-links are removed.
      type(neighbor_list) :: out
      logical, allocatable :: adj(:, :)
      integer :: i

      adj = nb_adjacency_matrix(nb)
      do i = 1, size(adj, 1)
         adj(i, i) = .false.
      end do
      out = matrix_to_nb(adj)
      out%self_included = .false.
   end function remove_self

   pure function nb_union(a, b) result(out)
      type(neighbor_list), intent(in) :: a !! First neighbor graph in the set union; must have the same number of regions as b.
      type(neighbor_list), intent(in) :: b !! Second neighbor graph in the set union; must have the same number of regions as a.
      type(neighbor_list) :: out
      logical, allocatable :: aa(:, :)
      logical, allocatable :: bb(:, :)

      if (a%size() /= b%size()) then
         allocate(out%neighbors(0))
         return
      end if
      aa = nb_adjacency_matrix(a)
      bb = nb_adjacency_matrix(b)
      out = matrix_to_nb(aa .or. bb)
      out%self_included = a%self_included .or. b%self_included
   end function nb_union

   pure function nb_intersection(a, b) result(out)
      type(neighbor_list), intent(in) :: a !! First neighbor graph in the set intersection; must match b in size.
      type(neighbor_list), intent(in) :: b !! Second neighbor graph in the set intersection; must match a in size.
      type(neighbor_list) :: out
      logical, allocatable :: aa(:, :)
      logical, allocatable :: bb(:, :)

      if (a%size() /= b%size()) then
         allocate(out%neighbors(0))
         return
      end if
      aa = nb_adjacency_matrix(a)
      bb = nb_adjacency_matrix(b)
      out = matrix_to_nb(aa .and. bb)
      out%self_included = a%self_included .and. b%self_included
   end function nb_intersection

   pure function nb_difference(a, b) result(out)
      type(neighbor_list), intent(in) :: a !! Neighbor graph from which links present in b are removed.
      type(neighbor_list), intent(in) :: b !! Neighbor graph supplying links to remove from a; must match a in size.
      type(neighbor_list) :: out
      logical, allocatable :: aa(:, :)
      logical, allocatable :: bb(:, :)

      if (a%size() /= b%size()) then
         allocate(out%neighbors(0))
         return
      end if
      aa = nb_adjacency_matrix(a)
      bb = nb_adjacency_matrix(b)
      out = matrix_to_nb(aa .and. .not. bb)
      out%self_included = a%self_included .and. .not. b%self_included
   end function nb_difference

   pure function droplinks(nb, from, to, symmetric) result(out)
      type(neighbor_list), intent(in) :: nb !! Original neighbor graph.
      integer, intent(in) :: from(:) !! Source-region indices of links to remove.
      integer, intent(in) :: to(:) !! Destination-region indices paired elementwise with from.
      logical, intent(in), optional :: symmetric !! If true, remove the reverse of each supplied link as well.
      type(neighbor_list) :: out
      logical, allocatable :: adj(:, :)
      logical :: use_sym
      integer :: i
      integer :: n

      n = nb%size()
      use_sym = .false.
      if (present(symmetric)) use_sym = symmetric
      adj = nb_adjacency_matrix(nb)
      do i = 1, min(size(from), size(to))
         if (from(i) < 1 .or. from(i) > n .or. to(i) < 1 .or. to(i) > n) cycle
         adj(from(i), to(i)) = .false.
         if (use_sym) adj(to(i), from(i)) = .false.
      end do
      out = matrix_to_nb(adj)
      out%self_included = any([(adj(i, i), i = 1, n)])
   end function droplinks

   pure function addlinks(nb, from, to, symmetric) result(out)
      type(neighbor_list), intent(in) :: nb !! Original neighbor graph.
      integer, intent(in) :: from(:) !! Source-region indices of links to add.
      integer, intent(in) :: to(:) !! Destination-region indices paired elementwise with from.
      logical, intent(in), optional :: symmetric !! If true, add the reverse of each supplied link as well.
      type(neighbor_list) :: out
      logical, allocatable :: adj(:, :)
      logical :: use_sym
      integer :: i
      integer :: n

      n = nb%size()
      use_sym = .false.
      if (present(symmetric)) use_sym = symmetric
      adj = nb_adjacency_matrix(nb)
      do i = 1, min(size(from), size(to))
         if (from(i) < 1 .or. from(i) > n .or. to(i) < 1 .or. to(i) > n) cycle
         adj(from(i), to(i)) = .true.
         if (use_sym) adj(to(i), from(i)) = .true.
      end do
      out = matrix_to_nb(adj)
      out%self_included = any([(adj(i, i), i = 1, n)])
   end function addlinks

   pure function nblag(nb, maxlag) result(lags)
      type(neighbor_list), intent(in) :: nb !! Neighbor graph from which exact graph-distance lag neighbor lists are generated.
      integer, intent(in) :: maxlag !! Largest positive graph distance to return.
      type(neighbor_list), allocatable :: lags(:)
      real(dp), allocatable :: dist(:, :)
      logical, allocatable :: adj(:, :)
      integer :: h
      integer :: i
      integer :: j
      integer :: n

      n = nb%size()
      if (maxlag < 1) then
         allocate(lags(0))
         return
      end if
      allocate(lags(maxlag))
      dist = graph_distance_matrix(nb)
      do h = 1, maxlag
         allocate(adj(n, n))
         adj = .false.
         do i = 1, n
            do j = 1, n
               if (nint(dist(i, j)) == h) adj(i, j) = .true.
            end do
         end do
         lags(h) = matrix_to_nb(adj)
         deallocate(adj)
      end do
   end function nblag

   pure function nb_adjacency_matrix(nb) result(adj)
      type(neighbor_list), intent(in) :: nb !! Neighbor list converted to a logical n-by-n adjacency matrix.
      logical, allocatable :: adj(:, :)
      integer :: i
      integer :: j
      integer :: k
      integer :: n

      n = nb%size()
      allocate(adj(n, n))
      adj = .false.
      do i = 1, n
         if (.not. allocated(nb%neighbors(i)%values)) cycle
         do k = 1, size(nb%neighbors(i)%values)
            j = nb%neighbors(i)%values(k)
            if (j >= 1 .and. j <= n) adj(i, j) = .true.
         end do
      end do
   end function nb_adjacency_matrix

   pure function graph_distance_matrix(nb) result(dist)
      type(neighbor_list), intent(in) :: nb !! Neighbor graph whose all-pairs unweighted shortest-path distances are returned.
      real(dp), allocatable :: dist(:, :)
      logical, allocatable :: adj(:, :)
      integer, allocatable :: queue(:)
      integer, allocatable :: level(:)
      integer :: n
      integer :: source
      integer :: head
      integer :: tail
      integer :: v
      integer :: j

      n = nb%size()
      allocate(dist(n, n), adj(n, n), queue(max(1, n)), level(n))
      adj = nb_adjacency_matrix(nb)
      adj = adj .or. transpose(adj)
      dist = huge(1.0_dp)
      do source = 1, n
         level = -1
         head = 1
         tail = 1
         queue(1) = source
         level(source) = 0
         do while (head <= tail)
            v = queue(head)
            head = head + 1
            do j = 1, n
               if (adj(v, j) .and. level(j) < 0) then
                  level(j) = level(v) + 1
                  tail = tail + 1
                  queue(tail) = j
               end if
            end do
         end do
         do j = 1, n
            if (level(j) >= 0) dist(source, j) = real(level(j), dp)
         end do
      end do
   end function graph_distance_matrix

   pure function matrix_to_nb(adj) result(nb)
      logical, intent(in) :: adj(:, :) !! Logical adjacency matrix converted to compact neighbor vectors.
      type(neighbor_list) :: nb
      integer :: n
      integer :: i
      integer :: j
      integer :: k
      integer :: count_i

      n = size(adj, 1)
      if (size(adj, 2) /= n) then
         allocate(nb%neighbors(0))
         return
      end if
      allocate(nb%neighbors(n))
      nb%self_included = .false.
      do i = 1, n
         count_i = count(adj(i, :))
         allocate(nb%neighbors(i)%values(count_i))
         k = 0
         do j = 1, n
            if (adj(i, j)) then
               k = k + 1
               nb%neighbors(i)%values(k) = j
               if (i == j) nb%self_included = .true.
            end if
         end do
      end do
   end function matrix_to_nb

end module spdep_graph
