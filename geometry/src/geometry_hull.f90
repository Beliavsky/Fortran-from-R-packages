module geometry_hull
   use geometry_kinds, only : dp
   use geometry_types, only : hull_result, intersection_result
   use geometry_linalg, only : cross3, hyperplane_normal, next_combination
   use geometry_linalg, only : combination_count, simplex_facet_measure, simplex_measure
   use geometry_linalg, only : sort_real_with_index, solve_linear, row_norm
   implicit none
   private

   public :: convhulln_points, halfspacen_vertices, inhulln_points
   public :: feasible_point_hulls, intersectn_points

contains

   pure subroutine convhulln_points(p, hull, tol)
      real(dp), intent(in) :: p(:, :) !! Input points, one point per row and coordinate per column.
      type(hull_result), intent(out) :: hull !! Convex-hull facets, normals, generalized area, and volume.
      real(dp), intent(in), optional :: tol !! Geometric tolerance; defaults relative to coordinate scale.
      real(dp) :: eps
      integer :: d, n

      n = size(p, 1)
      d = size(p, 2)
      hull%success = .false.
      allocate(hull%points(n, d))
      hull%points = p
      if (n < d + 1 .or. d < 1) then
         allocate(hull%facets(0, max(1, d)), hull%normals(0, d + 1))
         return
      end if
      if (present(tol)) then
         eps = max(tol, 0.0_dp)
      else
         eps = 1.0e-10_dp*max(1.0_dp, maxval(abs(p)))
      end if
      if (d == 1) then
         call hull_1d(p, hull)
      else if (d == 2) then
         call hull_2d(p, hull, eps)
      else
         call hull_nd(p, hull, eps)
      end if
   end subroutine convhulln_points

   pure subroutine hull_1d(p, hull)
      real(dp), intent(in) :: p(:, :) !! One-dimensional input points.
      type(hull_result), intent(inout) :: hull !! Hull result updated with endpoint facets and extent.
      integer :: imin, imax

      imin = minloc(p(:, 1), dim=1)
      imax = maxloc(p(:, 1), dim=1)
      allocate(hull%facets(2, 1), hull%normals(2, 2))
      hull%facets(:, 1) = [imin, imax]
      hull%normals(1, :) = [-1.0_dp, p(imin, 1)]
      hull%normals(2, :) = [1.0_dp, -p(imax, 1)]
      hull%area = 2.0_dp
      hull%volume = p(imax, 1) - p(imin, 1)
      hull%success = hull%volume >= 0.0_dp
   end subroutine hull_1d

   pure subroutine hull_2d(p, hull, tol)
      real(dp), intent(in) :: p(:, :) !! Two-dimensional input points.
      type(hull_result), intent(inout) :: hull !! Hull result updated with ordered polygon edges and normals.
      real(dp), intent(in) :: tol !! Geometric tolerance used for collinearity tests.
      integer, allocatable :: ord(:), stack(:), poly(:)
      real(dp) :: cr, dx, dy, len, signed_area
      integer :: i, m, n, top, j

      n = size(p, 1)
      allocate(ord(n), stack(2*n), poly(2*n))
      ord = [(i, i=1, n)]
      call sort_points_2d(p, ord)

      top = 0
      do i = 1, n
         do while (top >= 2)
            cr = cross2(p(stack(top), :) - p(stack(top - 1), :), p(ord(i), :) - p(stack(top), :))
            if (cr > tol) exit
            top = top - 1
         end do
         top = top + 1
         stack(top) = ord(i)
      end do
      m = top
      do i = n - 1, 1, -1
         do while (top > m)
            cr = cross2(p(stack(top), :) - p(stack(top - 1), :), p(ord(i), :) - p(stack(top), :))
            if (cr > tol) exit
            top = top - 1
         end do
         top = top + 1
         stack(top) = ord(i)
      end do
      if (top > 1) top = top - 1
      if (top < 3) then
         allocate(hull%facets(0, 2), hull%normals(0, 3))
         return
      end if
      poly(1:top) = stack(1:top)
      allocate(hull%facets(top, 2), hull%normals(top, 3))
      hull%area = 0.0_dp
      hull%volume = 0.0_dp
      signed_area = 0.0_dp
      do i = 1, top
         j = i + 1
         if (j > top) j = 1
         hull%facets(i, :) = [poly(i), poly(j)]
         dx = p(poly(j), 1) - p(poly(i), 1)
         dy = p(poly(j), 2) - p(poly(i), 2)
         len = sqrt(dx*dx + dy*dy)
         hull%normals(i, 1:2) = [dy/len, -dx/len]
         hull%normals(i, 3) = -dot_product(hull%normals(i, 1:2), p(poly(i), :))
         hull%area = hull%area + len
         signed_area = signed_area + p(poly(i), 1)*p(poly(j), 2) - p(poly(j), 1)*p(poly(i), 2)
      end do
      hull%volume = 0.5_dp*abs(signed_area)
      hull%success = .true.
   end subroutine hull_2d

   pure subroutine hull_nd(p, hull, tol)
      real(dp), intent(in) :: p(:, :) !! Input points in dimension three or greater.
      type(hull_result), intent(inout) :: hull !! Hull result updated with supporting simplicial facets.
      real(dp), intent(in) :: tol !! Geometric tolerance for supporting-plane and coplanarity tests.
      integer, parameter :: cap = 500001
      integer :: n, d, ncomb, i, j, np, nf, m, plane_id
      integer, allocatable :: c(:), plane_seed(:, :), facet_buf(:, :), facet_plane(:), vids(:), idx(:)
      real(dp), allocatable :: planes(:, :), verts(:, :), dist(:), angles(:)
      real(dp), allocatable :: centroid(:), e1(:), e2(:)
      real(dp) :: offset, maxd, mind, plane_tol, volume_piece
      real(dp) :: normal(size(p, 2))
      logical :: ok, is_new, ok_next

      n = size(p, 1)
      d = size(p, 2)
      ncomb = combination_count(n, d, cap)
      if (ncomb >= cap) then
         allocate(hull%facets(0, d), hull%normals(0, d + 1))
         return
      end if
      allocate(c(d), planes(max(1, ncomb), d + 1), plane_seed(max(1, ncomb), d))
      allocate(verts(d, d), dist(n))
      np = 0
      c = [(i, i=1, d)]
      do
         do i = 1, d
            verts(i, :) = p(c(i), :)
         end do
         call hyperplane_normal(verts, normal, offset, ok)
         if (ok) then
            dist = matmul(p, normal) + offset
            maxd = maxval(dist)
            mind = minval(dist)
            if (maxd <= tol .or. mind >= -tol) then
               if (mind >= -tol) then
                  normal = -normal
                  offset = -offset
                  dist = -dist
               end if
               is_new = .true.
               do j = 1, np
                  if (dot_product(normal, planes(j, 1:d)) > 1.0_dp - 1.0e-8_dp .and. &
                      abs(offset - planes(j, d + 1)) <= 10.0_dp*tol) then
                     is_new = .false.
                     exit
                  end if
               end do
               if (is_new) then
                  np = np + 1
                  planes(np, 1:d) = normal
                  planes(np, d + 1) = offset
                  plane_seed(np, :) = c
               end if
            end if
         end if
         call next_combination(c, n, ok_next)
         if (.not. ok_next) exit
      end do
      if (np == 0) then
         allocate(hull%facets(0, d), hull%normals(0, d + 1))
         return
      end if

      allocate(facet_buf(max(1, ncomb), d), facet_plane(max(1, ncomb)))
      nf = 0
      plane_tol = max(10.0_dp*tol, 1.0e-10_dp*max(1.0_dp, maxval(abs(p))))
      if (d == 3) then
         allocate(vids(n), idx(n), angles(n), centroid(3), e1(3), e2(3))
         do plane_id = 1, np
            dist = matmul(p, planes(plane_id, 1:d)) + planes(plane_id, d + 1)
            m = 0
            do i = 1, n
               if (abs(dist(i)) <= plane_tol) then
                  m = m + 1
                  vids(m) = i
               end if
            end do
            if (m < 3) cycle
            centroid = 0.0_dp
            do i = 1, m
               centroid = centroid + p(vids(i), :)
            end do
            centroid = centroid/real(m, dp)
            e1 = 0.0_dp
            do i = 1, m
               e1 = p(vids(i), :) - centroid
               if (row_norm(e1) > plane_tol) exit
            end do
            if (row_norm(e1) <= plane_tol) cycle
            e1 = e1/row_norm(e1)
            normal = planes(plane_id, 1:d)
            e2 = cross3(normal, e1)
            do i = 1, m
               angles(i) = atan2(dot_product(p(vids(i), :) - centroid, e2), &
                                  dot_product(p(vids(i), :) - centroid, e1))
               idx(i) = vids(i)
            end do
            call sort_real_with_index(angles(1:m), idx(1:m))
            do i = 2, m - 1
               nf = nf + 1
               facet_buf(nf, :) = [idx(1), idx(i), idx(i + 1)]
               facet_plane(nf) = plane_id
            end do
         end do
      else
         do plane_id = 1, np
            nf = nf + 1
            facet_buf(nf, :) = plane_seed(plane_id, :)
            facet_plane(nf) = plane_id
         end do
      end if

      allocate(hull%facets(nf, d), hull%normals(nf, d + 1))
      if (nf > 0) then
         hull%facets = facet_buf(1:nf, :)
         do i = 1, nf
            hull%normals(i, :) = planes(facet_plane(i), :)
         end do
      end if
      hull%area = 0.0_dp
      hull%volume = 0.0_dp
      if (nf > 0) then
         if (.not. allocated(centroid)) allocate(centroid(d))
         centroid = sum(p, dim=1)/real(n, dp)
         do i = 1, nf
            do j = 1, d
               verts(j, :) = p(hull%facets(i, j), :)
            end do
            hull%area = hull%area + simplex_facet_measure(verts)
            volume_piece = simplex_with_centroid_measure(verts, centroid)
            hull%volume = hull%volume + volume_piece
         end do
      end if
      hull%success = nf > 0
   end subroutine hull_nd

   pure subroutine halfspacen_vertices(halfspaces, feasible, vertices, success, tol)
      real(dp), intent(in) :: halfspaces(:, :) !! Rows of normal coefficients followed by hyperplane offsets.
      real(dp), intent(in) :: feasible(:) !! Point expected to satisfy every halfspace, used for validation.
      real(dp), allocatable, intent(out) :: vertices(:, :) !! Vertices formed by intersections of active hyperplanes.
      logical, intent(out) :: success !! True when at least one bounded-intersection vertex was found.
      real(dp), intent(in), optional :: tol !! Feasibility and duplicate tolerance.
      integer, parameter :: cap = 1000001
      integer :: m, d, ncomb, nv, i, j
      integer, allocatable :: c(:)
      real(dp), allocatable :: a(:, :), b(:), x(:), buf(:, :)
      real(dp) :: eps, dist, scale
      logical :: ok, duplicate, ok_next

      m = size(halfspaces, 1)
      d = size(halfspaces, 2) - 1
      success = .false.
      if (d < 1 .or. size(feasible) /= d .or. m < d) then
         allocate(vertices(0, max(1, d)))
         return
      end if
      scale = max(1.0_dp, maxval(abs(halfspaces)))
      if (present(tol)) then
         eps = max(tol, 0.0_dp)
      else
         eps = 1.0e-9_dp*scale
      end if
      if (any(matmul(halfspaces(:, 1:d), feasible) + halfspaces(:, d + 1) > 100.0_dp*eps)) then
         allocate(vertices(0, d))
         return
      end if
      ncomb = combination_count(m, d, cap)
      if (ncomb >= cap) then
         allocate(vertices(0, d))
         return
      end if
      allocate(c(d), a(d, d), b(d), x(d), buf(max(1, ncomb), d))
      nv = 0
      c = [(i, i=1, d)]
      do
         do i = 1, d
            a(i, :) = halfspaces(c(i), 1:d)
            b(i) = -halfspaces(c(i), d + 1)
         end do
         call solve_linear(a, b, x, ok)
         if (ok) then
            if (all(matmul(halfspaces(:, 1:d), x) + halfspaces(:, d + 1) <= 10.0_dp*eps)) then
               duplicate = .false.
               do j = 1, nv
                  dist = sqrt(sum((x - buf(j, :))**2))
                  if (dist <= 20.0_dp*eps) then
                     duplicate = .true.
                     exit
                  end if
               end do
               if (.not. duplicate) then
                  nv = nv + 1
                  buf(nv, :) = x
               end if
            end if
         end if
         call next_combination(c, m, ok_next)
         if (.not. ok_next) exit
      end do
      allocate(vertices(nv, d))
      if (nv > 0) vertices = buf(1:nv, :)
      success = nv > 0
   end subroutine halfspacen_vertices

   pure subroutine inhulln_points(hull, p, inside, tol)
      type(hull_result), intent(in) :: hull !! Convex hull containing outward facet-normal inequalities.
      real(dp), intent(in) :: p(:, :) !! Query points, one point per row.
      logical, allocatable, intent(out) :: inside(:) !! True only for points strictly inside all hull facets.
      real(dp), intent(in), optional :: tol !! Boundary tolerance; points on facets are treated as outside.
      real(dp) :: eps
      integer :: i, d

      allocate(inside(size(p, 1)))
      inside = .false.
      if (.not. hull%success) return
      d = size(p, 2)
      if (.not. allocated(hull%normals)) return
      if (size(hull%normals, 2) /= d + 1) return
      if (present(tol)) then
         eps = max(tol, 0.0_dp)
      else
         eps = 1.0e-12_dp*max(1.0_dp, maxval(abs(hull%points)))
      end if
      do i = 1, size(p, 1)
         inside(i) = all(matmul(hull%normals(:, 1:d), p(i, :)) + hull%normals(:, d + 1) < -eps)
      end do
   end subroutine inhulln_points

   pure subroutine feasible_point_hulls(ch1, ch2, tol, point, success)
      type(hull_result), intent(in) :: ch1 !! First convex hull with outward normal inequalities.
      type(hull_result), intent(in) :: ch2 !! Second convex hull with outward normal inequalities.
      real(dp), intent(in) :: tol !! Required inward margin from every facet.
      real(dp), allocatable, intent(out) :: point(:) !! Feasible point lying in both hulls when found.
      logical, intent(out) :: success !! True when a feasible point was found.
      real(dp), allocatable :: hs(:, :), verts(:, :), candidate(:)
      integer :: d, m1, m2
      logical :: hv

      success = .false.
      if (.not. ch1%success .or. .not. ch2%success) then
         allocate(point(0))
         return
      end if
      d = size(ch1%points, 2)
      if (size(ch2%points, 2) /= d) then
         allocate(point(0))
         return
      end if
      m1 = size(ch1%normals, 1)
      m2 = size(ch2%normals, 1)
      allocate(hs(m1 + m2, d + 1))
      hs(1:m1, :) = ch1%normals
      hs(m1 + 1:, :) = ch2%normals
      allocate(candidate(d))
      candidate = 0.5_dp*(sum(ch1%points, dim=1)/real(size(ch1%points, 1), dp) + &
                          sum(ch2%points, dim=1)/real(size(ch2%points, 1), dp))
      call project_feasible(hs, tol, candidate, success)
      if (.not. success) then
         call halfspacen_vertices(hs, candidate, verts, hv, max(tol, 1.0e-10_dp))
         if (hv) then
            candidate = sum(verts, dim=1)/real(size(verts, 1), dp)
            call project_feasible(hs, tol, candidate, success)
         end if
      end if
      allocate(point(d))
      point = candidate
   end subroutine feasible_point_hulls

   pure subroutine intersectn_points(ps1, ps2, tol, result)
      real(dp), intent(in) :: ps1(:, :) !! First point set whose convex hull participates in the intersection.
      real(dp), intent(in) :: ps2(:, :) !! Second point set whose convex hull participates in the intersection.
      real(dp), intent(in) :: tol !! Duplicate and feasibility tolerance for intersection vertices.
      type(intersection_result), intent(out) :: result !! Original hulls, intersection points, and intersection hull.
      real(dp), allocatable :: fp(:), hs(:, :), verts(:, :)
      logical :: ok, hv
      integer :: d, m1, m2

      result%success = .false.
      call convhulln_points(ps1, result%ch1)
      call convhulln_points(ps2, result%ch2)
      if (.not. result%ch1%success .or. .not. result%ch2%success) then
         allocate(result%points(0, size(ps1, 2)))
         return
      end if
      d = size(ps1, 2)
      if (size(ps2, 2) /= d) then
         allocate(result%points(0, d))
         return
      end if
      call feasible_point_hulls(result%ch1, result%ch2, tol, fp, ok)
      if (.not. ok) then
         allocate(result%points(0, d))
         result%ch%volume = 0.0_dp
         return
      end if
      m1 = size(result%ch1%normals, 1)
      m2 = size(result%ch2%normals, 1)
      allocate(hs(m1 + m2, d + 1))
      hs(1:m1, :) = result%ch1%normals
      hs(m1 + 1:, :) = result%ch2%normals
      call halfspacen_vertices(hs, fp, verts, hv, max(tol, 1.0e-10_dp))
      if (.not. hv) then
         allocate(result%points(0, d))
         result%ch%volume = 0.0_dp
         return
      end if
      allocate(result%points(size(verts, 1), d))
      result%points = verts
      call convhulln_points(verts, result%ch)
      result%success = result%ch%success
   end subroutine intersectn_points

   pure subroutine sort_points_2d(p, idx)
      real(dp), intent(in) :: p(:, :) !! Two-dimensional points used as lexicographic sort keys.
      integer, intent(inout) :: idx(:) !! Point indices sorted by x then y coordinates.
      integer :: i, j, key
      do i = 2, size(idx)
         key = idx(i)
         j = i - 1
         do while (j >= 1)
            if (p(idx(j), 1) < p(key, 1)) exit
            if (p(idx(j), 1) <= p(key, 1) .and. p(idx(j), 2) <= p(key, 2)) exit
            idx(j + 1) = idx(j)
            j = j - 1
         end do
         idx(j + 1) = key
      end do
   end subroutine sort_points_2d

   pure function cross2(a, b) result(value)
      real(dp), intent(in) :: a(2) !! First two-dimensional vector.
      real(dp), intent(in) :: b(2) !! Second two-dimensional vector.
      real(dp) :: value
      value = a(1)*b(2) - a(2)*b(1)
   end function cross2

   pure function simplex_with_centroid_measure(facet_vertices, centroid) result(value)
      real(dp), intent(in) :: facet_vertices(:, :) !! D vertices defining a boundary (D-1)-simplex.
      real(dp), intent(in) :: centroid(:) !! Interior reference point used to form a D-simplex.
      real(dp) :: value
      real(dp), allocatable :: simplex(:, :)
      integer :: d
      d = size(centroid)
      allocate(simplex(d + 1, d))
      simplex(1:d, :) = facet_vertices
      simplex(d + 1, :) = centroid
      value = simplex_measure(simplex)
   end function simplex_with_centroid_measure

   pure subroutine project_feasible(halfspaces, margin, point, success)
      real(dp), intent(in) :: halfspaces(:, :) !! Halfspace rows normal coefficients followed by offsets.
      real(dp), intent(in) :: margin !! Required inward margin for every inequality.
      real(dp), intent(inout) :: point(:) !! Candidate point projected iteratively into violated halfspaces.
      logical, intent(out) :: success !! True when the final point satisfies every inequality with margin.
      integer :: d, i, iter
      real(dp) :: violation, denom, worst

      d = size(point)
      success = .false.
      do iter = 1, 1000
         worst = -huge(1.0_dp)
         do i = 1, size(halfspaces, 1)
            violation = dot_product(halfspaces(i, 1:d), point) + halfspaces(i, d + 1) + margin
            worst = max(worst, violation)
            if (violation > 0.0_dp) then
               denom = sum(halfspaces(i, 1:d)**2)
               if (denom > 0.0_dp) point = point - violation*halfspaces(i, 1:d)/denom
            end if
         end do
         if (worst <= 1.0e-10_dp*max(1.0_dp, maxval(abs(point)))) then
            success = .true.
            return
         end if
      end do
   end subroutine project_feasible

end module geometry_hull
