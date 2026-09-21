module geometry_delaunay
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use geometry_kinds, only : dp
   use geometry_types, only : delaunay_result, search_result
   use geometry_linalg, only : hyperplane_normal, next_combination, combination_count
   use geometry_linalg, only : simplex_measure, solve_linear, sort_integer
   implicit none
   private

   public :: delaunayn_points, cart2bary_points, bary2cart_points
   public :: tsearch_points, tsearchn_points

contains

   pure subroutine delaunayn_points(p, result, tol)
      real(dp), intent(in) :: p(:, :) !! Input points, one point per row in D-dimensional space.
      type(delaunay_result), intent(out) :: result !! Simplicial Delaunay triangulation with areas and neighbours.
      real(dp), intent(in), optional :: tol !! Geometric tolerance used for lower-hull classification.
      integer, parameter :: cap = 1000001
      integer :: n, d, ncomb, nt, i, j
      integer, allocatable :: c(:), tri_buf(:, :)
      real(dp), allocatable :: lift(:, :), verts(:, :), dist(:)
      real(dp) :: normal(size(p, 2) + 1), offset, eps, scale, jitter
      logical :: ok, ok_next

      n = size(p, 1)
      d = size(p, 2)
      result%success = .false.
      allocate(result%points(n, d))
      result%points = p
      if (d < 1 .or. n < d + 1) then
         allocate(result%tri(0, d + 1), result%areas(0), result%neighbours(0, d + 1))
         return
      end if
      ncomb = combination_count(n, d + 1, cap)
      if (ncomb >= cap) then
         allocate(result%tri(0, d + 1), result%areas(0), result%neighbours(0, d + 1))
         return
      end if
      allocate(lift(n, d + 1), verts(d + 1, d + 1), dist(n))
      lift(:, 1:d) = p
      scale = max(1.0_dp, maxval(abs(p)), maxval(sum(p*p, dim=2)))
      jitter = 1.0e-7_dp*scale/max(1, n)
      do i = 1, n
         lift(i, d + 1) = sum(p(i, :)*p(i, :)) + jitter*real(i*i, dp)
      end do
      if (present(tol)) then
         eps = max(tol, 0.0_dp)
      else
         eps = 1.0e-12_dp*max(1.0_dp, maxval(abs(lift)))
      end if

      allocate(c(d + 1), tri_buf(max(1, ncomb), d + 1))
      c = [(i, i=1, d + 1)]
      nt = 0
      do
         do i = 1, d + 1
            verts(i, :) = lift(c(i), :)
         end do
         call hyperplane_normal(verts, normal, offset, ok)
         if (ok) then
            dist = matmul(lift, normal) + offset
            if (maxval(dist) <= eps .or. minval(dist) >= -eps) then
               if (minval(dist) >= -eps) then
                  normal = -normal
                  offset = -offset
               end if
               if (normal(d + 1) < -100.0_dp*epsilon(1.0_dp)) then
                  if (simplex_measure(p(c, :)) > 100.0_dp*epsilon(1.0_dp)*scale) then
                     nt = nt + 1
                     tri_buf(nt, :) = c
                  end if
               end if
            end if
         end if
         call next_combination(c, n, ok_next)
         if (.not. ok_next) exit
      end do

      allocate(result%tri(nt, d + 1), result%areas(nt), result%neighbours(nt, d + 1))
      if (nt > 0) result%tri = tri_buf(1:nt, :)
      do i = 1, nt
         result%areas(i) = simplex_measure(p(result%tri(i, :), :))
      end do
      result%neighbours = 0
      do i = 1, nt
         do j = i + 1, nt
            call register_neighbours(result%tri, i, j, result%neighbours)
         end do
      end do
      result%success = nt > 0
   end subroutine delaunayn_points

   pure subroutine cart2bary_points(simplex, points, beta, success)
      real(dp), intent(in) :: simplex(:, :) !! (D+1)-by-D reference simplex coordinates.
      real(dp), intent(in) :: points(:, :) !! M-by-D Cartesian query points.
      real(dp), allocatable, intent(out) :: beta(:, :) !! M-by-(D+1) barycentric coordinates.
      logical, intent(out) :: success !! True when the reference simplex is nondegenerate.
      real(dp), allocatable :: a(:, :), rhs(:), x(:)
      integer :: d, m, i
      logical :: ok

      d = size(simplex, 2)
      m = size(points, 1)
      success = .false.
      if (size(simplex, 1) /= d + 1 .or. size(points, 2) /= d) then
         allocate(beta(0, 0))
         return
      end if
      allocate(beta(m, d + 1), a(d, d), rhs(d), x(d))
      do i = 1, d
         a(i, :) = simplex(i, :) - simplex(d + 1, :)
      end do
      a = transpose(a)
      do i = 1, m
         rhs = points(i, :) - simplex(d + 1, :)
         call solve_linear(a, rhs, x, ok)
         if (.not. ok) then
            beta = ieee_value(0.0_dp, ieee_quiet_nan)
            return
         end if
         beta(i, 1:d) = x
         beta(i, d + 1) = 1.0_dp - sum(x)
      end do
      success = .true.
   end subroutine cart2bary_points

   pure subroutine bary2cart_points(simplex, beta, points, success)
      real(dp), intent(in) :: simplex(:, :) !! (D+1)-by-D reference simplex coordinates.
      real(dp), intent(in) :: beta(:, :) !! M-by-(D+1) barycentric coordinates.
      real(dp), allocatable, intent(out) :: points(:, :) !! M-by-D Cartesian coordinates.
      logical, intent(out) :: success !! True when dimensions are compatible.
      integer :: d

      d = size(simplex, 2)
      success = .false.
      if (size(simplex, 1) /= d + 1 .or. size(beta, 2) /= d + 1) then
         allocate(points(0, 0))
         return
      end if
      allocate(points(size(beta, 1), d))
      points = matmul(beta, simplex)
      success = .true.
   end subroutine bary2cart_points

   pure subroutine tsearchn_points(p, tri, query, result, tol)
      real(dp), intent(in) :: p(:, :) !! Triangulation vertices, one D-dimensional point per row.
      integer, intent(in) :: tri(:, :) !! Simplex vertex indices, D+1 entries per row.
      real(dp), intent(in) :: query(:, :) !! Query points whose enclosing simplex is requested.
      type(search_result), intent(out) :: result !! First enclosing simplex and barycentric coordinates per query.
      real(dp), intent(in), optional :: tol !! Nonnegative barycentric tolerance for point inclusion.
      real(dp), allocatable :: beta(:, :), one_point(:, :)
      real(dp) :: eps, nan_value
      logical :: ok
      integer :: i, q, d

      d = size(p, 2)
      nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (present(tol)) then
         eps = max(tol, 0.0_dp)
      else
         eps = 1.0e-12_dp
      end if
      allocate(result%idx(size(query, 1)), result%bary(size(query, 1), d + 1))
      result%idx = 0
      result%bary = nan_value
      if (size(tri, 2) /= d + 1 .or. size(query, 2) /= d) return
      allocate(one_point(1, d))
      do i = 1, size(tri, 1)
         if (any(tri(i, :) < 1) .or. any(tri(i, :) > size(p, 1))) cycle
         do q = 1, size(query, 1)
            if (result%idx(q) /= 0) cycle
            one_point(1, :) = query(q, :)
            call cart2bary_points(p(tri(i, :), :), one_point, beta, ok)
            if (.not. ok) cycle
            if (all(beta(1, :) >= -eps)) then
               result%idx(q) = i
               result%bary(q, :) = beta(1, :)
            end if
         end do
      end do
   end subroutine tsearchn_points

   pure subroutine tsearch_points(x, y, tri, xi, yi, result, tol)
      real(dp), intent(in) :: x(:) !! X-coordinates of triangulation vertices.
      real(dp), intent(in) :: y(:) !! Y-coordinates of triangulation vertices.
      integer, intent(in) :: tri(:, :) !! Triangles, with three one-based vertex indices per row.
      real(dp), intent(in) :: xi(:) !! X-coordinates of query points.
      real(dp), intent(in) :: yi(:) !! Y-coordinates of query points.
      type(search_result), intent(out) :: result !! Last enclosing triangle and barycentric coordinates per query.
      real(dp), intent(in), optional :: tol !! Nonnegative barycentric tolerance for boundary inclusion.
      real(dp), allocatable :: p(:, :), qmat(:, :), beta(:, :), one_point(:, :)
      real(dp) :: eps, nan_value
      logical :: ok
      integer :: i, q

      nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (present(tol)) then
         eps = max(tol, 0.0_dp)
      else
         eps = 1.0e-12_dp
      end if
      allocate(result%idx(size(xi)), result%bary(size(xi), 3))
      result%idx = 0
      result%bary = nan_value
      if (size(x) /= size(y) .or. size(xi) /= size(yi) .or. size(tri, 2) /= 3) return
      allocate(p(size(x), 2), qmat(size(xi), 2), one_point(1, 2))
      p(:, 1) = x
      p(:, 2) = y
      qmat(:, 1) = xi
      qmat(:, 2) = yi
      do i = 1, size(tri, 1)
         if (any(tri(i, :) < 1) .or. any(tri(i, :) > size(p, 1))) cycle
         do q = 1, size(qmat, 1)
            one_point(1, :) = qmat(q, :)
            call cart2bary_points(p(tri(i, :), :), one_point, beta, ok)
            if (.not. ok) cycle
            if (all(beta(1, :) >= -eps)) then
               result%idx(q) = i
               result%bary(q, :) = beta(1, :)
            end if
         end do
      end do
   end subroutine tsearch_points

   pure subroutine register_neighbours(tri, ia, ib, neighbours)
      integer, intent(in) :: tri(:, :) !! Complete simplex index table.
      integer, intent(in) :: ia !! First simplex row considered for adjacency.
      integer, intent(in) :: ib !! Second simplex row considered for adjacency.
      integer, intent(inout) :: neighbours(:, :) !! Adjacency table updated for a shared facet.
      integer, allocatable :: a(:), b(:)
      integer :: d, ka, kb, na, nb, i

      d = size(tri, 2) - 1
      do ka = 1, d + 1
         allocate(a(d))
         na = 0
         do i = 1, d + 1
            if (i == ka) cycle
            na = na + 1
            a(na) = tri(ia, i)
         end do
         call sort_integer(a)
         do kb = 1, d + 1
            allocate(b(d))
            nb = 0
            do i = 1, d + 1
               if (i == kb) cycle
               nb = nb + 1
               b(nb) = tri(ib, i)
            end do
            call sort_integer(b)
            if (all(a == b)) then
               neighbours(ia, ka) = ib
               neighbours(ib, kb) = ia
               deallocate(b)
               deallocate(a)
               return
            end if
            deallocate(b)
         end do
         deallocate(a)
      end do
   end subroutine register_neighbours

end module geometry_delaunay
