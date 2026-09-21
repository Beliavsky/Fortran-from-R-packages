module geometry_mesh
   use geometry_kinds, only : dp
   use geometry_types, only : delaunay_result
   use geometry_linalg, only : cross3, sort_integer, row_norm
   use geometry_delaunay, only : delaunayn_points
   implicit none
   private

   public :: surf_triangles, distmesh2d_circle, distmeshnd_sphere

contains

   pure subroutine surf_triangles(p, tetra, triangles)
      real(dp), intent(in) :: p(:, :) !! Three-dimensional mesh vertices, one point per row.
      integer, intent(in) :: tetra(:, :) !! Tetrahedra, four one-based vertex indices per row.
      integer, allocatable, intent(out) :: triangles(:, :) !! Oriented boundary triangles occurring on one tetrahedron.
      integer, allocatable :: faces(:, :), sorted_faces(:, :), opposite(:), keep(:), work(:)
      integer :: nt, nf, i, j, count, dup_count
      real(dp) :: v1(3), v2(3), v3(3), cp(3)

      if (size(p, 2) /= 3 .or. size(tetra, 2) /= 4) then
         allocate(triangles(0, 3))
         return
      end if
      nt = size(tetra, 1)
      nf = 4*nt
      allocate(faces(nf, 3), sorted_faces(nf, 3), opposite(nf), keep(nf), work(3))
      do i = 1, nt
         faces(4*i - 3, :) = tetra(i, [1, 2, 3])
         opposite(4*i - 3) = tetra(i, 4)
         faces(4*i - 2, :) = tetra(i, [1, 2, 4])
         opposite(4*i - 2) = tetra(i, 3)
         faces(4*i - 1, :) = tetra(i, [1, 3, 4])
         opposite(4*i - 1) = tetra(i, 2)
         faces(4*i, :) = tetra(i, [2, 3, 4])
         opposite(4*i) = tetra(i, 1)
      end do
      sorted_faces = faces
      do i = 1, nf
         work = sorted_faces(i, :)
         call sort_integer(work)
         sorted_faces(i, :) = work
      end do
      count = 0
      do i = 1, nf
         dup_count = 0
         do j = 1, nf
            if (all(sorted_faces(i, :) == sorted_faces(j, :))) dup_count = dup_count + 1
         end do
         if (dup_count == 1) then
            count = count + 1
            keep(count) = i
         end if
      end do
      allocate(triangles(count, 3))
      do j = 1, count
         i = keep(j)
         triangles(j, :) = faces(i, :)
         v1 = p(triangles(j, 2), :) - p(triangles(j, 1), :)
         v2 = p(triangles(j, 3), :) - p(triangles(j, 1), :)
         v3 = p(opposite(i), :) - p(triangles(j, 1), :)
         cp = cross3(v1, v2)
         if (dot_product(cp, v3) > 0.0_dp) then
            work(1) = triangles(j, 2)
            triangles(j, 2) = triangles(j, 3)
            triangles(j, 3) = work(1)
         end if
      end do
   end subroutine surf_triangles

   pure subroutine distmesh2d_circle(h0, bbox, radius, maxiter, points)
      real(dp), intent(in) :: h0 !! Initial nominal node spacing for the hexagonal grid.
      real(dp), intent(in) :: bbox(2, 2) !! Bounding box rows min/max and columns x/y.
      real(dp), intent(in) :: radius !! Radius of the circular signed-distance region.
      integer, intent(in) :: maxiter !! Maximum relaxation iterations; nonpositive returns the initial grid.
      real(dp), allocatable, intent(out) :: points(:, :) !! Relaxed mesh-node coordinates inside the circle.
      real(dp), allocatable :: raw(:, :)
      integer :: nx, ny, i, j, n
      real(dp) :: x, y, dy

      if (h0 <= 0.0_dp .or. radius <= 0.0_dp) then
         allocate(points(0, 2))
         return
      end if
      dy = h0*sqrt(3.0_dp)/2.0_dp
      nx = int((bbox(2, 1) - bbox(1, 1))/h0) + 2
      ny = int((bbox(2, 2) - bbox(1, 2))/dy) + 2
      allocate(raw(max(1, nx*ny), 2))
      n = 0
      do j = 0, ny - 1
         y = bbox(1, 2) + real(j, dp)*dy
         do i = 0, nx - 1
            x = bbox(1, 1) + real(i, dp)*h0
            if (mod(j, 2) == 1) x = x + 0.5_dp*h0
            if (x > bbox(2, 1) + h0 .or. y > bbox(2, 2) + dy) cycle
            if (sqrt(x*x + y*y) <= radius + 0.001_dp*h0) then
               n = n + 1
               raw(n, :) = [x, y]
            end if
         end do
      end do
      allocate(points(n, 2))
      if (n > 0) points = raw(1:n, :)
      if (n > 3 .and. maxiter > 0) call relax_sphere_mesh(points, radius, h0, maxiter)
   end subroutine distmesh2d_circle

   pure subroutine distmeshnd_sphere(start_points, radius, h0, maxiter, points)
      real(dp), intent(in) :: start_points(:, :) !! Initial nodes for the N-dimensional spherical mesh.
      real(dp), intent(in) :: radius !! Radius of the N-dimensional sphere.
      real(dp), intent(in) :: h0 !! Nominal edge length used to scale convergence and forces.
      integer, intent(in) :: maxiter !! Maximum number of Delaunay force-relaxation iterations.
      real(dp), allocatable, intent(out) :: points(:, :) !! Relaxed nodes projected into the sphere.
      integer :: i

      allocate(points(size(start_points, 1), size(start_points, 2)))
      points = start_points
      do i = 1, size(points, 1)
         if (row_norm(points(i, :)) > radius .and. row_norm(points(i, :)) > 0.0_dp) then
            points(i, :) = radius*points(i, :)/row_norm(points(i, :))
         end if
      end do
      if (size(points, 1) > size(points, 2) + 1 .and. maxiter > 0) then
         call relax_sphere_mesh(points, radius, h0, maxiter)
      end if
   end subroutine distmeshnd_sphere

   pure subroutine relax_sphere_mesh(points, radius, h0, maxiter)
      real(dp), intent(inout) :: points(:, :) !! Mesh nodes updated in place by Delaunay edge forces.
      real(dp), intent(in) :: radius !! Radius of the spherical region used for boundary projection.
      real(dp), intent(in) :: h0 !! Nominal edge length used in the convergence criterion.
      integer, intent(in) :: maxiter !! Maximum number of relaxation iterations.
      type(delaunay_result) :: dt
      integer, allocatable :: bars(:, :)
      real(dp), allocatable :: force(:, :), vec(:), lengths(:)
      real(dp) :: l0, len, f, step, nrm
      integer :: iter, i, nb, d, a, b

      d = size(points, 2)
      step = 0.2_dp
      do iter = 1, maxiter
         call delaunayn_points(points, dt)
         if (.not. dt%success) exit
         call simplex_edges(dt%tri, bars)
         nb = size(bars, 1)
         if (nb == 0) exit
         allocate(lengths(nb), force(size(points, 1), d), vec(d))
         force = 0.0_dp
         do i = 1, nb
            vec = points(bars(i, 1), :) - points(bars(i, 2), :)
            lengths(i) = row_norm(vec)
         end do
         l0 = 1.2_dp*sqrt(sum(lengths*lengths)/real(nb, dp))
         do i = 1, nb
            a = bars(i, 1)
            b = bars(i, 2)
            vec = points(a, :) - points(b, :)
            len = max(lengths(i), 100.0_dp*epsilon(1.0_dp))
            f = max(l0 - len, 0.0_dp)
            force(a, :) = force(a, :) + f*vec/len
            force(b, :) = force(b, :) - f*vec/len
         end do
         points = points + step*force
         do i = 1, size(points, 1)
            nrm = row_norm(points(i, :))
            if (nrm > radius .and. nrm > 0.0_dp) points(i, :) = radius*points(i, :)/nrm
         end do
         if (maxval(sqrt(sum((step*force)**2, dim=2)))/max(h0, epsilon(1.0_dp)) < 0.001_dp) then
            deallocate(lengths, force, vec, bars)
            exit
         end if
         deallocate(lengths, force, vec, bars)
      end do
   end subroutine relax_sphere_mesh

   pure subroutine simplex_edges(tri, bars)
      integer, intent(in) :: tri(:, :) !! Simplices from which all undirected edges are extracted.
      integer, allocatable, intent(out) :: bars(:, :) !! Sorted unique pairs of vertex indices.
      integer, allocatable :: buf(:, :)
      integer :: ns, d1, maxe, i, j, k, a, b, nb, q
      logical :: duplicate

      ns = size(tri, 1)
      d1 = size(tri, 2)
      maxe = ns*d1*(d1 - 1)/2
      allocate(buf(max(1, maxe), 2))
      nb = 0
      do i = 1, ns
         do j = 1, d1 - 1
            do k = j + 1, d1
               a = min(tri(i, j), tri(i, k))
               b = max(tri(i, j), tri(i, k))
               duplicate = .false.
               do q = 1, nb
                  if (buf(q, 1) == a .and. buf(q, 2) == b) then
                     duplicate = .true.
                     exit
                  end if
               end do
               if (.not. duplicate) then
                  nb = nb + 1
                  buf(nb, :) = [a, b]
               end if
            end do
         end do
      end do
      allocate(bars(nb, 2))
      if (nb > 0) bars = buf(1:nb, :)
   end subroutine simplex_edges

end module geometry_mesh
