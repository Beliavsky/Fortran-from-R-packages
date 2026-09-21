program test_geometry
   use geometry
   implicit none

   real(dp) :: square(4, 2), cube(8, 3), shifted(8, 3), fp3(3)
   real(dp) :: simplex(3, 2), q(2, 2)
   real(dp), allocatable :: vertices(:, :), beta(:, :), back(:, :), mesh(:, :)
   integer :: tri2(2, 3), tetra(1, 4)
   integer, allocatable :: surface(:, :)
   logical, allocatable :: inside(:)
   logical :: ok
   type(hull_result) :: hsquare, hcube
   type(delaunay_result) :: dt2, dt3
   type(search_result) :: sr
   type(intersection_result) :: inter

   square = reshape([0.0_dp, 0.0_dp, &
                     1.0_dp, 0.0_dp, &
                     1.0_dp, 1.0_dp, &
                     0.0_dp, 1.0_dp], [4, 2], order=[2, 1])
   call convhulln_points(square, hsquare)
   call assert_true(hsquare%success, "2D convex hull success")
   call assert_true(size(hsquare%facets, 1) == 4, "2D convex hull facet count")
   call assert_close(hsquare%area, 4.0_dp, 1.0e-12_dp, "2D hull perimeter")
   call assert_close(hsquare%volume, 1.0_dp, 1.0e-12_dp, "2D hull area")

   call delaunayn_points(square, dt2)
   call assert_true(dt2%success, "2D Delaunay success")
   call assert_true(size(dt2%tri, 1) == 2, "2D Delaunay triangle count")
   call assert_close(sum(dt2%areas), 1.0_dp, 1.0e-12_dp, "2D Delaunay area")

   simplex = square([1, 2, 4], :)
   q = reshape([0.25_dp, 0.25_dp, 0.75_dp, 0.10_dp], [2, 2], order=[2, 1])
   call cart2bary_points(simplex, q, beta, ok)
   call assert_true(ok, "cart2bary success")
   call bary2cart_points(simplex, beta, back, ok)
   call assert_true(ok, "bary2cart success")
   call assert_close(maxval(abs(back - q)), 0.0_dp, 1.0e-12_dp, "barycentric round trip")

   tri2 = reshape([1, 2, 3, 1, 3, 4], [2, 3], order=[2, 1])
   call tsearchn_points(square, tri2, q, sr)
   call assert_true(all(sr%idx > 0), "tsearchn finds interior points")
   call tsearch_points(square(:, 1), square(:, 2), tri2, q(:, 1), q(:, 2), sr)
   call assert_true(all(sr%idx > 0), "tsearch finds interior points")

   cube = reshape([-0.5_dp, -0.5_dp, -0.5_dp, &
                    0.5_dp, -0.5_dp, -0.5_dp, &
                    0.5_dp,  0.5_dp, -0.5_dp, &
                   -0.5_dp,  0.5_dp, -0.5_dp, &
                   -0.5_dp, -0.5_dp,  0.5_dp, &
                    0.5_dp, -0.5_dp,  0.5_dp, &
                    0.5_dp,  0.5_dp,  0.5_dp, &
                   -0.5_dp,  0.5_dp,  0.5_dp], [8, 3], order=[2, 1])
   call convhulln_points(cube, hcube)
   call assert_true(hcube%success, "3D convex hull success")
   call assert_true(size(hcube%facets, 1) == 12, "cube triangulated facet count")
   call assert_close(hcube%area, 6.0_dp, 1.0e-10_dp, "cube surface area")
   call assert_close(hcube%volume, 1.0_dp, 1.0e-10_dp, "cube volume")

   fp3 = 0.0_dp
   call halfspacen_vertices(hcube%normals, fp3, vertices, ok)
   call assert_true(ok, "halfspace intersection success")
   call assert_true(size(vertices, 1) == 8, "halfspace cube vertex count")
   q = reshape([0.5_dp, 0.5_dp, 1.2_dp, 0.5_dp], [2, 2], order=[2, 1])
   call inhulln_points(hsquare, q, inside)
   call assert_true(inside(1), "inhulln interior point")
   call assert_true(.not. inside(2), "inhulln exterior point")

   shifted = cube + 0.5_dp
   call intersectn_points(cube, shifted, 1.0e-10_dp, inter)
   call assert_true(inter%success, "intersectn overlapping cubes")
   call assert_close(inter%ch%volume, 0.125_dp, 1.0e-9_dp, "intersection cube volume")

   call delaunayn_points(cube, dt3)
   call assert_true(dt3%success, "3D Delaunay success")
   call assert_close(sum(dt3%areas), 1.0_dp, 1.0e-9_dp, "3D Delaunay volume")

   tetra(1, :) = [1, 2, 3, 5]
   call surf_triangles(cube, tetra, surface)
   call assert_true(size(surface, 1) == 4, "surf.tri tetrahedron boundary")

   call distmesh2d_circle(0.5_dp, reshape([-1.0_dp, -1.0_dp, 1.0_dp, 1.0_dp], [2, 2], order=[2, 1]), &
                          1.0_dp, 1, mesh)
   call assert_true(size(mesh, 1) > 3, "distmesh2d circle creates nodes")
   call assert_true(maxval(sqrt(sum(mesh*mesh, dim=2))) <= 1.0_dp + 1.0e-10_dp, &
                    "distmesh2d nodes remain inside circle")

   print '(a)', 'All geometry computational tests passed'

contains

   subroutine assert_close(actual, expected, tolerance, message)
      real(dp), intent(in) :: actual !! Computed scalar value under test.
      real(dp), intent(in) :: expected !! Expected scalar reference value.
      real(dp), intent(in) :: tolerance !! Maximum allowed absolute error.
      character(*), intent(in) :: message !! Test label printed when the assertion fails.
      if (abs(actual - expected) > tolerance) then
         write(*, '(a,2(1x,es24.16))') trim(message), actual, expected
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition !! Boolean condition that must evaluate true.
      character(*), intent(in) :: message !! Test label printed when the assertion fails.
      if (.not. condition) then
         write(*, '(a)') trim(message)
         error stop 1
      end if
   end subroutine assert_true

end program test_geometry
