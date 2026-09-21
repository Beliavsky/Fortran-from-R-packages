program geometry_example
   use geometry
   implicit none

   real(dp) :: cube(8, 3)
   type(hull_result) :: hull
   type(delaunay_result) :: dt

   cube = reshape([-0.5_dp, -0.5_dp, -0.5_dp, &
                    0.5_dp, -0.5_dp, -0.5_dp, &
                    0.5_dp,  0.5_dp, -0.5_dp, &
                   -0.5_dp,  0.5_dp, -0.5_dp, &
                   -0.5_dp, -0.5_dp,  0.5_dp, &
                    0.5_dp, -0.5_dp,  0.5_dp, &
                    0.5_dp,  0.5_dp,  0.5_dp, &
                   -0.5_dp,  0.5_dp,  0.5_dp], [8, 3], order=[2, 1])

   call convhulln_points(cube, hull)
   call delaunayn_points(cube, dt)

   write(*, '(a,i0)') 'triangulated hull facets: ', size(hull%facets, 1)
   write(*, '(a,f8.4)') 'surface area: ', hull%area
   write(*, '(a,f8.4)') 'volume: ', hull%volume
   write(*, '(a,i0)') 'Delaunay tetrahedra: ', size(dt%tri, 1)
   write(*, '(a,f8.4)') 'sum tetrahedron volumes: ', sum(dt%areas)
end program geometry_example
