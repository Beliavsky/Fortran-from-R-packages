module geometry
   use geometry_kinds, only : dp
   use geometry_types, only : hull_result, delaunay_result, search_result, intersection_result
   use geometry_basic, only : cart2pol_points, pol2cart_points, cart2sph_points, sph2cart_points
   use geometry_basic, only : dot_columns, dot_rows, extprod3d_rows, polyarea_vector
   use geometry_basic, only : matmax_rows, matmin_rows, matsort_rows, matorder_rows, unique_rows
   use geometry_basic, only : entry_value_real, entry_set_real, rbox_points
   use geometry_basic, only : mesh_dcircle, mesh_drectangle, mesh_dsphere, mesh_hunif
   use geometry_basic, only : mesh_diff_values, mesh_union_values, mesh_intersect_values
   use geometry_hull, only : convhulln_points, halfspacen_vertices, inhulln_points
   use geometry_hull, only : feasible_point_hulls, intersectn_points
   use geometry_delaunay, only : delaunayn_points, cart2bary_points, bary2cart_points
   use geometry_delaunay, only : tsearch_points, tsearchn_points
   use geometry_mesh, only : surf_triangles, distmesh2d_circle, distmeshnd_sphere
   implicit none
   private

   public :: dp, hull_result, delaunay_result, search_result, intersection_result
   public :: cart2pol_points, pol2cart_points, cart2sph_points, sph2cart_points
   public :: dot_columns, dot_rows, extprod3d_rows, polyarea_vector
   public :: matmax_rows, matmin_rows, matsort_rows, matorder_rows, unique_rows
   public :: entry_value_real, entry_set_real, rbox_points
   public :: mesh_dcircle, mesh_drectangle, mesh_dsphere, mesh_hunif
   public :: mesh_diff_values, mesh_union_values, mesh_intersect_values
   public :: convhulln_points, halfspacen_vertices, inhulln_points
   public :: feasible_point_hulls, intersectn_points
   public :: delaunayn_points, cart2bary_points, bary2cart_points
   public :: tsearch_points, tsearchn_points
   public :: surf_triangles, distmesh2d_circle, distmeshnd_sphere
end module geometry
