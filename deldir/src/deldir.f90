module deldir
use deldir_kinds, only: dp
use deldir_types, only: delaunay_segment, dirichlet_segment, point_summary, deldir_result, &
    voronoi_component, voronoi_tile, triangle3
use deldir_core, only: deldir_compute
use deldir_geometry, only: deldir_tiles, tile_area, polygon_area, polygon_signed_area, tile_perimeter, &
    tile_centroid, tile_centroids, deldir_triangles, deldir_triangle_matrix, deldir_get_neighbors, &
    mean_nearest_neighbor_distance, which_tile, inside_rect, inside_polygon, centroidal_voronoi
use deldir_utils, only: deldir_bin_sort, duplicated_xy, corner_indices, midpoint_inside, find_new_in_old
use deldir_analysis, only: integer_list, dividing_segment, law_summary_result, tile_info_entry, tile_info_result, &
    deldir_dividing_chain, deldir_law_summary, deldir_tile_info
use polyclip, only: poly_path, poly_set, make_path, make_set, clear_set, append_path
implicit none
private
public :: dp
public :: delaunay_segment, dirichlet_segment, point_summary, deldir_result
public :: voronoi_component, voronoi_tile, triangle3
public :: deldir_compute, deldir_tiles, tile_area, polygon_area, polygon_signed_area, tile_perimeter
public :: tile_centroid, tile_centroids
public :: deldir_triangles, deldir_triangle_matrix, deldir_get_neighbors, mean_nearest_neighbor_distance
public :: which_tile, inside_rect, inside_polygon, centroidal_voronoi
public :: integer_list, dividing_segment, law_summary_result, tile_info_entry, tile_info_result
public :: deldir_dividing_chain, deldir_law_summary, deldir_tile_info
public :: deldir_bin_sort, duplicated_xy, corner_indices, midpoint_inside, find_new_in_old
! Re-export the polygon container used by clipping APIs so callers need only `use deldir`.
public :: poly_path, poly_set, make_path, make_set, clear_set, append_path
end module deldir
