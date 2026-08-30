module polyclip
   use polyclip_kinds, only: dp
   use polyclip_types, only: poly_path, poly_set, make_path, make_set, clear_set, append_path, &
      clip_intersection, clip_union, clip_difference, clip_xor, &
      fill_evenodd, fill_nonzero, fill_positive, fill_negative, &
      join_square, join_round, join_miter, &
      end_closed_polygon, end_closed_line, end_open_butt, end_open_square, end_open_round
   use polyclip_api, only: polyclip_apply, polysimplify, polyoffset, polylineoffset, polyminkowski, pointinpolygon
   implicit none
   private
   public :: dp
   public :: poly_path, poly_set, make_path, make_set, clear_set, append_path
   public :: clip_intersection, clip_union, clip_difference, clip_xor
   public :: fill_evenodd, fill_nonzero, fill_positive, fill_negative
   public :: join_square, join_round, join_miter
   public :: end_closed_polygon, end_closed_line, end_open_butt, end_open_square, end_open_round
   public :: polyclip_apply, polysimplify, polyoffset, polylineoffset, polyminkowski, pointinpolygon
end module polyclip
