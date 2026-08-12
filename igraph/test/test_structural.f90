program test_structural
  use igraph
  implicit none
  type(graph_t) :: g
  integer :: e(2,5)
  integer, allocatable :: ap(:), br(:), lt(:)
  real(dp), allocatable :: tl(:)
  e = reshape([1,2, 2,3, 3,1, 3,4, 4,5],[2,5])
  g = make_graph(5,e,.false.)
  if (triangle_count(g) /= 1) error stop 'triangles'
  if (abs(transitivity_global(g)-0.5_dp)>1.0e-12_dp) error stop 'transitivity'
  lt=local_triangle_count(g); if (any(lt(1:3)/=1)) error stop 'local triangles'
  tl=transitivity_local(g); if (abs(tl(1)-1.0_dp)>1.0e-12_dp) error stop 'local transitivity'
  ap=articulation_points(g); if (size(ap)/=2 .or. any(ap/=[3,4])) error stop 'articulation'
  br=bridges(g); if (size(br)/=2) error stop 'bridges'
  print *, 'test_structural: PASS'
end program test_structural
