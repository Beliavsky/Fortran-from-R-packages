program test_graph_core
  use igraph
  implicit none
  type(graph_t) :: g, h
  integer :: e(2,5)
  integer, allocatable :: d(:)
  real(dp), allocatable :: s(:)
  e = reshape([1,2, 2,3, 3,1, 3,4, 4,4],[2,5])
  g = make_graph(4,e,.false.)
  if (vertex_count(g) /= 4 .or. edge_count(g) /= 5) error stop 'size'
  if (.not. are_adjacent(g,1,2)) error stop 'adjacency'
  d = degree(g,'all',.true.)
  if (any(d /= [2,2,3,3])) error stop 'degree'
  s = strength(g)
  if (maxval(abs(s-real(d,dp))) > 1.0e-12_dp) error stop 'strength'
  h = simplify_graph(g,.true.,.true.)
  if (h%m /= 4) error stop 'simplify'
  print *, 'test_graph_core: PASS'
end program test_graph_core
