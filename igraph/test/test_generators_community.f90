program test_generators_community
  use igraph
  implicit none
  type(graph_t) :: g
  integer, allocatable :: mem(:)
  real(dp) :: q
  g = ring_graph(6)
  if (g%m /= 6) error stop 'ring'
  g = erdos_renyi_gnm(10,12,seed=12345_i8)
  if (g%m /= 12) error stop 'gnm'
  g = make_graph(6,reshape([1,2,2,3,1,3,4,5,5,6,4,6],[2,6]),.false.)
  mem = label_propagation(g,seed=7_i8)
  if (maxval(mem) /= 2) error stop 'label propagation'
  q = modularity(g,[1,1,1,2,2,2])
  if (q < 0.45_dp) error stop 'modularity'
  print *, 'test_generators_community: PASS'
end program test_generators_community
