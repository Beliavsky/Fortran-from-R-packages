program test_centrality
  use igraph
  implicit none
  type(graph_t) :: g
  real(dp), allocatable :: c(:), b(:), pr(:), ev(:)
  real(dp) :: val
  integer :: e(2,4), it
  e = reshape([1,2, 2,3, 3,4, 4,5],[2,4])
  g = make_graph(5,e,.false.)
  c = closeness_centrality(g)
  if (maxloc(c,dim=1) /= 3) error stop 'closeness'
  b = betweenness_centrality(g)
  if (maxloc(b,dim=1) /= 3) error stop 'betweenness'
  call pagerank(g,pr,iterations=it)
  if (abs(sum(pr)-1.0_dp)>1.0e-12_dp) error stop 'pagerank'
  call eigenvector_centrality(g,ev,val)
  if (maxloc(ev,dim=1) /= 3) error stop 'eigenvector'
  print *, 'test_centrality: PASS'
end program test_centrality
