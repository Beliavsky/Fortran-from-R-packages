program basic_graph
  use igraph
  implicit none
  type(graph_t) :: g
  type(components_t) :: comp
  real(dp), allocatable :: pr(:), bc(:)
  integer :: e(2,7)
  e = reshape([1,2,2,3,3,1,3,4,4,5,5,6,6,4],[2,7])
  g = make_graph(6,e,.false.)
  comp = weak_components(g)
  bc = betweenness_centrality(g,.true.)
  call pagerank(g,pr)
  print '(a,i0,a,i0)', 'vertices=',g%n,' edges=',g%m
  print '(a,i0)', 'components=',comp%count
  print '(a,*(f8.4,1x))', 'pagerank: ',pr
  print '(a,*(f8.4,1x))', 'betweenness: ',bc
end program basic_graph
