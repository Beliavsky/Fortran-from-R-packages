program flow_matching
  use igraph
  implicit none
  type(graph_t) :: g
  type(flow_result_t) :: f
  integer :: e(2,5)
  e=reshape([1,2,1,3,2,3,2,4,3,4],[2,5])
  g=make_graph(4,e,.true.)
  f=maxflow(g,1,4,[3.0_dp,2.0_dp,1.0_dp,2.0_dp,4.0_dp])
  print '(a,f8.3)', 'maximum flow = ',f%value
end program flow_matching
