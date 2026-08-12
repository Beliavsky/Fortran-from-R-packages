program test_optimization
  use igraph
  implicit none
  type(graph_t) :: g, t, bg
  type(flow_result_t) :: fr
  type(matching_result_t) :: mr
  integer :: e(2,5), fe(2,5), be(2,5)
  logical :: types(6)
  e = reshape([1,2, 1,3, 2,3, 2,4, 3,4],[2,5])
  g = make_graph(4,e,.false.,[1.0_dp,4.0_dp,2.0_dp,5.0_dp,3.0_dp])
  t = minimum_spanning_tree(g)
  if (t%m /= 3 .or. abs(sum(t%weight)-6.0_dp)>1.0e-12_dp) error stop 'mst'

  fe = reshape([1,2, 1,3, 2,3, 2,4, 3,4],[2,5])
  g = make_graph(4,fe,.true.)
  fr = maxflow(g,1,4,[3.0_dp,2.0_dp,1.0_dp,2.0_dp,4.0_dp])
  if (abs(fr%value-5.0_dp)>1.0e-12_dp) error stop 'maxflow'

  be = reshape([1,4, 1,5, 2,4, 2,6, 3,5],[2,5])
  bg = make_graph(6,be,.false.)
  types=[.false.,.false.,.false.,.true.,.true.,.true.]
  mr=bipartite_matching(bg,types)
  if (mr%size /= 3) error stop 'matching'
  print *, 'test_optimization: PASS'
end program test_optimization
