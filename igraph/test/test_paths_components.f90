program test_paths_components
  use igraph
  implicit none
  type(graph_t) :: g, dg
  type(bfs_result_t) :: br
  type(components_t) :: c
  integer :: e(2,5), de(2,5)
  integer, allocatable :: p(:), du(:,:), topo(:)
  real(dp), allocatable :: dd(:,:), bf(:)
  integer, allocatable :: par(:)
  logical :: neg, dag_ok
  e = reshape([1,2, 2,3, 3,4, 4,5, 2,5],[2,5])
  g = make_graph(6,e,.false.)
  br = bfs(g,1)
  if (br%distance(5) /= 2 .or. br%distance(6) /= -1) error stop 'bfs'
  p = shortest_path(g,1,5)
  if (size(p) /= 3 .or. p(1) /= 1 .or. p(3) /= 5) error stop 'path'
  c = weak_components(g)
  if (c%count /= 2) error stop 'weak components'
  du = distances_unweighted(g,'all')
  if (du(1,5) /= 2) error stop 'distances'

  de = reshape([1,2, 1,3, 2,4, 3,4, 4,5],[2,5])
  dg = make_graph(5,de,.true.,[1.0_dp,4.0_dp,2.0_dp,1.0_dp,3.0_dp])
  dd = distances_dijkstra(dg)
  if (abs(dd(1,5)-6.0_dp) > 1.0e-12_dp) error stop 'dijkstra'
  call bellman_ford(dg,1,bf,par,neg)
  if (neg .or. abs(bf(5)-6.0_dp)>1.0e-12_dp) error stop 'bellman ford'
  topo = topological_sort(dg)
  dag_ok = is_dag(dg)
  if (any(topo==0)) error stop 'dag order'
  if (.not. dag_ok) error stop 'dag flag'
  print *, 'test_paths_components: PASS'
end program test_paths_components
