program advanced_algorithms
  use iso_fortran_env, only : int64
  use igraph
  implicit none
  type(graph_t) :: g
  type(community_result_t) :: comm
  type(vertex_set_list_t) :: mc
  type(isomorphism_result_t) :: iso
  integer, allocatable :: e(:,:)
  real(dp), allocatable :: w(:)
  integer :: i,j,k

  allocate(e(2,7),w(7));k=0
  do i=1,3
    do j=i+1,3;k=k+1;e(:,k)=[i,j];w(k)=1.0_dp;end do
  end do
  do i=4,6
    do j=i+1,6;k=k+1;e(:,k)=[i,j];w(k)=1.0_dp;end do
  end do
  k=k+1;e(:,k)=[3,4];w(k)=0.1_dp
  g=make_graph(6,e,.false.,w)

  comm=cluster_louvain(g,seed=2026_int64)
  print *, 'Louvain membership:',comm%membership
  print *, 'modularity:',comm%modularity

  mc=maximal_cliques(g)
  print *, 'maximal cliques:',mc%count

  iso=vf2_isomorphic(ring_graph(6),ring_graph(6))
  print *, 'ring isomorphic:',iso%isomorphic
end program advanced_algorithms
