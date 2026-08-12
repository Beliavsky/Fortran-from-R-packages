program v3_algorithms
  use igraph
  implicit none
  type(graph_t) :: g
  type(community_result_t) :: comm
  type(gomory_hu_result_t) :: gh
  type(canonical_result_t) :: canon
  integer :: e(2,7)
  real(dp) :: w(7)

  e(:,1)=[1,2];e(:,2)=[2,3];e(:,3)=[1,3]
  e(:,4)=[4,5];e(:,5)=[5,6];e(:,6)=[4,6];e(:,7)=[3,4]
  w=[1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,0.1_dp]
  g=make_graph(6,e,.false.,w)

  comm=cluster_leiden(g,seed=2026_i8)
  print '(a,*(i0,1x))','Leiden membership: ',comm%membership
  print '(a,f12.8)','modularity: ',comm%modularity

  gh=gomory_hu_tree(g,w)
  print '(a,*(f7.3,1x))','Gomory-Hu cut values: ',gh%cut_value(2:)

  canon=canonical_labeling(ring_graph(5))
  print '(a,*(i0,1x))','canonical C5 order: ',canon%permutation
  print '(a,i0)','C5 automorphisms: ',automorphism_count(ring_graph(5))
end program v3_algorithms
