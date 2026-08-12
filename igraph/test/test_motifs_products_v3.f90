program test_motifs_products_v3
  use igraph
  implicit none
  type(graph_t)::g,h,a,b
  type(motif_census_result_t)::m1,m2,m3
  real(dp)::cut4(4)
  integer::e1(2,2),e2(2,2)

  g=full_graph(5)
  m1=motif_census(g,4)
  if(m1%occurrences/=5_i8 .or. m1%nclasses/=1)error stop 'K5 K4 motif census'
  cut4=0.0_dp;m2=randesu_motif_census(g,4,cut4,123_i8)
  if(m2%occurrences/=m1%occurrences .or. any(m2%code/=m1%code) .or. any(m2%count/=m1%count))then
    error stop 'RAND-ESU zero-cut exactness'
  end if
  cut4=1.0_dp;m3=randesu_motif_census(g,4,cut4,123_i8)
  if(m3%occurrences/=0_i8)error stop 'RAND-ESU full cut'

  a=ring_graph(3);b=ring_graph(4)
  h=cartesian_product(a,b)
  if(h%n/=12 .or. h%m/=24)error stop 'cartesian product size'
  h=tensor_product(a,b)
  if(h%n/=12 .or. h%m/=24)error stop 'tensor product size'
  h=strong_product(a,b)
  if(h%n/=12 .or. h%m/=48)error stop 'strong product size'
  h=graph_join(a,b)
  if(h%n/=7 .or. h%m/=19)error stop 'graph join size'

  e1(:,1)=[1,2];e1(:,2)=[2,3]
  e2(:,1)=[2,3];e2(:,2)=[3,1]
  a=make_graph(3,e1,.false.);b=make_graph(3,e2,.false.)
  h=graph_union(a,b);if(h%m/=3)error stop 'graph union'
  h=graph_intersection(a,b);if(h%m/=1)error stop 'graph intersection'
  h=graph_difference(a,b);if(h%m/=1)error stop 'graph difference'
  print *, 'test_motifs_products_v3: PASS'
end program test_motifs_products_v3
