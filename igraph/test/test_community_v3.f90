program test_community_v3
  use igraph
  implicit none
  type(graph_t)::g
  type(community_result_t)::fg,ld,lv
  integer::e(2,7)
  real(dp)::w(7)

  e(:,1)=[1,2];e(:,2)=[2,3];e(:,3)=[1,3]
  e(:,4)=[4,5];e(:,5)=[5,6];e(:,6)=[4,6];e(:,7)=[3,4]
  w=[1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,0.1_dp]
  g=make_graph(6,e,.false.,w)
  fg=cluster_fast_greedy(g)
  ld=cluster_leiden(g,seed=42_i8)
  lv=cluster_louvain(g,seed=42_i8)
  call check_two(fg,'fast greedy')
  call check_two(ld,'leiden')
  if(abs(fg%modularity-0.48360655737704927_dp)>1.0e-14_dp)error stop 'fast greedy modularity'
  if(abs(ld%modularity-0.48360655737704927_dp)>1.0e-14_dp)error stop 'leiden modularity'
  if(ld%modularity+1.0e-14_dp<lv%modularity)error stop 'leiden below louvain regression'
  if(ld%levels<1 .or. ld%moves<1_i8)error stop 'leiden diagnostics'
  print *, 'test_community_v3: PASS'
contains
  subroutine check_two(r,name)
    type(community_result_t),intent(in)::r
    character(len=*),intent(in)::name
    if(maxval(r%membership)/=2)error stop name//': community count'
    if(any(r%membership(1:3)/=r%membership(1)))error stop name//': first triangle'
    if(any(r%membership(4:6)/=r%membership(4)))error stop name//': second triangle'
    if(r%membership(1)==r%membership(4))error stop name//': triangles not separated'
  end subroutine check_two
end program test_community_v3
