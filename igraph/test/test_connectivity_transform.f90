program test_connectivity_transform
  use igraph
  implicit none
  type(graph_t)::cyc,path,k4,disc,h,c,t,l,dg,tr
  integer,allocatable::e(:,:),v(:),ids(:)

  cyc=ring_graph(5)
  call check(edge_connectivity(cyc)==2,'cycle edge connectivity')
  call check(vertex_connectivity(cyc)==2,'cycle vertex connectivity')
  call check(is_biconnected(cyc),'cycle biconnected')

  allocate(e(2,4));e(:,1)=[1,2];e(:,2)=[2,3];e(:,3)=[3,4];e(:,4)=[4,5]
  path=make_graph(5,e,.false.)
  call check(edge_connectivity(path)==1,'path edge connectivity')
  call check(vertex_connectivity(path)==1,'path vertex connectivity')
  call check(.not.is_biconnected(path),'path not biconnected')

  k4=full_graph(4)
  call check(edge_connectivity(k4)==3,'K4 edge connectivity')
  call check(vertex_connectivity(k4)==3,'K4 vertex connectivity')

  deallocate(e);allocate(e(2,2));e(:,1)=[1,2];e(:,2)=[3,4]
  disc=make_graph(4,e,.false.)
  call check(edge_connectivity(disc)==0,'disconnected edge connectivity')
  call check(vertex_connectivity(disc)==0,'disconnected vertex connectivity')

  v=[1,2,3]
  h=induced_subgraph(k4,v)
  call check(h%n==3 .and. h%m==3,'induced subgraph')
  c=complement_graph(path)
  call check(c%n==5 .and. c%m==6,'complement edge count')
  l=line_graph(path)
  call check(l%n==4 .and. l%m==3,'line graph path')

  deallocate(e);allocate(e(2,3));e(:,1)=[1,2];e(:,2)=[2,3];e(:,3)=[1,3]
  dg=make_graph(3,e,.true.)
  tr=transpose_graph(dg)
  call check(are_adjacent(tr,2,1) .and. are_adjacent(tr,3,2),'transpose graph')

  ids=[1,3]
  t=edge_subgraph(dg,ids,.false.)
  call check(t%m==2 .and. t%n==3,'edge subgraph')
  h=disjoint_union(path,cyc)
  call check(h%n==10 .and. h%m==9,'disjoint union')

  print *, 'test_connectivity_transform: PASS'
contains
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok)then
      print *, 'FAIL: ',trim(msg)
      error stop 1
    end if
  end subroutine check
end program test_connectivity_transform
