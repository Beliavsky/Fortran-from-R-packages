program test_vf2_cliques
  use iso_fortran_env, only : int64
  use igraph
  implicit none
  type(graph_t) :: g1,g2,tri,k4,gc
  type(isomorphism_result_t) :: iso,sub
  type(vertex_set_list_t) :: mc,lc,allc
  integer, allocatable :: e(:,:), ep(:,:)
  integer :: i,j,k

  g1=ring_graph(6)
  allocate(e(2,6))
  e(:,1)=[3,5];e(:,2)=[5,2];e(:,3)=[2,6];e(:,4)=[6,1];e(:,5)=[1,4];e(:,6)=[4,3]
  g2=make_graph(6,e,.false.)
  iso=vf2_isomorphic(g1,g2)
  call check(iso%isomorphic,'ring isomorphism')
  call check(all(iso%map12>=1),'isomorphism mapping')

  allocate(ep(2,6));k=0
  do i=1,3
    do j=i+1,3;k=k+1;ep(:,k)=[i,j];end do
  end do
  do i=4,6
    do j=i+1,6;k=k+1;ep(:,k)=[i,j];end do
  end do
  gc=make_graph(6,ep,.false.)
  iso=vf2_isomorphic(g1,gc)
  call check(.not.iso%isomorphic,'non-isomorphic graphs')

  tri=full_graph(3)
  k4=full_graph(4)
  sub=vf2_subisomorphic(tri,k4)
  call check(sub%isomorphic,'triangle subisomorphism')

  deallocate(e);allocate(e(2,7));k=0
  do i=1,4
    do j=i+1,4;k=k+1;e(:,k)=[i,j];end do
  end do
  k=k+1;e(:,k)=[4,5]
  gc=make_graph(5,e,.false.)
  mc=maximal_cliques(gc)
  call check(mc%count==2,'maximal clique count')
  call check(clique_number(gc)==4,'clique number')
  lc=largest_cliques(gc)
  call check(lc%count==1 .and. size(lc%set(1)%v)==4,'largest clique')
  allc=cliques(gc,3,4)
  call check(allc%count==5,'all clique count sizes 3-4')

  print *, 'test_vf2_cliques: PASS'
contains
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok)then
      print *, 'FAIL: ',trim(msg)
      error stop 1
    end if
  end subroutine check
end program test_vf2_cliques
