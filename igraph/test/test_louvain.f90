program test_louvain
  use iso_fortran_env, only : int64
  use igraph
  implicit none
  type(graph_t)::g
  type(community_result_t)::r
  integer,allocatable::e(:,:),mem(:)
  real(dp),allocatable::w(:)
  integer::k,i,j
  real(dp)::q

  allocate(e(2,7),w(7));k=0
  do i=1,3
    do j=i+1,3;k=k+1;e(:,k)=[i,j];w(k)=1.0_dp;end do
  end do
  do i=4,6
    do j=i+1,6;k=k+1;e(:,k)=[i,j];w(k)=1.0_dp;end do
  end do
  k=k+1;e(:,k)=[3,4];w(k)=0.1_dp
  g=make_graph(6,e,.false.,w)
  mem=[1,1,1,2,2,2]
  q=modularity(g,mem)
  call check(q>0.45_dp .and. q<0.50_dp,'weighted modularity')

  r=cluster_louvain(g,seed=12345_int64)
  call check(maxval(r%membership)==2,'Louvain community count')
  call check(all(r%membership(1:3)==r%membership(1)),'first triangle together')
  call check(all(r%membership(4:6)==r%membership(4)),'second triangle together')
  call check(r%membership(1)/=r%membership(4),'triangles separated')
  call check(abs(r%modularity-q)<1.0e-12_dp,'Louvain modularity')
  call check(r%levels>=1,'Louvain levels')

  print *, 'test_louvain: PASS'
contains
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok)then
      print *, 'FAIL: ',trim(msg)
      error stop 1
    end if
  end subroutine check
end program test_louvain
