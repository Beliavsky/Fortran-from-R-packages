program test_motifs
  use igraph
  implicit none
  integer,parameter::rep(16)=[0,1,3,18,33,17,35,19,37,38,51,30,45,29,61,63]
  type(graph_t)::g
  integer(i8)::tc(16),dc(3)
  integer::i
  do i=1,16
    g=from_mask(rep(i))
    tc=triad_census(g)
    call check(tc(i)==1_i8 .and. sum(tc)==1_i8,'directed triad class')
  end do
  g=full_graph(3)
  tc=triad_census(g)
  call check(tc(16)==1_i8,'undirected K3 is 300')
  g=from_mask(3)
  dc=dyad_census(g)
  call check(all(dc==[1_i8,0_i8,2_i8]),'dyad census')
  print *, 'test_motifs: PASS'
contains
  function from_mask(mask) result(h)
    integer,intent(in)::mask
    type(graph_t)::h
    integer,allocatable::e(:,:)
    integer::b,k,u,v
    allocate(e(2,popcnt(mask)));k=0
    do b=0,5
      if(.not.btest(mask,b))cycle
      call pair(b,u,v);k=k+1;e(:,k)=[u,v]
    end do
    h=make_graph(3,e,.true.)
  end function from_mask
  subroutine pair(bit,u,v)
    integer,intent(in)::bit
    integer,intent(out)::u,v
    select case(bit)
    case(0);u=1;v=2
    case(1);u=2;v=1
    case(2);u=1;v=3
    case(3);u=3;v=1
    case(4);u=2;v=3
    case default;u=3;v=2
    end select
  end subroutine pair
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok)then;print *,'FAIL: ',trim(msg);error stop 1;end if
  end subroutine check
end program test_motifs
