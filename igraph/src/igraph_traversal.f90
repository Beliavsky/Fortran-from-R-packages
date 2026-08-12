module igraph_traversal
  use igraph_kinds, only : dp, igraph_inf
  use igraph_graph, only : graph_t
  implicit none
  private

  type, public :: bfs_result_t
    integer, allocatable :: order(:)
    integer, allocatable :: parent(:)
    integer, allocatable :: distance(:)
  end type bfs_result_t

  type, public :: dfs_result_t
    integer, allocatable :: order(:)
    integer, allocatable :: parent(:)
    integer, allocatable :: discover(:)
    integer, allocatable :: finish(:)
  end type dfs_result_t

  public :: bfs, dfs, topological_sort, is_dag
  public :: shortest_path, distances_unweighted, distances_dijkstra
  public :: bellman_ford, floyd_warshall

contains

  function bfs(g, root, mode) result(r)
    type(graph_t), intent(in) :: g
    integer, intent(in) :: root
    character(len=*), intent(in), optional :: mode
    type(bfs_result_t) :: r
    integer, allocatable :: q(:)
    logical, allocatable :: seen(:)
    integer :: head, tail, u, v, p, k
    character(len=8) :: md
    if (root<1 .or. root>g%n) error stop 'bfs: invalid root'
    md='out'; if (present(mode)) md=adjustl(mode)
    allocate(r%parent(g%n),r%distance(g%n),seen(g%n),q(g%n),r%order(g%n))
    r%parent=0; r%distance=-1; seen=.false.; head=1; tail=1; q(1)=root
    seen(root)=.true.; r%distance(root)=0; k=0
    do while(head<=tail)
      u=q(head); head=head+1; k=k+1; r%order(k)=u
      select case(trim(md))
      case('in')
        do p=g%in_ptr(u),g%in_ptr(u+1)-1
          v=g%in_nei(p)
          if (.not. seen(v)) then
            seen(v)=.true.; r%parent(v)=u; r%distance(v)=r%distance(u)+1
            tail=tail+1; q(tail)=v
          end if
        end do
      case('all')
        do p=g%out_ptr(u),g%out_ptr(u+1)-1
          v=g%out_nei(p)
          if (.not. seen(v)) then
            seen(v)=.true.; r%parent(v)=u; r%distance(v)=r%distance(u)+1
            tail=tail+1; q(tail)=v
          end if
        end do
        if (g%directed) then
          do p=g%in_ptr(u),g%in_ptr(u+1)-1
            v=g%in_nei(p)
            if (.not. seen(v)) then
              seen(v)=.true.; r%parent(v)=u; r%distance(v)=r%distance(u)+1
              tail=tail+1; q(tail)=v
            end if
          end do
        end if
      case default
        do p=g%out_ptr(u),g%out_ptr(u+1)-1
          v=g%out_nei(p)
          if (.not. seen(v)) then
            seen(v)=.true.; r%parent(v)=u; r%distance(v)=r%distance(u)+1
            tail=tail+1; q(tail)=v
          end if
        end do
      end select
    end do
    if (k<g%n) r%order(k+1:) = 0
  end function bfs

  function dfs(g, root, mode) result(r)
    type(graph_t), intent(in) :: g
    integer, intent(in) :: root
    character(len=*), intent(in), optional :: mode
    type(dfs_result_t) :: r
    logical, allocatable :: seen(:)
    integer :: timer, k
    character(len=8) :: md
    if (root<1 .or. root>g%n) error stop 'dfs: invalid root'
    md='out'; if (present(mode)) md=adjustl(mode)
    allocate(r%order(g%n),r%parent(g%n),r%discover(g%n),r%finish(g%n),seen(g%n))
    r%order=0; r%parent=0; r%discover=0; r%finish=0; seen=.false.; timer=0; k=0
    call visit(root)
  contains
    recursive subroutine visit(u)
      integer, intent(in) :: u
      integer :: p, v
      seen(u)=.true.; timer=timer+1; r%discover(u)=timer; k=k+1; r%order(k)=u
      select case(trim(md))
      case('in')
        do p=g%in_ptr(u),g%in_ptr(u+1)-1
          v=g%in_nei(p)
          if (.not. seen(v)) then; r%parent(v)=u; call visit(v); end if
        end do
      case('all')
        do p=g%out_ptr(u),g%out_ptr(u+1)-1
          v=g%out_nei(p)
          if (.not. seen(v)) then; r%parent(v)=u; call visit(v); end if
        end do
        if (g%directed) then
          do p=g%in_ptr(u),g%in_ptr(u+1)-1
            v=g%in_nei(p)
            if (.not. seen(v)) then; r%parent(v)=u; call visit(v); end if
          end do
        end if
      case default
        do p=g%out_ptr(u),g%out_ptr(u+1)-1
          v=g%out_nei(p)
          if (.not. seen(v)) then; r%parent(v)=u; call visit(v); end if
        end do
      end select
      timer=timer+1; r%finish(u)=timer
    end subroutine visit
  end function dfs

  function topological_sort(g) result(order)
    type(graph_t), intent(in) :: g
    integer, allocatable :: order(:)
    integer, allocatable :: indeg(:), q(:)
    integer :: i,p,u,v,head,tail,k
    if (.not. g%directed) error stop 'topological_sort: graph must be directed'
    allocate(order(g%n),indeg(g%n),q(g%n)); indeg=0
    do i=1,g%n; indeg(i)=g%in_ptr(i+1)-g%in_ptr(i); end do
    head=1; tail=0
    do i=1,g%n
      if (indeg(i)==0) then; tail=tail+1; q(tail)=i; end if
    end do
    k=0
    do while(head<=tail)
      u=q(head); head=head+1; k=k+1; order(k)=u
      do p=g%out_ptr(u),g%out_ptr(u+1)-1
        v=g%out_nei(p); indeg(v)=indeg(v)-1
        if (indeg(v)==0) then; tail=tail+1; q(tail)=v; end if
      end do
    end do
    if (k<g%n) order=0
  end function topological_sort

  logical function is_dag(g)
    type(graph_t), intent(in) :: g
    integer, allocatable :: o(:)
    if (.not. g%directed) then; is_dag=.false.; return; end if
    o=topological_sort(g); is_dag=all(o>0)
  end function is_dag

  function shortest_path(g, source, target, weighted) result(path)
    type(graph_t), intent(in) :: g
    integer, intent(in) :: source, target
    logical, intent(in), optional :: weighted
    integer, allocatable :: path(:)
    integer, allocatable :: parent(:), tmp(:)
    integer :: cur, len
    logical :: wt
    type(bfs_result_t) :: br
    real(dp), allocatable :: d(:)
    wt=.false.; if (present(weighted)) wt=weighted
    if (source<1 .or. source>g%n .or. target<1 .or. target>g%n) error stop 'shortest_path: invalid vertex'
    if (.not. wt) then
      br=bfs(g,source,'out'); parent=br%parent
      if (br%distance(target)<0) then; allocate(path(0)); return; end if
    else
      call dijkstra_one(g,source,d,parent)
      if (d(target)>=igraph_inf/2) then; allocate(path(0)); return; end if
    end if
    allocate(tmp(g%n)); len=0; cur=target
    do
      len=len+1; tmp(len)=cur
      if (cur==source) exit
      cur=parent(cur); if (cur==0) exit
    end do
    allocate(path(len)); path=tmp(len:1:-1)
  end function shortest_path

  function distances_unweighted(g, mode) result(d)
    type(graph_t), intent(in) :: g
    character(len=*), intent(in), optional :: mode
    integer, allocatable :: d(:,:)
    integer :: s
    type(bfs_result_t) :: br
    character(len=8) :: md
    md='out'; if (present(mode)) md=adjustl(mode)
    allocate(d(g%n,g%n))
    do s=1,g%n
      br=bfs(g,s,md); d(s,:)=br%distance
    end do
  end function distances_unweighted

  function distances_dijkstra(g, mode) result(d)
    type(graph_t), intent(in) :: g
    character(len=*), intent(in), optional :: mode
    real(dp), allocatable :: d(:,:)
    real(dp), allocatable :: ds(:)
    integer, allocatable :: parent(:)
    integer :: s
    character(len=8) :: md
    md='out'; if (present(mode)) md=adjustl(mode)
    if (trim(md)/='out') error stop 'distances_dijkstra: v0.1 supports mode=out only'
    allocate(d(g%n,g%n))
    do s=1,g%n
      call dijkstra_one(g,s,ds,parent); d(s,:)=ds
    end do
  end function distances_dijkstra

  subroutine dijkstra_one(g,source,d,parent)
    type(graph_t), intent(in) :: g
    integer, intent(in) :: source
    real(dp), allocatable, intent(out) :: d(:)
    integer, allocatable, intent(out) :: parent(:)
    logical, allocatable :: used(:)
    integer :: i,u,v,p,it
    real(dp) :: best, alt
    if (any(g%weight<0.0_dp)) error stop 'dijkstra: negative edge weight'
    allocate(d(g%n),parent(g%n),used(g%n)); d=igraph_inf; parent=0; used=.false.; d(source)=0.0_dp
    do it=1,g%n
      u=0; best=igraph_inf
      do i=1,g%n
        if (.not. used(i) .and. d(i)<best) then; best=d(i); u=i; end if
      end do
      if (u==0) exit
      used(u)=.true.
      do p=g%out_ptr(u),g%out_ptr(u+1)-1
        v=g%out_nei(p); alt=d(u)+g%out_w(p)
        if (alt<d(v)) then; d(v)=alt; parent(v)=u; end if
      end do
    end do
  end subroutine dijkstra_one

  subroutine bellman_ford(g,source,d,parent,negative_cycle)
    type(graph_t), intent(in) :: g
    integer, intent(in) :: source
    real(dp), allocatable, intent(out) :: d(:)
    integer, allocatable, intent(out) :: parent(:)
    logical, intent(out) :: negative_cycle
    integer :: it,k,u,v
    logical :: changed
    allocate(d(g%n),parent(g%n)); d=igraph_inf; parent=0; d(source)=0.0_dp
    do it=1,g%n-1
      changed=.false.
      do k=1,g%m
        u=g%edge(1,k); v=g%edge(2,k)
        if (d(u)<igraph_inf/2 .and. d(u)+g%weight(k)<d(v)) then
          d(v)=d(u)+g%weight(k); parent(v)=u; changed=.true.
        end if
        if (.not. g%directed) then
          if (d(v)<igraph_inf/2 .and. d(v)+g%weight(k)<d(u)) then
            d(u)=d(v)+g%weight(k); parent(u)=v; changed=.true.
          end if
        end if
      end do
      if (.not. changed) exit
    end do
    negative_cycle=.false.
    do k=1,g%m
      u=g%edge(1,k); v=g%edge(2,k)
      if (d(u)<igraph_inf/2 .and. d(u)+g%weight(k)<d(v)) negative_cycle=.true.
      if (.not. g%directed) then
        if (d(v)<igraph_inf/2 .and. d(v)+g%weight(k)<d(u)) negative_cycle=.true.
      end if
    end do
  end subroutine bellman_ford

  function floyd_warshall(g) result(d)
    type(graph_t), intent(in) :: g
    real(dp), allocatable :: d(:,:)
    integer :: i,j,k,e,u,v
    allocate(d(g%n,g%n)); d=igraph_inf
    do i=1,g%n; d(i,i)=0.0_dp; end do
    do e=1,g%m
      u=g%edge(1,e); v=g%edge(2,e); d(u,v)=min(d(u,v),g%weight(e))
      if (.not. g%directed) d(v,u)=min(d(v,u),g%weight(e))
    end do
    do k=1,g%n
      do i=1,g%n
        if (d(i,k)>=igraph_inf/2) cycle
        do j=1,g%n
          if (d(k,j)>=igraph_inf/2) cycle
          d(i,j)=min(d(i,j),d(i,k)+d(k,j))
        end do
      end do
    end do
  end function floyd_warshall

end module igraph_traversal
