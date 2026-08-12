module igraph_components
  use igraph_graph, only : graph_t
  implicit none
  private

  type, public :: components_t
    integer :: count = 0
    integer, allocatable :: membership(:)
    integer, allocatable :: csize(:)
  end type components_t

  public :: weak_components, strong_components, articulation_points, bridges

contains

  function weak_components(g) result(c)
    type(graph_t), intent(in) :: g
    type(components_t) :: c
    integer, allocatable :: q(:), sizes(:)
    logical, allocatable :: seen(:)
    integer :: s,head,tail,u,v,p,nc
    allocate(c%membership(g%n),seen(g%n),q(g%n),sizes(max(1,g%n)))
    c%membership=0; seen=.false.; sizes=0; nc=0
    do s=1,g%n
      if (seen(s)) cycle
      nc=nc+1; head=1; tail=1; q(1)=s; seen(s)=.true.; c%membership(s)=nc
      do while(head<=tail)
        u=q(head); head=head+1; sizes(nc)=sizes(nc)+1
        do p=g%out_ptr(u),g%out_ptr(u+1)-1
          v=g%out_nei(p)
          if (.not. seen(v)) then
            seen(v)=.true.; c%membership(v)=nc; tail=tail+1; q(tail)=v
          end if
        end do
        if (g%directed) then
          do p=g%in_ptr(u),g%in_ptr(u+1)-1
            v=g%in_nei(p)
            if (.not. seen(v)) then
              seen(v)=.true.; c%membership(v)=nc; tail=tail+1; q(tail)=v
            end if
          end do
        end if
      end do
    end do
    c%count=nc; allocate(c%csize(nc)); if (nc>0) c%csize=sizes(1:nc)
  end function weak_components

  function strong_components(g) result(c)
    type(graph_t), intent(in) :: g
    type(components_t) :: c
    logical, allocatable :: seen(:)
    integer, allocatable :: order(:), sizes(:)
    integer :: idx,k,v,nc
    if (.not. g%directed) then; c=weak_components(g); return; end if
    allocate(seen(g%n),order(g%n),c%membership(g%n),sizes(max(1,g%n)))
    seen=.false.; idx=0; c%membership=0; sizes=0
    do v=1,g%n
      if (.not. seen(v)) call dfs1(v)
    end do
    seen=.false.; nc=0
    do k=g%n,1,-1
      v=order(k)
      if (.not. seen(v)) then
        nc=nc+1; call dfs2(v,nc)
      end if
    end do
    c%count=nc; allocate(c%csize(nc)); if (nc>0) c%csize=sizes(1:nc)
  contains
    recursive subroutine dfs1(u)
      integer,intent(in)::u
      integer::p,w
      seen(u)=.true.
      do p=g%out_ptr(u),g%out_ptr(u+1)-1
        w=g%out_nei(p); if (.not. seen(w)) call dfs1(w)
      end do
      idx=idx+1; order(idx)=u
    end subroutine dfs1
    recursive subroutine dfs2(u,id)
      integer,intent(in)::u,id
      integer::p,w
      seen(u)=.true.; c%membership(u)=id; sizes(id)=sizes(id)+1
      do p=g%in_ptr(u),g%in_ptr(u+1)-1
        w=g%in_nei(p); if (.not. seen(w)) call dfs2(w,id)
      end do
    end subroutine dfs2
  end function strong_components

  function articulation_points(g) result(points)
    type(graph_t), intent(in) :: g
    integer, allocatable :: points(:)
    integer, allocatable :: disc(:), low(:), parent(:), tmp(:)
    logical, allocatable :: ap(:)
    integer :: time,s,k
    if (g%directed) error stop 'articulation_points: undirected graph required'
    allocate(disc(g%n),low(g%n),parent(g%n),ap(g%n),tmp(g%n))
    disc=0; low=0; parent=0; ap=.false.; time=0
    do s=1,g%n
      if (disc(s)==0) call visit(s)
    end do
    k=0
    do s=1,g%n
      if (ap(s)) then; k=k+1; tmp(k)=s; end if
    end do
    allocate(points(k)); if (k>0) points=tmp(1:k)
  contains
    recursive subroutine visit(u)
      integer,intent(in)::u
      integer::p,v,children
      children=0; time=time+1; disc(u)=time; low(u)=time
      do p=g%out_ptr(u),g%out_ptr(u+1)-1
        v=g%out_nei(p); if (v==u) cycle
        if (disc(v)==0) then
          children=children+1; parent(v)=u; call visit(v); low(u)=min(low(u),low(v))
          if (parent(u)==0 .and. children>1) ap(u)=.true.
          if (parent(u)/=0 .and. low(v)>=disc(u)) ap(u)=.true.
        else if (v/=parent(u)) then
          low(u)=min(low(u),disc(v))
        end if
      end do
    end subroutine visit
  end function articulation_points

  function bridges(g) result(eids)
    type(graph_t), intent(in) :: g
    integer, allocatable :: eids(:)
    integer, allocatable :: disc(:), low(:), parent_edge(:), tmp(:)
    integer :: time,s,k
    if (g%directed) error stop 'bridges: undirected graph required'
    allocate(disc(g%n),low(g%n),parent_edge(g%n),tmp(g%m)); disc=0; low=0; parent_edge=0; time=0; k=0
    do s=1,g%n
      if (disc(s)==0) call visit(s)
    end do
    allocate(eids(k)); if (k>0) eids=tmp(1:k)
  contains
    recursive subroutine visit(u)
      integer,intent(in)::u
      integer::p,v,e
      time=time+1; disc(u)=time; low(u)=time
      do p=g%out_ptr(u),g%out_ptr(u+1)-1
        v=g%out_nei(p); e=g%out_eid(p); if (e==parent_edge(u) .or. v==u) cycle
        if (disc(v)==0) then
          parent_edge(v)=e; call visit(v); low(u)=min(low(u),low(v))
          if (low(v)>disc(u)) then; k=k+1; tmp(k)=e; end if
        else
          low(u)=min(low(u),disc(v))
        end if
      end do
    end subroutine visit
  end function bridges

end module igraph_components
