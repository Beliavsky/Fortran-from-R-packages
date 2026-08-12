module igraph_connectivity
  use igraph_kinds, only : dp
  use igraph_graph, only : graph_t, make_graph, are_adjacent, degree
  use igraph_components, only : components_t, weak_components, strong_components, articulation_points
  use igraph_optimization, only : flow_result_t, maxflow
  implicit none
  private

  type, public :: cut_result_t
    real(dp) :: value = 0.0_dp
    logical, allocatable :: partition(:)
  end type cut_result_t

  public :: minimum_cut, edge_connectivity, vertex_connectivity, is_biconnected

contains

  function minimum_cut(g, source, sink, capacity) result(cut)
    type(graph_t), intent(in) :: g
    integer, intent(in) :: source, sink
    real(dp), intent(in), optional :: capacity(:)
    type(cut_result_t) :: cut
    type(flow_result_t) :: f
    if(present(capacity))then
      f=maxflow(g,source,sink,capacity)
    else
      f=maxflow(g,source,sink)
    end if
    cut%value=f%value
    allocate(cut%partition(g%n));cut%partition=f%source_side
  end function minimum_cut

  integer function edge_connectivity(g) result(lambda)
    type(graph_t), intent(in) :: g
    type(flow_result_t) :: f
    real(dp), allocatable :: cap(:)
    integer :: s,t,best
    if(g%n<=1)then;lambda=0;return;end if
    allocate(cap(g%m),source=1.0_dp)
    best=huge(1)
    if(g%directed)then
      do s=1,g%n
        do t=1,g%n
          if(s==t)cycle
          f=maxflow(g,s,t,cap)
          best=min(best,nint(f%value))
          if(best==0)then;lambda=0;return;end if
        end do
      end do
    else
      s=1
      do t=2,g%n
        f=maxflow(g,s,t,cap)
        best=min(best,nint(f%value))
        if(best==0)then;lambda=0;return;end if
      end do
    end if
    if(best==huge(1))best=0
    lambda=best
  end function edge_connectivity

  integer function vertex_connectivity(g) result(kappa)
    type(graph_t), intent(in) :: g
    type(components_t) :: comp
    integer, allocatable :: d(:)
    integer :: s,t,best
    logical :: has_pair
    if(g%n<=1)then;kappa=0;return;end if

    if(g%directed)then
      comp=strong_components(g)
    else
      comp=weak_components(g)
    end if
    if(comp%count>1)then;kappa=0;return;end if

    d=degree(g,'all',.false.)
    if(.not.g%directed)then
      if(all(d>=g%n-1))then;kappa=g%n-1;return;end if
    end if

    best=g%n-1;has_pair=.false.
    do s=1,g%n
      do t=1,g%n
        if(s==t)cycle
        if(.not.g%directed .and. t<s)cycle
        if(are_adjacent(g,s,t))then
          if(.not.g%directed)cycle
        end if
        has_pair=.true.
        best=min(best,local_vertex_connectivity(g,s,t))
        if(best==0)then;kappa=0;return;end if
      end do
    end do
    if(.not.has_pair)best=g%n-1
    kappa=best
  end function vertex_connectivity

  logical function is_biconnected(g) result(ok)
    type(graph_t), intent(in) :: g
    type(components_t) :: comp
    integer, allocatable :: aps(:)
    if(g%directed)error stop 'is_biconnected: undirected graph required'
    if(g%n<2)then;ok=.false.;return;end if
    comp=weak_components(g)
    if(comp%count/=1)then;ok=.false.;return;end if
    aps=articulation_points(g)
    ok=size(aps)==0
  end function is_biconnected

  integer function local_vertex_connectivity(g,s,t) result(val)
    type(graph_t), intent(in) :: g
    integer, intent(in) :: s,t
    type(graph_t) :: h
    type(flow_result_t) :: f
    integer, allocatable :: e(:,:)
    real(dp), allocatable :: cap(:)
    integer :: maxe,k,m,u,v,vin,vout
    real(dp) :: infcap

    infcap=real(g%n+1,dp)
    maxe=g%n+merge(g%m,2*g%m,g%directed)
    allocate(e(2,maxe),cap(maxe));m=0
    do u=1,g%n
      vin=2*u-1;vout=2*u;m=m+1;e(:,m)=[vin,vout]
      if(u==s .or. u==t)then;cap(m)=infcap;else;cap(m)=1.0_dp;end if
    end do
    do k=1,g%m
      u=g%edge(1,k);v=g%edge(2,k)
      m=m+1;e(:,m)=[2*u,2*v-1];cap(m)=infcap
      if(.not.g%directed)then
        m=m+1;e(:,m)=[2*v,2*u-1];cap(m)=infcap
      end if
    end do
    h=make_graph(2*g%n,e(:,:m),.true.,cap(:m))
    f=maxflow(h,2*s,2*t-1,cap(:m))
    val=min(g%n-1,nint(f%value))
  end function local_vertex_connectivity

end module igraph_connectivity
