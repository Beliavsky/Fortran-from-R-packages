module igraph_flow_ext
  use igraph_kinds, only : dp, igraph_inf
  use igraph_graph, only : graph_t, make_graph
  use igraph_optimization, only : flow_result_t, maxflow
  implicit none
  private

  type, public :: gomory_hu_result_t
    type(graph_t) :: tree
    integer, allocatable :: parent(:)
    real(dp), allocatable :: cut_value(:)
  end type gomory_hu_result_t

  type, public :: mincost_result_t
    real(dp) :: flow = 0.0_dp
    real(dp) :: cost = 0.0_dp
    real(dp), allocatable :: edge_flow(:)
    logical :: feasible = .true.
    integer :: augmentations = 0
  end type mincost_result_t

  public :: gomory_hu_tree, min_cost_flow

contains

  function gomory_hu_tree(g, capacity) result(res)
    type(graph_t), intent(in) :: g
    real(dp), intent(in), optional :: capacity(:)
    type(gomory_hu_result_t) :: res
    type(flow_result_t) :: f
    integer, allocatable :: e(:,:)
    real(dp), allocatable :: w(:), cap(:)
    integer :: s,t,i

    if (g%directed) error stop 'gomory_hu_tree: undirected graph required'
    if (present(capacity)) then
      if (size(capacity) /= g%m) error stop 'gomory_hu_tree: capacity length mismatch'
      if (any(capacity < 0.0_dp)) error stop 'gomory_hu_tree: negative capacity'
      allocate(cap(g%m));cap=capacity
    else
      allocate(cap(g%m));cap=g%weight
      if (any(cap < 0.0_dp)) error stop 'gomory_hu_tree: negative graph weight'
    end if
    allocate(res%parent(g%n),res%cut_value(g%n))
    res%parent=0;res%cut_value=0.0_dp
    if (g%n == 0) then
      allocate(e(2,0));res%tree=make_graph(0,e,.false.);return
    end if
    res%parent(1)=0
    if (g%n >= 2) res%parent(2:g%n)=1
    do s=2,g%n
      t=res%parent(s)
      f=maxflow(g,s,t,cap)
      res%cut_value(s)=f%value
      do i=s+1,g%n
        if(res%parent(i)==t .and. f%source_side(i))res%parent(i)=s
      end do
      if(t/=1)then
        if(f%source_side(res%parent(t)))then
          res%parent(s)=res%parent(t)
          res%parent(t)=s
          res%cut_value(s)=res%cut_value(t)
          res%cut_value(t)=f%value
        end if
      end if
    end do
    allocate(e(2,max(0,g%n-1)),w(max(0,g%n-1)))
    do s=2,g%n
      e(:,s-1)=[s,res%parent(s)]
      w(s-1)=res%cut_value(s)
    end do
    res%tree=make_graph(g%n,e,.false.,w)
  end function gomory_hu_tree

  function min_cost_flow(g,source,sink,demand,capacity,cost) result(res)
    type(graph_t),intent(in)::g
    integer,intent(in)::source,sink
    real(dp),intent(in),optional::demand
    real(dp),intent(in),optional::capacity(:),cost(:)
    type(mincost_result_t)::res
    integer::maxarc,ei,k,u,v,e,iter
    integer,allocatable::from(:),to(:),rev(:),orig(:),prev(:)
    real(dp),allocatable::cap(:),ecost(:),dist(:),cap0(:)
    real(dp)::c,cc,target,remain,aug,tol
    logical::changed,has_demand

    if(.not.g%directed)error stop 'min_cost_flow: directed graph required'
    if(source<1 .or. source>g%n .or. sink<1 .or. sink>g%n .or. source==sink)then
      error stop 'min_cost_flow: invalid source/sink'
    end if
    if(present(capacity))then
      if(size(capacity)/=g%m)error stop 'min_cost_flow: capacity length mismatch'
      if(any(capacity<0.0_dp))error stop 'min_cost_flow: negative capacity'
    end if
    if(present(cost))then
      if(size(cost)/=g%m)error stop 'min_cost_flow: cost length mismatch'
    end if
    has_demand=present(demand)
    if(has_demand)then
      if(demand<0.0_dp)error stop 'min_cost_flow: negative demand'
      target=demand
    else
      target=igraph_inf
    end if
    tol=1.0e-12_dp
    maxarc=max(2,2*g%m)
    allocate(from(maxarc),to(maxarc),rev(maxarc),orig(maxarc),cap(maxarc),ecost(maxarc),cap0(maxarc))
    ei=0
    do k=1,g%m
      u=g%edge(1,k);v=g%edge(2,k)
      if(present(capacity))then;c=capacity(k);else;c=1.0_dp;end if
      if(present(cost))then;cc=cost(k);else;cc=g%weight(k);end if
      call add_arc(u,v,c,cc,k)
    end do
    allocate(dist(g%n),prev(g%n),res%edge_flow(g%m));res%edge_flow=0.0_dp
    res%flow=0.0_dp;res%cost=0.0_dp;res%augmentations=0
    do
      if(has_demand)then
        remain=target-res%flow
        if(remain<=tol)exit
      else
        remain=igraph_inf
      end if
      dist=igraph_inf;prev=0;dist(source)=0.0_dp
      do iter=1,g%n-1
        changed=.false.
        do e=1,ei
          if(cap(e)<=tol)cycle
          u=from(e);v=to(e)
          if(dist(u)>=igraph_inf/2.0_dp)cycle
          if(dist(u)+ecost(e)<dist(v)-1.0e-14_dp)then
            dist(v)=dist(u)+ecost(e);prev(v)=e;changed=.true.
          end if
        end do
        if(.not.changed)exit
      end do
      if(prev(sink)==0)exit
      aug=remain;v=sink
      do while(v/=source)
        e=prev(v);if(e==0)exit
        aug=min(aug,cap(e));v=from(e)
      end do
      if(aug<=tol .or. aug>=igraph_inf/2.0_dp)then
        if(.not.has_demand .and. aug>=igraph_inf/2.0_dp)error stop 'min_cost_flow: unbounded capacity path'
        exit
      end if
      v=sink
      do while(v/=source)
        e=prev(v);cap(e)=cap(e)-aug;cap(rev(e))=cap(rev(e))+aug
        res%cost=res%cost+aug*ecost(e);v=from(e)
      end do
      res%flow=res%flow+aug;res%augmentations=res%augmentations+1
    end do
    do e=1,ei,2
      if(orig(e)>0)res%edge_flow(orig(e))=cap0(e)-cap(e)
    end do
    if(has_demand)then
      res%feasible=res%flow>=target-tol
    else
      res%feasible=.true.
    end if

  contains
    subroutine add_arc(a,b,costcap,unitcost,oeid)
      integer,intent(in)::a,b,oeid
      real(dp),intent(in)::costcap,unitcost
      integer::x,y
      x=ei+1;y=ei+2;ei=ei+2
      from(x)=a;to(x)=b;cap(x)=costcap;cap0(x)=costcap;ecost(x)=unitcost;rev(x)=y;orig(x)=oeid
      from(y)=b;to(y)=a;cap(y)=0.0_dp;cap0(y)=0.0_dp;ecost(y)=-unitcost;rev(y)=x;orig(y)=0
    end subroutine add_arc
  end function min_cost_flow

end module igraph_flow_ext
