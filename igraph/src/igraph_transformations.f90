module igraph_transformations
  use igraph_kinds, only : dp
  use igraph_graph, only : graph_t, make_graph, empty_graph, are_adjacent
  implicit none
  private
  public :: induced_subgraph, edge_subgraph, transpose_graph, complement_graph
  public :: disjoint_union, line_graph
  public :: cartesian_product, tensor_product, strong_product
  public :: graph_union, graph_intersection, graph_difference, graph_join
contains

  function induced_subgraph(g, vertices) result(h)
    type(graph_t), intent(in) :: g
    integer, intent(in) :: vertices(:)
    type(graph_t) :: h
    integer, allocatable :: map(:), e(:,:)
    real(dp), allocatable :: w(:)
    integer :: i, k, m, u, v

    allocate(map(g%n), source=0)
    do i = 1, size(vertices)
      if (vertices(i) < 1 .or. vertices(i) > g%n) error stop 'induced_subgraph: invalid vertex'
      if (map(vertices(i)) /= 0) error stop 'induced_subgraph: duplicate vertex'
      map(vertices(i)) = i
    end do
    m = 0
    do k = 1, g%m
      if (map(g%edge(1,k)) > 0 .and. map(g%edge(2,k)) > 0) m = m + 1
    end do
    allocate(e(2,m),w(m));m=0
    do k = 1, g%m
      u=map(g%edge(1,k));v=map(g%edge(2,k))
      if (u == 0 .or. v == 0) cycle
      m=m+1;e(:,m)=[u,v];w(m)=g%weight(k)
    end do
    h=make_graph(size(vertices),e,g%directed,w)
  end function induced_subgraph

  function edge_subgraph(g, edge_ids, delete_vertices) result(h)
    type(graph_t), intent(in) :: g
    integer, intent(in) :: edge_ids(:)
    logical, intent(in), optional :: delete_vertices
    type(graph_t) :: h
    logical :: del
    integer, allocatable :: e(:,:),map(:),used(:)
    real(dp), allocatable :: w(:)
    integer :: i,k,nnew
    del=.false.;if(present(delete_vertices))del=delete_vertices
    allocate(e(2,size(edge_ids)),w(size(edge_ids)))
    do i=1,size(edge_ids)
      k=edge_ids(i)
      if(k<1 .or. k>g%m)error stop 'edge_subgraph: invalid edge id'
      e(:,i)=g%edge(:,k);w(i)=g%weight(k)
    end do
    if(.not.del)then
      h=make_graph(g%n,e,g%directed,w)
      return
    end if
    allocate(used(g%n),source=0)
    do i=1,size(edge_ids)
      used(e(1,i))=1;used(e(2,i))=1
    end do
    allocate(map(g%n),source=0);nnew=0
    do i=1,g%n
      if(used(i)/=0)then;nnew=nnew+1;map(i)=nnew;end if
    end do
    do i=1,size(edge_ids)
      e(1,i)=map(e(1,i));e(2,i)=map(e(2,i))
    end do
    h=make_graph(nnew,e,g%directed,w)
  end function edge_subgraph

  function transpose_graph(g) result(h)
    type(graph_t), intent(in) :: g
    type(graph_t) :: h
    integer, allocatable :: e(:,:)
    if(.not.g%directed)then
      h=g
      return
    end if
    allocate(e(2,g%m))
    if(g%m>0)then
      e(1,:)=g%edge(2,:);e(2,:)=g%edge(1,:)
    end if
    h=make_graph(g%n,e,.true.,g%weight)
  end function transpose_graph

  function complement_graph(g, loops) result(h)
    type(graph_t), intent(in) :: g
    logical, intent(in), optional :: loops
    type(graph_t) :: h
    logical :: lp
    integer :: i,j,k,maxm
    integer, allocatable :: e(:,:)
    lp=.false.;if(present(loops))lp=loops
    if(g%directed)then
      maxm=g%n*g%n;if(.not.lp)maxm=g%n*(g%n-1)
    else
      maxm=g%n*(g%n-1)/2;if(lp)maxm=maxm+g%n
    end if
    allocate(e(2,maxm));k=0
    if(g%directed)then
      do i=1,g%n
        do j=1,g%n
          if(.not.lp .and. i==j)cycle
          if(.not.are_adjacent(g,i,j))then;k=k+1;e(:,k)=[i,j];end if
        end do
      end do
    else
      do i=1,g%n
        if(lp)then
          if(.not.are_adjacent(g,i,i))then;k=k+1;e(:,k)=[i,i];end if
        end if
        do j=i+1,g%n
          if(.not.are_adjacent(g,i,j))then;k=k+1;e(:,k)=[i,j];end if
        end do
      end do
    end if
    h=make_graph(g%n,e(:,:k),g%directed)
  end function complement_graph

  function disjoint_union(g1,g2) result(h)
    type(graph_t), intent(in) :: g1,g2
    type(graph_t) :: h
    integer, allocatable :: e(:,:)
    real(dp), allocatable :: w(:)
    if(g1%directed .neqv. g2%directed)error stop 'disjoint_union: directedness mismatch'
    allocate(e(2,g1%m+g2%m),w(g1%m+g2%m))
    if(g1%m>0)then;e(:,:g1%m)=g1%edge;w(:g1%m)=g1%weight;end if
    if(g2%m>0)then
      e(:,g1%m+1:)=g2%edge+g1%n
      w(g1%m+1:)=g2%weight
    end if
    h=make_graph(g1%n+g2%n,e,g1%directed,w)
  end function disjoint_union

  function line_graph(g) result(h)
    type(graph_t), intent(in) :: g
    type(graph_t) :: h
    integer, allocatable :: e(:,:)
    integer :: i,j,k,maxm
    logical :: touch
    if(g%directed)then
      maxm=g%m*max(0,g%m-1)
      allocate(e(2,maxm));k=0
      do i=1,g%m
        do j=1,g%m
          if(i==j)cycle
          if(g%edge(2,i)==g%edge(1,j))then;k=k+1;e(:,k)=[i,j];end if
        end do
      end do
      h=make_graph(g%m,e(:,:k),.true.)
    else
      maxm=g%m*max(0,g%m-1)/2
      allocate(e(2,maxm));k=0
      do i=1,g%m
        do j=i+1,g%m
          touch = g%edge(1,i)==g%edge(1,j) .or. g%edge(1,i)==g%edge(2,j) .or. &
                  g%edge(2,i)==g%edge(1,j) .or. g%edge(2,i)==g%edge(2,j)
          if(touch)then;k=k+1;e(:,k)=[i,j];end if
        end do
      end do
      h=make_graph(g%m,e(:,:k),.false.)
    end if
  end function line_graph

  function cartesian_product(g1,g2) result(h)
    type(graph_t),intent(in)::g1,g2
    type(graph_t)::h
    if(g1%directed .neqv. g2%directed)error stop 'cartesian_product: directedness mismatch'
    h=product_by_rule(g1,g2,1)
  end function cartesian_product

  function tensor_product(g1,g2) result(h)
    type(graph_t),intent(in)::g1,g2
    type(graph_t)::h
    if(g1%directed .neqv. g2%directed)error stop 'tensor_product: directedness mismatch'
    h=product_by_rule(g1,g2,2)
  end function tensor_product

  function strong_product(g1,g2) result(h)
    type(graph_t),intent(in)::g1,g2
    type(graph_t)::h
    if(g1%directed .neqv. g2%directed)error stop 'strong_product: directedness mismatch'
    h=product_by_rule(g1,g2,3)
  end function strong_product

  function product_by_rule(g1,g2,rule) result(h)
    type(graph_t),intent(in)::g1,g2
    integer,intent(in)::rule
    type(graph_t)::h
    integer,allocatable::e(:,:)
    integer::a,b,u1,u2,v1,v2,k,maxm,n
    logical::adj1,adj2,yes
    n=g1%n*g2%n
    if(g1%directed)then;maxm=n*n;else;maxm=n*max(0,n-1)/2+n;end if
    allocate(e(2,maxm));k=0
    if(g1%directed)then
      do a=1,n
        call decode(a,u1,u2)
        do b=1,n
          call decode(b,v1,v2)
          adj1=are_adjacent(g1,u1,v1);adj2=are_adjacent(g2,u2,v2)
          select case(rule)
          case(1);yes=(u1==v1 .and. adj2) .or. (u2==v2 .and. adj1)
          case(2);yes=adj1 .and. adj2
          case default;yes=(u1==v1 .and. adj2) .or. (u2==v2 .and. adj1) .or. (adj1 .and. adj2)
          end select
          if(yes)then;k=k+1;e(:,k)=[a,b];end if
        end do
      end do
    else
      do a=1,n
        call decode(a,u1,u2)
        do b=a,n
          call decode(b,v1,v2)
          adj1=are_adjacent(g1,u1,v1);adj2=are_adjacent(g2,u2,v2)
          select case(rule)
          case(1);yes=(u1==v1 .and. adj2) .or. (u2==v2 .and. adj1)
          case(2);yes=adj1 .and. adj2
          case default;yes=(u1==v1 .and. adj2) .or. (u2==v2 .and. adj1) .or. (adj1 .and. adj2)
          end select
          if(yes)then;k=k+1;e(:,k)=[a,b];end if
        end do
      end do
    end if
    h=make_graph(n,e(:,:k),g1%directed)
  contains
    subroutine decode(x,a1,a2)
      integer,intent(in)::x
      integer,intent(out)::a1,a2
      a1=(x-1)/g2%n+1;a2=mod(x-1,g2%n)+1
    end subroutine decode
  end function product_by_rule

  function graph_union(g1,g2) result(h)
    type(graph_t),intent(in)::g1,g2
    type(graph_t)::h
    h=set_operation(g1,g2,1)
  end function graph_union

  function graph_intersection(g1,g2) result(h)
    type(graph_t),intent(in)::g1,g2
    type(graph_t)::h
    h=set_operation(g1,g2,2)
  end function graph_intersection

  function graph_difference(g1,g2) result(h)
    type(graph_t),intent(in)::g1,g2
    type(graph_t)::h
    h=set_operation(g1,g2,3)
  end function graph_difference

  function set_operation(g1,g2,op) result(h)
    type(graph_t),intent(in)::g1,g2
    integer,intent(in)::op
    type(graph_t)::h
    integer,allocatable::e(:,:)
    integer::i,j,k,maxm
    logical::a,b,yes
    if(g1%n/=g2%n .or. g1%directed .neqv. g2%directed)error stop 'graph set operation: incompatible graphs'
    if(g1%directed)then;maxm=g1%n*g1%n;else;maxm=g1%n*(g1%n+1)/2;end if
    allocate(e(2,maxm));k=0
    if(g1%directed)then
      do i=1,g1%n;do j=1,g1%n
        a=are_adjacent(g1,i,j);b=are_adjacent(g2,i,j)
        select case(op);case(1);yes=a.or.b;case(2);yes=a.and.b;case default;yes=a.and.(.not.b);end select
        if(yes)then;k=k+1;e(:,k)=[i,j];end if
      end do;end do
    else
      do i=1,g1%n;do j=i,g1%n
        a=are_adjacent(g1,i,j);b=are_adjacent(g2,i,j)
        select case(op);case(1);yes=a.or.b;case(2);yes=a.and.b;case default;yes=a.and.(.not.b);end select
        if(yes)then;k=k+1;e(:,k)=[i,j];end if
      end do;end do
    end if
    h=make_graph(g1%n,e(:,:k),g1%directed)
  end function set_operation

  function graph_join(g1,g2) result(h)
    type(graph_t),intent(in)::g1,g2
    type(graph_t)::h
    integer,allocatable::e(:,:)
    real(dp),allocatable::w(:)
    integer::k,i,j,m
    if(g1%directed .or. g2%directed)error stop 'graph_join: undirected graphs required'
    m=g1%m+g2%m+g1%n*g2%n
    allocate(e(2,m),w(m));k=0
    do i=1,g1%m;k=k+1;e(:,k)=g1%edge(:,i);w(k)=g1%weight(i);end do
    do i=1,g2%m;k=k+1;e(:,k)=g2%edge(:,i)+g1%n;w(k)=g2%weight(i);end do
    do i=1,g1%n;do j=1,g2%n;k=k+1;e(:,k)=[i,g1%n+j];w(k)=1.0_dp;end do;end do
    h=make_graph(g1%n+g2%n,e,.false.,w)
  end function graph_join
end module igraph_transformations
