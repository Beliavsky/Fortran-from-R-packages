module igraph_graph
  use igraph_kinds, only : dp, i8
  implicit none
  private

  type, public :: graph_t
    integer :: n = 0
    integer :: m = 0
    logical :: directed = .false.
    integer, allocatable :: edge(:,:)      ! (2,m), 1-based vertices
    real(dp), allocatable :: weight(:)
    integer, allocatable :: out_ptr(:)
    integer, allocatable :: out_nei(:)
    integer, allocatable :: out_eid(:)
    real(dp), allocatable :: out_w(:)
    integer, allocatable :: in_ptr(:)
    integer, allocatable :: in_nei(:)
    integer, allocatable :: in_eid(:)
    real(dp), allocatable :: in_w(:)
  contains
    procedure :: finalize => graph_finalize
    procedure :: validate => graph_validate
  end type graph_t

  public :: make_graph, empty_graph, full_graph, ring_graph, star_graph
  public :: add_edges, simplify_graph, are_adjacent, edge_count, vertex_count
  public :: degree, strength, edgelist

contains

  function make_graph(n, edges, directed, weights) result(g)
    integer, intent(in) :: n
    integer, intent(in) :: edges(:,:)
    logical, intent(in), optional :: directed
    real(dp), intent(in), optional :: weights(:)
    type(graph_t) :: g
    integer :: m

    if (n < 0) error stop 'make_graph: n must be nonnegative'
    if (size(edges,1) /= 2) error stop 'make_graph: edges must have shape (2,m)'
    m = size(edges,2)
    g%n = n
    g%m = m
    g%directed = .false.
    if (present(directed)) g%directed = directed
    allocate(g%edge(2,m), g%weight(m))
    if (m > 0) g%edge = edges
    if (present(weights)) then
      if (size(weights) /= m) error stop 'make_graph: weights length mismatch'
      if (m > 0) g%weight = weights
    else
      if (m > 0) g%weight = 1.0_dp
    end if
    call g%validate()
    call g%finalize()
  end function make_graph

  function empty_graph(n, directed) result(g)
    integer, intent(in) :: n
    logical, intent(in), optional :: directed
    type(graph_t) :: g
    integer, allocatable :: e(:,:)
    allocate(e(2,0))
    g = make_graph(n, e, directed)
  end function empty_graph

  function full_graph(n, directed, loops) result(g)
    integer, intent(in) :: n
    logical, intent(in), optional :: directed, loops
    type(graph_t) :: g
    logical :: dir, lp
    integer :: i, j, k, m
    integer, allocatable :: e(:,:)
    dir = .false.; if (present(directed)) dir = directed
    lp = .false.; if (present(loops)) lp = loops
    if (n < 0) error stop 'full_graph: n must be nonnegative'
    if (dir) then
      m = n*n
      if (.not. lp) m = n*(n-1)
      allocate(e(2,m)); k = 0
      do i = 1, n
        do j = 1, n
          if (.not. lp .and. i == j) cycle
          k = k + 1; e(:,k) = [i,j]
        end do
      end do
    else
      m = n*(n-1)/2
      if (lp) m = m + n
      allocate(e(2,m)); k = 0
      do i = 1, n
        if (lp) then
          k = k + 1; e(:,k) = [i,i]
        end if
        do j = i+1, n
          k = k + 1; e(:,k) = [i,j]
        end do
      end do
    end if
    g = make_graph(n, e, dir)
  end function full_graph

  function ring_graph(n, directed, mutual, circular) result(g)
    integer, intent(in) :: n
    logical, intent(in), optional :: directed, mutual, circular
    type(graph_t) :: g
    logical :: dir, mut, cir
    integer :: m0, m, i, k
    integer, allocatable :: e(:,:)
    dir = .false.; if (present(directed)) dir = directed
    mut = .false.; if (present(mutual)) mut = mutual
    cir = .true.; if (present(circular)) cir = circular
    if (n <= 0) then
      g = empty_graph(max(n,0), dir); return
    end if
    m0 = max(0,n-1); if (cir .and. n > 1) m0 = n
    m = m0
    if (dir .and. mut) m = 2*m0
    allocate(e(2,m)); k = 0
    do i = 1, n-1
      k = k + 1; e(:,k) = [i,i+1]
      if (dir .and. mut) then
        k = k + 1; e(:,k) = [i+1,i]
      end if
    end do
    if (cir .and. n > 1) then
      k = k + 1; e(:,k) = [n,1]
      if (dir .and. mut) then
        k = k + 1; e(:,k) = [1,n]
      end if
    end if
    g = make_graph(n, e, dir)
  end function ring_graph

  function star_graph(n, center, directed, mode) result(g)
    integer, intent(in) :: n
    integer, intent(in), optional :: center
    logical, intent(in), optional :: directed
    character(len=*), intent(in), optional :: mode
    type(graph_t) :: g
    integer :: c, i, k, m
    logical :: dir
    character(len=8) :: md
    integer, allocatable :: e(:,:)
    if (n < 0) error stop 'star_graph: n must be nonnegative'
    c = 1; if (present(center)) c = center
    if (n > 0 .and. (c < 1 .or. c > n)) error stop 'star_graph: invalid center'
    dir = .false.; if (present(directed)) dir = directed
    md = 'out'; if (present(mode)) md = adjustl(mode)
    m = max(0,n-1); if (dir .and. trim(md) == 'mutual') m = 2*m
    allocate(e(2,m)); k = 0
    do i = 1, n
      if (i == c) cycle
      select case (trim(md))
      case ('in')
        k = k + 1; e(:,k) = [i,c]
      case ('mutual')
        if (dir) then
          k = k + 1; e(:,k) = [c,i]
          k = k + 1; e(:,k) = [i,c]
        else
          k = k + 1; e(:,k) = [c,i]
        end if
      case default
        k = k + 1; e(:,k) = [c,i]
      end select
    end do
    g = make_graph(n, e, dir)
  end function star_graph

  subroutine graph_validate(self)
    class(graph_t), intent(in) :: self
    integer :: k
    if (self%n < 0 .or. self%m < 0) error stop 'graph_validate: negative size'
    if (.not. allocated(self%edge)) return
    if (size(self%edge,1) /= 2 .or. size(self%edge,2) /= self%m) then
      error stop 'graph_validate: edge array mismatch'
    end if
    do k = 1, self%m
      if (self%edge(1,k) < 1 .or. self%edge(1,k) > self%n .or. &
          self%edge(2,k) < 1 .or. self%edge(2,k) > self%n) then
        error stop 'graph_validate: vertex id out of range'
      end if
    end do
  end subroutine graph_validate

  subroutine graph_finalize(self)
    class(graph_t), intent(inout) :: self
    integer :: i, k, u, v, total
    integer, allocatable :: cnt(:), pos(:)

    if (allocated(self%out_ptr)) deallocate(self%out_ptr, self%out_nei, self%out_eid, self%out_w)
    if (allocated(self%in_ptr)) deallocate(self%in_ptr, self%in_nei, self%in_eid, self%in_w)
    allocate(cnt(self%n), source=0)
    do k = 1, self%m
      u = self%edge(1,k); v = self%edge(2,k)
      cnt(u) = cnt(u) + 1
      if (.not. self%directed .and. v /= u) cnt(v) = cnt(v) + 1
    end do
    allocate(self%out_ptr(self%n+1)); self%out_ptr(1) = 1
    do i = 1, self%n
      self%out_ptr(i+1) = self%out_ptr(i) + cnt(i)
    end do
    total = self%out_ptr(self%n+1)-1
    allocate(self%out_nei(total), self%out_eid(total), self%out_w(total))
    allocate(pos(self%n)); pos = self%out_ptr(1:self%n)
    do k = 1, self%m
      u = self%edge(1,k); v = self%edge(2,k)
      i = pos(u); self%out_nei(i) = v; self%out_eid(i) = k; self%out_w(i) = self%weight(k); pos(u)=i+1
      if (.not. self%directed .and. v /= u) then
        i = pos(v); self%out_nei(i) = u; self%out_eid(i) = k; self%out_w(i) = self%weight(k); pos(v)=i+1
      end if
    end do
    deallocate(cnt, pos)

    if (self%directed) then
      allocate(cnt(self%n), source=0)
      do k = 1, self%m
        cnt(self%edge(2,k)) = cnt(self%edge(2,k)) + 1
      end do
      allocate(self%in_ptr(self%n+1)); self%in_ptr(1)=1
      do i = 1, self%n
        self%in_ptr(i+1)=self%in_ptr(i)+cnt(i)
      end do
      total = self%in_ptr(self%n+1)-1
      allocate(self%in_nei(total), self%in_eid(total), self%in_w(total), pos(self%n))
      pos = self%in_ptr(1:self%n)
      do k = 1, self%m
        u=self%edge(1,k); v=self%edge(2,k); i=pos(v)
        self%in_nei(i)=u; self%in_eid(i)=k; self%in_w(i)=self%weight(k); pos(v)=i+1
      end do
      deallocate(cnt,pos)
    else
      allocate(self%in_ptr(size(self%out_ptr)), self%in_nei(size(self%out_nei)), &
               self%in_eid(size(self%out_eid)), self%in_w(size(self%out_w)))
      self%in_ptr=self%out_ptr; self%in_nei=self%out_nei; self%in_eid=self%out_eid; self%in_w=self%out_w
    end if
  end subroutine graph_finalize

  subroutine add_edges(g, edges, weights)
    type(graph_t), intent(inout) :: g
    integer, intent(in) :: edges(:,:)
    real(dp), intent(in), optional :: weights(:)
    integer :: oldm, addm
    integer, allocatable :: e2(:,:)
    real(dp), allocatable :: w2(:)
    if (size(edges,1) /= 2) error stop 'add_edges: edges must be (2,m)'
    addm=size(edges,2); oldm=g%m
    if (present(weights)) then
      if (size(weights) /= addm) error stop 'add_edges: weights length mismatch'
    end if
    allocate(e2(2,oldm+addm), w2(oldm+addm))
    if (oldm>0) then; e2(:,1:oldm)=g%edge; w2(1:oldm)=g%weight; end if
    if (addm>0) then
      e2(:,oldm+1:)=edges
      if (present(weights)) then; w2(oldm+1:)=weights; else; w2(oldm+1:)=1.0_dp; end if
    end if
    call move_alloc(e2,g%edge); call move_alloc(w2,g%weight); g%m=oldm+addm
    call g%validate(); call g%finalize()
  end subroutine add_edges

  function simplify_graph(g, remove_loops, combine_weights) result(h)
    type(graph_t), intent(in) :: g
    logical, intent(in), optional :: remove_loops, combine_weights
    type(graph_t) :: h
    logical :: rl, cw, found
    integer :: k, j, keep, u, v, a, b
    integer, allocatable :: e(:,:)
    real(dp), allocatable :: w(:)
    rl=.true.; if (present(remove_loops)) rl=remove_loops
    cw=.true.; if (present(combine_weights)) cw=combine_weights
    allocate(e(2,g%m), w(g%m)); keep=0
    do k=1,g%m
      u=g%edge(1,k); v=g%edge(2,k)
      if (rl .and. u==v) cycle
      a=u; b=v
      if (.not. g%directed .and. a>b) then; j=a; a=b; b=j; end if
      found=.false.
      do j=1,keep
        if (e(1,j)==a .and. e(2,j)==b) then
          found=.true.; if (cw) w(j)=w(j)+g%weight(k); exit
        end if
      end do
      if (.not. found) then
        keep=keep+1; e(:,keep)=[a,b]; w(keep)=g%weight(k)
      end if
    end do
    h=make_graph(g%n,e(:,1:keep),g%directed,w(1:keep))
  end function simplify_graph

  logical function are_adjacent(g, u, v)
    type(graph_t), intent(in) :: g
    integer, intent(in) :: u, v
    integer :: p
    are_adjacent=.false.
    if (u<1 .or. u>g%n .or. v<1 .or. v>g%n) return
    do p=g%out_ptr(u),g%out_ptr(u+1)-1
      if (g%out_nei(p)==v) then; are_adjacent=.true.; return; end if
    end do
  end function are_adjacent

  integer function edge_count(g)
    type(graph_t), intent(in) :: g
    edge_count=g%m
  end function edge_count

  integer function vertex_count(g)
    type(graph_t), intent(in) :: g
    vertex_count=g%n
  end function vertex_count

  function edgelist(g) result(e)
    type(graph_t), intent(in) :: g
    integer, allocatable :: e(:,:)
    allocate(e(2,g%m)); if (g%m>0) e=g%edge
  end function edgelist

  function degree(g, mode, loops) result(d)
    type(graph_t), intent(in) :: g
    character(len=*), intent(in), optional :: mode
    logical, intent(in), optional :: loops
    integer, allocatable :: d(:)
    character(len=8) :: md
    logical :: lp
    integer :: i, k, u, v
    allocate(d(g%n),source=0)
    md='all'; if (present(mode)) md=adjustl(mode)
    lp=.true.; if (present(loops)) lp=loops
    if (.not. g%directed) then
      do k=1,g%m
        u=g%edge(1,k); v=g%edge(2,k)
        if (u==v) then
          if (lp) d(u)=d(u)+2
        else
          d(u)=d(u)+1; d(v)=d(v)+1
        end if
      end do
    else
      select case(trim(md))
      case('out')
        do i=1,g%n; d(i)=g%out_ptr(i+1)-g%out_ptr(i); end do
      case('in')
        do i=1,g%n; d(i)=g%in_ptr(i+1)-g%in_ptr(i); end do
      case default
        do k=1,g%m
          u=g%edge(1,k); v=g%edge(2,k)
          if (u==v) then
            if (lp) d(u)=d(u)+2
          else
            d(u)=d(u)+1; d(v)=d(v)+1
          end if
        end do
      end select
      if (.not. lp) then
        do k=1,g%m
          if (g%edge(1,k)==g%edge(2,k)) then
            select case(trim(md))
            case('out','in'); d(g%edge(1,k))=d(g%edge(1,k))-1
            case default; d(g%edge(1,k))=d(g%edge(1,k))-2
            end select
          end if
        end do
      end if
    end if
  end function degree

  function strength(g, mode, loops) result(s)
    type(graph_t), intent(in) :: g
    character(len=*), intent(in), optional :: mode
    logical, intent(in), optional :: loops
    real(dp), allocatable :: s(:)
    character(len=8) :: md
    logical :: lp
    integer :: k, u, v
    allocate(s(g%n),source=0.0_dp)
    md='all'; if (present(mode)) md=adjustl(mode)
    lp=.true.; if (present(loops)) lp=loops
    do k=1,g%m
      u=g%edge(1,k); v=g%edge(2,k)
      if (u==v .and. .not. lp) cycle
      if (.not. g%directed) then
        if (u==v) then; s(u)=s(u)+2.0_dp*g%weight(k)
        else; s(u)=s(u)+g%weight(k); s(v)=s(v)+g%weight(k); end if
      else
        select case(trim(md))
        case('out'); s(u)=s(u)+g%weight(k)
        case('in'); s(v)=s(v)+g%weight(k)
        case default
          if (u==v) then; s(u)=s(u)+2.0_dp*g%weight(k)
          else; s(u)=s(u)+g%weight(k); s(v)=s(v)+g%weight(k); end if
        end select
      end if
    end do
  end function strength

end module igraph_graph
