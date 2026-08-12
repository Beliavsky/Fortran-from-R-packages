module igraph_canonical
  use igraph_kinds, only : i8
  use igraph_graph, only : graph_t, make_graph, degree, are_adjacent
  implicit none
  private

  type, public :: canonical_result_t
    integer, allocatable :: permutation(:) ! canonical position -> original vertex
    integer, allocatable :: inverse(:)     ! original vertex -> canonical position
    integer, allocatable :: code(:)        ! flattened 0/1 adjacency in canonical order
    integer(i8) :: states = 0_i8
  end type canonical_result_t

  type, public :: automorphism_result_t
    integer(i8) :: count = 0_i8
    integer :: stored = 0
    logical :: truncated = .false.
    integer, allocatable :: maps(:,:) ! maps(original vertex, automorphism index)
  end type automorphism_result_t

  public :: canonical_labeling, canonical_code, canonical_form
  public :: automorphism_count, automorphisms

contains

  function canonical_labeling(g, colors) result(res)
    type(graph_t), intent(in) :: g
    integer, intent(in), optional :: colors(:)
    type(canonical_result_t) :: res
    integer, allocatable :: init(:), best_perm(:), best_code(:)
    integer :: i
    logical :: have_best

    if (present(colors)) then
      if (size(colors) /= g%n) error stop 'canonical_labeling: colors length mismatch'
      if (any(colors < 0)) error stop 'canonical_labeling: colors must be nonnegative'
    end if
    allocate(init(g%n))
    call initial_colors(g, colors, init)
    allocate(best_perm(g%n), best_code(g%n*g%n))
    best_perm = 0
    best_code = 0
    have_best = .false.
    call search(init)
    allocate(res%permutation(g%n), res%inverse(g%n), res%code(g%n*g%n))
    res%permutation = best_perm
    res%code = best_code
    res%inverse = 0
    do i = 1, g%n
      if (g%n > 0) res%inverse(best_perm(i)) = i
    end do

  contains

    recursive subroutine search(color_in)
      integer, intent(in) :: color_in(:)
      integer, allocatable :: color(:), perm(:), code(:), child(:), members(:)
      integer :: nc, c, sz, chosen, v, newc

      res%states = res%states + 1_i8
      allocate(color(size(color_in)))
      color = color_in
      call refine_colors(g, color)
      if (all_singletons(color)) then
        call permutation_from_colors(color, perm)
        call graph_code(g, perm, code)
        if (.not. have_best) then
          best_code = code
          best_perm = perm
          have_best = .true.
        else if (lex_less(code, best_code)) then
          best_code = code
          best_perm = perm
        end if
        return
      end if

      nc = maxval(color)
      chosen = 0
      do c = 1, nc
        sz = count(color == c)
        if (sz > 1) then
          if (chosen == 0) then
            chosen = c
          else if (sz < count(color == chosen)) then
            chosen = c
          end if
        end if
      end do
      allocate(members(count(color == chosen)))
      members = pack([(v, v=1,size(color))], color == chosen)
      newc = nc + 1
      do c = 1, size(members)
        allocate(child(size(color)))
        child = color
        child(members(c)) = newc
        call search(child)
        deallocate(child)
      end do
    end subroutine search

  end function canonical_labeling

  function canonical_code(g, colors) result(code)
    type(graph_t), intent(in) :: g
    integer, intent(in), optional :: colors(:)
    integer, allocatable :: code(:)
    type(canonical_result_t) :: r
    if (present(colors)) then
      r = canonical_labeling(g, colors)
    else
      r = canonical_labeling(g)
    end if
    allocate(code(size(r%code)))
    code = r%code
  end function canonical_code

  function canonical_form(g, colors) result(h)
    type(graph_t),intent(in)::g
    integer,intent(in),optional::colors(:)
    type(graph_t)::h
    type(canonical_result_t)::r
    integer,allocatable::e(:,:)
    integer::k
    if(present(colors))then;r=canonical_labeling(g,colors);else;r=canonical_labeling(g);end if
    allocate(e(2,g%m))
    do k=1,g%m
      e(1,k)=r%inverse(g%edge(1,k));e(2,k)=r%inverse(g%edge(2,k))
    end do
    h=make_graph(g%n,e,g%directed,g%weight)
  end function canonical_form

  function automorphisms(g,colors,max_maps) result(res)
    type(graph_t),intent(in)::g
    integer,intent(in),optional::colors(:)
    integer,intent(in),optional::max_maps
    type(automorphism_result_t)::res
    integer,allocatable::map(:),inv(:),order(:),dout(:),din(:),storage(:,:)
    integer::i,j,t,limit

    if(present(colors))then
      if(size(colors)/=g%n)error stop 'automorphisms: colors length mismatch'
    end if
    limit=10000;if(present(max_maps))limit=max(0,max_maps)
    allocate(map(g%n),inv(g%n),order(g%n),storage(g%n,limit));map=0;inv=0;storage=0
    order=[(i,i=1,g%n)];dout=degree(g,'out',.true.)
    if(g%directed)then;din=degree(g,'in',.true.);else;allocate(din(g%n),source=0);end if
    do i=2,g%n
      t=order(i);j=i-1
      do while(j>=1)
        if(dout(order(j))+din(order(j))>=dout(t)+din(t))exit
        order(j+1)=order(j);j=j-1
      end do
      order(j+1)=t
    end do
    call enum(1)
    allocate(res%maps(g%n,res%stored))
    if(res%stored>0)res%maps=storage(:,:res%stored)
    res%truncated=res%count>int(res%stored,i8)
  contains
    recursive subroutine enum(depth)
      integer,intent(in)::depth
      integer::u,v
      if(depth>g%n)then
        res%count=res%count+1_i8
        if(res%stored<limit)then;res%stored=res%stored+1;storage(:,res%stored)=map;end if
        return
      end if
      u=order(depth)
      do v=1,g%n
        if(inv(v)/=0)cycle
        if(.not.candidate(u,v))cycle
        map(u)=v;inv(v)=u;call enum(depth+1);map(u)=0;inv(v)=0
      end do
    end subroutine enum
    logical function candidate(a,b) result(ok)
      integer,intent(in)::a,b
      integer::x,y
      ok=.false.
      if(dout(a)/=dout(b))return
      if(g%directed)then;if(din(a)/=din(b))return;end if
      if(present(colors))then;if(colors(a)/=colors(b))return;end if
      if(are_adjacent(g,a,a) .neqv. are_adjacent(g,b,b))return
      do x=1,g%n
        y=map(x);if(y==0)cycle
        if(are_adjacent(g,a,x) .neqv. are_adjacent(g,b,y))return
        if(g%directed)then
          if(are_adjacent(g,x,a) .neqv. are_adjacent(g,y,b))return
        end if
      end do
      ok=.true.
    end function candidate
  end function automorphisms

  function automorphism_count(g, colors, max_count) result(count_aut)
    type(graph_t), intent(in) :: g
    integer, intent(in), optional :: colors(:)
    integer(i8), intent(in), optional :: max_count
    integer(i8) :: count_aut
    integer, allocatable :: map(:), inv(:), order(:), dout(:), din(:)
    integer :: i, j, t
    integer(i8) :: cap

    if (present(colors)) then
      if (size(colors) /= g%n) error stop 'automorphism_count: colors length mismatch'
    end if
    allocate(map(g%n), inv(g%n), order(g%n))
    map = 0
    inv = 0
    order = [(i, i=1,g%n)]
    dout = degree(g, 'out', .true.)
    if (g%directed) then
      din = degree(g, 'in', .true.)
    else
      allocate(din(g%n), source=0)
    end if
    do i = 2, g%n
      t = order(i)
      j = i - 1
      do while (j >= 1)
        if (score(order(j)) >= score(t)) exit
        order(j+1) = order(j)
        j = j - 1
      end do
      order(j+1) = t
    end do
    cap = huge(0_i8)
    if (present(max_count)) cap = max_count
    count_aut = 0_i8
    call enumerate(1)

  contains

    integer function score(v) result(s)
      integer, intent(in) :: v
      s = dout(v) + din(v)
    end function score

    recursive subroutine enumerate(depth)
      integer, intent(in) :: depth
      integer :: u, v
      if (count_aut >= cap) return
      if (depth > g%n) then
        count_aut = count_aut + 1_i8
        return
      end if
      u = order(depth)
      do v = 1, g%n
        if (inv(v) /= 0) cycle
        if (.not. auto_candidate(u,v)) cycle
        map(u) = v
        inv(v) = u
        call enumerate(depth+1)
        map(u) = 0
        inv(v) = 0
        if (count_aut >= cap) return
      end do
    end subroutine enumerate

    logical function auto_candidate(a,b) result(ok)
      integer, intent(in) :: a,b
      integer :: x, y
      ok = .false.
      if (dout(a) /= dout(b)) return
      if (g%directed) then
        if (din(a) /= din(b)) return
      end if
      if (present(colors)) then
        if (colors(a) /= colors(b)) return
      end if
      if (are_adjacent(g,a,a) .neqv. are_adjacent(g,b,b)) return
      do x = 1, g%n
        y = map(x)
        if (y == 0) cycle
        if (are_adjacent(g,a,x) .neqv. are_adjacent(g,b,y)) return
        if (g%directed) then
          if (are_adjacent(g,x,a) .neqv. are_adjacent(g,y,b)) return
        end if
      end do
      ok = .true.
    end function auto_candidate

  end function automorphism_count

  subroutine initial_colors(g, user, color)
    type(graph_t), intent(in) :: g
    integer, intent(in), optional :: user(:)
    integer, intent(out) :: color(:)
    integer, allocatable :: dout(:), din(:), raw(:,:), idx(:)
    integer :: i, j, t, nc

    dout = degree(g, 'out', .true.)
    if (g%directed) then
      din = degree(g, 'in', .true.)
    else
      allocate(din(g%n), source=0)
    end if
    allocate(raw(4,g%n), idx(g%n))
    do i = 1, g%n
      raw(1,i) = dout(i)
      raw(2,i) = din(i)
      raw(3,i) = merge(1,0,are_adjacent(g,i,i))
      raw(4,i) = 0
      if (present(user)) raw(4,i) = user(i)
      idx(i) = i
    end do
    do i = 2, g%n
      t = idx(i)
      j = i - 1
      do while (j >= 1)
        if (.not. sig_less(raw(:,t), raw(:,idx(j)))) exit
        idx(j+1) = idx(j)
        j = j - 1
      end do
      idx(j+1) = t
    end do
    nc = 0
    do i = 1, g%n
      if (i == 1) then
        nc = nc + 1
      else if (any(raw(:,idx(i)) /= raw(:,idx(i-1)))) then
        nc = nc + 1
      end if
      color(idx(i)) = nc
    end do
  end subroutine initial_colors

  subroutine refine_colors(g, color)
    type(graph_t), intent(in) :: g
    integer, intent(inout) :: color(:)
    integer, allocatable :: counts_out(:,:), counts_in(:,:), idx(:), newc(:), old(:)
    integer :: n, nc, i, j, t, p, v, u, nextc
    logical :: changed

    n = g%n
    if (n == 0) return
    do
      old = color
      nc = maxval(color)
      allocate(counts_out(nc,n), source=0)
      allocate(counts_in(nc,n), source=0)
      do u = 1, n
        do p = g%out_ptr(u), g%out_ptr(u+1)-1
          v = g%out_nei(p)
          counts_out(color(v),u) = counts_out(color(v),u) + 1
        end do
        if (g%directed) then
          do p = g%in_ptr(u), g%in_ptr(u+1)-1
            v = g%in_nei(p)
            counts_in(color(v),u) = counts_in(color(v),u) + 1
          end do
        end if
      end do
      allocate(idx(n), newc(n))
      idx = [(i, i=1,n)]
      do i = 2, n
        t = idx(i)
        j = i - 1
        do while (j >= 1)
          if (.not. vertex_sig_less(t, idx(j), color, counts_out, counts_in)) exit
          idx(j+1) = idx(j)
          j = j - 1
        end do
        idx(j+1) = t
      end do
      nextc = 0
      do i = 1, n
        if (i == 1) then
          nextc = 1
        else if (.not. same_vertex_sig(idx(i),idx(i-1),color,counts_out,counts_in)) then
          nextc = nextc + 1
        end if
        newc(idx(i)) = nextc
      end do
      color = newc
      changed = any(color /= old)
      deallocate(counts_out, counts_in, idx, newc)
      if (.not. changed) exit
    end do
  end subroutine refine_colors

  logical function vertex_sig_less(a,b,color,cout,cin) result(less)
    integer, intent(in) :: a,b,color(:),cout(:,:),cin(:,:)
    integer :: c
    less = .false.
    if (color(a) /= color(b)) then
      less = color(a) < color(b)
      return
    end if
    do c = 1, size(cout,1)
      if (cout(c,a) /= cout(c,b)) then
        less = cout(c,a) < cout(c,b)
        return
      end if
    end do
    do c = 1, size(cin,1)
      if (cin(c,a) /= cin(c,b)) then
        less = cin(c,a) < cin(c,b)
        return
      end if
    end do
  end function vertex_sig_less

  logical function same_vertex_sig(a,b,color,cout,cin) result(same)
    integer, intent(in) :: a,b,color(:),cout(:,:),cin(:,:)
    same = color(a) == color(b) .and. all(cout(:,a) == cout(:,b)) .and. all(cin(:,a) == cin(:,b))
  end function same_vertex_sig

  logical function all_singletons(color) result(ok)
    integer, intent(in) :: color(:)
    integer :: c
    ok = .true.
    if (size(color) == 0) return
    do c = 1, maxval(color)
      if (count(color == c) > 1) then
        ok = .false.
        return
      end if
    end do
  end function all_singletons

  subroutine permutation_from_colors(color, perm)
    integer, intent(in) :: color(:)
    integer, allocatable, intent(out) :: perm(:)
    integer :: i, c
    allocate(perm(size(color)))
    do i = 1, size(color)
      do c = 1, size(color)
        if (color(c) == i) then
          perm(i) = c
          exit
        end if
      end do
    end do
  end subroutine permutation_from_colors

  subroutine graph_code(g, perm, code)
    type(graph_t), intent(in) :: g
    integer, intent(in) :: perm(:)
    integer, allocatable, intent(out) :: code(:)
    integer :: i,j,k
    allocate(code(g%n*g%n))
    k = 0
    do i = 1, g%n
      do j = 1, g%n
        k = k + 1
        code(k) = merge(1,0,are_adjacent(g,perm(i),perm(j)))
      end do
    end do
  end subroutine graph_code

  logical function lex_less(a,b) result(less)
    integer, intent(in) :: a(:), b(:)
    integer :: i
    less = .false.
    do i = 1, min(size(a),size(b))
      if (a(i) < b(i)) then
        less = .true.
        return
      else if (a(i) > b(i)) then
        return
      end if
    end do
    less = size(a) < size(b)
  end function lex_less

  logical function sig_less(a,b) result(less)
    integer, intent(in) :: a(:), b(:)
    integer :: i
    less = .false.
    do i = 1, min(size(a),size(b))
      if (a(i) < b(i)) then
        less = .true.
        return
      else if (a(i) > b(i)) then
        return
      end if
    end do
  end function sig_less

end module igraph_canonical
