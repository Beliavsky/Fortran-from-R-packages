module igraph_cliques
  use igraph_graph, only : graph_t, are_adjacent
  implicit none
  private

  type, public :: vertex_set_t
    integer, allocatable :: v(:)
  end type vertex_set_t

  type, public :: vertex_set_list_t
    type(vertex_set_t), allocatable :: set(:)
    integer :: count = 0
  end type vertex_set_list_t

  public :: cliques, maximal_cliques, largest_cliques, clique_number

contains

  function cliques(g, min_size, max_size) result(out)
    type(graph_t), intent(in) :: g
    integer, intent(in), optional :: min_size, max_size
    type(vertex_set_list_t) :: out
    integer :: lo, hi, i
    integer, allocatable :: cur(:), cand(:)

    lo = 1
    if (present(min_size)) lo = max(0,min_size)
    hi = g%n
    if (present(max_size)) hi = min(g%n,max_size)
    if (hi < lo) then
      allocate(out%set(0))
      return
    end if
    allocate(out%set(0), cur(g%n), cand(g%n))
    cand = [(i,i=1,g%n)]
    if (lo == 0) call append_set(out, cur(:0))
    call enumerate(cur, 0, cand, g%n)

  contains
    recursive subroutine enumerate(r, nr, c, nc)
      integer, intent(inout) :: r(:)
      integer, intent(in) :: nr, c(:), nc
      integer :: pos, v, j, nn
      integer, allocatable :: next(:)
      if (nr >= hi) return
      do pos = 1, nc
        v = c(pos)
        r(nr+1) = v
        if (nr+1 >= lo .and. nr+1 <= hi) call append_set(out, r(:nr+1))
        if (nr+1 < hi .and. pos < nc) then
          allocate(next(nc-pos))
          nn = 0
          do j = pos+1, nc
            if (joined(g,v,c(j))) then
              nn = nn + 1
              next(nn) = c(j)
            end if
          end do
          if (nn > 0) call enumerate(r, nr+1, next, nn)
          deallocate(next)
        end if
      end do
    end subroutine enumerate
  end function cliques

  function maximal_cliques(g, min_size, max_size) result(out)
    type(graph_t), intent(in) :: g
    integer, intent(in), optional :: min_size, max_size
    type(vertex_set_list_t) :: out
    logical, allocatable :: p(:), x(:)
    integer, allocatable :: r(:)
    integer :: lo, hi

    lo = 1
    if (present(min_size)) lo = max(0,min_size)
    hi = g%n
    if (present(max_size)) hi = min(g%n,max_size)
    allocate(out%set(0), p(g%n), x(g%n), r(g%n))
    p = .true.
    x = .false.
    call bronk(r,0,p,x)

  contains
    recursive subroutine bronk(r, nr, pset, xset)
      integer, intent(inout) :: r(:)
      integer, intent(in) :: nr
      logical, intent(inout) :: pset(:), xset(:)
      logical, allocatable :: p2(:), x2(:), todo(:)
      integer :: v, u, best, countn, j

      if (.not. any(pset) .and. .not. any(xset)) then
        if (nr >= lo .and. nr <= hi) call append_set(out,r(:nr))
        return
      end if
      if (nr >= hi) return

      ! Pivot on P union X with the largest number of neighbors in P.
      u = 0
      best = -1
      do j = 1, g%n
        if (.not. pset(j) .and. .not. xset(j)) cycle
        countn = 0
        do v = 1, g%n
          if (pset(v)) then
            if (joined(g,j,v)) countn = countn + 1
          end if
        end do
        if (countn > best) then
          best = countn
          u = j
        end if
      end do
      allocate(todo(g%n))
      todo = pset
      if (u > 0) then
        do v = 1, g%n
          if (joined(g,u,v)) todo(v) = .false.
        end do
      end if

      do v = 1, g%n
        if (.not. todo(v)) cycle
        r(nr+1) = v
        allocate(p2(g%n),x2(g%n))
        p2 = .false.
        x2 = .false.
        do j = 1, g%n
          if (joined(g,v,j)) then
            p2(j) = pset(j)
            x2(j) = xset(j)
          end if
        end do
        call bronk(r,nr+1,p2,x2)
        deallocate(p2,x2)
        pset(v) = .false.
        xset(v) = .true.
      end do
      deallocate(todo)
    end subroutine bronk
  end function maximal_cliques

  function largest_cliques(g) result(out)
    type(graph_t), intent(in) :: g
    type(vertex_set_list_t) :: out
    type(vertex_set_list_t) :: all
    integer :: i, mx
    all = maximal_cliques(g)
    mx = 0
    do i = 1, all%count
      mx = max(mx,size(all%set(i)%v))
    end do
    allocate(out%set(0))
    do i = 1, all%count
      if (size(all%set(i)%v) == mx) call append_set(out,all%set(i)%v)
    end do
  end function largest_cliques

  integer function clique_number(g) result(omega)
    type(graph_t), intent(in) :: g
    type(vertex_set_list_t) :: all
    integer :: i
    all = maximal_cliques(g)
    omega = 0
    do i = 1, all%count
      omega = max(omega,size(all%set(i)%v))
    end do
  end function clique_number

  logical function joined(g,a,b) result(ok)
    type(graph_t), intent(in) :: g
    integer, intent(in) :: a,b
    ok = are_adjacent(g,a,b)
    if (g%directed .and. .not. ok) ok = are_adjacent(g,b,a)
  end function joined

  subroutine append_set(list, values)
    type(vertex_set_list_t), intent(inout) :: list
    integer, intent(in) :: values(:)
    type(vertex_set_t), allocatable :: tmp(:)
    integer :: n
    n = list%count
    allocate(tmp(n+1))
    if (n > 0) tmp(:n) = list%set
    allocate(tmp(n+1)%v(size(values)))
    if (size(values) > 0) tmp(n+1)%v = values
    call move_alloc(tmp,list%set)
    list%count = n+1
  end subroutine append_set

end module igraph_cliques
