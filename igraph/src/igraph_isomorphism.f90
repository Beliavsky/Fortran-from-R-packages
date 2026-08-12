module igraph_isomorphism
  use igraph_kinds, only : i8
  use igraph_graph, only : graph_t, degree, are_adjacent
  implicit none
  private

  type, public :: isomorphism_result_t
    logical :: isomorphic = .false.
    integer, allocatable :: map12(:)
    integer, allocatable :: map21(:)
    integer(i8) :: states = 0_i8
  end type isomorphism_result_t

  public :: vf2_isomorphic, vf2_subisomorphic

contains

  function vf2_isomorphic(g1, g2, color1, color2) result(res)
    type(graph_t), intent(in) :: g1, g2
    integer, intent(in), optional :: color1(:), color2(:)
    type(isomorphism_result_t) :: res

    if (g1%n /= g2%n .or. g1%m /= g2%m .or. g1%directed .neqv. g2%directed) then
      allocate(res%map12(g1%n), source=0)
      allocate(res%map21(g2%n), source=0)
      return
    end if
    res = vf2_core(g1, g2, .true., color1, color2)
  end function vf2_isomorphic

  function vf2_subisomorphic(pattern, target, color_pattern, color_target, induced) result(res)
    type(graph_t), intent(in) :: pattern, target
    integer, intent(in), optional :: color_pattern(:), color_target(:)
    logical, intent(in), optional :: induced
    type(isomorphism_result_t) :: res
    logical :: ind

    ind = .false.
    if (present(induced)) ind = induced
    if (pattern%n > target%n .or. pattern%m > target%m .or. pattern%directed .neqv. target%directed) then
      allocate(res%map12(pattern%n), source=0)
      allocate(res%map21(target%n), source=0)
      return
    end if
    res = vf2_core(pattern, target, ind, color_pattern, color_target)
  end function vf2_subisomorphic

  function vf2_core(g1, g2, induced, color1, color2) result(res)
    type(graph_t), intent(in) :: g1, g2
    logical, intent(in) :: induced
    integer, intent(in), optional :: color1(:), color2(:)
    type(isomorphism_result_t) :: res
    integer, allocatable :: order(:), map12(:), map21(:), d1(:), d2(:), di1(:), di2(:)
    integer :: i, j, tmp
    logical :: ok

    if (present(color1)) then
      if (size(color1) /= g1%n) error stop 'vf2: color1 length mismatch'
    end if
    if (present(color2)) then
      if (size(color2) /= g2%n) error stop 'vf2: color2 length mismatch'
    end if
    if (present(color1) .neqv. present(color2)) error stop 'vf2: both color arrays are required'

    allocate(order(g1%n), map12(g1%n), map21(g2%n))
    order = [(i, i=1,g1%n)]
    map12 = 0
    map21 = 0
    d1 = degree(g1, 'out', .true.)
    d2 = degree(g2, 'out', .true.)
    if (g1%directed) then
      di1 = degree(g1, 'in', .true.)
      di2 = degree(g2, 'in', .true.)
    else
      allocate(di1(g1%n), source=0)
      allocate(di2(g2%n), source=0)
    end if

    ! Highest degree first is a simple VF2-compatible ordering heuristic.
    do i = 2, g1%n
      tmp = order(i)
      j = i - 1
      do while (j >= 1)
        if (vertex_score(order(j)) >= vertex_score(tmp)) exit
        order(j+1) = order(j)
        j = j - 1
      end do
      order(j+1) = tmp
    end do

    ok = search(1)
    res%isomorphic = ok
    allocate(res%map12(g1%n), res%map21(g2%n))
    if (ok) then
      res%map12 = map12
      res%map21 = map21
    else
      res%map12 = 0
      res%map21 = 0
    end if

  contains

    integer function vertex_score(v) result(s)
      integer, intent(in) :: v
      s = d1(v) + di1(v)
    end function vertex_score

    recursive logical function search(depth) result(found)
      integer, intent(in) :: depth
      integer :: u, v
      found = .false.
      res%states = res%states + 1_i8
      if (depth > g1%n) then
        found = .true.
        return
      end if
      u = order(depth)
      do v = 1, g2%n
        if (map21(v) /= 0) cycle
        if (.not. candidate_ok(u,v)) cycle
        map12(u) = v
        map21(v) = u
        if (search(depth+1)) then
          found = .true.
          return
        end if
        map12(u) = 0
        map21(v) = 0
      end do
    end function search

    logical function candidate_ok(u,v) result(okc)
      integer, intent(in) :: u, v
      integer :: x, y
      logical :: a12, a21, b12, b21
      okc = .false.

      if (induced) then
        if (d1(u) /= d2(v)) return
        if (g1%directed) then
          if (di1(u) /= di2(v)) return
        end if
      else
        if (d1(u) > d2(v)) return
        if (g1%directed) then
          if (di1(u) > di2(v)) return
        end if
      end if

      if (present(color1)) then
        if (color1(u) /= color2(v)) return
      end if

      ! Loops are part of the adjacency relation.
      a12 = are_adjacent(g1,u,u)
      b12 = are_adjacent(g2,v,v)
      if (induced) then
        if (a12 .neqv. b12) return
      else
        if (a12 .and. .not. b12) return
      end if

      do x = 1, g1%n
        y = map12(x)
        if (y == 0) cycle
        a12 = are_adjacent(g1,u,x)
        b12 = are_adjacent(g2,v,y)
        if (induced) then
          if (a12 .neqv. b12) return
        else
          if (a12 .and. .not. b12) return
        end if
        if (g1%directed) then
          a21 = are_adjacent(g1,x,u)
          b21 = are_adjacent(g2,y,v)
          if (induced) then
            if (a21 .neqv. b21) return
          else
            if (a21 .and. .not. b21) return
          end if
        end if
      end do
      okc = .true.
    end function candidate_ok
  end function vf2_core

end module igraph_isomorphism
