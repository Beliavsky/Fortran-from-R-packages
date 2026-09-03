module grbase_igraph
  use igraph_graph, only : graph_t, make_graph
  use igraph_cliques, only : vertex_set_list_t, maximal_cliques
  use igraph_components, only : components_t, weak_components
  use grbase_sets, only : append_set
  use grbase_types, only : set_list_t
  implicit none
  private

  public :: maximal_cliques_adjacency
  public :: connected_components_adjacency

contains

  function maximal_cliques_adjacency(adj) result(out)
    integer, intent(in) :: adj(:, :) !! Symmetric zero-one adjacency matrix whose maximal cliques are requested.
    type(set_list_t) :: out
    type(graph_t) :: g
    type(vertex_set_list_t) :: raw
    integer, allocatable :: edges(:, :)
    integer :: i
    integer :: j
    integer :: k
    integer :: m
    integer :: n

    allocate(out%set(0))
    if (.not. valid_undirected(adj)) return
    n = size(adj, 1)
    m = 0
    do i = 1, n - 1
      do j = i + 1, n
        if (adj(i, j) /= 0) m = m + 1
      end do
    end do
    allocate(edges(2, m))
    k = 0
    do i = 1, n - 1
      do j = i + 1, n
        if (adj(i, j) /= 0) then
          k = k + 1
          edges(:, k) = [i, j]
        end if
      end do
    end do
    g = make_graph(n, edges, directed=.false.)
    raw = maximal_cliques(g)
    do i = 1, raw%count
      call append_set(out, raw%set(i)%v)
    end do
  end function maximal_cliques_adjacency

  subroutine connected_components_adjacency(adj, membership, sizes, count_components)
    integer, intent(in) :: adj(:, :) !! Symmetric zero-one adjacency matrix whose connected components are requested.
    integer, allocatable, intent(out) :: membership(:) !! Component label from 1..count_components for each vertex.
    integer, allocatable, intent(out) :: sizes(:) !! Number of vertices in each connected component.
    integer, intent(out) :: count_components !! Number of weak/undirected connected components.
    type(components_t) :: comp
    type(graph_t) :: g
    integer, allocatable :: edges(:, :)
    integer :: i
    integer :: j
    integer :: k
    integer :: m
    integer :: n

    if (.not. valid_undirected(adj)) then
      allocate(membership(0), sizes(0))
      count_components = 0
      return
    end if
    n = size(adj, 1)
    m = 0
    do i = 1, n - 1
      do j = i + 1, n
        if (adj(i, j) /= 0) m = m + 1
      end do
    end do
    allocate(edges(2, m))
    k = 0
    do i = 1, n - 1
      do j = i + 1, n
        if (adj(i, j) /= 0) then
          k = k + 1
          edges(:, k) = [i, j]
        end if
      end do
    end do
    g = make_graph(n, edges, directed=.false.)
    comp = weak_components(g)
    membership = comp%membership
    sizes = comp%csize
    count_components = comp%count
  end subroutine connected_components_adjacency

  pure logical function valid_undirected(adj) result(ok)
    integer, intent(in) :: adj(:, :) !! Candidate zero-one symmetric adjacency matrix.
    integer :: i

    ok = size(adj, 1) == size(adj, 2)
    if (.not. ok) return
    ok = all(adj == 0 .or. adj == 1)
    if (.not. ok) return
    ok = all(adj == transpose(adj))
    if (.not. ok) return
    do i = 1, size(adj, 1)
      if (adj(i, i) /= 0) then
        ok = .false.
        return
      end if
    end do
  end function valid_undirected

end module grbase_igraph
