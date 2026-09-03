module grbase_decompositions
  use grbase_graphs, only : is_symmetric_adjacency, maximum_cardinality_search
  use grbase_graphs, only : junction_tree_from_cliques, is_complete_set
  use grbase_igraph, only : maximal_cliques_adjacency
  use grbase_sets, only : append_set, set_intersection, set_union
  use grbase_types, only : set_list_t, rip_order_t
  implicit none
  private

  public :: rip_from_adjacency
  public :: maximal_prime_decomposition

contains

  function rip_from_adjacency(adj) result(tree)
    integer, intent(in) :: adj(:, :) !! Chordal undirected adjacency matrix whose RIP/junction-tree representation is requested.
    type(rip_order_t) :: tree
    type(set_list_t) :: cliques
    integer, allocatable :: order(:)
    logical :: perfect

    if (.not. is_symmetric_adjacency(adj)) then
      call empty_tree(tree, size(adj, 1))
      return
    end if
    call maximum_cardinality_search(adj, order, perfect)
    if (.not. perfect) then
      call empty_tree(tree, size(adj, 1))
      return
    end if
    cliques = maximal_cliques_adjacency(adj)
    tree = junction_tree_from_cliques(order, cliques)
  end function rip_from_adjacency

  function maximal_prime_decomposition(original, triangulated) result(out)
    integer, intent(in) :: original(:, :) !! Original undirected graph whose maximal prime subgraphs are requested.
    integer, intent(in) :: triangulated(:, :) !! Inclusion-minimal chordal supergraph used to construct the initial junction tree.
    type(set_list_t) :: out
    type(rip_order_t) :: tree
    type(set_list_t) :: work
    integer, allocatable :: parent(:)
    integer, allocatable :: sep(:)
    logical, allocatable :: alive(:)
    integer :: i
    integer :: j
    integer :: p

    allocate(out%set(0))
    if (.not. is_symmetric_adjacency(original)) return
    if (.not. is_symmetric_adjacency(triangulated)) return
    if (any((original /= 0) .and. (triangulated == 0))) return
    tree = rip_from_adjacency(triangulated)
    if (tree%cliques%count == 0) return
    work = tree%cliques
    parent = tree%parent
    allocate(alive(work%count))
    alive = .true.
    do i = work%count, 2, -1
      if (.not. alive(i)) cycle
      p = parent(i)
      if (p <= 0) cycle
      sep = set_intersection(work%set(i)%value, work%set(p)%value)
      if (.not. is_complete_set(original, sep)) then
        work%set(p)%value = set_union(work%set(p)%value, work%set(i)%value)
        alive(i) = .false.
        do j = i + 1, work%count
          if (alive(j) .and. parent(j) == i) parent(j) = p
        end do
      end if
    end do
    do i = 1, work%count
      if (alive(i)) call append_set(out, work%set(i)%value)
    end do
  end function maximal_prime_decomposition

  pure subroutine empty_tree(tree, nvertices)
    type(rip_order_t), intent(out) :: tree !! Tree object initialized to contain no cliques or separators.
    integer, value :: nvertices !! Number of graph vertices, used to size the empty host map.

    allocate(tree%cliques%set(0), tree%separators%set(0))
    allocate(tree%parent(0), tree%child(0), tree%host(max(0, nvertices)), tree%order(0))
    tree%cliques%count = 0
    tree%separators%count = 0
    tree%host = 0
  end subroutine empty_tree

end module grbase_decompositions
