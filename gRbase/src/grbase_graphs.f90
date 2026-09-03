module grbase_graphs
  use r_kinds, only : dp
  use grbase_sets, only : set_intersection, set_union, is_subset_of
  use grbase_types, only : integer_set_t, set_list_t, rip_order_t
  implicit none
  private

  public :: adjacency_from_edges
  public :: is_adjacency_matrix
  public :: is_symmetric_adjacency
  public :: topological_sort
  public :: is_dag
  public :: moralize_graph
  public :: maximum_cardinality_search
  public :: is_chordal
  public :: triangulate_elimination
  public :: triangulate_mcwh
  public :: minimal_triangulation
  public :: neighbors_of
  public :: parents_of
  public :: children_of
  public :: ancestors_of
  public :: descendants_of
  public :: is_complete_set
  public :: simplicial_nodes
  public :: separates_sets
  public :: junction_tree_from_cliques

contains

  pure function adjacency_from_edges(n, edges, directed) result(adj)
    integer, value :: n !! Number of graph vertices labeled 1..n.
    integer, intent(in) :: edges(:, :) !! Edge endpoints with shape (2,m), using one-based vertex labels.
    logical, value, optional :: directed !! True for directed edges; defaults to an undirected graph.
    integer, allocatable :: adj(:, :)
    logical :: is_directed
    integer :: e
    integer :: u
    integer :: v

    if (n < 0 .or. size(edges, 1) /= 2) then
      allocate(adj(0, 0))
      return
    end if
    allocate(adj(n, n))
    adj = 0
    is_directed = .false.
    if (present(directed)) is_directed = directed
    do e = 1, size(edges, 2)
      u = edges(1, e)
      v = edges(2, e)
      if (u < 1 .or. u > n .or. v < 1 .or. v > n) then
        deallocate(adj)
        allocate(adj(0, 0))
        return
      end if
      if (u == v) cycle
      adj(u, v) = 1
      if (.not. is_directed) adj(v, u) = 1
    end do
  end function adjacency_from_edges

  pure logical function is_adjacency_matrix(adj) result(ok)
    integer, intent(in) :: adj(:, :) !! Candidate zero-one adjacency matrix.
    integer :: i

    ok = size(adj, 1) == size(adj, 2)
    if (.not. ok) return
    ok = all(adj == 0 .or. adj == 1)
    if (.not. ok) return
    do i = 1, size(adj, 1)
      if (adj(i, i) /= 0) then
        ok = .false.
        return
      end if
    end do
  end function is_adjacency_matrix

  pure logical function is_symmetric_adjacency(adj) result(ok)
    integer, intent(in) :: adj(:, :) !! Candidate undirected zero-one adjacency matrix.

    ok = is_adjacency_matrix(adj)
    if (ok) ok = all(adj == transpose(adj))
  end function is_symmetric_adjacency

  pure subroutine topological_sort(adj, order, acyclic)
    integer, intent(in) :: adj(:, :) !! Directed adjacency matrix with row-to-column edge orientation.
    integer, allocatable, intent(out) :: order(:) !! One-based topological vertex order; empty when the input is not square.
    logical, intent(out) :: acyclic !! True when all vertices can be ordered without encountering a directed cycle.
    integer, allocatable :: indegree(:)
    logical, allocatable :: used(:)
    integer :: i
    integer :: j
    integer :: k
    integer :: n

    n = size(adj, 1)
    if (.not. is_adjacency_matrix(adj)) then
      allocate(order(0))
      acyclic = .false.
      return
    end if
    allocate(order(n), indegree(n), used(n))
    indegree = sum(adj, dim=1)
    used = .false.
    k = 0
    do
      i = 0
      do j = 1, n
        if (.not. used(j) .and. indegree(j) == 0) then
          i = j
          exit
        end if
      end do
      if (i == 0) exit
      k = k + 1
      order(k) = i
      used(i) = .true.
      do j = 1, n
        if (adj(i, j) /= 0) indegree(j) = indegree(j) - 1
      end do
    end do
    acyclic = k == n
    if (.not. acyclic) order = 0
  end subroutine topological_sort

  pure logical function is_dag(adj) result(ok)
    integer, intent(in) :: adj(:, :) !! Directed adjacency matrix to test for acyclicity.
    integer, allocatable :: order(:)

    call topological_sort(adj, order, ok)
  end function is_dag

  pure function moralize_graph(adj) result(out)
    integer, intent(in) :: adj(:, :) !! Directed adjacency matrix whose parents are married and directions removed.
    integer, allocatable :: out(:, :)
    integer :: child
    integer :: i
    integer :: j
    integer :: n

    if (.not. is_adjacency_matrix(adj)) then
      allocate(out(0, 0))
      return
    end if
    n = size(adj, 1)
    allocate(out(n, n))
    out = 0
    do i = 1, n
      do j = i + 1, n
        if (adj(i, j) /= 0 .or. adj(j, i) /= 0) then
          out(i, j) = 1
          out(j, i) = 1
        end if
      end do
    end do
    do child = 1, n
      do i = 1, n - 1
        if (adj(i, child) == 0) cycle
        do j = i + 1, n
          if (adj(j, child) /= 0) then
            out(i, j) = 1
            out(j, i) = 1
          end if
        end do
      end do
    end do
  end function moralize_graph

  pure subroutine maximum_cardinality_search(adj, order, perfect, preference)
    integer, intent(in) :: adj(:, :) !! Symmetric zero-one adjacency matrix for maximum-cardinality search.
    integer, allocatable, intent(out) :: order(:) !! Vertex order in the order vertices become passive.
    logical, intent(out) :: perfect !! True when each selected vertex's previously selected neighbors form a clique.
    integer, intent(in), optional :: preference(:) !! Optional permutation of 1..n used to resolve score ties at each position.
    integer, allocatable :: pref(:)
    integer, allocatable :: selected_neighbors(:)
    logical, allocatable :: selected(:)
    integer :: best
    integer :: count_selected
    integer :: i
    integer :: j
    integer :: k
    integer :: n
    integer :: score
    integer :: target
    integer :: v

    if (.not. is_symmetric_adjacency(adj)) then
      allocate(order(0))
      perfect = .false.
      return
    end if
    n = size(adj, 1)
    allocate(order(n), pref(n), selected(n), selected_neighbors(n))
    pref = [(i, i = 1, n)]
    if (present(preference)) then
      if (valid_permutation(preference, n)) pref = preference
    end if
    selected = .false.
    order = 0
    perfect = .true.
    do k = 1, n
      best = -1
      v = 0
      do i = 1, n
        if (selected(i)) cycle
        score = 0
        do j = 1, n
          if (selected(j) .and. adj(j, i) /= 0) score = score + 1
        end do
        if (score > best) then
          best = score
          v = i
        end if
      end do
      target = pref(k)
      if (.not. selected(target)) then
        score = 0
        do j = 1, n
          if (selected(j) .and. adj(j, target) /= 0) score = score + 1
        end do
        if (score == best) v = target
      end if
      order(k) = v
      count_selected = 0
      do i = 1, n
        if (selected(i) .and. adj(i, v) /= 0) then
          count_selected = count_selected + 1
          selected_neighbors(count_selected) = i
        end if
      end do
      if (count_selected > 1) then
        if (.not. is_complete_set(adj, selected_neighbors(:count_selected))) perfect = .false.
      end if
      selected(v) = .true.
    end do
  end subroutine maximum_cardinality_search

  pure logical function is_chordal(adj) result(ok)
    integer, intent(in) :: adj(:, :) !! Symmetric adjacency matrix tested for chordality.
    integer, allocatable :: order(:)

    call maximum_cardinality_search(adj, order, ok)
  end function is_chordal

  pure function triangulate_elimination(adj, elimination_order) result(out)
    integer, intent(in) :: adj(:, :) !! Symmetric adjacency matrix to triangulate by vertex elimination.
    integer, intent(in) :: elimination_order(:) !! Permutation of 1..n specifying the elimination order.
    integer, allocatable :: out(:, :)
    logical, allocatable :: active(:)
    integer, allocatable :: nbr(:)
    integer :: i
    integer :: j
    integer :: k
    integer :: n
    integer :: nn
    integer :: v

    if (.not. is_symmetric_adjacency(adj)) then
      allocate(out(0, 0))
      return
    end if
    n = size(adj, 1)
    if (.not. valid_permutation(elimination_order, n)) then
      allocate(out(0, 0))
      return
    end if
    out = adj
    allocate(active(n), nbr(n))
    active = .true.
    do k = 1, n
      v = elimination_order(k)
      nn = 0
      do i = 1, n
        if (active(i) .and. i /= v .and. out(i, v) /= 0) then
          nn = nn + 1
          nbr(nn) = i
        end if
      end do
      do i = 1, nn - 1
        do j = i + 1, nn
          out(nbr(i), nbr(j)) = 1
          out(nbr(j), nbr(i)) = 1
        end do
      end do
      active(v) = .false.
    end do
  end function triangulate_elimination

  pure function triangulate_mcwh(adj, nlevels) result(out)
    integer, intent(in) :: adj(:, :) !! Symmetric adjacency matrix triangulated using gRbase's minimum-closure-weight heuristic.
    real(dp), intent(in) :: nlevels(:) !! Positive per-vertex state-space weights used by the upstream heuristic.
    integer, allocatable :: out(:, :)
    logical, allocatable :: active(:)
    integer, allocatable :: nbr(:)
    integer :: i
    integer :: j
    integer :: n
    integer :: nn
    integer :: remaining
    integer :: v
    real(dp) :: best_weight
    real(dp) :: weight

    if (.not. is_symmetric_adjacency(adj) .or. size(nlevels) /= size(adj, 1)) then
      allocate(out(0, 0))
      return
    end if
    n = size(adj, 1)
    out = adj
    allocate(active(n), nbr(n))
    active = .true.
    remaining = n
    do while (remaining > 0)
      v = 0
      best_weight = huge(best_weight)
      do i = 1, n
        if (.not. active(i)) cycle
        weight = nlevels(i)
        do j = 1, n
          if (active(j) .and. out(j, i) /= 0) weight = weight + nlevels(j)
        end do
        if (weight < best_weight) then
          best_weight = weight
          v = i
        end if
      end do
      nn = 0
      do i = 1, n
        if (active(i) .and. i /= v .and. out(i, v) /= 0) then
          nn = nn + 1
          nbr(nn) = i
        end if
      end do
      do i = 1, nn - 1
        do j = i + 1, nn
          out(nbr(i), nbr(j)) = 1
          out(nbr(j), nbr(i)) = 1
        end do
      end do
      active(v) = .false.
      remaining = remaining - 1
    end do
  end function triangulate_mcwh

  pure function minimal_triangulation(adj, nlevels) result(out)
    integer, intent(in) :: adj(:, :) !! Symmetric adjacency matrix whose inclusion-minimal triangulation is requested.
    real(dp), intent(in) :: nlevels(:) !! Per-vertex heuristic weights passed to the initial MCWH triangulation.
    integer, allocatable :: out(:, :)
    integer, allocatable :: common(:)
    integer :: i
    integer :: j
    integer :: k
    integer :: n
    integer :: nc
    logical :: changed

    out = triangulate_mcwh(adj, nlevels)
    if (size(out, 1) == 0) return
    n = size(out, 1)
    allocate(common(n))
    do
      changed = .false.
      do i = 1, n - 1
        do j = i + 1, n
          if (adj(i, j) /= 0 .or. out(i, j) == 0) cycle
          nc = 0
          do k = 1, n
            if (k == i .or. k == j) cycle
            if (out(i, k) /= 0 .and. out(j, k) /= 0) then
              nc = nc + 1
              common(nc) = k
            end if
          end do
          if (is_complete_set(out, common(:nc))) then
            out(i, j) = 0
            out(j, i) = 0
            changed = .true.
          end if
        end do
      end do
      if (.not. changed) exit
    end do
  end function minimal_triangulation

  pure function neighbors_of(adj, vertices) result(out)
    integer, intent(in) :: adj(:, :) !! Graph adjacency matrix; an edge in either direction counts as adjacency.
    integer, intent(in) :: vertices(:) !! Vertex labels whose external neighborhood is requested.
    integer, allocatable :: out(:)
    logical, allocatable :: mark(:)
    integer :: i
    integer :: j
    integer :: n
    integer :: count_out

    if (.not. is_adjacency_matrix(adj)) then
      allocate(out(0))
      return
    end if
    n = size(adj, 1)
    allocate(mark(n))
    mark = .false.
    do i = 1, size(vertices)
      if (vertices(i) < 1 .or. vertices(i) > n) cycle
      do j = 1, n
        if (adj(vertices(i), j) /= 0 .or. adj(j, vertices(i)) /= 0) mark(j) = .true.
      end do
    end do
    do i = 1, size(vertices)
      if (vertices(i) >= 1 .and. vertices(i) <= n) mark(vertices(i)) = .false.
    end do
    count_out = count(mark)
    allocate(out(count_out))
    count_out = 0
    do i = 1, n
      if (mark(i)) then
        count_out = count_out + 1
        out(count_out) = i
      end if
    end do
  end function neighbors_of

  pure function parents_of(adj, vertices) result(out)
    integer, intent(in) :: adj(:, :) !! Directed adjacency matrix with parent-to-child orientation.
    integer, intent(in) :: vertices(:) !! Vertex labels whose external parents are requested.
    integer, allocatable :: out(:)
    logical, allocatable :: mark(:)
    integer :: i
    integer :: j
    integer :: n
    integer :: count_out

    if (.not. is_adjacency_matrix(adj)) then
      allocate(out(0))
      return
    end if
    n = size(adj, 1)
    allocate(mark(n))
    mark = .false.
    do i = 1, size(vertices)
      if (vertices(i) < 1 .or. vertices(i) > n) cycle
      do j = 1, n
        if (adj(j, vertices(i)) /= 0) mark(j) = .true.
      end do
    end do
    do i = 1, size(vertices)
      if (vertices(i) >= 1 .and. vertices(i) <= n) mark(vertices(i)) = .false.
    end do
    count_out = count(mark)
    allocate(out(count_out))
    count_out = 0
    do i = 1, n
      if (mark(i)) then
        count_out = count_out + 1
        out(count_out) = i
      end if
    end do
  end function parents_of

  pure function children_of(adj, vertices) result(out)
    integer, intent(in) :: adj(:, :) !! Directed adjacency matrix with parent-to-child orientation.
    integer, intent(in) :: vertices(:) !! Vertex labels whose external children are requested.
    integer, allocatable :: out(:)

    out = parents_of(transpose(adj), vertices)
  end function children_of

  pure function ancestors_of(adj, vertices) result(out)
    integer, intent(in) :: adj(:, :) !! Directed adjacency matrix with parent-to-child orientation.
    integer, intent(in) :: vertices(:) !! Starting vertices whose proper ancestors are requested.
    integer, allocatable :: out(:)

    out = closure_directed(adj, vertices, .true.)
  end function ancestors_of

  pure function descendants_of(adj, vertices) result(out)
    integer, intent(in) :: adj(:, :) !! Directed adjacency matrix with parent-to-child orientation.
    integer, intent(in) :: vertices(:) !! Starting vertices whose proper descendants are requested.
    integer, allocatable :: out(:)

    out = closure_directed(adj, vertices, .false.)
  end function descendants_of

  pure logical function is_complete_set(adj, vertices) result(ok)
    integer, intent(in) :: adj(:, :) !! Symmetric adjacency matrix containing the candidate clique.
    integer, intent(in) :: vertices(:) !! Vertex labels tested for pairwise adjacency.
    integer :: i
    integer :: j

    ok = is_symmetric_adjacency(adj)
    if (.not. ok) return
    do i = 1, size(vertices)
      if (vertices(i) < 1 .or. vertices(i) > size(adj, 1)) then
        ok = .false.
        return
      end if
    end do
    do i = 1, size(vertices) - 1
      do j = i + 1, size(vertices)
        if (adj(vertices(i), vertices(j)) == 0) then
          ok = .false.
          return
        end if
      end do
    end do
  end function is_complete_set

  pure function simplicial_nodes(adj) result(out)
    integer, intent(in) :: adj(:, :) !! Symmetric adjacency matrix whose simplicial vertices are requested.
    integer, allocatable :: out(:)
    integer, allocatable :: nbr(:)
    integer, allocatable :: tmp(:)
    integer :: i
    integer :: j
    integer :: n
    integer :: nn
    integer :: count_out

    if (.not. is_symmetric_adjacency(adj)) then
      allocate(out(0))
      return
    end if
    n = size(adj, 1)
    allocate(nbr(n), tmp(n))
    count_out = 0
    do i = 1, n
      nn = 0
      do j = 1, n
        if (adj(i, j) /= 0) then
          nn = nn + 1
          nbr(nn) = j
        end if
      end do
      if (is_complete_set(adj, nbr(:nn))) then
        count_out = count_out + 1
        tmp(count_out) = i
      end if
    end do
    allocate(out(count_out))
    if (count_out > 0) out = tmp(:count_out)
  end function simplicial_nodes

  pure logical function separates_sets(adj, a, b, separator) result(separates)
    integer, intent(in) :: adj(:, :) !! Symmetric adjacency matrix in which graph separation is tested.
    integer, intent(in) :: a(:) !! First vertex set.
    integer, intent(in) :: b(:) !! Second vertex set.
    integer, intent(in) :: separator(:) !! Vertices removed before testing reachability from `a` to `b`.
    logical, allocatable :: blocked(:)
    logical, allocatable :: seen(:)
    integer, allocatable :: queue(:)
    integer :: head
    integer :: i
    integer :: n
    integer :: tail
    integer :: u
    integer :: v

    separates = .false.
    if (.not. is_symmetric_adjacency(adj)) return
    n = size(adj, 1)
    allocate(blocked(n), seen(n), queue(n))
    blocked = .false.
    seen = .false.
    do i = 1, size(separator)
      if (separator(i) >= 1 .and. separator(i) <= n) blocked(separator(i)) = .true.
    end do
    head = 1
    tail = 0
    do i = 1, size(a)
      if (a(i) >= 1 .and. a(i) <= n) then
        if (.not. blocked(a(i)) .and. .not. seen(a(i))) then
          tail = tail + 1
          queue(tail) = a(i)
          seen(a(i)) = .true.
        end if
      end if
    end do
    do while (head <= tail)
      u = queue(head)
      head = head + 1
      do i = 1, size(b)
        if (b(i) == u) return
      end do
      do v = 1, n
        if (blocked(v) .or. seen(v)) cycle
        if (adj(u, v) /= 0) then
          seen(v) = .true.
          tail = tail + 1
          queue(tail) = v
        end if
      end do
    end do
    separates = .true.
  end function separates_sets

  pure function junction_tree_from_cliques(mcs_order, cliques) result(tree)
    integer, intent(in) :: mcs_order(:) !! Perfect maximum-cardinality vertex order used to rank the cliques.
    type(set_list_t), intent(in) :: cliques !! Maximal cliques of a chordal graph, expressed with one-based vertex labels.
    type(rip_order_t) :: tree
    integer, allocatable :: maxpos(:)
    integer, allocatable :: pos(:)
    integer, allocatable :: ord(:)
    integer, allocatable :: sep(:)
    integer, allocatable :: past(:)
    integer :: i
    integer :: j
    integer :: k
    integer :: n
    integer :: nc
    integer :: tmp

    n = size(mcs_order)
    nc = cliques%count
    allocate(tree%cliques%set(0), tree%separators%set(0))
    allocate(tree%parent(nc), tree%child(nc), tree%host(n), tree%order(nc))
    tree%parent = 0
    tree%child = 0
    tree%host = 0
    if (nc == 0) return
    allocate(pos(n), maxpos(nc), ord(nc))
    pos = 0
    do i = 1, n
      if (mcs_order(i) >= 1 .and. mcs_order(i) <= n) pos(mcs_order(i)) = i
    end do
    do i = 1, nc
      maxpos(i) = 0
      do k = 1, size(cliques%set(i)%value)
        j = cliques%set(i)%value(k)
        if (j >= 1 .and. j <= n) maxpos(i) = max(maxpos(i), pos(j))
      end do
      ord(i) = i
    end do
    do i = 2, nc
      tmp = ord(i)
      j = i - 1
      do while (j >= 1)
        if (maxpos(ord(j)) <= maxpos(tmp)) exit
        ord(j + 1) = ord(j)
        j = j - 1
      end do
      ord(j + 1) = tmp
    end do
    tree%order = ord
    do i = 1, nc
      call append_set_local(tree%cliques, cliques%set(ord(i))%value)
      call append_set_local(tree%separators, [integer ::])
      do k = 1, size(tree%cliques%set(i)%value)
        j = tree%cliques%set(i)%value(k)
        if (j >= 1 .and. j <= n) tree%host(j) = i
      end do
    end do
    past = tree%cliques%set(1)%value
    do i = 2, nc
      sep = set_intersection(past, tree%cliques%set(i)%value)
      tree%separators%set(i)%value = sep
      if (size(sep) > 0) then
        do j = i - 1, 1, -1
          if (is_subset_of(sep, tree%cliques%set(j)%value)) then
            tree%parent(i) = j
            tree%child(j) = i
            exit
          end if
        end do
      end if
      past = set_union(past, tree%cliques%set(i)%value)
    end do
  end function junction_tree_from_cliques

  pure function closure_directed(adj, vertices, reverse) result(out)
    integer, intent(in) :: adj(:, :) !! Directed adjacency matrix searched for transitive reachability.
    integer, intent(in) :: vertices(:) !! Starting vertices excluded from the returned proper closure.
    logical, value :: reverse !! True follows incoming edges for ancestors; false follows outgoing edges for descendants.
    integer, allocatable :: out(:)
    logical, allocatable :: start(:)
    logical, allocatable :: seen(:)
    integer, allocatable :: queue(:)
    integer :: head
    integer :: i
    integer :: n
    integer :: tail
    integer :: u
    integer :: v
    integer :: count_out

    if (.not. is_adjacency_matrix(adj)) then
      allocate(out(0))
      return
    end if
    n = size(adj, 1)
    allocate(start(n), seen(n), queue(n))
    start = .false.
    seen = .false.
    head = 1
    tail = 0
    do i = 1, size(vertices)
      if (vertices(i) < 1 .or. vertices(i) > n) cycle
      if (.not. seen(vertices(i))) then
        tail = tail + 1
        queue(tail) = vertices(i)
        seen(vertices(i)) = .true.
        start(vertices(i)) = .true.
      end if
    end do
    do while (head <= tail)
      u = queue(head)
      head = head + 1
      do v = 1, n
        if (seen(v)) cycle
        if (reverse) then
          if (adj(v, u) == 0) cycle
        else
          if (adj(u, v) == 0) cycle
        end if
        seen(v) = .true.
        tail = tail + 1
        queue(tail) = v
      end do
    end do
    count_out = count(seen .and. .not. start)
    allocate(out(count_out))
    count_out = 0
    do i = 1, n
      if (seen(i) .and. .not. start(i)) then
        count_out = count_out + 1
        out(count_out) = i
      end if
    end do
  end function closure_directed

  pure logical function valid_permutation(perm, n) result(ok)
    integer, intent(in) :: perm(:) !! Candidate permutation entries.
    integer, value :: n !! Required length and largest legal entry.
    integer :: i

    ok = size(perm) == n
    if (.not. ok) return
    do i = 1, n
      if (count(perm == i) /= 1) then
        ok = .false.
        return
      end if
    end do
  end function valid_permutation

  pure subroutine append_set_local(list, values)
    type(set_list_t), intent(inout) :: list !! Set list receiving one appended set without changing the set's ordering.
    integer, intent(in) :: values(:) !! Vertex labels stored in the appended set.
    type(integer_set_t), allocatable :: tmp_set(:)
    integer :: n

    n = list%count
    allocate(tmp_set(n + 1))
    if (n > 0) tmp_set(:n) = list%set
    tmp_set(n + 1)%value = values
    call move_alloc(tmp_set, list%set)
    list%count = n + 1
  end subroutine append_set_local

end module grbase_graphs
