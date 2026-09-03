program test_grbase
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use grbase
  implicit none

  integer :: failures

  failures = 0
  call test_arrays(failures)
  call test_sets(failures)
  call test_tables(failures)
  call test_reductions(failures)
  call test_graphs(failures)
  call test_igraph_bridge(failures)
  call test_decompositions(failures)
  call test_statistics(failures)

  if (failures /= 0) then
    write (*, '(a,i0)') 'gRbase tests failed: ', failures
    error stop 1
  end if
  write (*, '(a)') 'All gRbase deterministic tests passed.'

contains

  subroutine test_arrays(failures)
    integer, intent(inout) :: failures !! Running count of failed deterministic checks.
    integer, allocatable :: cell(:)
    integer, allocatable :: cmb(:, :)
    integer, allocatable :: entries(:)
    integer, allocatable :: map(:)
    integer :: mask(4)
    logical :: has_next

    call assert_integer(cell_to_entry([2, 2], [2, 3]), 4, 'cell_to_entry', failures)
    cell = entry_to_cell(6, [2, 3])
    call assert_integer_vector(cell, [2, 3], 'entry_to_cell', failures)
    call assert_integer_vector(make_plevels([2, 3, 4]), [1, 2, 6], 'make_plevels', failures)
    entries = slice_to_entries([2], [1], [2, 3])
    call assert_integer_vector(entries, [2, 4, 6], 'slice_to_entries', failures)
    map = perm_cell_entries([2, 3], [2, 1])
    call assert_integer_vector(map, [1, 3, 5, 2, 4, 6], 'perm_cell_entries', failures)
    call assert_integer(choose_integer(5, 2), 10, 'choose_integer', failures)
    cmb = combinations(4, 2)
    call assert_integer(size(cmb, 2), 6, 'combinations count', failures)
    call assert_integer_vector(cmb(:, 1), [1, 2], 'combinations first', failures)
    call assert_integer_vector(cmb(:, 6), [3, 4], 'combinations last', failures)
    mask = [1, 0, 0, 1]
    call next_combination_mask(mask, has_next)
    call assert_logical(has_next, .true., 'next_combination_mask flag', failures)
    call assert_integer_vector(mask, [0, 1, 1, 0], 'next_combination_mask value', failures)
  end subroutine test_arrays

  subroutine test_sets(failures)
    integer, intent(inout) :: failures !! Running count of failed deterministic checks.
    type(set_list_t) :: sets
    type(set_list_t) :: mx
    type(set_list_t) :: mn
    integer, allocatable :: pairs(:, :)

    call assert_integer_vector(unique_sorted([3, 1, 2, 1]), [1, 2, 3], 'unique_sorted', failures)
    call assert_integer_vector(set_union([3, 1], [2, 3]), [1, 2, 3], 'set_union', failures)
    call assert_integer_vector(set_intersection([1, 3, 4], [2, 3, 4]), [3, 4], 'set_intersection', failures)
    call assert_integer_vector(set_difference([1, 2, 3], [2, 4]), [1, 3], 'set_difference', failures)
    call assert_logical(is_subset_of([1, 3], [1, 2, 3]), .true., 'is_subset_of', failures)
    sets = all_subsets([1, 2, 3])
    call assert_integer(sets%count, 8, 'all_subsets count', failures)

    if (allocated(sets%set)) deallocate(sets%set)
    allocate(sets%set(0))
    sets%count = 0
    call append_set_for_test(sets, [1])
    call append_set_for_test(sets, [1, 2])
    call append_set_for_test(sets, [2, 3])
    mx = maximal_sets(sets)
    mn = minimal_sets(sets)
    call assert_integer(mx%count, 2, 'maximal_sets count', failures)
    call assert_integer(mn%count, 2, 'minimal_sets count', failures)
    pairs = all_pairs([1, 2, 3])
    call assert_integer(size(pairs, 2), 3, 'all_pairs count', failures)
    call assert_integer_vector(pairs(:, 1), [1, 2], 'all_pairs first', failures)
  end subroutine test_sets

  subroutine test_tables(failures)
    integer, intent(inout) :: failures !! Running count of failed deterministic checks.
    type(table_t) :: tab
    type(table_t) :: tab2
    type(table_t) :: out
    type(table_t) :: perm
    type(table_t) :: zero1
    type(table_t) :: zero2
    type(table_t) :: tabs(2)
    real(dp), parameter :: tol = 1.0e-12_dp

    tab = make_table([1, 2], [2, 3], [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp])
    call assert_logical(valid_table(tab), .true., 'valid_table', failures)
    out = table_margin(tab, [1])
    call assert_real_vector(out%value, [9.0_dp, 12.0_dp], tol, 'table_margin first', failures)
    out = table_margin(tab, [2])
    call assert_real_vector(out%value, [3.0_dp, 7.0_dp, 11.0_dp], tol, 'table_margin second', failures)

    perm = table_permute(tab, [2, 1])
    call assert_integer_vector(perm%var, [2, 1], 'table_permute vars', failures)
    call assert_integer_vector(perm%dim, [3, 2], 'table_permute dims', failures)
    call assert_real_vector(perm%value, [1.0_dp, 3.0_dp, 5.0_dp, 2.0_dp, 4.0_dp, 6.0_dp], &
      tol, 'table_permute values', failures)
    call assert_logical(table_equal(tab, perm), .true., 'table_equal reordered', failures)

    out = table_expand(tab, [2, 3], [3, 2], 0)
    call assert_integer_vector(out%var, [1, 2, 3], 'table_expand vars', failures)
    call assert_real_vector(out%value, &
      [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp], &
      tol, 'table_expand replicate', failures)
    out = table_expand(tab, [2, 3], [3, 2], 1)
    call assert_real_vector(out%value(:6), 0.5_dp * [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp], &
      tol, 'table_expand weighted', failures)
    out = table_expand(tab, [2, 3], [3, 2], 2)
    call assert_real_vector(out%value(7:), [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
      tol, 'table_expand zero fill', failures)

    tab2 = make_table([2], [3], [10.0_dp, 20.0_dp, 30.0_dp])
    tabs = [tab, tab2]
    out = table_list_multiply(tabs)
    call assert_real_vector(out%value, [10.0_dp, 20.0_dp, 60.0_dp, 80.0_dp, 150.0_dp, 180.0_dp], &
      tol, 'table_list_multiply', failures)
    out = table_list_add(tabs)
    call assert_real_vector(out%value, [11.0_dp, 12.0_dp, 23.0_dp, 24.0_dp, 35.0_dp, 36.0_dp], &
      tol, 'table_list_add', failures)

    out = table_multiply(tab, tab2)
    call assert_real_vector(out%value, [10.0_dp, 20.0_dp, 60.0_dp, 80.0_dp, 150.0_dp, 180.0_dp], &
      tol, 'table_multiply', failures)
    out = table_add(tab, tab2)
    call assert_real_vector(out%value, [11.0_dp, 12.0_dp, 23.0_dp, 24.0_dp, 35.0_dp, 36.0_dp], &
      tol, 'table_add', failures)

    out = table_normalize_first(make_table([1, 2], [2, 2], [1.0_dp, 3.0_dp, 2.0_dp, 2.0_dp]))
    call assert_real_vector(out%value, [0.25_dp, 0.75_dp, 0.5_dp, 0.5_dp], tol, 'table_normalize_first', failures)
    out = table_normalize_all(make_table([1], [4], [1.0_dp, 1.0_dp, 2.0_dp, 4.0_dp]))
    call assert_real_vector(out%value, [0.125_dp, 0.125_dp, 0.25_dp, 0.5_dp], tol, 'table_normalize_all', failures)
    out = table_slice(tab, [2], [2])
    call assert_integer_vector(out%var, [1], 'table_slice vars', failures)
    call assert_real_vector(out%value, [3.0_dp, 4.0_dp], tol, 'table_slice values', failures)

    out = table_align(perm, tab)
    call assert_integer_vector(out%var, [1, 2], 'table_align vars', failures)
    call assert_real_vector(out%value, tab%value, tol, 'table_align values', failures)
    call assert_integer_matrix(table_sample_from_uniforms(make_table([1, 2], [2, 2], &
      [1.0_dp, 1.0_dp, 2.0_dp, 4.0_dp]), [0.0_dp, 0.2_dp, 0.3_dp, 0.9_dp]), &
      reshape([1, 2, 1, 2, 1, 1, 2, 2], [4, 2]), 'table_sample_from_uniforms', failures)

    zero1 = make_table([1], [2], [0.0_dp, 1.0_dp])
    zero2 = make_table([1], [2], [0.0_dp, 0.0_dp])
    out = table_divide_zero(zero1, zero2)
    call assert_real_vector(out%value, [0.0_dp, 0.0_dp], tol, 'table_divide_zero', failures)
    call assert_logical(all(ieee_is_finite(out%value)), .true., 'table_divide_zero finite', failures)
  end subroutine test_tables

  subroutine test_reductions(failures)
    integer, intent(inout) :: failures !! Running count of failed deterministic checks.
    real(dp) :: x(2, 3)
    real(dp) :: y(2, 4)
    real(dp), parameter :: tol = 1.0e-12_dp

    x(1, :) = [1.0_dp, 2.0_dp, 3.0_dp]
    x(2, :) = [4.0_dp, 5.0_dp, 6.0_dp]
    call assert_real_vector(row_sums(x), [6.0_dp, 15.0_dp], tol, 'row_sums', failures)
    call assert_real_vector(column_sums(x), [5.0_dp, 7.0_dp, 9.0_dp], tol, 'column_sums', failures)

    y(:, 1) = [1.0_dp, 2.0_dp]
    y(:, 2) = [3.0_dp, 4.0_dp]
    y(:, 3) = [5.0_dp, 6.0_dp]
    y(:, 4) = [7.0_dp, 8.0_dp]
    call assert_real_matrix(columnwise_product([10.0_dp, -1.0_dp], y), &
      reshape([10.0_dp, 20.0_dp, -3.0_dp, -4.0_dp, 50.0_dp, 60.0_dp, -7.0_dp, -8.0_dp], [2, 4]), &
      tol, 'columnwise_product', failures)
    call assert_integer_matrix(matrix_nonzero_indices(reshape([0.0_dp, 2.0_dp, 3.0_dp, 0.0_dp], [2, 2])), &
      reshape([1, 2, 2, 1], [2, 2]), 'matrix_nonzero_indices', failures)
  end subroutine test_reductions

  subroutine test_graphs(failures)
    integer, intent(inout) :: failures !! Running count of failed deterministic checks.
    integer, allocatable :: dag(:, :)
    integer, allocatable :: moral(:, :)
    integer, allocatable :: order(:)
    integer, allocatable :: square(:, :)
    integer, allocatable :: tri(:, :)
    integer, allocatable :: vals(:)
    logical :: acyclic
    logical :: perfect
    real(dp) :: levels(4)

    dag = adjacency_from_edges(4, reshape([1, 2, 1, 3, 2, 4, 3, 4], [2, 4]), directed=.true.)
    call topological_sort(dag, order, acyclic)
    call assert_logical(acyclic, .true., 'topological_sort acyclic', failures)
    call assert_integer_vector(order, [1, 2, 3, 4], 'topological_sort order', failures)
    call assert_logical(is_dag(dag), .true., 'is_dag', failures)
    moral = moralize_graph(dag)
    call assert_integer(moral(2, 3), 1, 'moralize_graph co-parents', failures)

    square = adjacency_from_edges(4, reshape([1, 2, 3, 4, 2, 3, 4, 1], [2, 4]))
    call maximum_cardinality_search(square, order, perfect)
    call assert_logical(perfect, .false., 'MCS detects nonchordal cycle', failures)
    call assert_logical(is_chordal(square), .false., 'is_chordal cycle', failures)
    levels = 2.0_dp
    tri = triangulate_mcwh(square, levels)
    call assert_logical(is_chordal(tri), .true., 'triangulate_mcwh chordal', failures)
    call assert_integer((sum(tri) - sum(square)) / 2, 1, 'triangulate_mcwh fill count', failures)
    tri = minimal_triangulation(square, levels)
    call assert_logical(is_chordal(tri), .true., 'minimal_triangulation chordal', failures)
    call assert_integer((sum(tri) - sum(square)) / 2, 1, 'minimal_triangulation fill count', failures)

    tri = triangulate_elimination(square, [1, 2, 3, 4])
    call assert_logical(is_chordal(tri), .true., 'triangulate_elimination chordal', failures)
    vals = simplicial_nodes(tri)
    call assert_logical(size(vals) >= 2, .true., 'simplicial_nodes', failures)
    call assert_logical(is_complete_set(tri, [1, 2]), .true., 'is_complete_set', failures)
    call assert_logical(separates_sets(square, [1], [3], [2, 4]), .true., 'separates_sets', failures)
    call assert_integer_vector(parents_of(dag, [4]), [2, 3], 'parents_of', failures)
    call assert_integer_vector(children_of(dag, [1]), [2, 3], 'children_of', failures)
    call assert_integer_vector(ancestors_of(dag, [4]), [1, 2, 3], 'ancestors_of', failures)
    call assert_integer_vector(descendants_of(dag, [1]), [2, 3, 4], 'descendants_of', failures)
  end subroutine test_graphs

  subroutine test_igraph_bridge(failures)
    integer, intent(inout) :: failures !! Running count of failed deterministic checks.
    type(set_list_t) :: cliques
    integer, allocatable :: adj(:, :)
    integer, allocatable :: membership(:)
    integer, allocatable :: sizes(:)
    integer :: ncomp

    adj = adjacency_from_edges(4, reshape([1, 2, 2, 3], [2, 2]))
    call connected_components_adjacency(adj, membership, sizes, ncomp)
    call assert_integer(ncomp, 2, 'connected_components count', failures)
    call assert_integer_vector(sizes, [3, 1], 'connected_components sizes', failures)

    adj = adjacency_from_edges(4, reshape([1, 2, 1, 3, 2, 3, 2, 4, 3, 4], [2, 5]))
    cliques = maximal_cliques_adjacency(adj)
    call assert_integer(cliques%count, 2, 'maximal_cliques count', failures)
    call assert_logical(contains_set(cliques, [1, 2, 3]), .true., 'maximal_cliques first', failures)
    call assert_logical(contains_set(cliques, [2, 3, 4]), .true., 'maximal_cliques second', failures)
  end subroutine test_igraph_bridge

  subroutine test_decompositions(failures)
    integer, intent(inout) :: failures !! Running count of failed deterministic checks.
    type(rip_order_t) :: tree
    type(set_list_t) :: mpd
    integer, allocatable :: chordal(:, :)
    integer, allocatable :: cycle(:, :)
    integer, allocatable :: tri(:, :)
    real(dp) :: levels(4)

    chordal = adjacency_from_edges(4, reshape([1, 2, 1, 3, 2, 3, 2, 4, 3, 4], [2, 5]))
    tree = rip_from_adjacency(chordal)
    call assert_integer(tree%cliques%count, 2, 'RIP clique count', failures)
    call assert_integer_vector(tree%separators%set(2)%value, [2, 3], 'RIP separator', failures)
    call assert_integer(tree%parent(2), 1, 'RIP parent', failures)
    call assert_integer_vector(tree%host, [1, 2, 2, 2], 'RIP host', failures)

    cycle = adjacency_from_edges(4, reshape([1, 2, 3, 4, 2, 3, 4, 1], [2, 4]))
    levels = 2.0_dp
    tri = minimal_triangulation(cycle, levels)
    mpd = maximal_prime_decomposition(cycle, tri)
    call assert_integer(mpd%count, 1, 'MPD cycle count', failures)
    call assert_integer_vector(mpd%set(1)%value, [1, 2, 3, 4], 'MPD cycle vertices', failures)
  end subroutine test_decompositions

  subroutine test_statistics(failures)
    integer, intent(inout) :: failures !! Running count of failed deterministic checks.
    real(dp), allocatable :: inv(:, :)
    real(dp), allocatable :: pcor(:, :)
    real(dp) :: covariance(2, 2)
    real(dp) :: concentration(2, 2)
    integer :: info
    real(dp), parameter :: tol = 1.0e-12_dp

    covariance(1, :) = [4.0_dp, 2.0_dp]
    covariance(2, :) = [2.0_dp, 9.0_dp]
    call inverse_spd(covariance, inv, info)
    call assert_integer(info, 0, 'inverse_spd info', failures)
    call assert_real(inv(1, 1), 0.28125_dp, tol, 'inverse_spd 11', failures)
    call assert_real(inv(1, 2), -0.0625_dp, tol, 'inverse_spd 12', failures)
    call covariance_to_partial_correlation(covariance, pcor, info)
    call assert_integer(info, 0, 'covariance_to_partial_correlation info', failures)
    call assert_real(pcor(1, 2), 1.0_dp / 3.0_dp, tol, 'covariance partial correlation', failures)

    concentration(1, :) = [4.0_dp, -1.0_dp]
    concentration(2, :) = [-1.0_dp, 9.0_dp]
    call concentration_to_partial_correlation(concentration, pcor, info)
    call assert_integer(info, 0, 'concentration_to_partial_correlation info', failures)
    call assert_real(pcor(1, 2), 1.0_dp / 6.0_dp, tol, 'concentration partial correlation', failures)
    call assert_real(pcor(1, 1), 1.0_dp, tol, 'partial correlation diagonal', failures)
  end subroutine test_statistics

  subroutine append_set_for_test(list, values)
    type(set_list_t), intent(inout) :: list !! Test set list receiving one appended set.
    integer, intent(in) :: values(:) !! Integer values placed in the appended set.
    type(integer_set_t), allocatable :: tmp(:)
    integer :: n

    n = list%count
    allocate(tmp(n + 1))
    if (n > 0) tmp(:n) = list%set
    tmp(n + 1)%value = unique_sorted(values)
    call move_alloc(tmp, list%set)
    list%count = n + 1
  end subroutine append_set_for_test

  subroutine assert_integer(actual, expected, name, failures)
    integer, value :: actual !! Integer value produced by the routine under test.
    integer, value :: expected !! Exact expected integer value.
    character(len=*), intent(in) :: name !! Human-readable label identifying the deterministic check.
    integer, intent(inout) :: failures !! Running count incremented when the assertion fails.

    if (actual /= expected) then
      failures = failures + 1
      write (*, '(a,2(1x,i0))') 'FAIL '//trim(name)//':', actual, expected
    end if
  end subroutine assert_integer

  subroutine assert_integer_vector(actual, expected, name, failures)
    integer, intent(in) :: actual(:) !! Integer vector produced by the routine under test.
    integer, intent(in) :: expected(:) !! Exact expected integer vector.
    character(len=*), intent(in) :: name !! Human-readable label identifying the deterministic check.
    integer, intent(inout) :: failures !! Running count incremented when the assertion fails.

    if (size(actual) /= size(expected)) then
      failures = failures + 1
      write (*, '(a)') 'FAIL '//trim(name)//': size mismatch'
    else if (any(actual /= expected)) then
      failures = failures + 1
      write (*, '(a)') 'FAIL '//trim(name)//': value mismatch'
    end if
  end subroutine assert_integer_vector

  subroutine assert_real(actual, expected, tolerance, name, failures)
    real(dp), value :: actual !! Real value produced by the routine under test.
    real(dp), value :: expected !! Reference real value.
    real(dp), value :: tolerance !! Maximum allowed absolute difference.
    character(len=*), intent(in) :: name !! Human-readable label identifying the deterministic check.
    integer, intent(inout) :: failures !! Running count incremented when the assertion fails.

    if (abs(actual - expected) > tolerance) then
      failures = failures + 1
      write (*, '(a,2(1x,es16.8))') 'FAIL '//trim(name)//':', actual, expected
    end if
  end subroutine assert_real

  subroutine assert_real_vector(actual, expected, tolerance, name, failures)
    real(dp), intent(in) :: actual(:) !! Real vector produced by the routine under test.
    real(dp), intent(in) :: expected(:) !! Reference real vector.
    real(dp), value :: tolerance !! Maximum allowed elementwise absolute difference.
    character(len=*), intent(in) :: name !! Human-readable label identifying the deterministic check.
    integer, intent(inout) :: failures !! Running count incremented when the assertion fails.

    if (size(actual) /= size(expected)) then
      failures = failures + 1
      write (*, '(a)') 'FAIL '//trim(name)//': size mismatch'
    else if (any(abs(actual - expected) > tolerance)) then
      failures = failures + 1
      write (*, '(a)') 'FAIL '//trim(name)//': value mismatch'
    end if
  end subroutine assert_real_vector

  subroutine assert_integer_matrix(actual, expected, name, failures)
    integer, intent(in) :: actual(:, :) !! Integer matrix produced by the routine under test.
    integer, intent(in) :: expected(:, :) !! Exact expected integer matrix.
    character(len=*), intent(in) :: name !! Human-readable label identifying the deterministic check.
    integer, intent(inout) :: failures !! Running count incremented when the assertion fails.

    if (any(shape(actual) /= shape(expected))) then
      failures = failures + 1
      write (*, '(a)') 'FAIL '//trim(name)//': shape mismatch'
    else if (any(actual /= expected)) then
      failures = failures + 1
      write (*, '(a)') 'FAIL '//trim(name)//': value mismatch'
    end if
  end subroutine assert_integer_matrix

  subroutine assert_real_matrix(actual, expected, tolerance, name, failures)
    real(dp), intent(in) :: actual(:, :) !! Real matrix produced by the routine under test.
    real(dp), intent(in) :: expected(:, :) !! Reference real matrix.
    real(dp), value :: tolerance !! Maximum allowed elementwise absolute difference.
    character(len=*), intent(in) :: name !! Human-readable label identifying the deterministic check.
    integer, intent(inout) :: failures !! Running count incremented when the assertion fails.

    if (any(shape(actual) /= shape(expected))) then
      failures = failures + 1
      write (*, '(a)') 'FAIL '//trim(name)//': shape mismatch'
    else if (any(abs(actual - expected) > tolerance)) then
      failures = failures + 1
      write (*, '(a)') 'FAIL '//trim(name)//': value mismatch'
    end if
  end subroutine assert_real_matrix

  subroutine assert_logical(actual, expected, name, failures)
    logical, value :: actual !! Logical result produced by the routine under test.
    logical, value :: expected !! Expected logical result.
    character(len=*), intent(in) :: name !! Human-readable label identifying the deterministic check.
    integer, intent(inout) :: failures !! Running count incremented when the assertion fails.

    if (actual .neqv. expected) then
      failures = failures + 1
      write (*, '(a)') 'FAIL '//trim(name)
    end if
  end subroutine assert_logical

end program test_grbase
