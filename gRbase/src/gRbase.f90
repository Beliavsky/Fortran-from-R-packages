module grbase
  use r_kinds, only : dp
  use grbase_types, only : integer_set_t, set_list_t, table_t, rip_order_t
  use grbase_arrays, only : cell_to_entry, entry_to_cell, make_plevels
  use grbase_arrays, only : next_cell, next_cell_slice, slice_to_entries
  use grbase_arrays, only : cell_to_entry_perm, perm_cell_entries
  use grbase_arrays, only : choose_integer, combinations, next_combination_mask
  use grbase_sets, only : unique_sorted, set_union, set_intersection, set_difference
  use grbase_sets, only : is_subset_of, contains_set, all_subsets
  use grbase_sets, only : maximal_sets, minimal_sets, all_pairs
  use grbase_tables, only : make_table, valid_table, table_permute, table_expand
  use grbase_tables, only : table_align, table_margin, table_add, table_subtract
  use grbase_tables, only : table_multiply, table_divide, table_divide_zero
  use grbase_tables, only : table_equal, table_normalize_all, table_normalize_first
  use grbase_tables, only : table_slice, table_sample_from_uniforms
  use grbase_tables, only : table_list_add, table_list_multiply
  use grbase_graphs, only : adjacency_from_edges, is_adjacency_matrix
  use grbase_graphs, only : is_symmetric_adjacency, topological_sort, is_dag
  use grbase_graphs, only : moralize_graph, maximum_cardinality_search, is_chordal
  use grbase_graphs, only : triangulate_elimination, triangulate_mcwh
  use grbase_graphs, only : minimal_triangulation, neighbors_of, parents_of
  use grbase_graphs, only : children_of, ancestors_of, descendants_of
  use grbase_graphs, only : is_complete_set, simplicial_nodes, separates_sets
  use grbase_graphs, only : junction_tree_from_cliques
  use grbase_igraph, only : maximal_cliques_adjacency, connected_components_adjacency
  use grbase_decompositions, only : rip_from_adjacency, maximal_prime_decomposition
  use grbase_stats, only : inverse_spd, concentration_to_partial_correlation
  use grbase_stats, only : covariance_to_partial_correlation
  use grbase_reductions, only : row_sums, column_sums, columnwise_product
  use grbase_reductions, only : matrix_nonzero_indices
  implicit none
  public
end module grbase
