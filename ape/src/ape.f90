! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
!
! Public umbrella module for the modern Fortran translation of ape 5.8-1.
module ape
   use r_kinds, only : dp
   use ape_types, only : phylo_tree, make_phylo_tree, parent_vector, child_counts, edge_index_to_child
   use ape_reconstruction, only : nj, bionj, mvr, njs, bionjs, mvrs, additive_completion, ultrametric_completion
   use ape_triangle_reconstruction, only : triang_mtd, triang_mtds
   use ape_fastme, only : fastme_ols, fastme_bal
   use ape_tree_algorithms, only : node_depth_edgelength, node_depth_count, dist_nodes, mrca, node_path, &
      descendant_tip_counts, balance_counts, cherry_count, branching_times, pic, ace_pic, chrono_mpl, compute_brtime
   use ape_statistics, only : moran_result, yule_result, gamma_stat, yule_fit, coalescent_intervals, &
      ltt_coordinates, moran_i, minimum_spanning_tree, delta_plot_statistics
   use ape_misc_statistics, only : tree_count_result, chi_square_result, diversification_gof_result, &
      diversification_time_result, howmanytrees, gene_distance_matrix, splits_compatible, all_splits_compatible, &
      diversification_gof, diversification_time, slowinski_guyer_test, mcconway_sims_test, diversity_contrasts, &
      chi_square_survival
   use ape_topology, only : is_binary_tree, is_ultrametric_tree, phylogenetic_vcv, mrca_many, &
      is_monophyletic, tip_descendant_matrix, clade_tips, topological_distance_ph85, branch_score_distance
   use ape_tree_edit, only : has_singles, is_rooted_tree, collapse_singles, drop_tips, keep_tips, extract_clade, &
      reroot_node, root_outgroup, unroot_tree, di2multi, multi2di
   use ape_splits, only : split_collection, prop_part, bitsplits, count_bipartitions, consensus_tree, &
      prop_clades, tree_bipartitions
   use ape_ordination, only : pcoa_result, pcoa_none, pcoa_lingoes, pcoa_cailliez, pcoa
   use ape_ace, only : ace_continuous_result, ace_continuous_ml, ace_continuous_reml, ace_continuous_gls
   use ape_discrete_ace, only : ace_discrete_result, ace_discrete_fit, ace_discrete_likelihood, ace_rate_index_matrix
   use ape_pgls, only : pgls_result, cor_brownian, cor_martins, grafen_tree, cor_grafen, cor_pagel, &
      cor_blomberg, pgls_fit, pgls_fit_model
   use ape_chronopl, only : chronopl_result, chronopl_fit, chronopl_objective, chronos_result, chronos_fit, &
      chronos_objective, chronos_clock_result, chronos_clock_fit
   use ape_compar_ou, only : compar_ou_result, compar_ou_fit, compar_ou_likelihood
   use ape_compar_lynch, only : compar_lynch_result, compar_lynch_fit
   use ape_corphylo, only : corphylo_result, corphylo_fit, corphylo_objective
   use ape_binary_pglmm, only : binary_pglmm_result, binary_pglmm_fit, binary_pglmm_reml_objective
   use ape_reconstruct, only : reconstruct_result, reconstruct_fit, reconstruct_gls_bm, reconstruct_gls_abm, &
      reconstruct_gls_ou_stationary, reconstruct_gls_ou
   use ape_skyline, only : skyline_result, collapsed_intervals, skyline_from_intervals, skyline_tree, &
      find_skyline_epsilon
   use ape_birthdeath, only : birthdeath_result, birthdeath_fit, birthdeath_from_times, birthdeath_deviance
   use ape_birthdeath_extended, only : birthdeath_extended_result, birthdeath_extended_fit, &
      birthdeath_extended_from_data, birthdeath_extended_deviance
   use ape_dna, only : dna_unknown, dna_a, dna_c, dna_g, dna_t, dna_gap, dna_r, dna_m, dna_w, dna_s, dna_k, dna_y, &
      dna_v, dna_h, dna_d, dna_b, dna_n, dna_distance, dna_distance_matrix, &
      dna_bh87_matrix, dna_distance_with_variance, dna_distance_matrix_with_variance, dna_base_frequencies, &
      dna_base_proportions, dna_leading_trailing_gaps_to_n, dna_global_deletion_mask, dna_segregating_sites, &
      dna_contingency_table, dna_pattern_positions, translate_dna
   implicit none
   public
end module ape
