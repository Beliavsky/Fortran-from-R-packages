! SPDX-License-Identifier: GPL-2.0-or-later
! Computational translation derived from R package spdep 1.4-2.
! Translated/modified 2026-08-30 for modern free-form Fortran.
module spdep
   use spdep_kinds, only : dp
   use spdep_types, only : int_vector, real_vector, neighbor_list, spatial_weights, &
      knn_result, weights_constants, spatial_test_result, local_stat_result, &
      eb_result, mst_result, spatial_delta_result
   use spdep_graph, only : cell2nb, dnearneigh, knearneigh, knn2nb, nbdists, &
      gabrielneigh, relative_neighborhood, tri2nb, connected_components, &
      is_symmetric_nb, make_symmetric_nb, include_self, remove_self, nb_union, &
      nb_intersection, nb_difference, droplinks, addlinks, nblag, &
      nb_adjacency_matrix, graph_distance_matrix, euclidean_distance, &
      great_circle_distance
   use spdep_weights, only : card, nb2listw, nb2mat, mat2listw, listw2mat, &
      lag_listw, lag_listw_matrix, spweights_constants, szero, listw2U
   use spdep_statistics, only : moran, moran_test, moran_mc, moran_bv, &
      local_moran, local_moran_bv, geary, geary_test, geary_mc, &
      global_g_test, local_g, lee, local_lee, losh, joincount_test, &
      ebest, eblocal, choynowski, ebi_moran
   use spdep_clustering, only : nbcosts, ssw, mstree, prunecost, skater_groups
   use spdep_additional, only : rotation, complement_nb, nblag_cumul, nb2blocknb, &
      nb2listwdist, autocov_dist, local_geary
   use spdep_delta, only : spatialdelta, linearised_diffusive_weights, &
      metropolis_hastings_weights, iterative_proportional_fitting_weights, &
      graph_distance_weights, localdelta, cornish_fisher
   implicit none
   private

   public :: dp
   public :: int_vector, real_vector, neighbor_list, spatial_weights
   public :: knn_result, weights_constants, spatial_test_result, local_stat_result
   public :: eb_result, mst_result, spatial_delta_result
   public :: cell2nb, dnearneigh, knearneigh, knn2nb, nbdists
   public :: gabrielneigh, relative_neighborhood, tri2nb, connected_components
   public :: is_symmetric_nb, make_symmetric_nb, include_self, remove_self
   public :: nb_union, nb_intersection, nb_difference, droplinks, addlinks, nblag
   public :: nb_adjacency_matrix, graph_distance_matrix
   public :: euclidean_distance, great_circle_distance
   public :: card, nb2listw, nb2mat, mat2listw, listw2mat
   public :: lag_listw, lag_listw_matrix, spweights_constants, szero, listw2U
   public :: moran, moran_test, moran_mc, moran_bv, local_moran, local_moran_bv
   public :: geary, geary_test, geary_mc, global_g_test, local_g
   public :: lee, local_lee, losh, joincount_test, ebest, eblocal, choynowski
   public :: ebi_moran
   public :: nbcosts, ssw, mstree, prunecost, skater_groups
   public :: rotation, complement_nb, nblag_cumul, nb2blocknb
   public :: nb2listwdist, autocov_dist, local_geary
   public :: spatialdelta, linearised_diffusive_weights, metropolis_hastings_weights
   public :: iterative_proportional_fitting_weights, graph_distance_weights
   public :: localdelta, cornish_fisher

end module spdep
