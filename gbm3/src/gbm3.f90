! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of gbm3 3.0.3; see NOTICE.md and upstream/.
module gbm3
   use gbm3_kinds, only : dp
   use gbm3_constants
   use gbm3_types, only : gbm_options, gbm_node, gbm_tree, gbm_model, gbm_cv_result
   use gbm3_core, only : gbm_fit, gbm_fit_vector, gbm_fit_cox, gbm_continue, gbm_continue_vector, &
                         gbm_continue_cox, gbm_predict, gbm_predict_response, gbm_predict_trees, &
                         gbm_predict_staged, gbm_predict_response_staged, gbm_partial_dependence, &
                         gbm_relative_influence, gbm_set_seed
   use gbm3_cox, only : cox_baseline_hazard
   use gbm3_cv, only : gbm_cross_validate, gbm_cross_validate_vector, gbm_cross_validate_cox, gbm_cv_best_iteration
   use gbm3_diagnostics, only : gbm_best_iteration, gbm_permutation_importance, &
                                gbm_permutation_importance_vector, gbm_permutation_importance_cox, &
                                gbm_get_tree, gbm_get_node, gbm_num_nodes, gbm_interaction_strength
   implicit none
   public
end module gbm3
