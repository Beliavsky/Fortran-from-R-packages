! SPDX-License-Identifier: GPL-2.0-or-later
module randomforest
   use r_kinds, only : dp
   use rf_types, only : rf_options, rf_tree, rf_classification_forest, rf_regression_forest
   use rf_classification, only : fit_classification, predict_classification, fit_unsupervised, make_synthetic_class
   use rf_regression, only : fit_regression, predict_regression
   use rf_utilities, only : classification_margin, outlier_scores, roughfix_numeric, class_centers
   use rf_utilities, only : tree_sizes, variable_usage, partial_dependence_regression
   use rf_utilities, only : partial_dependence_classification
   use rf_workflows, only : rf_impute_classification, rf_impute_regression
   use rf_workflows, only : tune_classification_mtry, tune_regression_mtry
   use rf_workflows, only : rfcv_classification, rfcv_regression, mds_coordinates
   implicit none
   private

   public :: dp
   public :: rf_options, rf_tree, rf_classification_forest, rf_regression_forest
   public :: fit_classification, predict_classification, fit_unsupervised, make_synthetic_class
   public :: fit_regression, predict_regression
   public :: classification_margin, outlier_scores, roughfix_numeric, class_centers
   public :: tree_sizes, variable_usage, partial_dependence_regression, partial_dependence_classification
   public :: rf_impute_classification, rf_impute_regression
   public :: tune_classification_mtry, tune_regression_mtry
   public :: rfcv_classification, rfcv_regression, mds_coordinates

end module randomforest
