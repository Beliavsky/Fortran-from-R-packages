! SPDX-License-Identifier: GPL-3.0-only
! Modern Fortran translation of ranger 0.18.0 computational code.
! Upstream C++ core: Copyright (c) 2014-2018 Marvin N. Wright; MIT licensed.
! The upstream R package is GPL-3. See NOTICE.md for full provenance.
module ranger
   use r_kinds, only : dp
   use ranger_types, only : ranger_options, ranger_tree
   use ranger_types, only : ranger_classification_forest, ranger_regression_forest
   use ranger_types, only : ranger_probability_forest, ranger_survival_forest
   use ranger_types, only : RANGER_IMPORTANCE_NONE, RANGER_IMPORTANCE_IMPURITY
   use ranger_types, only : RANGER_IMPORTANCE_PERMUTATION, RANGER_IMPORTANCE_IMPURITY_CORRECTED
   use ranger_types, only : RANGER_IMPORTANCE_PERMUTATION_CASEWISE
   use ranger_types, only : RANGER_SPLIT_STANDARD, RANGER_SPLIT_AUC, RANGER_SPLIT_AUC_IGNORE_TIES
   use ranger_types, only : RANGER_SPLIT_MAXSTAT, RANGER_SPLIT_EXTRATREES, RANGER_SPLIT_BETA
   use ranger_types, only : RANGER_SPLIT_HELLINGER, RANGER_SPLIT_POISSON
   use ranger_types, only : RANGER_UNORDERED_IGNORE, RANGER_UNORDERED_ORDER, RANGER_UNORDERED_PARTITION
   use ranger_classification, only : fit_ranger_classification, predict_ranger_classification
   use ranger_classification, only : fit_ranger_probability, predict_ranger_probability
   use ranger_regression, only : fit_ranger_regression, predict_ranger_regression, predict_ranger_quantiles
   use ranger_regression, only : predict_ranger_quantiles_oob
   use ranger_survival, only : fit_ranger_survival, predict_ranger_survival
   use ranger_utilities, only : importance_pvalues_janitza, importance_pvalues_altmann, infinitesimal_jackknife
   use ranger_utilities, only : hierarchical_shrink_regression, hierarchical_shrink_probability
   use ranger_utilities, only : case_specific_weights, tree_sizes, variable_usage
   implicit none
   private

   public :: dp
   public :: ranger_options, ranger_tree
   public :: ranger_classification_forest, ranger_regression_forest, ranger_probability_forest, ranger_survival_forest
   public :: RANGER_IMPORTANCE_NONE, RANGER_IMPORTANCE_IMPURITY, RANGER_IMPORTANCE_PERMUTATION
   public :: RANGER_IMPORTANCE_IMPURITY_CORRECTED, RANGER_IMPORTANCE_PERMUTATION_CASEWISE
   public :: RANGER_SPLIT_STANDARD, RANGER_SPLIT_AUC, RANGER_SPLIT_AUC_IGNORE_TIES, RANGER_SPLIT_MAXSTAT
   public :: RANGER_SPLIT_EXTRATREES, RANGER_SPLIT_BETA, RANGER_SPLIT_HELLINGER, RANGER_SPLIT_POISSON
   public :: RANGER_UNORDERED_IGNORE, RANGER_UNORDERED_ORDER, RANGER_UNORDERED_PARTITION
   public :: fit_ranger_classification, predict_ranger_classification
   public :: fit_ranger_probability, predict_ranger_probability
   public :: fit_ranger_regression, predict_ranger_regression, predict_ranger_quantiles, predict_ranger_quantiles_oob
   public :: fit_ranger_survival, predict_ranger_survival
   public :: importance_pvalues_janitza, importance_pvalues_altmann, infinitesimal_jackknife
   public :: hierarchical_shrink_regression, hierarchical_shrink_probability
   public :: case_specific_weights, tree_sizes, variable_usage

end module ranger
