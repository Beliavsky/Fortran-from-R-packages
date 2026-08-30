! SPDX-License-Identifier: GPL-3.0-only
! Modern Fortran translation of ranger 0.18.0 computational code.
! Upstream C++ core: Copyright (c) 2014-2018 Marvin N. Wright; MIT licensed.
! The upstream R package is GPL-3. See NOTICE.md for full provenance.
module ranger_types
   use r_kinds, only : dp, i64
   implicit none
   private

   integer, parameter, public :: RANGER_TERMINAL = -1
   integer, parameter, public :: RANGER_INTERIOR = 1

   integer, parameter, public :: RANGER_IMPORTANCE_NONE = 0
   integer, parameter, public :: RANGER_IMPORTANCE_IMPURITY = 1
   integer, parameter, public :: RANGER_IMPORTANCE_PERMUTATION = 3
   integer, parameter, public :: RANGER_IMPORTANCE_IMPURITY_CORRECTED = 5
   integer, parameter, public :: RANGER_IMPORTANCE_PERMUTATION_CASEWISE = 6

   integer, parameter, public :: RANGER_SPLIT_STANDARD = 1
   integer, parameter, public :: RANGER_SPLIT_AUC = 2
   integer, parameter, public :: RANGER_SPLIT_AUC_IGNORE_TIES = 3
   integer, parameter, public :: RANGER_SPLIT_MAXSTAT = 4
   integer, parameter, public :: RANGER_SPLIT_EXTRATREES = 5
   integer, parameter, public :: RANGER_SPLIT_BETA = 6
   integer, parameter, public :: RANGER_SPLIT_HELLINGER = 7
   integer, parameter, public :: RANGER_SPLIT_POISSON = 8

   integer, parameter, public :: RANGER_UNORDERED_IGNORE = 0
   integer, parameter, public :: RANGER_UNORDERED_ORDER = 1
   integer, parameter, public :: RANGER_UNORDERED_PARTITION = 2

   type, public :: ranger_options
      integer :: num_trees = 500
      integer :: mtry = 0
      integer :: min_node_size = 0
      integer :: min_bucket = 0
      integer :: max_depth = 0
      integer :: split_rule = RANGER_SPLIT_STANDARD
      integer :: num_random_splits = 1
      integer :: importance_mode = RANGER_IMPORTANCE_NONE
      integer :: respect_unordered_factors = RANGER_UNORDERED_IGNORE
      integer :: n_perm = 1
      real(dp) :: sample_fraction = -1.0_dp
      real(dp) :: alpha = 0.5_dp
      real(dp) :: minprop = 0.1_dp
      real(dp) :: poisson_tau = 1.0_dp
      logical :: replace = .true.
      logical :: scale_permutation_importance = .false.
      logical :: local_importance = .false.
      logical :: keep_inbag = .false.
      logical :: holdout = .false.
      logical :: quantreg = .false.
      logical :: oob_error = .true.
      logical :: node_stats = .false.
      logical :: regularization_usedepth = .false.
      logical :: na_learn = .true.
      integer(i64) :: seed = 1976_i64
      real(dp), allocatable :: split_select_weights(:)
      real(dp), allocatable :: split_select_weights_by_tree(:,:)
      real(dp), allocatable :: regularization_factor(:)
      integer, allocatable :: always_split_variables(:)
   end type ranger_options

   type, public :: ranger_tree
      integer :: n_nodes = 0
      integer :: nclass = 0
      integer :: ntime = 0
      integer :: maxcat = 1
      integer, allocatable :: left(:)
      integer, allocatable :: right(:)
      integer, allocatable :: split_var(:)
      integer, allocatable :: status(:)
      integer, allocatable :: node_class(:)
      integer, allocatable :: node_n(:)
      real(dp), allocatable :: split_value(:)
      real(dp), allocatable :: node_mean(:)
      real(dp), allocatable :: split_stat(:)
      logical, allocatable :: nan_go_right(:)
      logical, allocatable :: cat_left(:,:)
      real(dp), allocatable :: class_prob(:,:)
      real(dp), allocatable :: chf(:,:)
      real(dp), allocatable :: impurity_decrease(:)
   end type ranger_tree

   type, public :: ranger_classification_forest
      integer :: nclass = 0
      integer :: nvar = 0
      integer :: nobs = 0
      integer, allocatable :: ncat(:)
      integer, allocatable :: category_map(:,:)
      type(ranger_tree), allocatable :: trees(:)
      integer, allocatable :: oob_prediction(:)
      integer, allocatable :: oob_count(:)
      real(dp), allocatable :: oob_votes(:,:)
      real(dp), allocatable :: prediction_error(:)
      real(dp), allocatable :: variable_importance(:)
      real(dp), allocatable :: variable_importance_sd(:)
      real(dp), allocatable :: local_importance(:,:)
      integer, allocatable :: inbag(:,:)
   end type ranger_classification_forest

   type, public :: ranger_regression_forest
      integer :: nvar = 0
      integer :: nobs = 0
      integer, allocatable :: ncat(:)
      integer, allocatable :: category_map(:,:)
      type(ranger_tree), allocatable :: trees(:)
      real(dp), allocatable :: oob_prediction(:)
      integer, allocatable :: oob_count(:)
      real(dp), allocatable :: prediction_error(:)
      real(dp), allocatable :: variable_importance(:)
      real(dp), allocatable :: variable_importance_sd(:)
      real(dp), allocatable :: local_importance(:,:)
      integer, allocatable :: inbag(:,:)
      real(dp), allocatable :: training_y(:)
      integer, allocatable :: training_terminal(:,:)
      real(dp), allocatable :: random_node_value(:,:)
      logical, allocatable :: random_node_has_value(:,:)
      real(dp), allocatable :: random_node_value_oob(:,:)
      integer :: quantile_oob_count = 0
   end type ranger_regression_forest

   type, public :: ranger_probability_forest
      integer :: nclass = 0
      integer :: nvar = 0
      integer :: nobs = 0
      integer, allocatable :: ncat(:)
      integer, allocatable :: category_map(:,:)
      type(ranger_tree), allocatable :: trees(:)
      real(dp), allocatable :: oob_probability(:,:)
      integer, allocatable :: oob_count(:)
      real(dp), allocatable :: prediction_error(:)
      real(dp), allocatable :: variable_importance(:)
      real(dp), allocatable :: variable_importance_sd(:)
      real(dp), allocatable :: local_importance(:,:)
      integer, allocatable :: inbag(:,:)
   end type ranger_probability_forest

   type, public :: ranger_survival_forest
      integer :: nvar = 0
      integer :: nobs = 0
      integer :: ntime = 0
      integer, allocatable :: ncat(:)
      integer, allocatable :: category_map(:,:)
      type(ranger_tree), allocatable :: trees(:)
      real(dp), allocatable :: unique_timepoints(:)
      real(dp), allocatable :: oob_chf(:,:)
      integer, allocatable :: oob_count(:)
      real(dp), allocatable :: prediction_error(:)
      real(dp), allocatable :: variable_importance(:)
      real(dp), allocatable :: variable_importance_sd(:)
      real(dp), allocatable :: local_importance(:,:)
      integer, allocatable :: inbag(:,:)
   end type ranger_survival_forest

end module ranger_types
