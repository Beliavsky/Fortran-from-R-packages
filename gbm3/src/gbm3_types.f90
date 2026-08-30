! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of gbm3 3.0.3; see NOTICE.md and upstream/.
module gbm3_types
   use gbm3_kinds, only : dp
   use gbm3_constants
   implicit none
   private

   type, public :: gbm_options
      integer :: distribution = GBM_GAUSSIAN
      integer :: num_trees = 100
      integer :: interaction_depth = 1
      integer :: min_num_obs_in_node = 10
      real(dp) :: shrinkage = 0.001_dp
      real(dp) :: bag_fraction = 0.5_dp
      integer :: num_features = 0
      integer :: num_train = 0
      real(dp) :: quantile_alpha = 0.25_dp
      real(dp) :: t_df = 4.0_dp
      real(dp) :: tweedie_power = 1.5_dp
      integer :: cox_ties = GBM_TIES_EFRON
      real(dp) :: cox_prior_node_coeff_var = 1000.0_dp
      integer :: pairwise_metric = GBM_METRIC_NDCG
      integer :: pairwise_max_rank = 0
   end type gbm_options

   type, public :: gbm_node
      logical :: is_terminal = .true.
      integer :: split_var = 0
      integer :: left = 0
      integer :: right = 0
      integer :: missing = 0
      real(dp) :: split_value = 0.0_dp
      real(dp) :: improvement = 0.0_dp
      real(dp) :: prediction = 0.0_dp
      real(dp) :: total_weight = 0.0_dp
      integer :: num_obs = 0
      integer, allocatable :: left_categories(:)
   end type gbm_node

   type, public :: gbm_tree
      type(gbm_node), allocatable :: nodes(:)
      integer :: n_nodes = 0
      real(dp) :: shrinkage = 0.0_dp
   end type gbm_tree

   type, public :: gbm_model
      type(gbm_options) :: options
      integer :: n_features = 0
      integer :: n_rows = 0
      integer :: n_train = 0
      integer :: n_trees = 0
      real(dp) :: init_f = 0.0_dp
      type(gbm_tree), allocatable :: trees(:)
      real(dp), allocatable :: train_error(:)
      real(dp), allocatable :: validation_error(:)
      real(dp), allocatable :: oob_improvement(:)
      real(dp), allocatable :: fitted(:)
      integer, allocatable :: var_classes(:)
      integer, allocatable :: monotone(:)
   end type gbm_model

   type, public :: gbm_cv_result
      integer :: n_folds = 0
      integer :: best_iteration = 0
      real(dp), allocatable :: error(:)
      real(dp), allocatable :: fitted(:)
      integer, allocatable :: fold_id(:)
   end type gbm_cv_result

   type :: node_stat
      integer :: n = 0
      real(dp) :: wr = 0.0_dp
      real(dp) :: w = 0.0_dp
   contains
      procedure :: prediction => node_stat_prediction
   end type node_stat

   type :: split_candidate
      logical :: valid = .false.
      integer :: node_id = 0
      integer :: split_var = 0
      integer :: split_class = 0
      real(dp) :: split_value = 0.0_dp
      real(dp) :: improvement = -huge(1.0_dp)
      type(node_stat) :: left
      type(node_stat) :: right
      type(node_stat) :: missing
      integer, allocatable :: ordering(:)
      integer :: bias = huge(0)
   end type split_candidate

   public :: node_stat, split_candidate

contains
   pure real(dp) function node_stat_prediction(self) result(v)
      class(node_stat), intent(in) :: self
      if (abs(self%w) <= tiny(1.0_dp)) then
         v = 0.0_dp
      else
         v = self%wr / self%w
      end if
   end function node_stat_prediction
end module gbm3_types
