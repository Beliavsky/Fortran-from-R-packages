! SPDX-License-Identifier: GPL-2.0-or-later
module rf_types
   use r_kinds, only : dp, i64
   implicit none
   private

   integer, parameter, public :: RF_TERMINAL = -1
   integer, parameter, public :: RF_INTERIOR = 1

   type, public :: rf_options
      integer :: ntree = 500
      integer :: mtry = 0
      integer :: nodesize = 0
      integer :: maxnodes = 0
      integer :: sample_size = 0
      integer :: n_perm = 1
      logical :: replace = .true.
      logical :: importance = .false.
      logical :: proximity = .false.
      logical :: oob_proximity = .false.
      logical :: keep_inbag = .false.
      logical :: bias_correct = .false.
      integer(i64) :: seed = 1976_i64
   end type rf_options

   type, public :: rf_tree
      integer :: n_nodes = 0
      integer :: maxcat = 1
      integer, allocatable :: left(:)
      integer, allocatable :: right(:)
      integer, allocatable :: split_var(:)
      integer, allocatable :: status(:)
      integer, allocatable :: node_class(:)
      real(dp), allocatable :: split_value(:)
      real(dp), allocatable :: node_mean(:)
      logical, allocatable :: cat_left(:,:)
      real(dp), allocatable :: impurity_decrease(:)
   end type rf_tree

   type, public :: rf_classification_forest
      integer :: nclass = 0
      integer :: nvar = 0
      integer :: nobs = 0
      integer, allocatable :: ncat(:)
      type(rf_tree), allocatable :: trees(:)
      real(dp), allocatable :: cutoff(:)
      integer, allocatable :: oob_votes(:,:)
      integer, allocatable :: oob_count(:)
      integer, allocatable :: oob_prediction(:)
      real(dp), allocatable :: error_curve(:,:)
      real(dp), allocatable :: importance_accuracy(:,:)
      real(dp), allocatable :: importance_sd(:,:)
      real(dp), allocatable :: importance_gini(:)
      real(dp), allocatable :: proximity(:,:)
      integer, allocatable :: inbag(:,:)
   end type rf_classification_forest

   type, public :: rf_regression_forest
      integer :: nvar = 0
      integer :: nobs = 0
      integer, allocatable :: ncat(:)
      type(rf_tree), allocatable :: trees(:)
      real(dp), allocatable :: oob_prediction(:)
      integer, allocatable :: oob_count(:)
      real(dp), allocatable :: mse_curve(:)
      real(dp), allocatable :: importance_accuracy(:)
      real(dp), allocatable :: importance_sd(:)
      real(dp), allocatable :: importance_gini(:)
      real(dp), allocatable :: proximity(:,:)
      integer, allocatable :: inbag(:,:)
      real(dp) :: bias_intercept = 0.0_dp
      real(dp) :: bias_slope = 1.0_dp
   end type rf_regression_forest

end module rf_types
