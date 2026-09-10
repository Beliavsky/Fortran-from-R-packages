module grf_types
   use grf_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: grf_regression = 1
   integer, parameter, public :: grf_multi_regression = 2
   integer, parameter, public :: grf_probability = 3
   integer, parameter, public :: grf_quantile = 4
   integer, parameter, public :: grf_causal = 5
   integer, parameter, public :: grf_instrumental = 6
   integer, parameter, public :: grf_survival = 7
   integer, parameter, public :: grf_causal_survival = 8
   integer, parameter, public :: grf_lm = 9

   type, public :: grf_options
      integer :: num_trees = 200
      integer :: ci_group_size = 1
      real(dp) :: sample_fraction = 0.5_dp
      integer :: mtry = 0
      integer :: min_node_size = 5
      logical :: honesty = .true.
      real(dp) :: honesty_fraction = 0.5_dp
      logical :: honesty_prune_leaves = .true.
      real(dp) :: alpha = 0.05_dp
      real(dp) :: imbalance_penalty = 0.0_dp
      logical :: stabilize_splits = .true.
      logical :: fast_logrank = .false.
      integer :: seed = 42
      logical :: tune_parameters = .false.
      integer :: tune_num_trees = 10
      integer :: tune_num_reps = 100
      integer :: tune_num_draws = 1000
      integer, allocatable :: clusters(:)
      logical :: equalize_cluster_weights = .false.
   end type grf_options

   type, public :: grf_tree
      integer :: n_nodes = 0
      integer :: root = 1
      integer, allocatable :: left(:)
      integer, allocatable :: right(:)
      integer, allocatable :: split_var(:)
      real(dp), allocatable :: split_value(:)
      logical, allocatable :: missing_left(:)
      integer, allocatable :: depth(:)
      integer, allocatable :: leaf_start(:)
      integer, allocatable :: leaf_count(:)
      integer, allocatable :: leaf_samples(:)
      logical, allocatable :: inbag(:)
   end type grf_tree

   type, public :: grf_forest
      integer :: family = 0
      integer :: n_train = 0
      integer :: p = 0
      integer :: n_outputs = 0
      integer :: n_treatments = 0
      integer :: n_classes = 0
      integer :: n_times = 0
      integer :: n_clusters = 0
      integer :: samples_per_cluster = 0
      real(dp) :: horizon = 0.0_dp
      logical :: survival_probability_target = .false.
      real(dp) :: reduced_form_weight = 0.0_dp
      type(grf_options) :: options
      type(grf_tree), allocatable :: trees(:)
      real(dp), allocatable :: x(:,:)
      real(dp), allocatable :: y(:,:)
      real(dp), allocatable :: w(:,:)
      real(dp), allocatable :: z(:)
      real(dp), allocatable :: sample_weights(:)
      integer, allocatable :: clusters(:)
      logical :: equalize_cluster_weights = .false.
      real(dp), allocatable :: y_hat(:,:)
      real(dp), allocatable :: w_hat(:,:)
      real(dp), allocatable :: arm_propensity_hat(:,:)
      logical :: multi_arm_categorical = .false.
      real(dp), allocatable :: z_hat(:)
      real(dp), allocatable :: causal_survival_numerator(:)
      real(dp), allocatable :: causal_survival_denominator(:)
      real(dp), allocatable :: censor_survival_at_y(:)
      real(dp), allocatable :: treatment_variance_hat(:)
      real(dp), allocatable :: compliance_hat(:)
      logical :: ll_split_enabled = .false.
      real(dp) :: ll_split_lambda = 0.1_dp
      logical :: ll_split_weight_penalty = .false.
      integer :: ll_split_cutoff = 0
      integer, allocatable :: ll_split_variables(:)
      real(dp), allocatable :: ll_split_overall_beta(:)
      real(dp), allocatable :: failure_times(:)
      integer, allocatable :: survival_time_index(:)
      real(dp), allocatable :: quantiles(:)
      integer, allocatable :: classes(:)
      integer, allocatable :: event(:)
   end type grf_forest

   type, public :: grf_boosted_forest
      integer :: n_stages = 0
      real(dp) :: intercept = 0.0_dp
      real(dp) :: learning_rate = 0.1_dp
      type(grf_forest), allocatable :: stages(:)
      real(dp), allocatable :: training_prediction(:)
      real(dp), allocatable :: stage_error(:)
   end type grf_boosted_forest

   type, public :: grf_ate_result
      real(dp), allocatable :: estimate(:)
      real(dp), allocatable :: std_err(:)
   end type grf_ate_result

   type, public :: grf_linear_result
      real(dp), allocatable :: coefficient(:)
      real(dp), allocatable :: std_err(:)
      integer :: info = 0
   end type grf_linear_result

   type, public :: grf_rate_result
      real(dp), allocatable :: q(:)
      real(dp), allocatable :: toc(:)
      real(dp), allocatable :: toc_std_err(:)
      real(dp) :: estimate = 0.0_dp
      real(dp) :: std_err = 0.0_dp
   end type grf_rate_result

end module grf_types
