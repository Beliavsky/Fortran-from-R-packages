! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of computational routines from wavethresh 4.7.3.
module wavethresh_types
   use r_kinds, only : dp
   implicit none
   private

   public :: dp
   public :: wt_vector_t, wt_matrix_t, wt_complex_vector_t, wt_filter_t, multiple_filter_t
   public :: wd_t, cwd_t, wp_t, mwd_t, imwd_level_t, imwd_t, wd3d_level_t, wd3d_t, interval_wavelet_t
   public :: signal_set_t, chirp_t, test_signal_t, support_t
   public :: rss_result_t, wavelet_cv_result_t, ewspec_result_t
   public :: complex_threshold_parameters_t, complex_threshold_result_t
   public :: first_last_t, grid_data_t, irregular_wavelet_t, scaling_function_t
   public :: density_projection_t, density_grid_t, density_wavelet_t
   public :: basis_selection_t, wpst_matrix_t, wpst_regression_t
   public :: lda_model_t, wpst_discrimination_t, wpst_classification_t
   public :: node_level_t, node_vector_t

   type :: wt_vector_t
      real(dp), allocatable :: values(:)
   end type wt_vector_t

   type :: wt_matrix_t
      real(dp), allocatable :: values(:,:)
   end type wt_matrix_t

   type :: wt_complex_vector_t
      complex(dp), allocatable :: values(:)
   end type wt_complex_vector_t

   type :: wt_filter_t
      character(len=24) :: family = ""
      real(dp) :: filter_number = 0.0_dp
      character(len=40) :: name = ""
      logical :: is_complex = .false.
      logical :: ok = .false.
      character(len=160) :: message = ""
      real(dp), allocatable :: low(:)
      real(dp), allocatable :: high(:)
      complex(dp), allocatable :: low_complex(:)
      complex(dp), allocatable :: high_complex(:)
   end type wt_filter_t

   type :: multiple_filter_t
      character(len=24) :: filter_type = ""
      character(len=48) :: name = ""
      integer :: nphi = 0
      integer :: npsi = 0
      integer :: nh = 0
      integer :: ndecim = 0
      real(dp), allocatable :: h(:)
      real(dp), allocatable :: g(:)
      logical :: ok = .false.
      character(len=160) :: message = ""
   end type multiple_filter_t

   type :: wd_t
      integer :: n_original = 0
      integer :: nlevels = 0
      character(len=16) :: transform_type = "wavelet"
      character(len=16) :: boundary = "periodic"
      type(wt_filter_t) :: filter
      type(wt_vector_t), allocatable :: detail(:)  ! 0:nlevels-1, coarse to fine.
      type(wt_vector_t), allocatable :: scaling(:) ! 0:nlevels, coarse to original.
      logical :: ok = .false.
      character(len=160) :: message = ""
   end type wd_t

   type :: cwd_t
      integer :: n_original = 0
      integer :: nlevels = 0
      character(len=16) :: transform_type = "wavelet"
      character(len=16) :: boundary = "periodic"
      type(wt_filter_t) :: filter
      type(wt_complex_vector_t), allocatable :: detail(:)
      type(wt_complex_vector_t), allocatable :: scaling(:)
      logical :: ok = .false.
      character(len=160) :: message = ""
   end type cwd_t

   type :: wp_t
      integer :: n_original = 0
      integer :: nlevels = 0
      logical :: stationary = .false.
      type(wt_filter_t) :: filter
      type(wt_vector_t), allocatable :: packet(:,:) ! (0:nlevels,0:2**nlevels-1)
      logical :: ok = .false.
      character(len=160) :: message = ""
   end type wp_t

   type :: mwd_t
      integer :: n_original = 0
      integer :: nlevels = 0
      character(len=16) :: transform_type = "wavelet"
      character(len=16) :: boundary = "periodic"
      character(len=16) :: prefilter = "default"
      type(multiple_filter_t) :: filter
      type(wt_matrix_t), allocatable :: detail(:) !! Multiple-wavelet detail matrices, indexed 0:nlevels-1.
      type(wt_matrix_t), allocatable :: scaling(:) !! Multiple-wavelet scaling matrices, indexed 0:nlevels.
      logical :: ok = .false.
      character(len=160) :: message = ""
   end type mwd_t

   type :: imwd_level_t
      real(dp), allocatable :: smooth(:,:) !! Scaling coefficients retained at this decomposition level.
      real(dp), allocatable :: lh(:,:)
      real(dp), allocatable :: hl(:,:)
      real(dp), allocatable :: hh(:,:)
   end type imwd_level_t

   type :: imwd_t
      integer :: nrow_original = 0
      integer :: ncol_original = 0
      integer :: nlevels = 0
      logical :: stationary = .false.
      type(wt_filter_t) :: filter
      type(imwd_level_t), allocatable :: level(:) ! 0:nlevels-1, coarse to fine.
      real(dp), allocatable :: smooth(:,:)
      logical :: ok = .false.
      character(len=160) :: message = ""
   end type imwd_t

   type :: wd3d_level_t
      real(dp), allocatable :: band(:,:,:,:) ! Seven detail bands.
   end type wd3d_level_t

   type :: wd3d_t
      integer :: n1_original = 0
      integer :: n2_original = 0
      integer :: n3_original = 0
      integer :: nlevels = 0
      type(wt_filter_t) :: filter
      type(wd3d_level_t), allocatable :: level(:) ! 0:nlevels-1, coarse to fine.
      real(dp), allocatable :: smooth(:,:,:)
      logical :: ok = .false.
      character(len=160) :: message = ""
   end type wd3d_t

   type :: interval_wavelet_t
      real(dp), allocatable :: transformed(:) !! Packed interval-wavelet scaling and detail coefficients.
      integer, allocatable :: filters_used(:) !! Boundary-filter orders, ordered from finest to coarsest step.
      integer :: current_scale = 0 !! Coarsest scale retained by the decomposition.
      logical :: preconditioned = .false. !! Whether boundary preconditioning was applied.
      logical :: ok = .false. !! True when the transform completed successfully.
      character(len=160) :: message = "" !! Validation or transform status.
   end type interval_wavelet_t

   type :: complex_threshold_parameters_t
      real(dp), allocatable :: nonzero_probability(:) !! Estimated mixture probability by R-style level.
      real(dp), allocatable :: signal_covariance(:,:,:) !! Estimated nonzero covariance, indexed (level,2,2).
      real(dp), allocatable :: noise_covariance(:,:,:) !! Wavelet-noise covariance, indexed (level,2,2).
      logical :: ok = .false. !! True when all requested level fits completed.
      character(len=160) :: message = "" !! Validation or optimizer status.
   end type complex_threshold_parameters_t

   type :: complex_threshold_result_t
      type(cwd_t) :: transform !! Original complex wavelet decomposition.
      type(cwd_t) :: thresholded !! Thresholded complex wavelet decomposition.
      type(complex_threshold_parameters_t) :: parameters !! Empirical-Bayes parameters when requested.
      complex(dp), allocatable :: estimate(:) !! Reconstructed complex-valued estimate.
      real(dp) :: noise_variance = 0.0_dp !! Estimated input noise variance.
      logical :: ok = .false. !! True when thresholding and reconstruction completed.
      character(len=160) :: message = "" !! Validation or threshold status.
   end type complex_threshold_result_t

   type :: signal_set_t
      real(dp), allocatable :: blocks(:)
      real(dp), allocatable :: bumps(:)
      real(dp), allocatable :: heavi_sine(:)
      real(dp), allocatable :: doppler(:)
   end type signal_set_t

   type :: chirp_t
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: y(:)
   end type chirp_t


   type :: test_signal_t
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: signal(:)
      real(dp), allocatable :: noisy(:)
      character(len=24) :: signal_type = ""
      logical :: ok = .false.
      character(len=160) :: message = ""
   end type test_signal_t

   type :: support_t
      real(dp) :: left = 0.0_dp
      real(dp) :: right = 0.0_dp
      real(dp) :: psi_left = 0.0_dp
      real(dp) :: psi_right = 0.0_dp
      real(dp) :: phi_left = 0.0_dp
      real(dp) :: phi_right = 0.0_dp
      logical :: ok = .false.
      character(len=160) :: message = ""
   end type support_t

   type :: rss_result_t
      real(dp) :: ssq = 0.0_dp
      real(dp) :: threshold = 0.0_dp
      integer :: df = 0
      logical :: ok = .false.
      character(len=160) :: message = ""
   end type rss_result_t


   type :: first_last_t
      integer, allocatable :: scaling(:,:) !! Scaling-coefficient first, last, and packed-offset columns.
      integer, allocatable :: detail(:,:) !! Detail-coefficient first, last, and packed-offset columns.
      integer :: ntotal = 0
      integer :: ntotal_detail = 0
      logical :: ok = .false.
      character(len=160) :: message = ""
   end type first_last_t

   type :: grid_data_t
      real(dp), allocatable :: grid_t(:) !! Regular grid locations in (0, 1).
      real(dp), allocatable :: grid_y(:) !! Interpolated observations on the regular grid.
      real(dp), allocatable :: weight_left(:) !! Weight on the left sorted observation for each grid point.
      integer, allocatable :: left_index(:) !! One-based index in the sorted observations for each grid point.
      integer :: n_observations = 0 !! Number of irregular observations represented by the interpolation map.
      logical :: ok = .false.
      character(len=160) :: message = ""
   end type grid_data_t

   type :: irregular_wavelet_t
      type(wd_t) :: transform !! Wavelet transform of the interpolated regular-grid observations.
      type(wt_vector_t), allocatable :: coefficient_variance(:) !! Noise variance multipliers by R-style detail level.
      logical :: ok = .false. !! True when decomposition and variance propagation succeeded.
      character(len=160) :: message = "" !! Validation or transform status.
   end type irregular_wavelet_t

   type :: scaling_function_t
      real(dp), allocatable :: x(:) !! Sample locations of the scaling function.
      real(dp), allocatable :: y(:) !! Scaling-function values at the sample locations.
      logical :: ok = .false.
      character(len=160) :: message = ""
   end type scaling_function_t

   type :: density_projection_t
      real(dp), allocatable :: coefficients(:) !! Empirical high-resolution scaling coefficients for k_min:k_max.
      real(dp), allocatable :: covariance(:,:) !! Upper covariance bands; column one is the main diagonal.
      type(wt_filter_t) :: filter !! Scaling-function filter used for the projection.
      real(dp) :: primary_resolution = 0.0_dp !! Product tau*2**resolution_level.
      real(dp) :: tau = 0.0_dp !! Fine-tuning multiplier applied to the dyadic resolution.
      integer :: resolution_level = 0 !! Requested high-resolution level J.
      integer :: iterations = 0 !! Daubechies-Lagarias binary-product iterations.
      integer :: sample_size = 0 !! Number of observations used to estimate the coefficients.
      integer :: k_min = 0 !! Smallest integer translate represented by coefficients.
      integer :: k_max = -1 !! Largest integer translate represented by coefficients.
      logical :: has_covariance = .false. !! True when covariance bands were requested and computed.
      logical :: ok = .false.
      character(len=160) :: message = ""
   end type density_projection_t

   type :: density_grid_t
      real(dp), allocatable :: x(:) !! Locations at which the projected density was evaluated.
      real(dp), allocatable :: y(:) !! Projected density values at the corresponding locations.
      integer, allocatable :: scaling_indices(:) !! Integer scaling-function translations used by CWavDE.
      integer, allocatable :: wavelet_k_min(:) !! Smallest wavelet translation at each CWavDE level.
      integer, allocatable :: wavelet_k_max(:) !! Largest wavelet translation at each CWavDE level.
      logical :: ok = .false.
      character(len=160) :: message = ""
   end type density_grid_t

   type :: density_wavelet_t
      type(wt_filter_t) :: filter !! Wavelet filter used for the zero-boundary decomposition.
      type(wt_vector_t), allocatable :: scaling(:) !! Scaling coefficients indexed from level zero through nlevels.
      type(wt_vector_t), allocatable :: detail(:) !! Detail coefficients indexed from level zero through nlevels-1.
      integer, allocatable :: scaling_first(:) !! First integer scaling-coefficient index at each level.
      integer, allocatable :: scaling_last(:) !! Last integer scaling-coefficient index at each level.
      integer, allocatable :: detail_first(:) !! First integer detail-coefficient index at each level.
      integer, allocatable :: detail_last(:) !! Last integer detail-coefficient index at each level.
      integer :: nlevels = 0 !! Number of zero-boundary decomposition steps.
      logical :: ok = .false.
      character(len=160) :: message = ""
   end type density_wavelet_t

   type :: basis_selection_t
      integer, allocatable :: index(:) !! One-based selected basis-column indices in descending score order.
      real(dp), allocatable :: score(:) !! Absolute correlation scores corresponding to selected columns.
      logical :: ok = .false.
      character(len=160) :: message = ""
   end type basis_selection_t

   type :: wpst_matrix_t
      real(dp), allocatable :: matrix(:,:) !! Stationary packet features with observations in rows and packets in columns.
      integer, allocatable :: level(:) !! Upstream R-style packet level associated with each feature column.
      integer, allocatable :: packet_index(:) !! Zero-based packet index associated with each feature column.
      integer, allocatable :: groups(:) !! Optional integer group labels associated with the observation rows.
      integer :: nlevels = 0
      logical :: ok = .false.
      character(len=160) :: message = ""
   end type wpst_matrix_t

   type :: wpst_regression_t
      real(dp), allocatable :: response(:) !! Training response values paired with stationary-packet feature rows.
      real(dp), allocatable :: matrix(:,:) !! Selected training basis matrix, with observations in rows.
      integer, allocatable :: original_index(:) !! One-based column indices in the complete wpst2m basis matrix.
      integer, allocatable :: level(:) !! R-style stationary packet level for each selected basis column.
      integer, allocatable :: packet_index(:) !! Zero-based stationary packet index for each selected basis column.
      real(dp), allocatable :: score(:) !! Signed response correlation for each selected basis column.
      integer :: nlevels = 0
      type(wt_filter_t) :: filter
      character(len=24) :: transform = "logabs"
      logical :: ok = .false.
      character(len=160) :: message = ""
   end type wpst_regression_t

   type :: lda_model_t
      integer, allocatable :: classes(:) !! Sorted integer class labels represented by the fitted discriminant model.
      real(dp), allocatable :: prior(:) !! Empirical class priors, one probability per class.
      real(dp), allocatable :: means(:,:) !! Class means with classes in rows and selected features in columns.
      real(dp), allocatable :: inverse_covariance(:,:) !! Inverse pooled within-class covariance matrix.
      logical :: ok = .false.
      character(len=160) :: message = ""
   end type lda_model_t

   type :: wpst_discrimination_t
      real(dp), allocatable :: matrix(:,:) !! Selected training discrimination features with observations in rows.
      integer, allocatable :: groups(:) !! Integer training class labels, one per feature row.
      integer, allocatable :: level(:) !! R-style stationary packet level for each selected feature.
      integer, allocatable :: packet_index(:) !! Zero-based stationary packet index for each selected feature.
      real(dp), allocatable :: basiscoef(:) !! Absolute group-correlation score for each selected feature.
      type(wt_filter_t) :: filter !! Wavelet filter used to construct the stationary packet features.
      type(lda_model_t) :: lda !! Pooled-covariance linear discriminant model fitted to the selected features.
      logical :: ok = .false.
      character(len=160) :: message = ""
   end type wpst_discrimination_t

   type :: wpst_classification_t
      real(dp), allocatable :: basis_matrix(:,:) !! Selected packet features extracted from the new time series.
      real(dp), allocatable :: discriminant_score(:,:) !! LDA class scores with observations in rows and classes in columns.
      integer, allocatable :: predicted_group(:) !! Predicted integer class label for each observation.
      logical :: ok = .false.
      character(len=160) :: message = ""
   end type wpst_classification_t

   type :: node_level_t
      character(len=1), allocatable :: upperctrl(:) !! S/L/R stationary, T/B packet, or L/R path controls for this node level.
      real(dp), allocatable :: upperl(:) !! Minimum-entropy node scores or zero-filled path metadata, depending on producer.
   end type node_level_t

   type :: node_vector_t
      type(node_level_t), allocatable :: node(:) !! Node-control levels corresponding to the upstream node.list entries.
      integer :: nlevels = 0
      logical :: ok = .false.
      character(len=160) :: message = ""
   end type node_vector_t

   type :: ewspec_result_t
      type(wd_t) :: spectrum !! Corrected evolutionary wavelet-spectrum coefficients stored in a stationary wavelet object.
      type(wd_t) :: periodogram !! Raw squared stationary-wavelet periodogram before autocorrelation-wavelet correction.
      real(dp), allocatable :: rm(:,:) !! Autocorrelation-wavelet inner-product matrix A.
      real(dp), allocatable :: irm(:,:) !! Inverse autocorrelation-wavelet inner-product matrix used for correction.
      logical :: ok = .false.
      character(len=160) :: message = ""
   end type ewspec_result_t

   type :: wavelet_cv_result_t
      real(dp) :: cv_threshold = 0.0_dp
      real(dp), allocatable :: cv_thresholds(:) !! Optional level-specific CV thresholds followed by any final top-level threshold.
      real(dp) :: universal_threshold = 0.0_dp
      integer :: cv_df = 0
      integer :: universal_df = 0
      real(dp), allocatable :: cv_reconstruction(:)
      real(dp), allocatable :: universal_reconstruction(:)
      real(dp), allocatable :: trial_thresholds(:)
      real(dp), allocatable :: trial_errors(:)
      logical :: ok = .false.
      character(len=160) :: message = ""
   end type wavelet_cv_result_t

end module wavethresh_types
