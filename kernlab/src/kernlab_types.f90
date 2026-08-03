! SPDX-License-Identifier: GPL-2.0-only
module kernlab_types
  use kernlab_kinds, only: dp, KL_SUCCESS
  implicit none
  private

  integer, parameter, public :: KERNEL_RBF = 1
  integer, parameter, public :: KERNEL_LAPLACE = 2
  integer, parameter, public :: KERNEL_BESSEL = 3
  integer, parameter, public :: KERNEL_POLY = 4
  integer, parameter, public :: KERNEL_TANH = 5
  integer, parameter, public :: KERNEL_LINEAR = 6
  integer, parameter, public :: KERNEL_ANOVA = 7
  integer, parameter, public :: KERNEL_SPLINE = 8
  integer, parameter, public :: KERNEL_STRING = 9
  integer, parameter, public :: KERNEL_FOURIER = 10

  integer, parameter, public :: MODEL_REGRESSION = 1
  integer, parameter, public :: MODEL_CLASSIFICATION = 2
  integer, parameter, public :: MODEL_NOVELTY = 3

  type, public :: kernel_spec
    integer :: kind = KERNEL_RBF
    real(dp) :: sigma = 1.0_dp
    real(dp) :: scale = 1.0_dp
    real(dp) :: offset = 1.0_dp
    real(dp) :: lambda = 1.0_dp
    integer :: degree = 1
    integer :: order = 1
    integer :: string_length = 4
    logical :: normalized = .true.
  end type kernel_spec

  type, public :: inchol_result
    real(dp), allocatable :: factor(:,:)
    integer, allocatable :: pivots(:)
    real(dp), allocatable :: diag_residues(:)
    real(dp), allocatable :: max_residuals(:)
    integer :: rank = 0
    integer :: status = KL_SUCCESS
  end type inchol_result

  type, public :: ipop_result
    real(dp), allocatable :: primal(:)
    real(dp), allocatable :: dual(:)
    real(dp) :: objective = huge(1.0_dp)
    real(dp) :: primal_infeasibility = huge(1.0_dp)
    integer :: iterations = 0
    integer :: status = KL_SUCCESS
  end type ipop_result

  type, public :: kpca_result
    type(kernel_spec) :: kernel
    real(dp), allocatable :: train(:,:)
    real(dp), allocatable :: coefficients(:,:)
    real(dp), allocatable :: rotated(:,:)
    real(dp), allocatable :: eigenvalues(:)
    real(dp), allocatable :: kernel_row_mean(:)
    real(dp) :: kernel_grand_mean = 0.0_dp
    integer :: status = KL_SUCCESS
  end type kpca_result

  type, public :: kcca_result
    type(kernel_spec) :: kernel
    real(dp), allocatable :: xtrain(:,:), ytrain(:,:)
    real(dp), allocatable :: xcoef(:,:), ycoef(:,:), correlations(:)
    integer :: status = KL_SUCCESS
  end type kcca_result

  type, public :: cluster_result
    integer, allocatable :: labels(:)
    real(dp), allocatable :: embedding(:,:)
    real(dp), allocatable :: centers(:,:)
    real(dp), allocatable :: withinss(:)
    integer :: iterations = 0
    integer :: status = KL_SUCCESS
  end type cluster_result

  type, public :: kernel_model
    type(kernel_spec) :: kernel
    integer :: model_type = MODEL_REGRESSION
    real(dp), allocatable :: train(:,:)
    real(dp), allocatable :: coefficients(:,:)
    real(dp), allocatable :: bias(:)
    integer, allocatable :: class_labels(:)
    real(dp), allocatable :: auxiliary(:,:)
    real(dp) :: noise = 0.0_dp
    integer :: iterations = 0
    integer :: status = KL_SUCCESS
  end type kernel_model

  type, public :: mmd_result
    real(dp) :: mmd1 = 0.0_dp
    real(dp) :: mmd3 = 0.0_dp
    real(dp) :: rademacher_bound = 0.0_dp
    real(dp) :: bootstrap_bound = 0.0_dp
    logical :: reject_rademacher = .false.
    logical :: reject_bootstrap = .false.
    integer :: status = KL_SUCCESS
  end type mmd_result

  type, public :: ranking_result
    real(dp), allocatable :: score(:)
    real(dp), allocatable :: rank(:)
    real(dp), allocatable :: convergence(:,:)
    real(dp), allocatable :: edge_graph(:,:)
    integer :: status = KL_SUCCESS
  end type ranking_result

  type, public :: csi_result
    real(dp), allocatable :: g(:,:), q(:,:), r(:,:)
    integer, allocatable :: pivots(:)
    real(dp), allocatable :: predicted_gain(:), true_gain(:)
    integer :: rank = 0
    integer :: status = KL_SUCCESS
  end type csi_result

end module kernlab_types
