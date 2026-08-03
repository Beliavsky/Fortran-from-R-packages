! SPDX-License-Identifier: GPL-3.0-only
module imputefin_types
  use imputefin_kinds, only : dp
  implicit none
  private

  integer, parameter, public :: impute_ok = 0
  integer, parameter, public :: impute_invalid_input = 1
  integer, parameter, public :: impute_insufficient_data = 2
  integer, parameter, public :: impute_singular = 3
  integer, parameter, public :: impute_not_converged = 4

  type, public :: ar1_options
    logical :: random_walk = .false.
    logical :: zero_mean = .false.
    logical :: remove_outliers = .false.
    logical :: fast_and_heuristic = .true.
    real(dp) :: outlier_prob_th = 1.0e-3_dp
    real(dp) :: tol = 1.0e-8_dp
    integer :: maxiter = 100
    integer :: n_chain = 10
    integer :: n_thin = 1
    integer :: saem_burn = 30
  end type ar1_options

  type, public :: imputation_options
    integer :: n_samples = 1
    integer :: n_burn = 100
    integer :: n_thin = 1
    integer :: rolling_window = 252
    integer(kind=8) :: seed = 5489_8
  end type imputation_options

  type, public :: ar1_fit_result
    real(dp) :: phi0 = 0.0_dp
    real(dp) :: phi1 = 0.0_dp
    real(dp) :: sigma2 = 0.0_dp
    real(dp) :: nu = huge(1.0_dp)
    logical :: converged = .false.
    integer :: iterations = 0
    integer :: status = impute_ok
    character(len=160) :: message = ''
    integer, allocatable :: index_miss(:)
    integer, allocatable :: index_outliers(:)
    real(dp), allocatable :: phi0_iterates(:)
    real(dp), allocatable :: phi1_iterates(:)
    real(dp), allocatable :: sigma2_iterates(:)
    real(dp), allocatable :: nu_iterates(:)
    real(dp), allocatable :: cond_mean(:)
    real(dp), allocatable :: cond_cov(:,:)
  end type ar1_fit_result

  type, public :: imputation_result
    real(dp), allocatable :: values(:,:,:)  ! time x series x sample
    type(ar1_fit_result), allocatable :: fits(:)
    integer :: status = impute_ok
    character(len=160) :: message = ''
  end type imputation_result

  type, public :: var_t_options
    integer :: p = 1
    logical :: omit_missing = .false.
    real(dp) :: tol = 1.0e-4_dp
    integer :: maxiter = 100
    integer :: saem_burn = 30
  end type var_t_options

  type, public :: var_t_result
    real(dp) :: nu = 6.0_dp
    real(dp), allocatable :: phi0(:)
    real(dp), allocatable :: phi(:,:,:)
    real(dp), allocatable :: scatter(:,:)
    real(dp), allocatable :: completed(:,:)
    logical :: converged = .false.
    integer :: iterations = 0
    integer :: n_used = 0
    integer :: status = impute_ok
    character(len=160) :: message = ''
  end type var_t_result

end module imputefin_types
