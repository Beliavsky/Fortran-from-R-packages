! SPDX-License-Identifier: GPL-3.0-or-later
module rpeif_types
  use rpeif_kinds, only : dp
  implicit none
  private

  integer, parameter, public :: rpeif_success = 0
  integer, parameter, public :: rpeif_invalid_argument = 1
  integer, parameter, public :: rpeif_numerical_failure = 2
  integer, parameter, public :: rpeif_unknown_estimator = 3

  type, public :: nuisance_parameters
    real(dp) :: mu = 0.01_dp
    real(dp) :: sd = 0.05_dp
    real(dp) :: c = 0.0_dp
    real(dp) :: alpha = 0.1_dp
    real(dp) :: beta = 0.1_dp
    real(dp) :: q_alpha = 0.0_dp
    real(dp) :: es_alpha = 0.0_dp
    real(dp) :: lpm1 = 0.0_dp
    real(dp) :: lpm2 = 0.0_dp
    real(dp) :: semisd = 0.0_dp
    real(dp) :: semimean = 0.0_dp
    real(dp) :: fq_alpha = 0.0_dp
    real(dp) :: dsr = 0.0_dp
    real(dp) :: es_ratio = 0.0_dp
    real(dp) :: upm1 = 0.0_dp
    real(dp) :: omega = 0.0_dp
    real(dp) :: q_beta = 0.0_dp
    real(dp) :: eg_beta = 0.0_dp
    real(dp) :: rachev_ratio = 0.0_dp
    real(dp) :: psi_prime = 0.9038_dp
    real(dp) :: sor_c = 0.0_dp
    real(dp) :: sor_mu = 0.0_dp
    real(dp) :: sr = 0.0_dp
    real(dp) :: var_ratio = 0.0_dp
  end type nuisance_parameters

  type, public :: rpeif_options
    real(dp) :: alpha = -1.0_dp
    real(dp) :: beta = 0.1_dp
    real(dp) :: risk_free = 0.0_dp
    real(dp) :: threshold_constant = 0.0_dp
    integer :: lpm_order = 1
    character(len=8) :: sortino_threshold = 'const'
    character(len=16) :: robust_family = 'mopt'
    real(dp) :: efficiency = -1.0_dp
    logical :: clean_outliers = .false.
    logical :: prewhiten = .false.
    integer :: ar_order = 1
    logical :: source_compatibility = .true.
  end type rpeif_options

  type, public :: influence_result
    real(dp), allocatable :: x(:)
    real(dp), allocatable :: values(:)
    integer :: status = rpeif_success
    character(len=32) :: estimator = ''
    character(len=160) :: message = ''
  end type influence_result
end module rpeif_types
